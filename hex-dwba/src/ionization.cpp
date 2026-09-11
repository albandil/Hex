//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2026, Jakub Benda, Charles University in Prague                    //
//                                                                                   //
// MIT License:                                                                      //
//                                                                                   //
//  Permission is hereby granted, free of charge, to any person obtaining a          //
// copy of this software and associated documentation files (the "Software"),        //
// to deal in the Software without restriction, including without limitation         //
// the rights to use, copy, modify, merge, publish, distribute, sublicense,          //
// and/or sell copies of the Software, and to permit persons to whom the             //
// Software is furnished to do so, subject to the following conditions:              //
//                                                                                   //
//  The above copyright notice and this permission notice shall be included          //
// in all copies or substantial portions of the Software.                            //
//                                                                                   //
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS          //
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       //
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE       //
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, //
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF         //
// OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.  //
//                                                                                   //
//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //

#include <algorithm>
#include <cmath>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>

#include "hex-arrays.h"
#include "hex-hydrogen.h"
#include "hex-misc.h"
#include "hex-special.h"

#include "ionization.h"

using special::constant::pi;

// --------------------------------------------------------------------------------- //

//
// The direct first-order T-matrix is assembled exactly as in "pwba" of hex-born, with
// the bound final state replaced by a continuum state of the ejected electron,
//
//     P_{kappa,la}(r) = sqrt(2/pi) F_la(-1/kappa, kappa r) ,
//
// normalized to delta(kappa - kappa'), so that the sum over the final states reads
// int dkappa sum_{la m}. Because the initial state is spherical (Li = 0) the
// Clebsch-Gordan coefficient of the entrance channel collapses to delta(li,L), and the
// Wigner symbol hidden in "computef" forces the single multipole lambda = la. What
// remains for a given L is
//
//     T(la,lf,Mf) = (4 pi)^2 / (ki kf) i^(L-lf) sqrt((2L+1)/(4 pi))
//                   * <la Mf, lf -Mf | L 0> f(la,la,lf,0,L,L) I(la,lf) ,
//
//     I(la,lf) = int jhat_L(ki r2) jhat_lf(kf r2) V_la(r2) dr2 ,
//
//     V_la(r2) = int P_{kappa,la}(r1) ( r_<^la / r_>^(la+1) - delta_la0 / r2 )
//                P_1s(r1) dr1 .
//
// The magnetic sublevels of the ejected electron are summed over, and the
// Clebsch-Gordan coefficients are normalized, sum_{Mf} <...|L 0>^2 = 1. Without
// exchange the singlet and triplet T-matrices coincide, so the spin sum of the cross
// section formula used throughout hex-dwba, kf/ki (sigma_S + 3 sigma_T) / (16 pi^2),
// contributes a plain factor of four, and everything collapses to
//
//     sigma_L = 16 pi (2L+1) / ki^3 int dkappa (1/kf)
//               sum_{la,lf} f(la,la,lf,0,L,L)^2 I(la,lf)^2 .
//
// This has been checked term by term against the momentum-transfer representation of
// hex-fullborn: for a fixed kappa and la, the sum over L of the above reproduces the
// contribution of that la to dsigma/dkappa to six significant digits.
//
// Numerically the work is arranged so that no radial function is ever evaluated twice.
// The inner integral V_la depends on kappa and la alone, not on L or lf, so it is
// tabulated once per (kappa, la) on a common quadrature grid; the products of the two
// Riccati-Bessel functions are tabulated as well, the one of the projectile once for
// the whole run. Every (L, lf) pair is then a dot product over the grid.
//
// V_la is built by accumulating the two cumulative integrals node to node, which costs
// one small Gauss-Legendre rule per interval and gives both branches of the multipole
// in a single sweep. Beyond the radius where P_1s has died out the inner integral is
// the exact power law C_la r^(-la-1), and that region needs no radial function at all,
// so it is simply carried on the same tabulated grid, out to a radius where what is
// left beyond it no longer matters.
//
// An earlier version stopped the grid short and picked the rest up with an oscillatory
// quadrature over the nodes of the beating Riccati-Bessel product. That remainder is
// worth a few parts in a million of the partial cross sections, but its parcels alternate
// and cancel almost completely, so the sum of them is far smaller than any one of them
// and no term-against-sum stopping rule ever fires; it ran to its parcel cap every time
// and, at the low impact energies where everything else is cheap, came to dominate the
// whole calculation. Reaching further out on the grid buys the same accuracy for a small
// fraction of the cost.
//

// --------------------------------------------------------------------------------- //

/// Radius beyond which the 1s function, and with it the inner integrand, has died out.
static const double RBOUND = 40.;

/// Extent of the tabulated region beyond RBOUND, in wavelengths of the projectile.
/// Past RBOUND the inner integral is the exact power law C r^(-la-1) and the integrand
/// costs nothing but the two Riccati-Bessel tables, so the region is made long enough
/// for the remainder to be negligible. Measuring it in wavelengths rather than in bohrs
/// keeps the number of quadrature nodes -- and hence the cost -- the same at every
/// impact energy, while automatically reaching furthest out at the low energies, where
/// the grid is coarse and the reach is cheap.
static const double FAR_WAVELENGTHS = 432.;

/// Bounds on the resulting radius, to keep both extremes sane.
static const double RMAX_MIN = 400.;
static const double RMAX_MAX = 3000.;

/// Gauss-Legendre panels per wavelength of the fastest oscillation, and their order.
/// Beyond RBOUND the inner integral is a smooth power law and only the two projectile
/// waves are left to resolve, so half the density is enough there; it was checked to
/// give bit-identical partial cross sections.
static const double PANELS_PER_WAVE     = 2.;
static const double FAR_PANELS_PER_WAVE = 1.;
static const int    GL_ORDER        = 16;

/// Order of the small rule that advances the cumulative inner integrals node to node.
static const int    SUB_ORDER       = 8;

/// Panels x GL_ORDER nodes are used for the integral over the ejected momentum.
static const int    K_PANELS        = 4;

/// Highest angular momentum of the ejected electron that will ever be tried.
static const int    LA_LIMIT        = 60;

// --------------------------------------------------------------------------------- //

namespace
{

/// Composite Gauss-Legendre nodes and weights, appended to the supplied arrays.
void append_nodes
(
    double a, double b, int panels, int order,
    std::vector<double> & x, std::vector<double> & w
)
{
    gsl_integration_glfixed_table * tab = gsl_integration_glfixed_table_alloc(order);
    double h = (b - a) / panels;

    for (int p = 0; p < panels; p++)
    {
        for (int i = 0; i < order; i++)
        {
            double xi, wi;
            gsl_integration_glfixed_point(a + p*h, a + (p+1)*h, i, &xi, &wi, tab);
            x.push_back(xi); w.push_back(wi);
        }
    }

    gsl_integration_glfixed_table_free(tab);
}

/// Number of panels needed to resolve the wavenumber "k" over a stretch of length "len".
int panels_for (double len, double k, double per_wave = PANELS_PER_WAVE)
{
    return std::max(1, (int)std::ceil(len * std::max(k,1e-3) / (2*pi) * per_wave));
}

} // anonymous namespace

// --------------------------------------------------------------------------------- //

void pwba_ionization
(
    double Ei,
    int Lmin, int Lmax,
    rArray & sigma,
    int & lamax,
    double eps
)
{
    sigma = rArray(Lmax - Lmin + 1, 0.);
    lamax = 0;

    if (Ei <= 1.)
        return;

    double ki = std::sqrt(Ei);
    double kappa_max = std::sqrt(Ei - 1.);

    //
    // common radial grid: resolved on [0,RBOUND], where the inner integrand lives and
    // oscillates with the ejected momentum as well, and coarser on [RBOUND,RMAX],
    // where only the two projectile waves still oscillate
    //

    double RMAX = std::min(RMAX_MAX, std::max(RMAX_MIN, RBOUND + FAR_WAVELENGTHS * 2*pi / ki));

    std::vector<double> r, w;
    append_nodes(0., RBOUND, panels_for(RBOUND, kappa_max + 2*ki), GL_ORDER, r, w);
    std::size_t Nb = r.size();
    append_nodes(RBOUND, RMAX, panels_for(RMAX - RBOUND, 2*ki, FAR_PANELS_PER_WAVE), GL_ORDER, r, w);
    std::size_t N = r.size();

    // the 1s function and the projectile waves do not depend on anything iterated below
    std::vector<double> P1s (Nb);
    for (std::size_t i = 0; i < Nb; i++)
        P1s[i] = Hydrogen::P(1,0,r[i]);

    // Riccati-Bessel functions of the incoming projectile, jhat_L(ki r) = ki r j_L(ki r).
    // The whole column of orders is taken in one sweep: evaluating them one by one would
    // start a separate continued fraction for each order and cost O(Lmax^2) per radius.
    std::vector<double> jL ((std::size_t)(Lmax+1) * N);
    {
        std::vector<double> col (Lmax+1);
        for (std::size_t i = 0; i < N; i++)
        {
            special::sph_jv(Lmax, ki*r[i], col.data());

            for (int L = 0; L <= Lmax; L++)
                jL[(std::size_t)L*N + i] = ki * r[i] * col[L];
        }
    }

    // nodes of the integral over the ejected momentum, kappa = kappa_max sin(phi)
    std::vector<double> phi, wphi;
    append_nodes(0., 0.5*pi, K_PANELS, GL_ORDER, phi, wphi);

    std::vector<double> total (Lmax - Lmin + 1, 0.);
    int lamax_used = 0;

    # pragma omp parallel
    {
        std::vector<double> mine (Lmax - Lmin + 1, 0.);
        int mine_lamax = 0;

        # pragma omp for schedule (dynamic)
        for (std::size_t ik = 0; ik < phi.size(); ik++)
        {
            double kappa = kappa_max * std::sin(phi[ik]);
            double jac   = wphi[ik] * kappa_max * std::cos(phi[ik]);
            double Ef    = Ei - 1. - kappa*kappa;

            if (kappa <= 0. or Ef <= 0.)
                continue;

            double kf = std::sqrt(Ef);

            // Riccati-Bessel functions of the scattered projectile. The multipole loop
            // below stops on its own accord, usually far below LA_LIMIT, and only orders
            // up to Lmax + la are ever read; the table is therefore grown to match as the
            // multipole rises, instead of being filled once for the worst case.
            std::vector<double> jf;
            int lf_top = -1;

            auto grow_jf = [&](int need)
            {
                if (need <= lf_top)
                    return;

                // Growing costs one sweep over the whole grid, and the array routine
                // returns every order below the top anyway, so the rows already held
                // cannot be kept and extended cheaply. Grow geometrically instead, to
                // keep the number of sweeps down to a couple per ejected momentum.
                lf_top = std::min(std::max(need + 12, 2*lf_top + 1), Lmax + LA_LIMIT);
                jf.assign((std::size_t)(lf_top+1) * N, 0.);

                std::vector<double> col (lf_top+1);
                for (std::size_t i = 0; i < N; i++)
                {
                    special::sph_jv(lf_top, kf*r[i], col.data());

                    for (int lf = 0; lf <= lf_top; lf++)
                        jf[(std::size_t)lf*N + i] = kf * r[i] * col[lf];
                }
            };

            std::vector<double> V (N), dA (Nb), dB (Nb), dB0 (Nb), A (Nb), B (Nb), B0 (Nb);
            gsl_integration_glfixed_table * tab = gsl_integration_glfixed_table_alloc(SUB_ORDER);

            double running = 0;

            for (int la = 0; la <= LA_LIMIT; la++)
            {
                //
                // inner multipole integral V_la, accumulated node to node
                //

                for (std::size_t i = 0; i < Nb; i++)
                {
                    double a = (i == 0 ? 0. : r[i-1]), b = r[i];
                    dA[i] = dB[i] = dB0[i] = 0.;

                    for (int s = 0; s < SUB_ORDER; s++)
                    {
                        double x, ww;
                        gsl_integration_glfixed_point(a, b, s, &x, &ww, tab);

                        double g = kappa * Hydrogen::F(kappa,la,x) * Hydrogen::P(1,0,x);

                        dA[i]  += ww * g * std::pow(x, (double)la);      // r1 < r2 branch
                        dB[i]  += ww * g * std::pow(x, -la-1.);          // r1 > r2 branch
                        dB0[i] += ww * g;                                // screening of the monopole
                    }
                }

                double acc = 0;
                for (std::size_t i = 0; i < Nb; i++) { acc += dA[i]; A[i] = acc; }

                // dB spans (r_{i-1}, r_i], so the running sum has to be taken before the
                // store for B[i] to be the integral from r_i upwards, not from r_{i-1}
                acc = 0;
                for (std::size_t i = Nb; i-- > 0; ) { B[i]  = acc; acc += dB[i];  }
                acc = 0;
                for (std::size_t i = Nb; i-- > 0; ) { B0[i] = acc; acc += dB0[i]; }

                double C = A[Nb-1];    // int_0^inf P_kappa,la P_1s r^la dr

                for (std::size_t i = 0; i < Nb; i++)
                {
                    V[i] = (la == 0 ? B[i] - B0[i]/r[i]
                                    : A[i]*std::pow(r[i],-la-1.) + B[i]*std::pow(r[i],(double)la));
                }
                for (std::size_t i = Nb; i < N; i++)
                {
                    // beyond RBOUND only the exact power-law tail is left
                    V[i] = (la == 0 ? 0. : C * std::pow(r[i], -la-1.));
                }

                //
                // outer integral and the cross section
                //

                grow_jf(Lmax + la);

                double added = 0;

                for (int L = Lmin; L <= Lmax; L++)
                {
                    double partial = 0;

                    for (int lf = std::abs(L - la); lf <= L + la; lf++)
                    {
                        if ((L + la + lf) % 2 != 0)
                            continue;

                        double f = special::computef(la, la, lf, 0, L, L);

                        if (f == 0. or not std::isfinite(f))
                            continue;

                        const double * jLp = &jL[(std::size_t)L*N];
                        const double * jfp = &jf[(std::size_t)lf*N];

                        double I = 0;
                        for (std::size_t i = 0; i < N; i++)
                            I += w[i] * jLp[i] * jfp[i] * V[i];

                        partial += f * f * I * I;
                    }

                    double contrib = jac * 16.*pi*(2*L + 1) / (ki*ki*ki) / kf * partial;

                    mine[L - Lmin] += contrib;
                    added += contrib;
                }

                running += added;
                mine_lamax = std::max(mine_lamax, la);

                // the terms are non-negative, so the neglected rest is bounded by the last one
                if (la > 2 and added < eps * running)
                    break;
            }

            gsl_integration_glfixed_table_free(tab);
        }

        # pragma omp critical
        {
            for (std::size_t i = 0; i < mine.size(); i++)
                total[i] += mine[i];
            lamax_used = std::max(lamax_used, mine_lamax);
        }
    }

    for (std::size_t i = 0; i < total.size(); i++)
        sigma[i] = total[i];

    lamax = lamax_used;
}
