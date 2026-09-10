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
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>

#include "hex-arrays.h"
#include "hex-gausskronrod.h"
#include "hex-hydrogen.h"
#include "hex-misc.h"
#include "hex-special.h"

// --------------------------------------------------------------------------------- //

//
// Total plane-wave Born cross section (no exchange) for the transition
//
//     e(k_i) + H(1s)  ->  e(k_f) + H(nl) ,
//
// summed over the final magnetic sublevels "m". Both the excitations and the elastic
// scattering (n = 1, l = 0) are covered. Rydberg units are used throughout, as
// everywhere else in Hex: energies in Ry, momenta in a.u. (E = k^2), lengths in a_0
// and the resulting cross section in a_0^2.
//
// The direct first-order amplitude for the interaction  V = 1/|R - r| - 1/R  is
//
//     f(q) = -(2/q^2) [ <nlm| exp(i q.r) |1s> - <nlm|1s> ] ,     q = k_i - k_f ,
//
// the second term being the screening of the projectile-electron repulsion by the
// attraction to the nucleus. Expanding the plane wave into partial waves,
//
//     exp(i q.r) = 4 pi sum_{L M} i^L j_L(qr) Y*_{LM}(q) Y_{LM}(r) ,
//
// and using the sphericity of the initial state, only the single multipole L = l
// survives the angular integration. Because j_0(0) = 1, the overlap is precisely the
// q = 0 value of that multipole, so both terms can be kept together,
//
//     <nlm| exp(i q.r) |1s> - <nlm|1s> = sqrt(4 pi) i^l Y*_{lm}(q) M_l(q) ,
//
// with the screened multipole Born form factor
//
//     M_l(q) = int_0^infty P_nl(r) [ j_l(qr) - delta_{l0} ] P_1s(r) dr .    (*)
//
// The overlap is nonzero only for the elastic channel, where it is the contribution
// of the nucleus, without which the cross section would diverge; every excitation has
// <ns|1s> = 0 and is left untouched.
//
// The sum over "m" is then trivial, sum_m |Y_lm|^2 = (2l+1)/(4 pi), so that
//
//     dsigma/dOmega = (k_f/k_i) (4/q^4) (2l+1) M_l(q)^2 .
//
// Finally, dOmega = 2 pi q dq / (k_i k_f) turns the angular integration into an
// integration over the momentum transfer,
//
//     sigma = 8 pi (2l+1) / k_i^2  int_{|k_i - k_f|}^{k_i + k_f} M_l(q)^2 dq / q^3 .
//
// For the elastic channel the lower bound is q = 0, where the integrand vanishes as
// q^1 because M_0(q) = O(q^2) there.
//

// --------------------------------------------------------------------------------- //

/**
 * @brief Coefficients of the closed-form moment of the spherical Bessel function.
 * 
 * The integral
 * @f[
 *     \int_0^\infty r^m \mathrm{e}^{-cr} j_l(qr) \mathrm{d}r \,, \qquad m \ge l + 1 \,,
 * @f]
 * is known in closed form. The lowest power has the elementary value
 * @f[
 *     \int_0^\infty r^{l+1} \mathrm{e}^{-cr} j_l(qr) \mathrm{d}r
 *     = \frac{2^l l! q^l}{(c^2 + q^2)^{l+1}} \,,
 * @f]
 * and every further power of @f$ r @f$ follows from
 * @f$ r^{m+1} \mathrm{e}^{-cr} = -\partial_c\, r^m \mathrm{e}^{-cr} @f$. The set of
 * functions
 * @f[
 *     \sum_i a_i \frac{c^{e_i}}{(c^2 + q^2)^{l+1+i}}
 * @f]
 * is closed under that differentiation, and the exponent @f$ e_i = 2i - D @f$ is fixed
 * by the index and by the number @f$ D = m - l - 1 @f$ of derivatives taken, so the
 * whole expansion is carried in the single array of coefficients returned here.
 * 
 * @param m Power of the radial variable.
 * @param l Order of the spherical Bessel function.
 */
std::vector<double> bessel_moment_coefficients (int m, int l)
{
    // number of derivatives with respect to "c" that separate "m" from the base case
    int D = m - l - 1;

    if (D < 0)
        HexException("The moment r^%d of j_%d is not implemented (needs m >= l + 1).", m, l);

    std::vector<double> a (D + 1, 0.);
    a[0] = 1.;

    for (int d = 0; d < D; d++)
    {
        std::vector<double> b (D + 1, 0.);

        for (int i = 0; i <= d; i++)
        {
            //  -d/dc [ a c^e (c^2+q^2)^-p ] = -a e c^(e-1) (c^2+q^2)^-p + 2 a p c^(e+1) (c^2+q^2)^(-p-1)
            int p = l + 1 + i, e = 2*i - d;

            if (e > 0)
                b[i] -= a[i] * e;

            b[i + 1] += a[i] * 2. * p;
        }

        a = std::move(b);
    }

    return a;
}

/**
 * @brief Moment of the spherical Bessel function.
 * 
 * Evaluates @f$ \int_0^\infty r^m \mathrm{e}^{-cr} j_l(qr) \mathrm{d}r @f$ by summing
 * the expansion set up in @ref bessel_moment_coefficients.
 * 
 * @param m Power of the radial variable.
 * @param l Order of the spherical Bessel function.
 * @param c Exponential decay constant.
 * @param q Momentum transfer.
 */
double bessel_moment (int m, int l, double c, double q)
{
    std::vector<double> a = bessel_moment_coefficients(m, l);

    int D = m - l - 1;
    double u = c * c + q * q, sum = 0;

    for (int i = 0; i <= D; i++)
        sum += a[i] * gsl_sf_pow_int(c, 2*i - D) / gsl_sf_pow_int(u, l + 1 + i);

    return gsl_sf_pow_int(2.,l) * gsl_sf_fact(l) * gsl_sf_pow_int(q,l) * sum;
}

/**
 * @brief Screened moment of the spherical Bessel function.
 * 
 * Evaluates the monopole moment with its @f$ q = 0 @f$ value removed,
 * @f[
 *     \int_0^\infty r^m \mathrm{e}^{-cr} \left[ j_0(qr) - 1 \right] \mathrm{d}r \,.
 * @f]
 * Forming the difference term by term would be catastrophic at small @f$ q @f$, where
 * the two moments agree to many digits; that matters, because @f$ q = 0 @f$ is the very
 * end of the integration range of the elastic channel. Each term is therefore
 * subtracted analytically using
 * @f[
 *     \frac{1}{u^p} - \frac{1}{v^p} = -q^2 \sum_{s=0}^{p-1} \frac{1}{u^{p-s} v^{s+1}} \,,
 *     \qquad u = c^2 + q^2 \,, \quad v = c^2 \,,
 * @f]
 * which is exact and carries the factor @f$ q^2 @f$ explicitly, so that no cancellation
 * takes place at all.
 * 
 * @param m Power of the radial variable.
 * @param c Exponential decay constant.
 * @param q Momentum transfer.
 */
double bessel_moment_screened (int m, double c, double q)
{
    std::vector<double> a = bessel_moment_coefficients(m, 0);

    int D = m - 1;
    double u = c * c + q * q, v = c * c, sum = 0;

    for (int i = 0; i <= D; i++)
    {
        if (a[i] == 0.)
            continue;

        int p = i + 1;

        double geom = 0;
        for (int s = 0; s < p; s++)
            geom += 1. / (gsl_sf_pow_int(u, p - s) * gsl_sf_pow_int(v, s + 1));

        sum += a[i] * gsl_sf_pow_int(c, 2*i - D) * geom;
    }

    return -q * q * sum;
}

/**
 * @brief Screened multipole Born form factor for the 1s -> nl transition.
 * 
 * Evaluates the radial integral (*) above in closed form. The bound radial functions
 * @f[
 *     P_{nl}(r) = \sqrt{\left(\frac{2}{n}\right)^3 \frac{(n-l-1)!}{2n(n+l)!}}
 *                 \left(\frac{2r}{n}\right)^l \mathrm{e}^{-r/n}
 *                 L_{n-l-1}^{2l+1}\!\left(\frac{2r}{n}\right) r \,,
 *     \qquad P_{1s}(r) = 2 r \mathrm{e}^{-r} \,,
 * @f]
 * multiply to a polynomial times @f$ \exp(-cr) @f$ with @f$ c = 1 + 1/n @f$, so the
 * integral reduces to a combination of the moments evaluated above.
 * 
 * The monopole is taken in whichever of its two equivalent forms is better
 * conditioned, see the comment in the loop below. Neither the elastic channel, where
 * the screening by the nucleus is physically required, nor the excitations, where the
 * bare monopole would cancel to zero as @f$ q \to 0 @f$ by the orthogonality of the two
 * bound states, then suffer any loss of accuracy.
 * 
 * @param n Principal quantum number of the final state.
 * @param l Orbital quantum number of the final state.
 * @param q Momentum transfer.
 */
double formfactor (int n, int l, double q)
{
    if (n <= l)
        HexException("Invalid final state (%d,%d).", n, l);

    // normalization of the final bound state
    double norm = std::sqrt(gsl_sf_pow_int(2./n,3) * gsl_sf_fact(n-l-1) / (2. * n * gsl_sf_fact(n+l)));

    // decay constant of the product of the two bound states
    double c = 1. + 1./n;

    // sum over the terms of the generalized Laguerre polynomial
    double sum = 0;
    for (int j = 0; j <= n - l - 1; j++)
    {
        // coefficient of x^j in L_{n-l-1}^{2l+1}(x), rescaled from "x" to "2r/n"
        double cj = (j % 2 == 0 ? +1. : -1.) * gsl_sf_choose(n+l, n-l-1-j) / gsl_sf_fact(j);

        // The monopole is evaluated with the overlap <nl|1s> subtracted whenever that
        // is the better conditioned of the two forms. For the elastic channel the
        // subtraction is mandatory anyway. For an excitation the overlap vanishes, so
        // both forms are the same integral, and the subtracted one is used only below
        // q ~ c, where the bare terms would cancel against each other; above that the
        // roles are reversed and the bare form is the accurate one.
        double moment = (l == 0 and (n == 1 or q < c) ? bessel_moment_screened(j + 2, c, q)
                                                      : bessel_moment(j + l + 2, l, c, q));

        sum += cj * std::pow(2./n, j) * moment;
    }

    // the remaining factors are the "2" of P_1s and the angular factor of P_nl
    return 2. * norm * std::pow(2./n, l) * sum;
}

/**
 * @brief Run the Gauss-Kronrod quadrature and check that it converged.
 */
template <class Functor> double integrate (Functor f, double a, double b)
{
    GaussKronrod<Functor> Q (f);
    Q.integrate(a, b);

    if (not Q.ok())
        HexException("Failed to integrate over the momentum transfer - %s.", Q.status().c_str());

    return Q.result();
}

/**
 * @brief Total Born cross section for the 1s -> nl transition.
 * 
 * Integrates the form factor over the momentum transfer as sketched above.
 * 
 * For the excitations the interval @f$ [\,k_i - k_f, k_i + k_f\,] @f$ spans many decades
 * at high impact energies, so the substitution @f$ q = \mathrm{e}^t @f$ is used; the
 * integrand is smooth there and, for the dipole transitions, nearly constant. The
 * elastic channel has @f$ k_i = k_f @f$ and starts right at @f$ q = 0 @f$, where there
 * is no dynamic range to compress: the plain variable is used instead, in which the
 * integrand rises linearly from the origin. (Taking the logarithm anyway would leave a
 * semi-infinite interval carrying an integrand that is flat and negligible over almost
 * all of it, which exhausts the subdivision limit of the quadrature.)
 * 
 * @param n Principal quantum number of the final state.
 * @param l Orbital quantum number of the final state.
 * @param Ei Impact energy (Ry).
 */
double cross_section (int n, int l, double Ei)
{
    // final energy from the conservation law
    double Ef = Ei - 1. + 1./(n*n);

    // nothing to do below the excitation threshold (never the case for n = 1)
    if (Ef <= 0)
        return 0;

    double ki = std::sqrt(Ei), kf = std::sqrt(Ef);
    double qmin = std::abs(ki - kf), qmax = ki + kf;

    // integrand of the momentum-transfer integral
    auto dsigma = [n,l](double q) -> double
    {
        double M = formfactor(n, l, q);

        return M * M / (q * q * q);
    };

    // integrate, in the logarithmic variable whenever the interval avoids the origin
    double integral;
    if (qmin > 0)
    {
        auto logarithmic = [dsigma](double t) -> double { return dsigma(std::exp(t)) * std::exp(t); };

        integral = integrate(logarithmic, std::log(qmin), std::log(qmax));
    }
    else
    {
        integral = integrate(dsigma, 0., qmax);
    }

    return 8. * special::constant::pi * (2*l + 1) / (ki * ki) * integral;
}

// --------------------------------------------------------------------------------- //

//
// Total plane-wave Born cross section (no exchange) for the ionization
//
//     e(k_i) + H(1s)  ->  e(k_f) + p + e(kappa) ,
//
// obtained from the very same first-order amplitude as the excitations above, with
// the bound final state replaced by a continuum state of the ejected electron. The
// continuum state is normalized to the delta-function in the momentum,
//
//     P_{kappa,l}(r) = sqrt(2/pi) F_l(-1/kappa, kappa r) ,
//
//     int_0^infty P_{kappa,l}(r) P_{kappa',l}(r) dr = delta(kappa - kappa') ,
//
// so that the sum over the final states is int dkappa sum_{l m}. This is the same
// function as Hydrogen::F up to the factor 1/kappa that the latter carries as part
// of the normalization of the full three-dimensional scattering state.
//
// Nothing else changes. The screening term of the potential V = 1/|R - r| - 1/R
// again enters only through the overlap <f|1s>, which now vanishes identically by
// the orthogonality of the continuum to the bound spectrum, exactly as it does for
// every excitation; the initial state is still spherical, so a final partial wave
// "l" is again fed by the single multipole "l" alone. Therefore
//
//     M_l(q,kappa) = int_0^infty P_{kappa,l}(r) j_l(qr) P_1s(r) dr
//
// plays the role of the bound form factor (*), and the derivation carries over term
// by term, leaving the singly differential cross section
//
//     dsigma/dkappa = 8 pi / k_i^2 sum_l (2l+1)
//                     int_{|k_i-k_f|}^{k_i+k_f} M_l(q,kappa)^2 dq / q^3
//
// with the energy conservation  E_i = 1 + kappa^2 + k_f^2  fixing k_f, and finally
//
//     sigma = int_0^{sqrt(E_i - 1)} dsigma/dkappa dkappa .
//
// The ejected electron is the atomic one and the scattered electron is the
// projectile: they are distinguishable at this order, the full range of kappa is
// integrated over, and each ionization event is counted once.
//
// Unlike the bound form factor, M_l(q,kappa) has no elementary closed form, so it is
// integrated numerically. Doing that adaptively for every momentum transfer anew is
// ruinously slow, because each evaluation of the integrand costs a Coulomb function.
// Instead the radial integral is discretized once per (kappa, l) on a composite
// Gauss-Legendre grid whose panels resolve the fastest oscillation present, and the
// form factor is then a plain dot product for every q. The integrand is the product
// of the continuum function with P_1s, so it is damped by exp(-r) and the grid can
// stop at a modest radius. Against the closed forms above, evaluated for the bound
// transitions, the same grid reproduces the form factor to thirteen digits.
//

// --------------------------------------------------------------------------------- //

/// Radial extent of the quadrature grid; the integrand is damped by exp(-r) from P_1s.
static const double ION_RMAX = 30.;

/// Gauss-Legendre panels per wavelength of the fastest oscillation, and their order.
static const double ION_PANELS_PER_WAVE = 2.;
static const int    ION_GL_ORDER        = 16;

/// Panels x order used for the momentum-transfer and the ejected-momentum integrals.
static const int    ION_Q_PANELS        = 4;
static const int    ION_K_PANELS        = 4;

/**
 * @brief Composite Gauss-Legendre nodes and weights on an interval.
 */
void gauss_legendre_nodes
(
    double a, double b, int panels, int order,
    std::vector<double> & x, std::vector<double> & w
)
{
    x.clear(); w.clear();
    x.reserve(panels*order); w.reserve(panels*order);

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

/**
 * @brief Regular Coulomb functions F_0 ... F_lmax at a single radius.
 *
 * The array routine of GSL degrades as its highest order grows: the recursion, which
 * descends from that order, fails outright for a part of the arguments needed here and
 * loses accuracy in the rest. The remedy is the same as for the Bessel functions --
 * ask only for the orders that can carry any weight. The classical turning point of
 * the partial wave @f$ l @f$ lies at
 * @f[
 *     r_t(l) = \frac{\sqrt{1 + k^2 l(l+1)} - 1}{k^2}
 *     \qquad \Longleftrightarrow \qquad
 *     l(l+1) = 2r + k^2 r^2 \,,
 * @f]
 * so at a given radius every order above @c lcut below sits under the centrifugal
 * barrier with a wide margin and is set to zero. Should the recursion still fail on
 * the orders that were asked for, the scalar routine takes over; it carries its own
 * asymptotic and WKB branches, but only within the barrier are those reliable, which
 * is why it is never called above @c lcut either.
 *
 * The functions are returned in the plain Coulomb normalization, i.e. without the
 * factor sqrt(2/pi) of the momentum-normalized radial function.
 *
 * @param k Momentum of the continuum electron.
 * @param lmax Highest order requested.
 * @param r Radius.
 * @param F Output array of lmax + 1 elements.
 */
void coulomb_F_array (double k, int lmax, double r, double * F)
{
    for (int l = 0; l <= lmax; l++)
        F[l] = 0.;

    if (r <= 0. or k <= 0.)
        return;

    // highest order that is not hopelessly under the centrifugal barrier
    double lturn = std::sqrt(2.*r + k*k*r*r);
    int lcut = std::min<int>(lmax, (int)std::ceil(lturn) + 20);

    double expF = 0.;
    int err = gsl_sf_coulomb_wave_F_array(0., lcut, -1./k, k*r, F, &expF);

    bool ok = ((err == GSL_SUCCESS or err == GSL_EUNDRFLW) and std::isfinite(expF));

    if (ok)
    {
        double scale = std::exp(expF);
        for (int l = 0; l <= lcut and ok; l++)
        {
            F[l] *= scale;
            ok = std::isfinite(F[l]);
        }
    }

    if (not ok)
    {
        for (int l = 0; l <= lcut; l++)
        {
            // Hydrogen::F carries the sqrt(2/pi)/k normalization of the free state
            double v = k * Hydrogen::F(k,l,r) / special::constant::sqrt_two * special::constant::sqrt_pi;
            F[l] = (std::isfinite(v) ? v : 0.);
        }
    }

    for (int l = lcut + 1; l <= lmax; l++)
        F[l] = 0.;
}

/**
 * @brief Singly differential Born cross section for ionization.
 *
 * Returns @f$ \mathrm{d}\sigma/\mathrm{d}\kappa @f$ at the ejected momentum "kappa",
 * summed over the angular momenta of the ejected electron. The sum is truncated where
 * the highest orders no longer contribute; the estimate of the truncation point is
 * kinematic, and is raised until the last few orders are negligible, so that the
 * result does not depend on it.
 *
 * @param Ei Impact energy (Ry).
 * @param kappa Momentum of the ejected electron (a.u.).
 */
double dsigma_dkappa (double Ei, double kappa)
{
    double ki = std::sqrt(Ei);
    double Ef = Ei - 1. - kappa * kappa;

    // no room left for the scattered electron (or none for the ejected one)
    if (Ef <= 0. or kappa <= 0.)
        return 0;

    double kf = std::sqrt(Ef);
    double qmin = std::abs(ki - kf), qmax = ki + kf;

    if (qmin <= 0. or qmax <= qmin)
        return 0;

    // radial grid, resolving the fastest oscillation of the integrand
    std::vector<double> r, wr;
    int rpanels = std::max(1, (int)std::ceil(ION_RMAX * (kappa + qmax)
                                             / (2 * special::constant::pi) * ION_PANELS_PER_WAVE));
    gauss_legendre_nodes(0., ION_RMAX, rpanels, ION_GL_ORDER, r, wr);
    std::size_t N = r.size();

    // momentum transfer grid, in the logarithmic variable (dq/q^3 -> exp(-2t) dt)
    std::vector<double> t, wt;
    gauss_legendre_nodes(std::log(qmin), std::log(qmax), ION_Q_PANELS, ION_GL_ORDER, t, wt);

    // the 1s function is independent of everything that is iterated below
    std::vector<double> P1s(N);
    for (std::size_t i = 0; i < N; i++)
        P1s[i] = Hydrogen::P(1,0,r[i]);

    for (int lmax = (int)std::ceil(4. * (kappa + qmax)) + 12; ; lmax += 32)
    {
        // A[l][i] = w_i P_{kappa,l}(r_i) P_1s(r_i)
        std::vector<double> A ((lmax+1) * N), F (lmax+1);
        for (std::size_t i = 0; i < N; i++)
        {
            coulomb_F_array(kappa, lmax, r[i], F.data());

            for (int l = 0; l <= lmax; l++)
                A[l*N+i] = wr[i] * special::constant::sqrt_two / special::constant::sqrt_pi * F[l] * P1s[i];
        }

        // S[l] = int M_l(q,kappa)^2 dq / q^3
        std::vector<double> S (lmax+1, 0.), M (lmax+1), j (lmax+1);
        for (std::size_t n = 0; n < t.size(); n++)
        {
            double q = std::exp(t[n]);

            std::fill(M.begin(), M.end(), 0.);
            for (std::size_t i = 0; i < N; i++)
            {
                special::sph_jv(lmax, q * r[i], j.data());

                for (int l = 0; l <= lmax; l++)
                    M[l] += A[l*N+i] * j[l];
            }

            for (int l = 0; l <= lmax; l++)
                S[l] += wt[n] * M[l] * M[l] / (q * q);
        }

        // every term is non-negative, so the truncation error is bounded by the tail
        double sum = 0, tail = 0;
        for (int l = 0; l <= lmax; l++)
        {
            double c = (2*l + 1) * S[l];
            sum += c;
            if (l >= lmax - 2)
                tail += c;
        }

        if (sum == 0. or tail < 1e-10 * sum or lmax > 400)
            return 8. * special::constant::pi / (ki * ki) * sum;
    }
}

/**
 * @brief Total Born cross section for ionization.
 *
 * Integrates @ref dsigma_dkappa over the momentum of the ejected electron. The
 * substitution @f$ \kappa = \kappa_{\max} \sin\varphi @f$ is used: at the upper end
 * the interval of the momentum transfer closes as
 * @f$ q_{\max} - q_{\min} = 2k_f = 2\sqrt{\kappa_{\max}^2 - \kappa^2} @f$, which leaves
 * a square-root edge in @f$ \kappa @f$ that would cripple the convergence of the
 * quadrature; in @f$ \varphi @f$ the integrand is smooth at both ends and a few dozen
 * nodes exhaust it.
 *
 * @param Ei Impact energy (Ry).
 */
double cross_section_ionization (double Ei)
{
    // nothing to do below the ionization threshold
    if (Ei <= 1.)
        return 0;

    double kmax = std::sqrt(Ei - 1.);

    std::vector<double> phi, wphi;
    gauss_legendre_nodes(0., 0.5 * special::constant::pi, ION_K_PANELS, ION_GL_ORDER, phi, wphi);

    double sigma = 0;

    # pragma omp parallel for schedule (dynamic) reduction (+:sigma)
    for (std::size_t i = 0; i < phi.size(); i++)
        sigma += wphi[i] * kmax * std::cos(phi[i]) * dsigma_dkappa(Ei, kmax * std::sin(phi[i]));

    return sigma;
}

// --------------------------------------------------------------------------------- //

const std::string usage =
    "\nUsage:\n\n"
    "  hex-fullborn [--form-factor] <n> <l> [<x> ...]\n"
    "  hex-fullborn --ionization [<x> ...]\n\n"
    "Writes the total plane-wave Born cross section (no exchange) of the 1s -> nl\n"
    "transition in hydrogen, summed over the final magnetic sublevels, for every\n"
    "impact energy <x> given. The elastic channel is obtained with n = 1, l = 0.\n"
    "With \"--form-factor\" the screened multipole Born form factor M_l(q) is written\n"
    "instead, for every momentum transfer <x> given. With \"--ionization\" the total\n"
    "cross section of 1s -> continuum is written, integrated over the energy and\n"
    "summed over the angular momentum of the ejected electron; <n> and <l> are then\n"
    "not used. When no <x> is given on the command line, the values are read from\n"
    "the standard input.\n\n"
    "Energies are in Rydbergs, momenta in atomic units, cross sections in a_0^2.\n\n";

int main (int argc, char * argv[])
{
    // do not let GSL abort on underflows; they are legitimate here
    gsl_set_error_handler_off();

    bool ff = false, ion = false;
    std::vector<const char*> args;

    for (int iarg = 1; iarg < argc; iarg++)
    {
        if (std::strcmp(argv[iarg],"--form-factor") == 0)
            ff = true;
        else if (std::strcmp(argv[iarg],"--ionization") == 0)
            ion = true;
        else
            args.push_back(argv[iarg]);
    }

    if (ff and ion)
    {
        std::cerr << "The options \"--form-factor\" and \"--ionization\" are mutually exclusive." << std::endl;
        return EXIT_FAILURE;
    }

    // the ionization channel is not labelled by a final bound state
    std::size_t first = (ion ? 0 : 2);

    if (args.size() < first)
    {
        std::cout << usage;
        return EXIT_FAILURE;
    }

    int n = 0, l = 0;

    if (not ion)
    {
        n = std::atoi(args[0]);
        l = std::atoi(args[1]);

        if (n < 1 or l < 0 or l >= n)
        {
            std::cerr << "Invalid final state (" << n << "," << l << ")." << std::endl;
            return EXIT_FAILURE;
        }
    }

    std::cout << std::scientific << std::setprecision(10);

    if (ion)
    {
        std::cout << "# Plane-wave Born (no exchange) cross section, 1s -> ionization" << std::endl;
        std::cout << "# ionization threshold: " << 1. << " Ry" << std::endl;
        std::cout << "# Ei [Ry]\tsigma [a_0^2]" << std::endl;
    }
    else
    {
        std::cout << "# Plane-wave Born (no exchange) cross section, 1s -> " << Hydrogen::stateName(n,l) << std::endl;

        if (n == 1)
            std::cout << "# elastic channel, screened by the attraction to the nucleus" << std::endl;
        else
            std::cout << "# excitation threshold: " << 1. - 1./(n*n) << " Ry" << std::endl;
        std::cout << (ff ? "# q [a.u.]\tM_l(q)" : "# Ei [Ry]\tsigma [a_0^2]") << std::endl;
    }

    // collect the abscissae from the command line, or from the standard input when there are none
    rArray xs;
    if (args.size() > first)
    {
        for (std::size_t i = first; i < args.size(); i++)
            xs.push_back(std::atof(args[i]));
    }
    else
    {
        for (double x; std::cin >> x; )
            xs.push_back(x);
    }

    for (double x : xs)
    {
        double y = (ion ? cross_section_ionization(x)
                        : (ff ? formfactor(n,l,x) : cross_section(n,l,x)));

        std::cout << x << "\t" << y << std::endl;
    }

    return EXIT_SUCCESS;
}
