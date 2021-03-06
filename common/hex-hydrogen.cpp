//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2017, Jakub Benda, Charles University in Prague                    //
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

#include <sstream>

#include <cstdlib>
#include <cmath>
#include <cstddef>

#include <gsl/gsl_sf.h>

#include "hex-hydrogen.h"
#include "hex-misc.h"
#include "hex-special.h"

namespace Hydrogen
{

CarthesianBoundWaveFunction::CarthesianBoundWaveFunction (int N, int L, int M)
    : N_(N), L_(L), M_(M)
{
    // compute product of norms (radial × angular)
    norm_ = std::sqrt
    (
        std::pow(2./N,3) * gsl_sf_fact(N-L-1) / (2 * N * gsl_sf_fact(N+L))
        * (2*L+1) / (4 * special::constant::pi)
        * gsl_sf_fact(L - std::abs(M)) / gsl_sf_fact(L + std::abs(M))
    ) * std::pow(1./N,L);

    // for all terms of the Laguerre polynomial
    for (int i = 0; i <= N - L - 1; i++)
    {
        Real Nfactor = (i % 2 == 0 ? 1. : -1.) * std::pow(2./N,i) * gsl_sf_choose(N+L,N-L-1-i) / gsl_sf_fact(i);

        // for all terms of the z-dependent angular factor
        for (int k = 0; k <= (L - std::abs(M))/2; k++)
        {
            Real Lfactor = (k % 2 == 0 ? 1. : -1.) * gsl_sf_choose(L,k) * gsl_sf_choose(2*L-2*k,L)
                                * gsl_sf_fact(L-2*k) / gsl_sf_fact(L-2*k-M);

            // for all terms of the x,y-dependent angular factor
            for (int p = 0; p <= std::abs(M); p++)
            {
                Complex Mfactor = gsl_sf_choose(std::abs(M),p);

                if (M > 0) // (-1)^m * (Am + i*Bm)
                {
                    Mfactor *= (std::abs(M) % 2 == 0 ? 1.0_r : -1.0_r) * std::pow(Complex(0.,1.),std::abs(M)-p);
                }

                if (M < 0) // Am - i * Bm
                {
                    Mfactor *= std::pow(Complex(0.,1.),p-std::abs(M));
                }

                std::cout << "i = " << i << ", k = " << k << ", p = " << p << std::endl;

                // append new term to the list of terms
                Term term;
                term.c = Nfactor * Lfactor * Mfactor;
                term.n = i + 2*k;
                term.u = p;
                term.v = std::abs(M) - p;
                term.w = L - 2*k - std::abs(M);
                terms_.push_back(term);
            }
        }
    }
}

CarthesianBoundWaveFunction::~CarthesianBoundWaveFunction ()
{

}

Complex CarthesianBoundWaveFunction::operator() (double x, double y, double z) const
{
    Complex result = 0;
    Real r = std::sqrt(x*x+y*y+z*z);

    for (Term const & T : terms_)
    {
        Complex term = T.c;

        if (T.u != 0) term *= std::pow(x,T.u);
        if (T.v != 0) term *= std::pow(y,T.v);
        if (T.w != 0) term *= std::pow(z,T.w);
        if (T.n != 0) term *= std::pow(r,T.n);

        result += term;
    }

    return result * std::exp(-r / N_);
}

std::string stateName (int n, int l)
{
    static const std::string ang = "spdfghiklmno";

    std::ostringstream oss;

    if (n == 0)
        oss << "ion";
    else if (l < (int)ang.size())
        oss << n << ang[l];
    else
        oss << n << '[' << l << ']';

    return oss.str();
}

std::string stateName (int n, int l, int m)
{
    std::string name = stateName(n,l);

    if (n != 0)
        name += "(" + std::to_string(m) + ")";

    return name;
}

std::string stateRName (int n, int l, int two_j, int two_m)
{
    static const char ang[] = { 's', 'p', 'd', 'f', 'g', 'h', 'i', 'k', 'l', 'm', 'n', 'o' };

    assert(l < (int)sizeof(ang));

    std::ostringstream oss;
    if (n != 0)
        oss << n << ang[l] << two_j << "/2(" << two_m << "/2)";
    else
        oss << "ion";

    return oss.str();
}

double lastZeroBound (int n, int l)
{
    // all Laguerre(r) roots are less or equal to
    return (n + l + (n - l - 2.) * sqrt(n + l));
}

double getBoundFar (int n, int l, double eps, int max_steps)
{
    // all Laguerre(2r/n) roots are less or equal to
    double last_zero_ubound = lastZeroBound(n, l) * 0.5 * n;

    // hunt for low value
    double far_r = std::max(last_zero_ubound, 1.);
    double far_val, der;
    do
    {
        // move further
        far_r *= 2;

        // evaluate function in 'far'
        far_val = fabs(P(n,l,far_r));

        // compute forward h-diference in 'far'
        der = far_val - std::fabs(P(n,l,far_r+1e-5));
    }
    while (der > 0 or far_val > eps);

    // bisect for exact value
    int steps = 0;
    double near_r = last_zero_ubound;
    double near_val = std::fabs(P(n,l,near_r));
    while (steps < max_steps)
    {
        double middle = 0.5 * (near_r + far_r);
        double middle_val = std::fabs(P(n,l,middle));
        if (middle_val > eps)
        {
            near_r = middle;
            near_val = middle_val;
        }
        else
        {
            far_r = middle;
            far_val = middle_val;
        }
        steps++;
        if (near_val == far_val)
            break;
    }

    // return bisection
    return 0.5 * (near_r + far_r);
}

double Norm (int n, int l)
{
    return std::sqrt(gsl_sf_pow_int(2./n,3) * gsl_sf_fact(n-l-1) / (2*n*gsl_sf_fact(n+l))) * gsl_sf_pow_int(2./n,l);
}

double getSturmFar (int n, int l, double lambda, double eps, int max_steps)
{
    // all Laguerre(2r) roots are less or equal to
    double last_zero_ubound = lastZeroBound(n,l) * 0.5;

    // hunt for low value
    double far_r = std::max(last_zero_ubound, 1.);
    double far_val, der;
    do
    {
        // move further
        far_r *= 2;

        // evaluate function in 'far'
        far_val = std::fabs(S(n,l,far_r,lambda));

        // compute forward h-diference in 'far'
        der = std::fabs(S(n,l,far_r+0.001,lambda)) - far_val;
    }
    while (der > 0 or far_val > eps);

    // bisect for exact value
    int steps = 0;
    double near_r = last_zero_ubound;
    double near_val = std::fabs(S(n,l,near_r,lambda));
    while (steps < max_steps)
    {
        double middle = 0.5 * (near_r + far_r);
        double middle_val = std::fabs(S(n,l,middle,lambda));
        if (middle_val > eps)
        {
            near_r = middle;
            near_val = middle_val;
        }
        else
        {
            far_r = middle;
            far_val = middle_val;
        }
        steps++;
        if (near_val == far_val)
            break;
    }

    // return bisection
    return 0.5 * (near_r + far_r);
}

double P (unsigned n, unsigned l, double r, double Z)
{
    // bound orbital is a regular solution
    if (r == 0.)
        return 0.;

    // compute the radial function by GSL
    gsl_sf_result R;
    int err = gsl_sf_hydrogenicR_e(n, l, Z, r, &R);

    // silent underflows -> zero
    if (err == GSL_EUNDRFLW)
    {
        return 0.;
    }

    // abort on other errors
    if (err != GSL_SUCCESS)
    {
        HexException
        (
            "Unable to evaluate the hydrogen radial function for Z = %g, n = %d, l = %d, r = %g (\"%s\").",
            Z, n, l, r, gsl_strerror(err)
        );
    }

    // return the radial orbital (multiply by radius)
    return r * R.val;
}

double F (double k, int l, double r, double sigma)
{
    // Coulomb wave is a regular solution
    if (r == 0.)
        return 0.;

    // classical turning point
    double rt = (std::sqrt(1 + k*k*l*(l+1))-1)/(k*k);

    // normalization
    double norm = special::constant::sqrt_two / (k * special::constant::sqrt_pi);

    // some local variables
    double F, exp_F, Fp;
    int err;

    // evaluate the function
    if ((err = gsl_sf_coulomb_wave_F_array(l, 0, -1./k, k*r, &F, &exp_F)) == GSL_SUCCESS)
        return norm * F;

    // high distances => use asymptotic form
    if (r > rt)
        return norm * special::coul_F_asy(l, k, r, (std::isfinite(sigma) ? sigma : special::coul_F_sigma(l,k)));

    // low energies => use uniform WKB approximation by Michel
    if ((err = special::coul_F_michel(l, k, r, F, Fp)) == GSL_SUCCESS)
        return norm * F;

    // some other problem
    HexException
    (
        "Evaluation of hydrogen free state failed for l = %d, k = %g, r = %g\nError: %d %s\n",
        l, k, r, err, gsl_strerror(err)
    );
}

double evalFreeStatePhase (double k, int l, double sigma)
{
    return l * 0.5 * special::constant::pi + (std::isfinite(sigma) ? sigma : special::coul_F_sigma(l,k));
}

double S (int n, int l, double r, double lambda)
{
    return n * n * std::pow(lambda,l+1) * std::sqrt(std::pow(2./n,3)*gsl_sf_fact(n-l-1.)/(2.*n*gsl_sf_fact(n+l))) 
        * std::exp(-lambda*r) * std::pow(2.*r,l) * gsl_sf_laguerre_n(n-l-1,2*l+1,2*r);
}

double evalFreeState_asy (double k, int l, double r, double sigma)
{
    return special::coul_F_asy(l,k,r,sigma);
}

double getFreeAsyZero (double k, int l, double Sigma, double eps, int max_steps, int nzero)
{
    // find solution (r) of
    //    n*pi = k*r - l*pi/2 + log(2*k*r)/k + Sigma
    // using trivial Banach contraction

    double kr = (2*nzero+l)*special::constant::pi_half - Sigma;
    if (kr < 0)
        return kr;

    while (max_steps-- > 0 and std::abs(kr) > eps)
        kr = (2*nzero+l)*special::constant::pi_half - Sigma - std::log(2*kr)/k;

    return kr / k;
}

double getFreeAsyTop (double k, int l, double Sigma, double eps, int max_steps, int ntop)
{
    // find solution (r) of
    //    (n+1/4)*2*pi = k*r - l*pi/2 + log(2*k*r)/k + Sigma
    // using trivial Banach contraction

    double kr = (2*ntop+0.5*(l+1))*special::constant::pi - Sigma;
    if (kr < 0)
        return kr;

    while (max_steps-- > 0 and std::abs(kr) > eps)
        kr = (2*ntop+0.5*(l+1))*special::constant::pi - Sigma - std::log(2*kr)/k;

    return kr / k;
}

double getFreeFar (double k, int l, double Sigma, double eps, int max_steps)
{
    // precompute Coulomb shift
    Sigma = std::isfinite(Sigma) ? Sigma : special::coul_F_sigma(l,k);

    //
    // hunt phase
    //

    int idx;
    double rzero = 0., eval = 0.;

    // loop over sine zeros
    int limit = max_steps;
    for (idx = 1; limit > 0; idx *= 2, limit--)
    {
        // get sine zero for current index
        rzero = getFreeAsyZero(k,l,Sigma,eps/100,max_steps,idx);
        if (rzero < 0)
            continue;

        // evaluate Coulomb wave
        eval = std::abs(F(k,l,rzero,Sigma));

        // terminate if requested precision met
        if (eval < eps)
            break;
    }

    // if a boundary was chosen, exit
    if (idx == 1 or limit == 0)
        return std::max(0., rzero);

    //
    // bisect phase
    //

    int idx_left = idx/2;
    int idx_right = idx;
    int idx_mid = (idx_right - idx_left) / 2;
    while (idx_right - idx_mid > 2)
    {
        // get sine zero for current index
        rzero = getFreeAsyZero(k,l,Sigma,eps/100,max_steps,idx_mid);

        // evaluate Coulomb wave
        eval = std::abs(F(k,l,rzero,Sigma));

        // move bisection guardians
        if (eval > eps)
            idx_left = idx_mid;
        else
            idx_right = idx_mid;

        idx_mid = (idx_right + idx_left) / 2;
    }

    return idx_mid * M_PI;
}

} // endof namespace Hydrogen

double HydrogenFunction::operator() (double r) const
{
    // if too near, return simplified value
//     if (r < 1e-5)
//         return Hydrogen::getBoundN(n,l) * gsl_sf_pow_int(r,l+1);

    // else compute precise value
    return Hydrogen::P(n_, l_, r);
}

double HydrogenFunction::getTurningPoint () const
{
    // bound state
//     if (n_ != 0)
        return n_ * (n_ - std::sqrt(n_*n_ - l_*(l_+1.)));

    // free state
//         return sqrt(l*(l+1))/k;
}
