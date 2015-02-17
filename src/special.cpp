//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2015, Jakub Benda, Charles University in Prague                    //
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

#include <cmath>
#include <iostream>

#include <gsl/gsl_sf.h>

#include "arrays.h"
#include "special.h"

#ifndef NO_LAPACK
/**
 * @brief Lapack routine: Eigenvalues, eigenvectors of a symmetric tridiagonal matrix.
 * 
 * Eigenvalues, eigenvectors of a symmetric tridiagonal matrix.
 * 
 * DSTEV computes the eigenvalues and, optionally, the left and/or right
 * eigenvectors of other matrices.
 * 
 * @param jobz 'N' for eigenvalues only, 'V' also eigenvectors.
 * @param n Order of the matrix.
 * @param d On entry, the N diagonal elements of the tridiagonal matrix A.
 *          On exit, if INFO = 0, the eigenvalues in ascending order.
 * @param e On entry, the (n-1) subdiagonal elements of the tridiagonal
 *          matrix A, stored in elements 1 to N-1 of E. On exit, the contents
 *          of E are destroyed.
 * @param z If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
 *          eigenvectors of the matrix A, with the i-th column of Z
 *          holding the eigenvector associated with D(i).
 *          If JOBZ = 'N', then Z is not referenced.
 * @param ldz The leading dimension of the array Z.  LDZ >= 1, and if
 *            JOBZ = 'V', LDZ >= max(1,N).
 * @param work If JOBZ = 'N', WORK is not referenced.
 * @param info = 0:  successful exit; 
 *             < 0:  if INFO = -i, the i-th argument had an illegal value; 
 *             > 0:  if INFO = i, the algorithm failed to converge; i
 *                   off-diagonal elements of E did not converge to zero.
 */
extern "C" void dstev_
(
    char * jobz,
    int * n,
    double * d,
    double * e,
    double * z,
    int * ldz,
    double * work,
    int * info
);

int special::coulomb_zeros (double eta, int L, int nzeros, double * zeros, double epsrel)
{
    // auxiliary variables
    int info, ldz = 1;
    char jobz = 'N';
    
    // memory array
    double oldzeros[nzeros];
    std::memset(oldzeros, 0, sizeof(oldzeros));
    
    // for different sizes of the tridiagonal matrix
    for (int n = nzeros; ; n *= 2)
    {
        // compose diagonal
        double d[n];
        for (int i = 1; i <= n; i++)
        {
            // compute diagonal element
            d[i-1] = -eta / ((L+i)*(L+i+1));
            
            // if the new element is zero, terminate and assume convergence
            if (not std::isfinite(d[i-1]) or d[i-1] == 0.)
            {
                std::cerr << "Warning: The Coulomb nodes are approximate." << std::endl;
                return n / 2;
            }
        }
        
        // compose subdiagonal
        double e[n-1];
        for (int i = 1; i <= n-1; i++)
            e[i-1] = std::sqrt((L+i+1)*(L+i+1)+eta*eta) / ((L+i+1)*std::sqrt((2*(L+i)+1)*(2*(L+i)+3)));
        
        // calculate eigenvalues
        dstev_(&jobz, &n, &d[0], &e[0], nullptr, &ldz, nullptr, &info);
        
        // check status information
        if (info < 0)
            HexException("Illegal value to DSTEV in coulomb_zeros (argument %d).", -info);
        if (info > 0)
            HexException("DSTEV failed to converge (%d offdiagonal elements).", info);
        
        // compute new zeros
        for (int i = 0; i < nzeros; i++)
        {
            oldzeros[i] = zeros[i];
            zeros[i] = 1./d[n-1-i];
        }
        
        // check that the last zero's shift is within tolerance
        if (std::abs(zeros[nzeros-1]-oldzeros[nzeros-1]) < epsrel * std::abs(zeros[nzeros-1]))
            return n;
    }
}
#endif /* NO_LAPACK */

Complex special::cfgamma (Complex s, Complex z)
{
    unsigned n = 0, k = 0;
    double epsrel = 1e-15;
    unsigned maxiter = 10000;
    
    // n = 0
    Complex App = 0, Bpp = 1;
    n++;
    
    // n = 1
    Complex Ap = 1, Bp = z;
    Complex cfp = Ap / Bp;
    n++; k++;
    
    // n = 2, 3, 4, ...
    Complex a, b, A, B, cf;
    do
    {
        // n = 2k
        a = Complex(k) - s;
        b = Complex(1);
        A = b * Ap + a * App;
        B = b * Bp + a * Bpp;
        App = Ap; Ap = A;
        Bpp = Bp; Bp = B;
        n++;
        
        // n = 2k + 1
        a = Complex(k);
        b = z;
        A = b * Ap + a * App;
        B = b * Bp + a * Bpp;
        cfp = cf; cf = A / B;
        App = Ap; Ap = A;
        Bpp = Bp; Bp = B;
        n++; k++;
    }
    while (n < maxiter and std::abs(cf - cfp) > epsrel * std::abs(cf));
    
    // scale by z^s e^(-z)
    return cf * std::pow(z,s) * std::exp(-z);
}

cArray special::ric_jv (int lmax, Complex z)
{
    // results
    cArray eval(lmax+1);
    
    // use library routine for pure real arguments
    if (z.imag() == 0.)
    {
        rArray ev(lmax+1);
        int err = gsl_sf_bessel_jl_steed_array(lmax, z.real(), &ev[0]);
        
        // check that all evaluations are finite
        bool all_finite = std::all_of
        (
            ev.begin(),
            ev.end(),
            [](double x) -> bool { return std::isfinite(x); }
        );
        
        // stop at failure
        if (err != GSL_SUCCESS or not all_finite)
        {
            if (z.real() < lmax)
            {
                // most probably due underflow in the classically forbidden region
                // -- occurs only for high angular momenta and small radii
                // => use zero-asymptotic form DLMF §10.52.1 jn(z) ~ z^n / (2n + 1)!!
                ev[0] = 1;
                for (int n = 1; n <= lmax; n++)
                {
                    ev[n] = ev[n-1] * z.real() / (2*n+1);
                }
            }
            else
            {
                HexException("Error %d while evaluating j[l≤%d](%g+%gi).", err, lmax, z.real(), z.imag());
            }
        }
        
        // Bessel -> Riccati-Bessel function
        for (int i = 0; i <= lmax; i++)
            eval[i] = z * Complex(ev[i]);
        
        return eval;
    }
    
    // shorthand
    Complex inv_z = 1. / z;
    
    // evaluate all angular momenta up to lmax
    for (int l = 0; l <= lmax; l++)
    {
        if (l == 0)
            eval[l] = std::sin(z);
        else if (l == 1)
            eval[l] = std::sin(z) * inv_z - std::cos(z);
        else
            eval[l] = Complex(2.*l - 1.) * eval[l-1] * inv_z - eval[l-2];
    }
    
    return eval;
}

Complex special::ric_j (int l, Complex z)
{
    return special::ric_jv(l, z).back();
}

cArray special::dric_jv (int lmax, Complex z)
{
    // evaluate first the Riccati-Bessel functions themselves
    cArray eval = ric_jv(lmax, z);
    
    // results
    cArray deval(lmax+1);
    
    // shorthand
    Complex inv_z = Complex(1.)/z;
    
    // evaluate all angular momenta up to lmax
    for (int l = 0; l <= lmax; l++)
    {
        if (l == 0)
            deval[l] = std::cos(z);
        else if (l == 1)
            deval[l] = inv_z * ( std::cos(z) - std::sin(z) * inv_z ) + std::sin(z);
        else
            deval[l] = -Complex(2.*l - 1) * (eval[l-1] * inv_z - deval[l-1] ) * inv_z - deval[l-2];
    }
    
    return deval;
}

Complex special::dric_j (int l, Complex z)
{
    return special::dric_jv(l, z).back();
}

// Slater-type-orbital data for hydrogen

static double a10[] = {2.};
static double a20[] = {.7071067811865475,-.3535533905932737};
static double a21[] = {.2041241452319315};
static double a30[] = {.3849001794597506,-.25660011963983370,.02851112440442597};
static double a31[] = {.12096245643373720,-.02016040940562287};
static double a32[] = {.00901600917703977};
static double a40[] = {.25,-.1875,.03125,-.001302083333333333};
static double a41[] = {.08068715304598784,-.02017178826149696,.001008589413074848};
static double a42[] = {.006987712429686843,-5.823093691405702E-4};
static double a43[] = {2.2009225383555117E-4};
static double a50[] = {.1788854381999831,-.1431083505599865,.02862167011199729,-.001908111340799819,3.8162226815996353E-5};
static double a51[] = {.05842373946721772,-.01752712184016532,.001402169747213225,-3.115932771584945E-5};
static double a52[] = {.005354624169818084,-7.139498893090778E-4,2.0398568265973652E-5};
static double a53[] = {2.039856826597365E-4,-1.0199284132986826E-5};
static double a54[] = {3.3997613776622754E-6};

static double* ak[6][5] = {
    {   0,  0,  0,   0,   0   },  // n = 0 gives only zeros
    { a10,  0,  0,   0,   0   },
    { a20, a21, 0,   0,   0   },
    { a30, a31, a32, 0,   0   },
    { a40, a41, a42, a43, 0   },
    { a50, a51, a52, a53, a54 }
};

static unsigned max_table_n = 5;

Complex hydro_P_table (unsigned n, unsigned l, Complex z)
{
    // slater-type poly term count
    int terms = n - l;
    
    // get the coefficients
    const double* const a = ak[n][l];
    
    // compute the sum
    Complex sum = 0;
    for (int i = terms - 1; i >= 0; i--)
        sum = a[i] + z * sum;
    
    // return the result
    return sum * std::pow(z, l + 1) * std::exp(-z/double(n));
}

Complex dhydro_P_table (unsigned n, unsigned l, Complex z)
{
    // slater-type poly term count
    int terms = n - l;
    
    // get the coefficients
    const double* a = ak[n][l];
    
    // compute the sum
    Complex sum = 0;
    for (int i = terms - 1; i >= 0; i--)
        sum = (l+i+1)*a[i] + z * sum;
    
    // return the result
    return sum * std::pow(z, l) * std::exp(-z/double(n)) - hydro_P_table(n,l,z) / double(n);
}

/*
 * Laguerre polynomial
 * Laguerre(k,s,x) := sum((-1)^j * (k!)^2 * x^(j-s) / ( (k-j)! * j! * (j-s)! ), j, s, k);
 */
Complex associated_laguerre_poly (int k, int s, Complex z)
{
    // value of the polynomial to be returned
    Complex val = 0;
    
    // begin with highest order
    val += std::pow(-1,k) / (fac(k) * fac(k-s));
    
    // continue with other orders
    for (int j = k - 1; j >= s; j--)
        val = z * val + std::pow(-1,j) / (fac(k-j) * fac(j) * fac(j-s));
    
    return val * std::pow(fac(k),2);
}

/*
 * Derivative of Laguerre polynomial
 * DLaguerre(k,s,x) := sum((-1)^j * (k!)^2 * x^(j-s-1) / ( (k-j)! * j! * (j-s-1)! ), j, s+1, k);
 */
Complex der_associated_laguerre_poly (int k, int s, Complex z)
{
    // value of the polynomial to be returned
    Complex val = 0;
    
    // begin with highest order
    val += std::pow(-1,k) / (fac(k) * fac(k-s-1));
    
    // continue with other orders
    for (int j = k - 1; j >= s + 1; j--)
        val = z * val + std::pow(-1,j) / (fac(k-j) * fac(j) * fac(j-s-1));
    
    return val * std::pow(fac(k),2);
}

/*
 * Hydrogen radial function normalization factor
 * sqrt((2/n)^3 * (n-l-1)! / (2*n*((n+l)!)^3));
 */
double hydrogen_wfn_normalization (int n, int l)
{
    return std::sqrt(std::pow(2./n,3) * fac(n-l-1) / (2*n*std::pow(fac(n+l),3)));
}

/*
 * Hydrogen radial function.
 * HydrogenP(n,l,r) := r * HydrogenN(n,l) * (2*r/n)^l * Laguerre(n+l,2*l+1,2*r/n) * exp(-r/n);
 */
Complex special::hydro_P (unsigned n, unsigned l, Complex z)
{
    // this is faster
    if (n <= max_table_n)
        return hydro_P_table(n, l, z);
    
    // this is general
    double Norm = hydrogen_wfn_normalization(n, l);
    Complex Lag = associated_laguerre_poly(n + l, 2 * l + 1, z);
    return z * Norm * std::pow(2.*z/double(n), l) * Lag * std::exp(-z/double(n));
}

/*
 * Derivative of hydrogen radial function.
 * DP(n,l,r) := HydrogenN(n,l) * (2*r/n)^l * (
 *   (l+1-z/n) * Laguerre(n+l,2*l+1,2*r/n) +
 *    z * DLaguerre(n+l,2*l+1,2*r/n)
 *    ) * exp(-r/n)
 */
Complex special::dhydro_P (unsigned n, unsigned l, Complex z)
{
    // this is faster
    if (n <= max_table_n)
        return dhydro_P_table(n, l, z);
    
    // this is general
    double Norm = hydrogen_wfn_normalization(n, l);
    Complex Lag = associated_laguerre_poly(n + l, 2 * l + 1, z);
    Complex DLag= der_associated_laguerre_poly(n + 1, 2 * l + 1, z);
    return Norm * std::pow(2.*z/double(n), l) * ((l + 1. - z/double(n)) * Lag + z * DLag) * std::exp(-z/double(n));
}


Complex special::sphY (int l, int m, double theta, double phi)
{
    if (l < std::abs(m))
        return 0.;
    
    return gsl_sf_legendre_sphPlm(l,std::abs(m),std::cos(theta))
           * Complex(std::cos(m*phi),std::sin(m*phi));
}

void clipang (double & theta, double & phi)
{
    // clip theta to (0,inf)
    if (theta < 0.)
    {
        // rotate plane
        phi += special::constant::pi;
        theta = -theta;
    }
    
    // clip theta to (0,2π)
    theta = std::fmod(theta, special::constant::two_pi);
    
    // clip theta to (0,π)
    if (theta > special::constant::pi)
    {
        // rotate plane
        phi += special::constant::pi;
        theta = 2 * special::constant::pi - theta;
    }
    
    // clip phi to (0,2π) FIXME (won't work for phi < 0)
    phi = std::fmod(phi, special::constant::two_pi);
}

Complex special::sphBiY (int l1, int l2, int L, int M, double theta1, double phi1, double theta2, double phi2)
{
    // NOTE This is very strange... To get the expected results one has to clip the
    //      angles to the 0..π and 0..2π intervals, respectively, and to include
    //      a bogus (-1)^m phase factor in the summation. Why the hell is that?
    
    // clip angles to the definition domain (to avoid possible uncontrolled phase factors)
    // theta = 0 .. π
    // phi = 0 .. 2π
    clipang(theta1,phi1);
    clipang(theta2,phi2);
    
    // evaluate the bi-polar spherical harmonic function
    Complex YY = 0;
    for (int m = -l1; m <= l1; m++)
    {
        Complex Y1 = special::sphY(l1, m, theta1, phi1);
        Complex Y2 = special::sphY(l2, M-m, theta2, phi2);
        double CG = special::ClebschGordan(l1,m,l2,M-m,L,M);
        YY += gsl_sf_pow_int(-1,m) * CG * Y1 * Y2;
    }
    return YY;
}

int special::coul_F_michel (int l, double k, double r, double& F, double& Fp)
{
    // initialize parameters
    double eta = -1/k;
    double rho_t = eta + std::sqrt(eta*eta + l*(l+1));
    double x = (k*r-rho_t)/rho_t;
    double a = 1 - 2*eta/rho_t;
    double phi, phip;
    
    // evaluate phi-function
    if (x < 0)
    {
        phi = -std::pow(1.5*(-std::sqrt(-x*(1+a+x)) + (1-a)/2*std::acos(1+2*x/(1+a))
              + 2*std::sqrt(a)*std::atanh(std::sqrt(-a*x/(1+a+x)))),2./3);
    }
    else
    {
        phi = std::pow(1.5*((1-a)*(std::log(std::sqrt(1+a)/(std::sqrt(x)+std::sqrt(1+x+a)))) + std::sqrt(x*(1+a+x))
              -2*std::sqrt(a)*std::atan(std::sqrt(a*x/(1+a+x)))), 2./3);
    }
    
    // evaluate derivative of the phi-function
    phip = std::sqrt((x/(1+x) + a*x/std::pow(x+1,2))/phi);
    
    // check the result and use asymptotics if the full turned unstable
    #define ASYEPS 1e-5 // TODO tune?
    //	if (std::abs(x) < ASYEPS and ( not finite(phi) or not finite(phip) ))
    if (std::abs(x) < ASYEPS or not std::isfinite(phi) or not std::isfinite(phip))
    {
        // use x ⟶ 0 asymptotic formulas
        phip = std::pow(1 + a, 1./3);
        phi = phip * x;
    }
    
    // evaluate the second derivative
    double phipp = ((x+a+1-a*x)/gsl_sf_pow_int(x+1,3) - gsl_sf_pow_int(phip,3)) / (2*phi*phip);
    
    // evaluate Airy function and its derivative
    gsl_sf_result ai, aip;
    double ai_arg = -pow(rho_t,2./3) * phi;
    int err = gsl_sf_airy_Ai_e(ai_arg, GSL_PREC_DOUBLE, &ai);
    int errp = gsl_sf_airy_Ai_deriv_e(ai_arg, GSL_PREC_DOUBLE, &aip);
    
    // evaluate the Coulomb wave function and its derivative
    F = special::constant::sqrt_pi * std::pow(rho_t,1./6)/std::sqrt(phip) * ai.val;
    Fp = special::constant::sqrt_pi * std::pow(rho_t,-5./6) * (
            0.5*std::pow(phip,-1.5)*phipp * ai.val
            - std::sqrt(phip)*aip.val*std::pow(rho_t,2./3)
    );
    
    // return corresponding GSL error
    return GSL_ERROR_SELECT_2(err, errp);
}

int special::coul_F (int l, double k, double r, double& F, double& Fp)
{
    if (r < 0.)
        return GSL_EDOM;
    
    gsl_sf_result f,g,fp,gp;
    double ef,eg;
    double eta = -1/k;
    int err;
    
    // evaluate non-S wave in origin (= zero)
    if (r == 0. and l != 0)
    {
        F = Fp = 0.;
        return GSL_SUCCESS;
    }
    
    // evaluate S wave in origin (Abramovitz & Stegun 14.6.2)
    if (r == 0. and l == 0)
    {
        gsl_sf_result C0;
        err = gsl_sf_coulomb_CL_e(l, eta, &C0);
        
        F = 0.;
        Fp = C0.val;
        
        return err;
    }
    
    err = gsl_sf_coulomb_wave_FG_e (eta, k*r, l, 0, &f, &fp, &g, &gp, &ef, &eg);
    
    // if the results are reliable, use them
    if (std::isfinite(f.val) and std::isfinite(fp.val))
    {
        F = f.val;
        Fp = fp.val;
        return GSL_SUCCESS;
    }
    
    // if the precision is insufficent, use uniform approximation
    //  if (err == GSL_ELOSS)
    {
        err = coul_F_michel(l, k, r, F, Fp);
        return err;
    }
    
    // otherwise pass the error up
    return err;
}

double special::coul_F_sigma (int l, double k)
{
    // return arg(gamma(Complex(l+1,-1./k)));
    
    gsl_sf_result lnr, arg;
    int err = gsl_sf_lngamma_complex_e(l+1, -1/k, &lnr, &arg);
    
    if (err != GSL_SUCCESS)
        HexException("Error while evaluating Coulomb phaseshift.");
    
    return arg.val;
}

double special::coul_F_asy (int l, double k, double r, double sigma)
{
    sigma = (std::isfinite(sigma) ? sigma : coul_F_sigma(l,k));
    
    return std::sin(k*r - l*special::constant::pi_half + std::log(2.*k*r)/k + sigma);
}

bool special::makes_triangle (int two_j1, int two_j2, int two_j3)
{
    return std::abs(two_j1 - two_j2) <= two_j3 and two_j3 <= two_j1 + two_j2
       and std::abs(two_j2 - two_j3) <= two_j1 and two_j1 <= two_j2 + two_j3
       and std::abs(two_j3 - two_j1) <= two_j2 and two_j2 <= two_j3 + two_j1;
}

double special::logdelta (int two_j1, int two_j2, int two_j3)
{
    return .5 * (
        std::lgamma(1 + (  two_j1 + two_j2 - two_j3) / 2)
      + std::lgamma(1 + (  two_j1 - two_j2 + two_j3) / 2)
      + std::lgamma(1 + (- two_j1 + two_j2 + two_j3) / 2)
      - std::lgamma(2 + (  two_j1 + two_j2 + two_j3) / 2)
    );
}

double special::Wigner3j_2 (int two_j1, int two_j2, int two_j3, int two_m1, int two_m2, int two_m3)
{
    // check conservation of the projection
    if (two_m1 + two_m2 + two_m3 != 0)
        return 0.;
    
    // check conservation of the magnitude
    // - sum of magnitudes has to be integer (=> twice the sum has to be even number)
    if ((two_j1 + two_j2 + two_j3) % 2 != 0)
        return 0.;
    // - in case of all zero projections the sum of magnitudes has to be even number
    if (two_m1 == 0 and two_m2 == 0 and two_m3 == 0 and ((two_j1 + two_j2 + two_j3) / 2) % 2 != 0)
        return 0.;
    
    // check projection magnitudes
    if (std::abs(two_m1) > two_j1 or std::abs(two_m2) > two_j2 or std::abs(two_m2) > two_j2)
        return 0.;
    
    // check triangle inequalities
    if (not special::makes_triangle (two_j1, two_j2, two_j3))
        return 0.;
    
    // determine sign
    int Sign = (((two_j1 - two_j2 - two_m3) / 2) % 2 == 0) ? 1 : -1;
    
    // evaluate Delta-coefficient
    double lnDelta = special::logdelta (two_j1, two_j2, two_j3);
    
    // precompute square root terms
    double lnSqrt = .5*(
        std::lgamma(1 + (two_j1 + two_m1) / 2)
      + std::lgamma(1 + (two_j1 - two_m1) / 2)
      + std::lgamma(1 + (two_j2 + two_m2) / 2)
      + std::lgamma(1 + (two_j2 - two_m2) / 2)
      + std::lgamma(1 + (two_j3 + two_m3) / 2)
      + std::lgamma(1 + (two_j3 - two_m3) / 2)
    );
    
    // compute the sum bounds
    int kmax = mmin (
        two_j1 + two_j2 - two_j3,
        two_j1 - two_m1,
        two_j2 + two_m2
    ) / 2;
    int kmin = mmax (
       -two_j3 + two_j2 - two_m1,
       -two_j3 + two_j1 + two_m2,
        0
    ) / 2;
    
    // compute the sum
    double Sum = 0.;
    for (int k = kmin; k <= kmax; k++)
    {
        // evaluate term
        double Term = std::exp (
              lnDelta
            + lnSqrt
            - std::lgamma(1 + k)
            - std::lgamma(1 + (two_j1 + two_j2 - two_j3) / 2 - k)
            - std::lgamma(1 + (two_j1 - two_m1) / 2 - k)
            - std::lgamma(1 + (two_j2 + two_m2) / 2 - k)
            - std::lgamma(1 + (two_j3 - two_j2 + two_m1) / 2 + k)
            - std::lgamma(1 + (two_j3 - two_j1 - two_m2) / 2 + k)
        );
        
        // update sum
        Sum = (k % 2 == 0) ? (Sum + Term) : (Sum - Term);
    }
    
    // return result
    return Sign * Sum;
}

double special::Wigner6j_2 (int two_j1, int two_j2, int two_j3, int two_j4, int two_j5, int two_j6)
{
    // check that triads have integer sum (i.e. twice the sum is even)
    if ((two_j1 + two_j2 + two_j3) % 2 != 0 or
        (two_j3 + two_j4 + two_j5) % 2 != 0 or
        (two_j1 + two_j5 + two_j6) % 2 != 0 or
        (two_j2 + two_j4 + two_j6) % 2 != 0)
        return 0.;
    
    // check triangle inequalities
    if (not makes_triangle (two_j1, two_j2, two_j3) or
        not makes_triangle (two_j3, two_j4, two_j5) or
        not makes_triangle (two_j1, two_j5, two_j6) or
        not makes_triangle (two_j2, two_j4, two_j6))
        return 0.;
    
    // compute Delta functions
    double d1 = logdelta (two_j1, two_j2, two_j3);
    double d2 = logdelta (two_j3, two_j4, two_j5);
    double d3 = logdelta (two_j1, two_j5, two_j6);
    double d4 = logdelta (two_j2, two_j4, two_j6);
    
    // determine sum bounds
    int kmin = mmax
    (
        two_j1 + two_j2 + two_j3,
        two_j1 + two_j5 + two_j6,
        two_j2 + two_j4 + two_j6,
        two_j3 + two_j4 + two_j5
    ) / 2;
    int kmax = mmin
    (
        two_j1 + two_j2 + two_j4 + two_j5,
        two_j1 + two_j3 + two_j4 + two_j6,
        two_j2 + two_j3 + two_j5 + two_j6
    ) / 2;
    
    // compute the sum
    double Sum = 0.;
    for (int k = kmin; k <= kmax; k++)
    {
        // compute the new term
        double Term = std::exp
        (
              d1 + d2 + d3 + d4
            + std::lgamma (2 + k)
            - std::lgamma (1 + k + (- two_j1 - two_j2 - two_j3) / 2)
            - std::lgamma (1 + k + (- two_j1 - two_j5 - two_j6) / 2)
            - std::lgamma (1 + k + (- two_j2 - two_j4 - two_j6) / 2)
            - std::lgamma (1 + k + (- two_j3 - two_j4 - two_j5) / 2)
            - std::lgamma (1 + (two_j1 + two_j2 + two_j4 + two_j5) / 2 - k)
            - std::lgamma (1 + (two_j1 + two_j3 + two_j4 + two_j6) / 2 - k)
            - std::lgamma (1 + (two_j2 + two_j3 + two_j5 + two_j6) / 2 - k)
        );
        
        // update the sum
        Sum = (k % 2 == 0) ? Sum + Term : Sum - Term;
    }
    
    // return the result
    return Sum;
}

double special::Wigner_d (int two_j, int two_ma, int two_mb, double beta)
{
    if (std::abs(two_ma) > two_j or (two_j + two_ma) % 2 != 0 or
        std::abs(two_mb) > two_j or (two_j + two_mb) % 2 != 0)
        return 0;
    
    int j_plus_ma = (two_j + two_ma) / 2;
    int j_plus_mb = (two_j + two_mb) / 2;
    int j_minus_ma = (two_j - two_ma) / 2;
    int j_minus_mb = (two_j - two_mb) / 2;
    int mb_minus_ma = (two_mb - two_ma) / 2;
    
    double cos_beta_half = std::cos(0.5 * beta);
    double sin_beta_half = std::sin(0.5 * beta);
    double suma = 0;
    
    for (int s = std::max(0,-mb_minus_ma); s <= std::min(j_plus_ma,j_minus_mb); s++)
    {
        suma += ((mb_minus_ma + s) % 2 == 0 ? 1. : -1.) / (
            gsl_sf_fact(j_plus_ma - s) * gsl_sf_fact(s) *
            gsl_sf_fact(mb_minus_ma + s) * gsl_sf_fact(j_minus_mb - s)
        ) * gsl_sf_pow_int(cos_beta_half,two_j-mb_minus_ma-2*s) *
        gsl_sf_pow_int(sin_beta_half,mb_minus_ma+2*s);
    }
    
    return suma * std::sqrt(
        gsl_sf_fact(j_plus_mb) * gsl_sf_fact(j_minus_mb) *
        gsl_sf_fact(j_plus_ma) * gsl_sf_fact(j_minus_ma)
    );
}

double special::computef (int lambda, int l1, int l2, int l1p, int l2p, int L)
{
    double A = Wigner6j(l1, l2, L, l2p, l1p, lambda);
    double B = Wigner3j(l1, lambda, l1p, 0, 0, 0);
    double C = Wigner3j(l2, lambda, l2p, 0, 0, 0);
    
    if (not std::isfinite(A))
        HexException("Wigner6j(%d,%d,%d,%d,%d,%d) not finite.", l1, l2, L, l2p, l1p, lambda);
    if (not std::isfinite(B))
        HexException("Wigner3j(%d,%d,%d,0,0,0) not finite.", l1, lambda, l1p);
    if (not std::isfinite(C))
        HexException("Wigner3j(%d,%d,%d,0,0,0) not finite.", l2, lambda, l2p);
    
    return pow(-1, L + l2 + l2p) * sqrt((2*l1 + 1) * (2*l2 + 1) * (2*l1p + 1) * (2*l2p + 1)) * A * B * C;
}

inline long double dfact (long double x)
{
    if (x < 0)
        return 0.;
    
    long double prod = 1.;
    
    while (x >= 0.0001) // = 0 + rounding errors
    {
        prod *= x;
        x -= 1.;
    }
    
    return prod;
}

double special::ClebschGordan (int in_j1, int in_m1, int in_j2, int in_m2, int in_J, int in_M)
{
    if((in_m1 + in_m2) != in_M) return 0.;
    if(abs(in_m1) > in_j1) return 0;
    if(abs(in_m2) > in_j2) return 0;
           
    // convert to pure integers (each 2*spin)
    int j1 = (int)(2.*in_j1);
    int m1 = (int)(2.*in_m1);
    int j2 = (int)(2.*in_j2);
    int m2 = (int)(2.*in_m2);
    int J = (int)(2.*in_J);
    int M = (int)(2.*in_M);
    
    long double n0,n1,n2,n3,n4,n5,d0,d1,d2,d3,d4,A,exp;
    int nu = 0;
    
    long double sum = 0;
    while (((d3=(j1-j2-M)/2+nu) < 0)||((n2=(j1-m1)/2+nu) < 0 ))
        nu++;
    while (((d1=(J-j1+j2)/2-nu) >= 0) && ((d2=(J+M)/2-nu) >= 0) && ((n1=(j2+J+m1)/2-nu) >= 0 ))
    {
        d3=((j1-j2-M)/2+nu);
        n2=((j1-m1)/2+nu);
        d0=dfact((double) nu);
        exp=nu+(j2+m2)/2;
        n0 = (double) std::pow(-1.,exp);
        sum += ((n0*dfact(n1)*dfact(n2))/(d0*dfact(d1)*dfact(d2)*dfact(d3)));
        nu++;
    }
    
    if (sum == 0)
        return 0;
    
    n0 = J+1;
    n1 = dfact((double) (J+j1-j2)/2);
    n2 = dfact((double) (J-j1+j2)/2);
    n3 = dfact((double) (j1+j2-J)/2);
    n4 = dfact((double) (J+M)/2);
    n5 = dfact((J-M)/2);
    
    d0 = dfact((double) (j1+j2+J)/2+1);
    d1 = dfact((double) (j1-m1)/2);
    d2 = dfact((double) (j1+m1)/2);
    d3 = dfact((double) (j2-m2)/2);
    d4 = dfact((double) (j2+m2)/2);
    
    A = ((long double) (n0*n1*n2*n3*n4*n5))/((long double) (d0*d1*d2*d3*d4));
    
    return std::sqrt(A)*sum;
}

double special::Gaunt (int l1, int m1, int l2, int m2, int l, int m)
{
    typedef std::tuple<int,int,int,int,int,int> TKey;
    
    // dictionary
    static std::map<TKey,double> dict;
    
    // dictionary key
    TKey key = std::make_tuple(l1,m1,l2,m2,l,m);
    
    // try to find this Gaunt's coefficient in the dictionary
    if (dict.find(key) != dict.end())
        return dict[key];
    
    // compute the value and store it to the dictionary
    return dict[key] = std::sqrt((2*l1+1)*(2*l2+1) / (4*special::constant::pi*(2*l+1))) *
        special::ClebschGordan(l1, m1, l2, m2, l, m) * special::ClebschGordan(l1,  0, l2,  0, l, 0);
}

int triangle_count (int L, int maxl)
{
    int n = 0;
    
    for (int l1 = 0; l1 <= maxl; l1++)
    for (int l2 = 0; l2 <= maxl; l2++)
    if (std::abs(l1-l2) <= L and l1+l2 >= L)
        n++;
    
    return n;
}

std::vector<std::vector<int>> special::FdB_partition (int n)
{
    // all conformant partitionings
    std::vector<std::vector<int>> partgs;
    
    // current partitioning (possibly violating the condition)
    std::vector<int> partg(n, 0);
    
    do
    {
        // increment the partitioning
        for (int digit = 0; digit < n; digit++)
        {
            if ((digit + 1) * (++partg[digit]) <= n)
                break;
            else
                partg[digit] = 0;
        }
        
        // compute Faa di Bruno sum
        int suma = 0;
        for (int i = 0; i < n; i++)
            suma += (i + 1) * partg[i];
        
        // check the condition
        if (suma == n)
            partgs.push_back(partg);
    }
    while (partg.back() == 0); // stop at the partitioning (0,0,...,1)
    
    return partgs;
}
