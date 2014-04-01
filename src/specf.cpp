/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2014                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cmath>
#include <iostream>

#include <gsl/gsl_sf.h>

#include "arrays.h"
#include "specf.h"
#include "complex.h"

// ----------------------------------------------------------------------- //
//  Special functions                                                      //
// ----------------------------------------------------------------------- //

cArray ric_jv(int lmax, Complex z)
{
    // results
    cArray eval(lmax+1);
    
    // use library routine for pure real arguments
    if (z.imag() == 0.)
    {
        rArray ev(lmax+1);
        int err = gsl_sf_bessel_jl_steed_array(lmax, z.real(), &ev[0]);
        
        // check that all evaluations are finite
        bool all_finite = std::all_of (
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
                    ev[n] = ev[n-1] * z.real() / (2*n+1);
            }
            else
            {
                throw exception ("Error %d while evaluating j[l≤%d](%g+%gi).", err, lmax, z.real(), z.imag());
            }
        }
        
        // Bessel -> Riccati-Bessel function
        for (int i = 0; i <= lmax; i++)
            eval[i] = z * Complex(ev[i]);
        
        return eval;
    }
    
    // shorthand
    Complex inv_z = Complex(1.)/z;
    
    // evaluate all angular momenta up to lmax
    for (int l = 0; l <= lmax; l++)
    {
        if (l == 0)
            eval[l] = sin(z);
        else if (l == 1)
            eval[l] = sin(z) * inv_z - cos(z);
        else
            eval[l] = Complex(2.*l - 1.) * eval[l-1] * inv_z - eval[l-2];
    }
    
    return eval;
}

Complex ric_j(int l, Complex z)
{
    return ric_jv(l, z).back();
}

cArray dric_jv(int lmax, Complex z)
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
            deval[l] = cos(z);
        else if (l == 1)
            deval[l] = inv_z * ( cos(z) - sin(z) * inv_z ) + sin(z);
        else
            deval[l] = -Complex(2.*l - 1) * (eval[l-1] * inv_z - deval[l-1] ) * inv_z - deval[l-2];
    }
    
    return deval;
}

Complex dric_j(int l, Complex z)
{
    return dric_jv(l, z).back();
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

Complex hydro_P_table(unsigned n, unsigned l, Complex z)
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
    return sum * pow(z, l + 1) * exp(-z/double(n));
}

Complex dhydro_P_table(unsigned n, unsigned l, Complex z)
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
    return sum * pow(z, l) * exp(-z/double(n)) - hydro_P_table(n,l,z) / double(n);
}

/*
 * Laguerre polynomial
 * Laguerre(k,s,x) := sum((-1)^j * (k!)^2 * x^(j-s) / ( (k-j)! * j! * (j-s)! ), j, s, k);
 */
Complex associated_laguerre_poly(int k, int s, Complex z)
{
    // value of the polynomial to be returned
    Complex val = 0;
    
    // begin with highest order
    val += pow(-1,k) / (fac(k) * fac(k-s));
    
    // continue with other orders
    for (int j = k - 1; j >= s; j--)
        val = z * val + pow(-1,j) / (fac(k-j) * fac(j) * fac(j-s));
    
    return val * pow(fac(k),2);
}

/*
 * Derivative of Laguerre polynomial
 * DLaguerre(k,s,x) := sum((-1)^j * (k!)^2 * x^(j-s-1) / ( (k-j)! * j! * (j-s-1)! ), j, s+1, k);
 */
Complex der_associated_laguerre_poly(int k, int s, Complex z)
{
    // value of the polynomial to be returned
    Complex val = 0;
    
    // begin with highest order
    val += pow(-1,k) / (fac(k) * fac(k-s-1));
    
    // continue with other orders
    for (int j = k - 1; j >= s + 1; j--)
        val = z * val + pow(-1,j) / (fac(k-j) * fac(j) * fac(j-s-1));
    
    return val * pow(fac(k),2);
}

/*
 * Hydrogen radial function normalization factor
 * sqrt((2/n)^3 * (n-l-1)! / (2*n*((n+l)!)^3));
 */
double hydrogen_wfn_normalization(int n, int l)
{
    return sqrt(pow(2./n,3) * fac(n-l-1) / (2*n*pow(fac(n+l),3)));
}

/*
 * Hydrogen radial function.
 * HydrogenP(n,l,r) := r * HydrogenN(n,l) * (2*r/n)^l * Laguerre(n+l,2*l+1,2*r/n) * exp(-r/n);
 */
Complex hydro_P(unsigned n, unsigned l, Complex z)
{
    // this is faster
    if (n <= max_table_n)
        return hydro_P_table(n, l, z);
    
    // this is general
    double Norm = hydrogen_wfn_normalization(n, l);
    Complex Lag = associated_laguerre_poly(n + l, 2 * l + 1, z);
    return z * Norm * pow(2.*z/double(n), l) * Lag * exp(-z/double(n));
}

/*
 * Derivative of hydrogen radial function.
 * DP(n,l,r) := HydrogenN(n,l) * (2*r/n)^l * (
 * 					(l+1-z/n) * Laguerre(n+l,2*l+1,2*r/n) +
 * 					z * DLaguerre(n+l,2*l+1,2*r/n)
 *              ) * exp(-r/n)
 */
Complex dhydro_P(unsigned n, unsigned l, Complex z)
{
    // this is faster
    if (n <= max_table_n)
        return dhydro_P_table(n, l, z);
    
    // this is general
    double Norm = hydrogen_wfn_normalization(n, l);
    Complex Lag = associated_laguerre_poly(n + l, 2 * l + 1, z);
    Complex DLag= der_associated_laguerre_poly(n + 1, 2 * l + 1, z);
    return Norm * pow(2.*z/double(n), l) * ((l + 1. - z/double(n)) * Lag + z * DLag) * exp(-z/double(n));
}


Complex sphY(int l, int m, double theta, double phi)
{
    if (l < abs(m))
        return 0.;
    return gsl_sf_legendre_sphPlm(l,abs(m),cos(theta)) * Complex(cos(m*phi),sin(m*phi));
}

void clipang (double & theta, double & phi)
{
    // clip theta to (0,inf)
    if (theta < 0.)
    {
        // rotate plane
        phi += M_PI;
        theta = -theta;
    }
    
    // clip theta to (0,2π)
    theta = fmod(theta, 2 * M_PI);
    
    // clip theta to (0,π)
    if (theta > M_PI)
    {
        // rotate plane
        phi += M_PI;
        theta = 2 * M_PI - theta;
    }
    
    // clip phi to (0,2π) FIXME (won't work for phi < 0)
    phi = fmod(phi, 2*M_PI);
}

Complex sphBiY(int l1, int l2, int L, int M, double theta1, double phi1, double theta2, double phi2)
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
        YY += gsl_sf_pow_int(-1,m) * ClebschGordan(l1,m,l2,M-m,L,M) * sphY(l1, m, theta1, phi1) * sphY(l2, M-m, theta2, phi2);
    return YY;
}

int coul_F_michel(int l, double k, double r, double& F, double& Fp)
{
    // initialize parameters
    double eta = -1/k;
    double rho_t = eta + sqrt(eta*eta + l*(l+1));
    double x = (k*r-rho_t)/rho_t;
    double a = 1 - 2*eta/rho_t;
    double phi, phip;
    
    // evaluate phi-function
    if (x < 0)
        phi = -pow(1.5*(-sqrt(-x*(1+a+x)) + (1-a)/2*acos(1+2*x/(1+a)) + 2*sqrt(a)*atanh(sqrt(-a*x/(1+a+x)))),2./3);
    else
        phi = pow(1.5*((1-a)*(log(sqrt(1+a)/(sqrt(x)+sqrt(1+x+a)))) + sqrt(x*(1+a+x)) -2*sqrt(a)*atan(sqrt(a*x/(1+a+x)))), 2./3);
    
    // evaluate derivative of the phi-function
    phip = sqrt((x/(1+x) + a*x/pow(x+1,2))/phi);
    
    // check the result and use asymptotics if the full turned unstable
    #define ASYEPS 1e-5 // TODO tune?
    //	if (std::abs(x) < ASYEPS and ( not finite(phi) or not finite(phip) ))
    if (std::abs(x) < ASYEPS or not std::isfinite(phi) or not std::isfinite(phip))
    {
        // use x ⟶ 0 asymptotic formulas
        phip = pow(1 + a, 1./3);
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
    F = sqrt(M_PI)*pow(rho_t,1./6)/sqrt(phip) * ai.val;
    Fp = sqrt(M_PI)*pow(rho_t,-5./6)*(0.5*pow(phip,-1.5)*phipp * ai.val - sqrt(phip)*aip.val*pow(rho_t,2./3));
    
    // return corresponding GSL error
    return GSL_ERROR_SELECT_2(err, errp);
}

int coul_F(int l, double k, double r, double& F, double& Fp)
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
    // 	if (err == GSL_ELOSS)
    {
        err = coul_F_michel(l, k, r, F, Fp);
        return err;
    }
    
    // otherwise pass the error up
    return err;
}

double coul_F_sigma(int l, double k)
{
    // 	return arg(gamma(Complex(l+1,-1./k)));
    
    gsl_sf_result lnr, arg;
    int err = gsl_sf_lngamma_complex_e(l+1, -1/k, &lnr, &arg);
    
    if (err != GSL_SUCCESS)
        throw exception ("Error while evaluating Coulomb phaseshift.");
    
    return arg.val;
}

double coul_F_asy(int l, double k, double r, double sigma)
{
    if (std::isfinite(sigma))
        return sqrt(M_2_PI)/k * sin(k*r - 0.5*l*M_PI + log(2*k*r)/k + sigma);
    else
        return sqrt(M_2_PI)/k * sin(k*r - 0.5*l*M_PI + log(2*k*r)/k + coul_F_sigma(l,k));
}

bool makes_triangle (int two_j1, int two_j2, int two_j3)
{
    return std::abs(two_j1 - two_j2) <= two_j3 and two_j3 <= two_j1 + two_j2
       and std::abs(two_j2 - two_j3) <= two_j1 and two_j1 <= two_j2 + two_j3
       and std::abs(two_j3 - two_j1) <= two_j2 and two_j2 <= two_j3 + two_j1;
}

double logdelta (int two_j1, int two_j2, int two_j3)
{
    return .5 * (
        lgamma(1 + (  two_j1 + two_j2 - two_j3) / 2)
      + lgamma(1 + (  two_j1 - two_j2 + two_j3) / 2)
      + lgamma(1 + (- two_j1 + two_j2 + two_j3) / 2)
      - lgamma(2 + (  two_j1 + two_j2 + two_j3) / 2)
    );
}

double Wigner3j_2 (int two_j1, int two_j2, int two_j3, int two_m1, int two_m2, int two_m3)
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
    if (not makes_triangle (two_j1, two_j2, two_j3))
        return 0.;
    
    // determine sign
    int Sign = (((two_j1 - two_j2 - two_m3) / 2) % 2 == 0) ? 1 : -1;
    
    // evaluate Delta-coefficient
    double lnDelta = logdelta (two_j1, two_j2, two_j3);
    
    // precompute square root terms
    double lnSqrt = .5*(
        lgamma(1 + (two_j1 + two_m1) / 2)
      + lgamma(1 + (two_j1 - two_m1) / 2)
      + lgamma(1 + (two_j2 + two_m2) / 2)
      + lgamma(1 + (two_j2 - two_m2) / 2)
      + lgamma(1 + (two_j3 + two_m3) / 2)
      + lgamma(1 + (two_j3 - two_m3) / 2)
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
        double Term = exp (
              lnDelta
            + lnSqrt
            - lgamma(1 + k)
            - lgamma(1 + (two_j1 + two_j2 - two_j3) / 2 - k)
            - lgamma(1 + (two_j1 - two_m1) / 2 - k)
            - lgamma(1 + (two_j2 + two_m2) / 2 - k)
            - lgamma(1 + (two_j3 - two_j2 + two_m1) / 2 + k)
            - lgamma(1 + (two_j3 - two_j1 - two_m2) / 2 + k)
        );
        
        // update sum
        Sum = (k % 2 == 0) ? (Sum + Term) : (Sum - Term);
    }
    
    // return result
    return Sign * Sum;
}

double Wigner6j_2 (int two_j1, int two_j2, int two_j3, int two_j4, int two_j5, int two_j6)
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
    int kmin = mmax (
        two_j1 + two_j2 + two_j3,
        two_j1 + two_j5 + two_j6,
        two_j2 + two_j4 + two_j6,
        two_j3 + two_j4 + two_j5
    ) / 2;
    int kmax = mmin (
        two_j1 + two_j2 + two_j4 + two_j5,
        two_j1 + two_j3 + two_j4 + two_j6,
        two_j2 + two_j3 + two_j5 + two_j6
    ) / 2;
    
    // compute the sum
    double Sum = 0.;
    for (int k = kmin; k <= kmax; k++)
    {
        // compute the new term
        double Term = exp (
              d1 + d2 + d3 + d4
            + lgamma (2 + k)
            - lgamma (1 + k + (- two_j1 - two_j2 - two_j3) / 2)
            - lgamma (1 + k + (- two_j1 - two_j5 - two_j6) / 2)
            - lgamma (1 + k + (- two_j2 - two_j4 - two_j6) / 2)
            - lgamma (1 + k + (- two_j3 - two_j4 - two_j5) / 2)
            - lgamma (1 + (two_j1 + two_j2 + two_j4 + two_j5) / 2 - k)
            - lgamma (1 + (two_j1 + two_j3 + two_j4 + two_j6) / 2 - k)
            - lgamma (1 + (two_j2 + two_j3 + two_j5 + two_j6) / 2 - k)
        );
        
        // update the sum
        Sum = (k % 2 == 0) ? Sum + Term : Sum - Term;
    }
    
    // return the result
    return Sum;
}

double computef (int lambda, int l1, int l2, int l1p, int l2p, int L)
{
    double A = Wigner6j(l1, l2, L, l2p, l1p, lambda);
    double B = Wigner3j(l1, lambda, l1p, 0, 0, 0);
    double C = Wigner3j(l2, lambda, l2p, 0, 0, 0);
    
    if (not std::isfinite(A))
        throw exception ("Wigner6j(%d,%d,%d,%d,%d,%d) not finite.", l1, l2, L, l2p, l1p, lambda);
    if (not std::isfinite(B))
        throw exception ("Wigner3j(%d,%d,%d,0,0,0) not finite.", l1, lambda, l1p);
    if (not std::isfinite(C))
        throw exception ("Wigner3j(%d,%d,%d,0,0,0) not finite.", l2, lambda, l2p);
    
    return pow(-1, L + l2 + l2p) * sqrt((2*l1 + 1) * (2*l2 + 1) * (2*l1p + 1) * (2*l2p + 1)) * A * B * C;
}

inline long double dfact(long double x)
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

double ClebschGordan(int __j1, int __m1, int __j2, int __m2, int __J, int __M)
{
    if((__m1 + __m2) != __M) return 0.;
    if(abs(__m1) > __j1) return 0;
    if(abs(__m2) > __j2) return 0;
           
    // convert to pure integers (each 2*spin)
    int j1 = (int)(2.*__j1);
    int m1 = (int)(2.*__m1);
    int j2 = (int)(2.*__j2);
    int m2 = (int)(2.*__m2);
    int J = (int)(2.*__J);
    int M = (int)(2.*__M);
    
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
        n0 = (double) pow(-1.,exp);
        sum += ((n0*dfact(n1)*dfact(n2))/(d0*dfact(d1)*dfact(d2)*dfact(d3)));
        nu++;
    }
    
    if (sum == 0) return 0;
    
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
    
    return sqrtl(A)*sum;           
}

double Gaunt(int l1, int m1, int l2, int m2, int l, int m)
{
    // dictionary
    static std::map<std::tuple<int,int,int,int,int,int>,double> dict;
    
    // dictionary key
    std::tuple<int,int,int,int,int,int> key = std::make_tuple(l1,m1,l2,m2,l,m);
    
    // try to find this Gaunt's coefficient in the dictionary
    if (dict.find(key) != dict.end())
        return dict[key];
    
    // compute the value and store it to the dictionary
    dict[key] = sqrt((2*l1+1)*(2*l2+1) / (4*M_PI*(2*l+1))) *
    ClebschGordan(l1, m1, l2, m2, l, m) *
    ClebschGordan(l1,  0, l2,  0, l, 0);
    return dict[key];
}

int triangle_count(int L, int maxl)
{
    int n = 0;
    
    for (int l1 = 0; l1 <= maxl; l1++)
    for (int l2 = 0; l2 <= maxl; l2++)
    if (std::abs(l1-l2) <= L and l1+l2 >= L)
        n++;
    
    return n;
}

double Hyper2F1 (double a, double b, double c, double x)
{
    if (x < -1.0)
    {
        double w = 1 / (1 - x);
        return std::pow(1-x,-a) * std::exp(gsl_sf_lngamma(c) + gsl_sf_lngamma(b-a) - gsl_sf_lngamma(b) - gsl_sf_lngamma(c-a)) * gsl_sf_hyperg_2F1(a,c-b,a-b+1,w)
             + std::pow(1-x,-b) * std::exp(gsl_sf_lngamma(c) + gsl_sf_lngamma(a-b) - gsl_sf_lngamma(a) - gsl_sf_lngamma(c-b)) * gsl_sf_hyperg_2F1(b,c-a,b-a+1,w);
    }
    else if (x < 0.0)
    {
        double w = x / (1 - x);
        return std::pow(1-x,-a) * gsl_sf_hyperg_2F1(a,c-b,c,w);
    }
    else if (x <= 0.5)
    {
        return gsl_sf_hyperg_2F1(a,b,c,x);
    }
    else if (x <= 1.0)
    {
        double w = 1 - x;
        return std::exp(gsl_sf_lngamma(c) + gsl_sf_lngamma(c-a-b) - gsl_sf_lngamma(c-a) - gsl_sf_lngamma(c-b)) * gsl_sf_hyperg_2F1(a,b,a+b+1-c,w)
             + std::pow(1-x,c-a-b) * std::exp(gsl_sf_lngamma(c) + gsl_sf_lngamma(a+b-c) - gsl_sf_lngamma(a) - gsl_sf_lngamma(b)) * gsl_sf_hyperg_2F1(c-a,c-b,c-a-b+1,w);
    }
    else if (x <= 2.0)
    {
        double w = 1 - 1/x;
        return std::pow(x,-a) * std::exp(gsl_sf_lngamma(c) + gsl_sf_lngamma(c-a-b) - gsl_sf_lngamma(c-a) - gsl_sf_lngamma(c-b)) * gsl_sf_hyperg_2F1(a,a-c+1,a+b-c+1,w)
             + std::pow(x,a-c) * std::pow(1-x,c-a-b) * std::exp(gsl_sf_lngamma(c) + gsl_sf_lngamma(a+b-c) - gsl_sf_lngamma(a) - gsl_sf_lngamma(b)) * gsl_sf_hyperg_2F1(c-a,1-a,c-a-b+1,w);
    }
    else /* x > 2 */
    {
        double w = 1/x;
        return std::pow(-x,-a)*std::exp(gsl_sf_lngamma(c) + gsl_sf_lngamma(b-a) - gsl_sf_lngamma(b) - gsl_sf_lngamma(c-b)) * gsl_sf_hyperg_2F1(a,a-c+1,a-b+1,w)
             + std::pow(-x,-b)*std::exp(gsl_sf_lngamma(c) + gsl_sf_lngamma(a-b) - gsl_sf_lngamma(a) - gsl_sf_lngamma(c-b)) * gsl_sf_hyperg_2F1(b-c+1,b,b-a+1,w);
    }
}
