/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2013                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <complex>
#include <cmath>
#include <cstdlib>

#include <gsl/gsl_sf.h>
#define fac(n) gsl_sf_fact(n)

#include "arrays.h"

// ----------------------------------------------------------------------- //
//  Special functions                                                      //
// ----------------------------------------------------------------------- //

LComplex ric_j(int n, LComplex z)
{
	if (z.imag() == 0.)
	{
		// use library routine for pure real arguments
		gsl_sf_result j;
		int err = gsl_sf_bessel_jl_e(n, z.real(), &j);
		if (err != GSL_SUCCESS)
			throw exception("Error %d while evaluating j[%d](%d+%di).", err, n, z.real(), z.imag());
		z *= j.val;
		return z;
	}
	
	if (n == 0)
		return sin(z);
	else if (n == 1)
		return sin(z)/z - cos(z);
	else
		return LComplex(2.*n - 1.) * ric_j(n-1,z) / z - ric_j(n-2,z);
}

LComplex dric_j(int n, LComplex z)
{
	if (n == 0)
		return cos(z);
	else if (n == 1)
		return cos(z)/z - sin(z)/(z*z) + sin(z);
	else
		return -LComplex(2.*n - 1) * ric_j(n-1,z) / (z*z) + LComplex(2.*n - 1.) * dric_j(n-1,z) / z - dric_j(n-2,z);
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
