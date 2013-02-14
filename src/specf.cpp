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

#include <cmath>

#include <gsl/gsl_sf.h>

#include "specf.h"
#include "complex.h"

Complex sphY(int l, int m, double theta, double phi)
{
	if (l < abs(m))
		return 0.;
	return gsl_sf_legendre_sphPlm(l,abs(m),cos(theta)) * Complex(cos(m*phi),sin(m*phi));
}

static const int g = 7; 
static const double pi = 3.1415926535897932384626433832795028841972;
static const double p[g+2] = {
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
};

Complex gamma(Complex z)
{
        if (real(z) < 0.5)
        {
                 return pi / (sin(pi*z) * gamma(1. - z));
        }
         z -= 1.;
        Complex x = p[0];
        for (int i = 1; i < g + 2; i++)
        {
                x += p[i] / (z + Complex(i,0));
        }
        Complex t = z + (g + 0.5);
        return sqrt(2*pi) * pow(t,z+0.5) * exp(-t) * x;
}

double F_sigma(int l, double k)
{
	return arg(gamma(Complex(l+1,-1./k)));
}

double F_asy(int l, double k, double r, double sigma)
{
	if (finite(sigma))
		return sqrt(2./M_PI)/k * sin(k*r - l*M_PI/2. + log(2*k*r)/k + sigma);
	else
		return sqrt(2./M_PI)/k * sin(k*r - l*M_PI/2. + log(2*k*r)/k + F_sigma(l,k));
}

Complex cgamma_cfrac (double a, Complex z, int max_iter, int* iter, double eps)
{
	Complex CF = 0, cf;
	int n;
	
	// for all continued fraction depths
	for (n = 0; n <= max_iter; n++)
	{
		// evaluate the continued fraction
		cf = 0.;
		for (int m = n; m > 0; m--)
			cf = m * (m - a) / (z + 2.*m + 1. - a - cf);
		cf = 1. / (z + 1. - a - cf);
		
		// check convergence
		if (n > 0 and abs(cf-CF) / abs(cf) < eps)
		{
			CF = cf;
			break;
		}
		
		// store this iteration
		CF = cf;
	}
	
	// store iteration count
	if (iter != 0)
		*iter = n;
	
	return CF * pow(z,a) * exp(-z);
}

Complex cgamma_series (double a, Complex z, int max_iter, int* iter, double eps)
{
	// return zero in origin
	if (z == 0.)
		return 0.;
	
	// auxiliary variables
	Complex sum = 1./a, term = 1./a;
	int n;
	
	// for all terms of the series
	for (n = 1; n <= max_iter; n++)
	{
		// evaluate new term using previous term
		term *= z / (a + n);
		
		// update sum
		sum += term;
		
		// check convergence
		if (n > 0 and abs(term) / abs(sum) < eps)
			break;
	}
	
	// store iteration count
	if (iter != 0)
		*iter = n;
	
	return sum * pow(z,a) * exp(-z);
}

Complex cgamma (double a, Complex z, int max_iter, int* iter, double eps)
{
	// choose evaluation method
	
	if (abs(z) <= a + 1)
		// use infinite series
		return tgamma(a) - cgamma_series(a,z,max_iter,iter,eps);
	
	else
		// use infinite continued fraction
		return cgamma_cfrac(a,z,max_iter,iter,eps);
}
