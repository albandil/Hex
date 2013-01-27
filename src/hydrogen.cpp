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
#include <cstddef>

#include <o2scl/exception.h>
#include <gsl/gsl_sf.h>

#include "hydrogen.h"

double Hydrogen::lastZeroBound(int n, int l)
{
	// all Laguerre(r) roots are less or equal to
	return (n + l + (n - l - 2.) * sqrt(n + l));
}

double Hydrogen::getFarRadius(int n, int l, double eps, int max_steps)
{
	// all Laguerre(2r/n) roots are less or equal to
	double last_zero_ubound = lastZeroBound(n, l) * 0.5 * n;
	
	// hunt for low value
	double far = std::max(last_zero_ubound, 1.);
	double far_val, der;
	do
	{
		// move further
		far *= 2;
		
		// evaluate function in 'far'
		far_val = fabs(evalBoundState(n,l,far));
		
		// compute forward h-diference in 'far'
		der = far_val - fabs(evalBoundState(n,l,far+1e-5));
	}
	while (der > 0 or far_val > eps);
	
	// bisect for exact value
	int steps = 0;
	double near = last_zero_ubound;
	double near_val = fabs(evalBoundState(n,l,near));
	while (steps < max_steps)
	{
		double middle = 0.5 * (near + far);
		double middle_val = fabs(evalBoundState(n,l,middle));
		if (middle_val > eps)
		{
			near = middle;
			near_val = middle_val;
		}
		else
		{
			far = middle;
			far_val = middle_val;
		}
		steps++;
		if (near_val == far_val)
			break;
	}

	// return bisection
	return 0.5 * (near + far);
}

double Hydrogen::getFarRadiusS(int n, int l, double lambda, double eps, int max_steps)
{
	// all Laguerre(2r) roots are less or equal to
	double last_zero_ubound = lastZeroBound(n,l) * 0.5;
	
	// hunt for low value
	double far = std::max(last_zero_ubound, 1.);
	double far_val, der;
	do
	{
		// move further
		far *= 2;
		
		// evaluate function in 'far'
		far_val = fabs(evalSturmian(n,l,far,lambda));
		
		// compute forward h-diference in 'far'
		der = fabs(evalSturmian(n,l,far+0.001,lambda)) - far_val;
	}
	while (der > 0 or far_val > eps);
	
	// bisect for exact value
	int steps = 0;
	double near = last_zero_ubound;
	double near_val = fabs(evalSturmian(n,l,near,lambda));
	while (steps < max_steps)
	{
		double middle = 0.5 * (near + far);
		double middle_val = fabs(evalSturmian(n,l,middle,lambda));
		if (middle_val > eps)
		{
			near = middle;
			near_val = middle_val;
		}
		else
		{
			far = middle;
			far_val = middle_val;
		}
		steps++;
		if (near_val == far_val)
			break;
	}
	
	// return bisection
	return 0.5 * (near + far);
}

double Hydrogen::evalBoundState(int n, int l, double r)
{
	double psi;
	
	try {
		psi = r*gsl_sf_hydrogenicR(n, l, 1., r);
	} catch (o2scl::exc_range_error e) {
		return 0.;
	}
	
	return psi;
}

double Hydrogen::evalFreeState(double k, int l, double r)
{
	// Coulomb wave is a regular solution
	if (r == 0.)
		return 0.;
	
	//normalization
	double norm = sqrt(2./M_PI) / k;
	
	// some local variables
	double F, exp_F;
		
	// evaluate the function
	try {
		
		int err = gsl_sf_coulomb_wave_F_array(l, 0, -1./k, k*r, &F, &exp_F);
		
		if (err != GSL_SUCCESS)
		{
			fprintf(stderr, "Error (%s) while evaluating free hydrogen function.\n", gsl_strerror(err));
			abort();
		}
		
		return norm * F;
		
	} catch (o2scl::exc_runtime_error e) {
		
		fprintf(stderr, "Error (%s) while evaluating F[%d](%g,%g).\n", e.what(), l, k, r);
		abort();
		
	} catch (o2scl::exc_range_error e) {
		
		fprintf(stderr, "Error (%s) while evaluating F[%d](%g,%g).\n", e.what(), l, k, r);
		abort();
		
	}
}

double Hydrogen::evalSturmian(int n, int l, double r, double lambda)
{
	double S;
	
	try {
		S = n * n * pow(lambda,l+1) * sqrt(pow(2./n,3)*gsl_sf_fact(n-l-1.)/(2.*n*gsl_sf_fact(n+l))) 
			* exp(-lambda*r) * pow(2.*r,l) * gsl_sf_laguerre_n(n-l-1,2*l+1,2*r);
	} catch (o2scl::exc_range_error e) {
		return 0.;
	}
	
	return S;
}

double Hydrogen::evalFreeState_asy(double k, int l, double r)
{
	return F_asy(l,k,r);
}
