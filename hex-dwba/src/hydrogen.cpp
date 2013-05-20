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

#include <cstdlib>
#include <cmath>
#include <cstddef>

#include <o2scl/exception.h>
#include <gsl/gsl_sf.h>

#include "hydrogen.h"
#include "interpolate.h"
#include "misc.h"
#include "specf.h"

namespace Hydrogen
{

double lastZeroBound(int n, int l)
{
	// all Laguerre(r) roots are less or equal to
	return (n + l + (n - l - 2.) * sqrt(n + l));
}

double getBoundFar(int n, int l, double eps, int max_steps)
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

double getBoundN(int n, int l)
{
	return sqrt(gsl_sf_pow_int(2./n,3) * gsl_sf_fact(n-l-1) / (2*n*gsl_sf_fact(n+l))) * gsl_sf_pow_int(2./n,l);
}

double getSturmFar(int n, int l, double lambda, double eps, int max_steps)
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

double evalBoundState(int n, int l, double r)
{
	double psi;
	
	try {
		
		psi = r*gsl_sf_hydrogenicR(n, l, 1., r);
		
	} catch (o2scl::exc_range_error e) {
		
		// probably due to the underflow error
		return 0.;
		
	}
	
	return psi;
}

double evalFreeState(double k, int l, double r, double sigma)
{
	// Coulomb wave is a regular solution
	if (r == 0.)
		return 0.;
	
	//normalization
	double norm = sqrt(2./M_PI)/k;
	
	// some local variables
	double F, exp_F;
		
	// evaluate the function
	try {
		
  		gsl_sf_coulomb_wave_F_array(l, 0, -1./k, k*r, &F, &exp_F);
		return norm * F;
		
	} catch (std::exception e) {
		
		// evaluation failed
		
		if (k * r > 1)
		{
			// probably due to "iteration process out of control" for large radii
			return F_asy(l, k, r, (finite(sigma) ? sigma : F_sigma(l,k)));
		}
		else
		{
			// some other problem
			throw exception("Evaluation of hydrogen free state failed for l = %d, k = %g, r = %g\nMessage: %s\n", l, k, r, e.what());
		}
		
	}
}

double evalFreeStatePhase(double k, int l, double sigma)
{
	return l * 0.5 * M_PI + (finite(sigma) ? sigma : F_sigma(l,k));
}

double evalSturmian(int n, int l, double r, double lambda)
{
	double S;
	
	try {
		
		S = n * n * pow(lambda,l+1) * sqrt(pow(2./n,3)*gsl_sf_fact(n-l-1.)/(2.*n*gsl_sf_fact(n+l))) 
			* exp(-lambda*r) * pow(2.*r,l) * gsl_sf_laguerre_n(n-l-1,2*l+1,2*r);
		
	} catch (o2scl::exc_range_error e) {
		
		// probably due to the underflow error
		return 0.;
		
	}
	
	return S;
}

double evalFreeState_asy(double k, int l, double r, double sigma)
{
	return F_asy(l,k,r,sigma);
}

double getFreeAsyZero(double k, int l, double Sigma, double eps, int max_steps, int nzero)
{
	// find solution (r) of
	//    n*pi = k*r - l*pi/2 + log(2*k*r)/k + Sigma
	// using trivial Banach contraction
	
	double kr = (2*nzero+l)*M_PI/2 - Sigma;
	if (kr < 0)
		return kr;
	
	while (max_steps-- > 0 and std::abs(kr) > eps)
		kr = (2*nzero+l)*M_PI/2 - Sigma - log(2*kr)/k;
	
	return kr / k;
}

double getFreeAsyTop(double k, int l, double Sigma, double eps, int max_steps, int ntop)
{
	// find solution (r) of
	//    (n+1/4)*2*pi = k*r - l*pi/2 + log(2*k*r)/k + Sigma
	// using trivial Banach contraction
	
	double kr = (2*ntop+0.5*(l+1))*M_PI - Sigma;
	if (kr < 0)
		return kr;
	
	while (max_steps-- > 0 and std::abs(kr) > eps)
		kr = (2*ntop+0.5*(l+1))*M_PI - Sigma - log(2*kr)/k;
	
	return kr / k;
}

double getFreeFar(double k, int l, double Sigma, double eps, int max_steps)
{
	// precompute Coulomb shift
	Sigma = finite(Sigma) ? Sigma : F_sigma(l,k);
	
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
		eval = std::abs(evalFreeState(k,l,rzero,Sigma));
		
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
		eval = std::abs(evalFreeState(k,l,rzero,Sigma));
		
		// move bisection guardians
		if (eval > eps)
			idx_left = idx_mid;
		else
			idx_right = idx_mid;
		
		idx_mid = (idx_right + idx_left) / 2;
	}
	
	return idx_mid * M_PI;
}

}; // endof namespace Hydrogen

double HydrogenFunction::operator()(double r) const
{
	if (n != 0)
	{
		// this state is bound
		// -------------------
		
		// if too far, return zero
		if (r > Far)
			return 0.;
		
		// if too near, return simplified value
// 		if (r < 1e-5)
// 			return Hydrogen::getBoundN(n,l) * gsl_sf_pow_int(r,l+1);
		
		// else compute precise value
		return Hydrogen::evalBoundState(n, l, r);
	}
	else
	{
		// this state is free
		// ------------------
		
		// if near enough, return precise value
		if (Far == 0 or r < Far)
			return Hydrogen::evalFreeState(k, l, r, Sigma == 0 ? Nan : Sigma);
		
		// else evaluate using the scale factor
 		return Hydrogen::evalFreeState_asy(k, l, r, Sigma) / (1. + 0.5/(r*k*k));
	}
}

Complex HydrogenFunction::operator()(Complex r) const
{
	if (n == 0)
		throw exception("[HydrogenFunction] Do not use operator() for evaluating complex free states!");
	
	if (n == 1 and l == 0)
		return 2. * r * exp(-r);
		
	throw exception("[HydrogenFunction] Complex argument for bound state P{%d,%d}(r) not implemented.", n, l);
}
