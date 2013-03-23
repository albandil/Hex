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

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <o2scl/ode_funct.h>
#include <o2scl/ode_iv_solve.h>

#include <gsl/gsl_sf.h>

#include "arrays.h"
#include "complex.h"
#include "wave_hyperb.h"

// -------------------------------------------------------------------------- //
//                                                                            //
//  HyperbolicWave class members                                              //
//                                                                            //
// -------------------------------------------------------------------------- //

HyperbolicWave::HyperbolicWave(double _kn, int _ln, DistortingPotential const & _U)
	: Evaluations(0), U(_U), kn(_kn), ln(_ln), Scaled(false)
{
	// get far coordinate
	double r = U.getFarRadius();
	
	// determine discretization
	samples = 10001;		// with both boundaries
	h = r / (samples - 1);	// average grid step
	
	// create grid
	grid.resize(samples);
	for (int i = 0; i < samples; i++)
		grid[i] = i * h;
	
	// look for relevant data
	char etaname[50];
	sprintf(etaname, "dwh-N%d-K%g-l%d-k%g.arr", U.n, U.k, ln, kn);
	if (not load_array(array, etaname))
	{
		// data arrays
		o2scl::ovector xg(samples-1);
		o2scl::omatrix yg(samples-1,2), ypg(samples-1,2), yerrg(samples-1,2);
		
		// derivatives callback
		o2scl::ode_funct_mfptr<HyperbolicWave> derivs(this, & HyperbolicWave::derivs);
		
		// faster adaptive Cash-Karp Runge-Kutta stepper
		o2scl::gsl_astep<decltype(derivs)> adapt_stepper;
		o2scl::gsl_rkck_fast<2,decltype(derivs)> stepper;
		adapt_stepper.set_step(stepper);
		
		// create grid
		for (int i = 0; i < samples - 1; i++)
			xg[i] = h * (i + 1);
		
		// set (left) boundary conditions from asymptotics of Riccati-Bessel function
		yg[0][0] = h * 1e-5;				// Θ (h)
		yg[0][1] = (ln + 1 - kn*h) * 1e-5;	// Θ'(h)
		
		solve2(xg, samples, h, yg, ypg, yerrg, adapt_stepper, derivs, NORMALIZE_ON_OVERFLOW);
		
		// save with correct normalization
		array = rArray(samples);
		array[0] = 0.;
		double norm = ric_i_scaled(ln,kn*r) / yg[samples-3][0];
		for (int i = 1; i < samples; i++)
			array[i] = yg[i-1][0] * norm;
		
		save_array(array, etaname);
	}
	
	// setup the interpolator
	this->interpolator.set(grid, array, array.size());
}

HyperbolicWave HyperbolicWave::operator=(HyperbolicWave const & W)
{
	kn = W.kn;
	ln = W.ln;
	U = W.U;
	grid = W.grid;
	array = W.array;
	samples = W.samples;
	h = W.h;
	
	interpolator.set(grid, array, array.size());
	
	return *this;
}

int HyperbolicWave::derivs(double x, size_t nv, const o2scl::ovector_base& y, o2scl::ovector_base& dydx)
{
	dydx[0] = y[1];
	dydx[1] = (x == 0.) ? 0. : (2.*U(x) + ln*(ln+1)/(x*x)) * y[0] - 2*kn*y[1];
	return 0;
}

double HyperbolicWave::operator()(double x) const
{
	if (x > grid.back())
	
		// extrapolate
		return Scaled ? ric_i_scaled(ln,kn*x) : ric_i(ln,kn*x);
	
	else
		
		// interpolate
		return Scaled ? interpolator.interp(x) : interpolator.interp(x) * exp(kn*x);
}

void HyperbolicWave::toFile(const char* filename) const
{
	write_array(grid, array, filename);
}

void HyperbolicWave::scale(bool s)
{
	Scaled = s;
}
