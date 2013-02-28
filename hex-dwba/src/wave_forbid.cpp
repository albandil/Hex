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
#include "wave_forbid.h"

// -------------------------------------------------------------------------- //
//                                                                            //
//  ForbiddenWave class members                                               //
//                                                                            //
// -------------------------------------------------------------------------- //

ForbiddenWave::ForbiddenWave(double _kn, int _ln, DistortingPotential const & _U)
	: Evaluations(0), U(_U), kn(_kn), ln(_ln)
{
	// get far coordinate
	double r = U.getFarRadius();
	
	// determine discretization
	samples = 10001;		// with both boundaries
	h = r / (samples - 1);	// grid step
	
	// create grid
	grid.resize(samples);
	for (int i = 0; i < samples; i++)
		grid[i] = i * h;
	
	// look for relevant data
	char filename[50];
	sprintf(filename, "dwf-N%d-K%g-l%d-k%g.arr", U.n, U.k, ln, kn);
	if (not load_array(array, filename))
	{
		// prepare derivative callback MFPTR-wrapper
		o2scl::ode_funct_mfptr<ForbiddenWave> derivs(this, &ForbiddenWave::derivs);
		
		// prepare faster adaptive Cash-Karp Runge-Kutta stepper
		o2scl::gsl_astep<decltype(derivs)> adapt_stepper;
		o2scl::gsl_rkck_fast<2,decltype(derivs)> stepper;
		adapt_stepper.set_step(stepper);
		
		// data arrays
		o2scl::ovector xg(samples);
		o2scl::omatrix yg(samples,2), ypg(samples,2), yerrg(samples,2);
		
		// create grid
		for (int i = 0; i < samples; i++)
			xg[i] = h * (samples - i - 1);
		
		// set FAR RIGHT boundary condition
		yg[0][0] =  ric_k_scaled(ln,kn*r);
		yg[0][1] = dric_k_scaled(ln,kn*r) * kn;
		
		// solve
		int nx = solve2(xg, samples, -h, yg, ypg, yerrg, adapt_stepper, derivs, RETURN_ON_OVERFLOW);
		
		// match modified irregular Riccati-Bessel function of the second kind
		if (nx != samples)
		{
			// matching grid point for ζ
			double rx = xg[nx];
			
			// evaluate "k" for the purpose of determining the scaling factor
			double eval_k = ric_k_scaled(ln,kn*rx);
			
			// scaling factor
			double scale = yg[nx][0] / eval_k;
			
			// scale "k" for the rest of the solutions
			for (int i = nx + 1; i < samples - 1; i++)
			{
				double r = xg[i];
				eval_k = ric_k_scaled(ln,kn*r);
				
				// scale "k"
				yg[i][0] = scale * eval_k;
			}
		}
		
		// copy array to class storage
		array = rArray(samples);
		array[0] = 0.;
		for (int i = 1; i < samples; i++)
			array[i] = yg[samples-i-1][0];
		
		// save for recyclation
		save_array(array, filename);
	}
	
	// setup the interpolator
	this->interpolator.set(grid, array, array.size());
}

ForbiddenWave ForbiddenWave::operator=(ForbiddenWave const & W)
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

int ForbiddenWave::derivs(double x, size_t nv, const o2scl::ovector_base& y, o2scl::ovector_base& dydx)
{
	dydx[0] = y[1];
	dydx[1] = (x == 0.) ? 0. : (2.*U(x) + ln*(ln+1)/(x*x)) * y[0] + 2*kn*y[1];
	return 0;
}

double ForbiddenWave::operator()(double x) const
{
	if (x > grid.back())
		
		// extrapolate
		return ric_k(ln,kn*x);
	
	else
		
		// interpolate
		return interpolator.interp(x) * exp(-kn*x);
}

void ForbiddenWave::toFile(const char* filename) const
{
	write_array(grid, array, filename);
}
