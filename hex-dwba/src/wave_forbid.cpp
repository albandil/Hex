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
{
	kn = _kn;
	ln = _ln;
	U = _U;
	
	this->kn = _kn;
	this->ln = _ln;
	this->U = _U;
	
	Evaluations = 0;
	
	// get far coordinate
	double r = U.getFarRadius();
	
	// determine discretization
	samples = 10001;		// with both boundaries
	h = r / (samples - 1);	// grid step
	
	// create grid
	grid = rArray(samples);
	for (int i = 0; i < samples; i++)
		grid[i] = h * i;
	
	// create auxiliary complex arrry
	cArray carray(samples);
	
	// look for relevant data
	char filename[50];
	sprintf(filename, "dwf-N%d-K%g-l%d-k%g.arr", U.n, U.k, ln, kn);
	if (not load_array(carray, filename))
	{
		// prepare derivative callback MFPTR-wrapper
		o2scl::ode_funct_mfptr<ForbiddenWave> derivs(this, &ForbiddenWave::derivs);
		
		// prepare faster adaptive Cash-Karp Runge-Kutta stepper
		o2scl::gsl_astep<decltype(derivs)> adapt_stepper;
		o2scl::gsl_rkck_fast<2,decltype(derivs)> stepper;
		adapt_stepper.set_step(stepper);
		
		// data arrays
		o2scl::ovector xg(samples-1);
		o2scl::omatrix yg(samples-1,2), ypg(samples-1,2), yerrg(samples-1,2);
		
		// create grid
		for (int i = 0; i < samples - 1; i++)
			xg[i] = h * (i + 1);
		
		// set RIGHT boundary condition
		yg[0][0] = h;	// FIXME
		yg[0][1] = h;	// FIXME
		
		// solve
		solve2(xg, samples, h, yg, ypg, yerrg, adapt_stepper, derivs, NORMALIZE_ON_OVERFLOW);
		
		// save for recyclation
		save_array(carray, filename);
	}
	else
	{
		// TODO
	}
	
	// setup the interpolator
	this->interpolator_re.set(grid, array_re, array_re.size());
	this->interpolator_im.set(grid, array_im, array_im.size());
}

ForbiddenWave ForbiddenWave::operator=(ForbiddenWave const & W)
{
	kn = W.kn;
	ln = W.ln;
	U = W.U;
	grid = W.grid;
	array_re = W.array_re;
	array_im = W.array_im;
	samples = W.samples;
	h = W.h;
	
	interpolator_re.set(grid, array_re, array_re.size());
	interpolator_im.set(grid, array_im, array_im.size());
	
	return *this;
}

int ForbiddenWave::derivs(double x, size_t nv, const o2scl::ovector_base& y, o2scl::ovector_base& dydx)
{
	dydx[0] = y[1];
	dydx[1] = (x == 0.) ? 0. : (kn*kn + 2.*U(x) + ln*(ln+1)/(x*x)) * y[0];
	return 0;
}

void ForbiddenWave::toFile(const char* filename) const
{
	int N = array_re.size();
	cArray array(N);
	for (int i = 0; i < N; i++)
		array[i] = Complex(array_re[i], array_im[i]);
	
	write_array(grid, array, filename);
}

Complex ForbiddenWave::operator()(double x) const
{
	if (x > grid.back())
		
		// extrapolate
		return pow(Complex(0.,1.),ln+1) * ric_i(ln,kn*x) + pow(Complex(0.,1.),-ln) * 2. / M_PI * ric_k(ln,kn*x);
	
	else
		
		// interpolate
		return Complex(interpolator_re.interp(x), interpolator_im.interp(x));
}

