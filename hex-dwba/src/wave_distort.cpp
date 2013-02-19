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
#include "ode.h"
#include "potential.h"
#include "wave_distort.h"

// -------------------------------------------------------------------------- //
//                                                                            //
//  DistortedWave class members                                               //
//                                                                            //
// -------------------------------------------------------------------------- //

DistortedWave DistortedWave::operator= (DistortedWave const& W)
{
	U = W.U;
	kn = W.kn;
	ln = W.ln;
	samples = W.samples;
	h = W.h;
	phase = W.phase;
	array = W.array;
	
	interpolator.set(W.grid, W.array, W.array.size());
	
	return *this;
}

int DistortedWave::derivs(double x, size_t nv, const o2scl::ovector_base& y, o2scl::ovector_base& dydx)
{
	dydx[0] = y[1];
	dydx[1] = (x == 0.) ? 0. : (-kn*kn + 2.*U(x) + ln*(ln+1)/(x*x)) * y[0];
	return 0;
}

DistortedWave::DistortedWave(double _kn, int _ln, DistortingPotential const & _U)
	: Evaluations(0), U(_U), kn(_kn), ln(_ln)
{
	// get far coordinate
	double r = U.getFarRadius();
	
	// determine discretization, use at least 1000 samples
	int N = 1000;				// N samples per wave length
	this->h = std::min (		// grid step
		2*M_PI/(N*kn),			//  -> N samples per wave length and
		r/1000					//  -> at least 1000 samples totally
	);
	this->samples = r/h + 1;	// with both boundaries
	
	// create grid
	grid = rArray(samples);
	for (int i = 0; i < samples; i++)
		grid[i] = h * i;
	
	// look for relevant data
	char filename[50];
	sprintf(filename, "dwr-N%d-K%g-l%d-k%g.arr", U.n, U.k, ln, kn);
	if (not load_array(this->array, filename, &this->phase))
	{
		// prepare derivative callback MFPTR-wrapper
		o2scl::ode_funct_mfptr<DistortedWave> derivs(this, &DistortedWave::derivs);
		
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
		
		// set (left) boundary conditions from asymptotics of Riccati-Bessel function
		yg[0][0] = h * 1e-5;	// η (h)
		yg[0][1] = (ln + 1) * 1e-5;	// η'(h)
		
		// solve
		solve2(xg, samples, h, yg, ypg, yerrg, adapt_stepper, derivs, NORMALIZE_ON_OVERFLOW);
		
		// compute phase shift; use last point
		double R = xg[samples-2];
		double val =  yg[samples-2][0];
		double der = ypg[samples-2][0];
		double D = der / val;
		double pilhalf = M_PI * ln / 2.;
		double tg_delta = (kn * cos(kn*R - pilhalf) - D * sin(kn*R - pilhalf)) /
							(kn * sin(kn*R - pilhalf) + D * cos(kn*R - pilhalf));
		phase = atan(tg_delta);
		
		// normalize ( far away should be η ≈ exp(iδ)sin(kr-πl/2+δ) )
		array = rArray(samples);
		double inverse_norm = 1./sqrt(val*val + der*der/(kn*kn));
		for (int ir = 0; ir < samples; ir++)
		{
			if (ir == 0)
				array[ir] = 0.;
			else
				array[ir] = yg[ir-1][0] * inverse_norm;
		}
		
		// save for recyclation
		save_array(array, filename, &phase);
	}
	
	// setup the interpolator
	this->interpolator.set(grid, array, array.size());
}

double DistortedWave::operator() (double x) const
{
	Evaluations++;
	
	// extrapolate using asymptotic form if far away
	if (x > grid.back())
		return sin(kn*x - ln*M_PI*0.5 + phase);
	
	// otherwise interpolate using cspline interpolator
	else
		return interpolator.interp(x);
}

double DistortedWave::farRadius() const
{
	return grid.back();
}

size_t DistortedWave::sampleCount() const
{
	return grid.size();
}

double DistortedWave::getPhase() const
{
	return phase;
}

Complex DistortedWave::getPhasef() const
{
	return Complex(cos(phase),sin(phase));
}

void DistortedWave::toFile(const char* filename) const
{
	write_array(grid, array, filename);
}
