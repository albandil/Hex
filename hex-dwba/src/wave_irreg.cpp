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
#include "wave_irreg.h"

// -------------------------------------------------------------------------- //
//                                                                            //
//  IrregularWave class members                                               //
//                                                                            //
// -------------------------------------------------------------------------- //

IrregularWave::IrregularWave(double _kn, int _ln, DistortingPotential const & _U)
	: Evaluations(0), U(_U), kn(_kn), ln(_ln)
{
	// get far coordinate
	double r = U.getFarRadius();
	
	// determine discretization
	int N = 1000;				// N samples per wave length
	this->h = std::min (		// grid step
		2*M_PI/(N*kn),			//  -> N samples per wave length and
		r/1000					//  -> at least 1000 samples total
	);
	samples = r/h + 1;			// with both boundaries
	
	cArray array; // auxiliary array for HDF I/O
	
	// create grid
	grid = rArray(samples);
	for (int i = 0; i < samples; i++)
		grid[i] = h * i;
	
	char etaname[50];
	sprintf(etaname, "dwi-N%d-K%g-l%d-k%g.dwf", U.n, U.k, ln, kn);
	if (not array.hdfload(etaname))
	{
		// data arrays
		o2scl::ovector xg(samples);
		o2scl::omatrix y1g(samples,2), y1pg(samples,2), y1errg(samples,2);	// Re η
		o2scl::omatrix y2g(samples,2), y2pg(samples,2), y2errg(samples,2);	// Im η
		
		// derivatives callback
		o2scl::ode_funct_mfptr<IrregularWave> derivs(this, & IrregularWave::derivs);
		
		// faster adaptive Cash-Karp Runge-Kutta stepper
		o2scl::gsl_astep<decltype(derivs)> adapt_stepper;
		o2scl::gsl_rkck_fast<2,decltype(derivs)> stepper;
		adapt_stepper.set_step(stepper);
		
		// create grid
		for (int i = 0; i < samples; i++)
			xg[i] = h * (samples - i - 1);
		
		// set (right) boundary conditions from asymptotics of Riccati-Bessel function
		y1g[0][0] = -ric_n(ln,kn*xg[0]);		// Re η(R)
		y1g[0][1] = -kn*dric_n(ln,kn*xg[0]);	// Re η'(R)
		y2g[0][0] = ric_j(ln,kn*xg[0]);			// Im η(R)
		y2g[0][1] = kn*dric_j(ln,kn*xg[0]);		// Im η'(R)
		
		// The solution diverges for x -> 0 with the same asymptotics
		// as the irregular Coulomb wave. In the (approximate) interval < 0 , R₁ >
		// the function behaves as x^{-ℓ}. In the interval < R₁ , ∞ > it has an
		// oscillatory nature. To speed up the computation, the full differential
		// equation is solved until it overflows some hard-coded limit. From that
		// point it continues as a appropriatelly scaled irregular Coulomb wave.
			
		int n1 = solve2(xg, samples, -h, y1g, y1pg, y1errg, adapt_stepper, derivs, RETURN_ON_OVERFLOW);
		int n2 = solve2(xg, samples, -h, y2g, y2pg, y2errg, adapt_stepper, derivs, RETURN_ON_OVERFLOW);
		
		int nx = std::min(n1,n2);		// matching grid point index
		double rx;						// matching grid point for Re η	
		
		// match irregular Coulomb function
		if (nx < samples - 1)
		{
			rx = xg[nx];
			
			// evaluate G so that we can scale missing values from it
			gsl_sf_result Fx, Fpx, Gx, Gpx, F, Fp, G, Gp;
			double exp_Fx, exp_Gx, exp_F, exp_G;
			int err = gsl_sf_coulomb_wave_FG_e(-1./kn, kn*rx, ln, 0, &Fx, &Fpx, &Gx, &Gpx, &exp_Fx, &exp_Gx);
			
			// only use the value if valid
			if (err == GSL_SUCCESS)
			{
				// scale the rest of the solutions
				for (int i = nx + 1; i < samples - 1; i++)
				{
					if (gsl_sf_coulomb_wave_FG_e(-1./kn, kn*xg[i], ln, 0, &F, &Fp, &G, &Gp, &exp_F, &exp_G) == GSL_SUCCESS)
					{
						double scale = G.val / Gx.val;
						
						// scale G if needed
						if (i > n1)
							y1g[i][0] = y1g[nx][0] * scale;
						if (i > n2)
							y2g[i][0] = y2g[nx][0] * scale;
					}
					else
					{
						// unsuccessfull evaluation
						nx = i - 1;
						break;
					}
				}
			}
		}
		
		// if still not done ...
		if (nx < samples - 1)
		{
			// ... match the "origin" solution 1/r^ℓ
			for (int i = nx + 1; i < samples - 1; i++)
			{
				double scale =  pow(xg[nx]/xg[i], ln);
				
				if (i > n1)
					y1g[i][0] = y1g[nx][0] * scale;
				if (i > n2)
					y2g[i][0] = y2g[nx][0] * scale;
			}
		}
		
		// combine real and imaginary part
		array.resize(samples);
		array_re.resize(samples);
		array_im.resize(samples);
		array[0] = array_im[0] = array_re[0] = 0.;
		for (int i = 1; i < samples; i++)
		{
			array_re[i] = y1g[samples-i-1][0];
			array_im[i] = y2g[samples-i-1][0];
			array[i] = Complex(array_re[i], array_im[i]);
		}
		
		array.hdfsave(etaname);
	}
	else
	{
		array_re = rArray(samples);
		array_im = rArray(samples);
		for (int i = 0; i < samples; i++)
		{
			array_re[i] = array[i].real();
			array_im[i] = array[i].imag();
		}
	}
	
	// setup the interpolator
	this->interpolator_re.set(grid, array_re, array_re.size());
	this->interpolator_im.set(grid, array_im, array_im.size());
}

IrregularWave IrregularWave::operator=(const IrregularWave& W)
{
	U = W.U;
	kn = W.kn;
	ln = W.ln;
	samples = W.samples;
	h = W.h;
	array_im = W.array_im;
	array_re = W.array_re;
	
	interpolator_im.set(W.grid, W.array_im, W.array_im.size());
	interpolator_re.set(W.grid, W.array_re, W.array_re.size());
	
	return *this;
}

int IrregularWave::derivs(double x, size_t nv, const o2scl::ovector_base& y, o2scl::ovector_base& dydx)
{
	dydx[0] = y[1];
	dydx[1] = (x == 0.) ? 0. : (-kn*kn + 2.*U(x) + ln*(ln+1)/(x*x)) * y[0];
	return 0;
}

Complex IrregularWave::operator()(double x) const
{
	Evaluations++;
	
	// extrapolate using asymptotic form if far away
	if (x > grid.back())
		return Complex(-ric_n(ln,kn*x),ric_j(ln,kn*x));
	
	// otherwise interpolate using cspline interpolator
	else
		return Complex(interpolator_re.interp(x), interpolator_im.interp(x));
}

void IrregularWave::toFile(const char* filename) const
{
	int N = array_re.size();
	cArray array(N);
	for (int i = 0; i < N; i++)
		array[i] = Complex(array_re[i], array_im[i]);
	
	write_array(grid, array, filename);
}
