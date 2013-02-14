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
{
	kn = _kn;
	ln = _ln;
	U = _U;
	
	// get far coordinate
	double r = U.getFarRadius();
	
	// determine discretization
	samples = 10001;		// with both boundaries
	h = r / (samples - 1);	// average grid step
	
	// allocate memory
	array = rArray(samples);
	
	// create logarithmic grid (densely spaced around origin)
	grid = logspace(1., r + 1., samples);
	for (auto g : grid)
		g -= 1.;	// shift back
	
	char etaname[50];
	sprintf(etaname, "dwi-%d-%g-%d-%g.arr", U.n, U.k, ln, kn);
	if (not load_array(array, etaname))
	{
		// data arrays
		o2scl::ovector xg(samples);
		o2scl::omatrix yg(samples,2), ypg(samples,2), yerrg(samples,2);
		
		// derivatives callback
		o2scl::ode_funct_mfptr<HyperbolicWave> derivs(this, & HyperbolicWave::derivs);
		
		// faster adaptive Cash-Karp Runge-Kutta stepper
		o2scl::gsl_astep<decltype(derivs)> adapt_stepper;
		o2scl::gsl_rkck_fast<2,decltype(derivs)> stepper;
		adapt_stepper.set_step(stepper);
		
		// create grid
		for (int i = 0; i < samples; i++)
		{
			xg[i] = h * (samples - i - 1);
			yg[i][0] = yg[i][1] = 0.;
			ypg[i][0] = ypg[i][1] = 0.;
			yerrg[i][0] = yerrg[i][1] = 0.;
		}
		
		// set (right) boundary conditions
		yg[0][0] = exp(-kn*xg[0]);
		yg[0][1] = -kn*exp(-kn*xg[0]);
		
		// The solution diverges for x -> 0 with the same asymptotics
		// as the irregular Coulomb wave. In the (approximate) interval < 0 , R₁ >
		// the function behaves as x^{-ℓ}. In the interval < R₁ , ∞ > it has an
		// oscillatory nature. To speed up the computation, the full differential
		// equation is solved until it overflows some hard-coded limit. From that
		// point it continues as a appropriatelly scaled irregular Coulomb wave.
			
		int n1 = solve2(xg, samples, -h, yg, ypg, yerrg, adapt_stepper, derivs, RETURN_ON_OVERFLOW);
		
		// match irregular Coulomb function
		gsl_sf_result F, Fp, G, Gp;
		double exp_F, exp_G;
		if (n1 != samples)
		{
			int nx = n1;					// matching grid point index
			double rx = xg[nx];				// matching grid point for η	
			bool warn = true;				// warn if something goes wrong
			
			// evaluate G for the purpose of determining the scaling factor
			int err = gsl_sf_coulomb_wave_FG_e(-1./kn, kn*rx, ln, 0, &F, &Fp, &G, &Gp, &exp_F, &exp_G);
			if (err != GSL_SUCCESS and err != GSL_ELOSS)
			{
				// if an error was encountered (other than loss of accuracy)
				fprintf (
					stderr,
					"[distort_irregular_allowed] Error %d (\"%s\") while evaluating G[%d](%g,%g).\n",
					err, gsl_strerror(err), ln, kn, rx
				);
				abort();
			}
			else
			{
				// scale G for the rest of the solutions
				double const_re = yg[nx][0] / G.val;	// scaling factor for real part
				for (int i = nx + 1; i < samples - 1; i++)
				{
					double r = xg[i];
					int err = gsl_sf_coulomb_wave_FG_e(-1./kn, kn*r, ln, 0, &F, &Fp, &G, &Gp, &exp_F, &exp_G);
					if (err == GSL_SUCCESS or err == GSL_ELOSS)
					{
						// take both precise and inaccurate solution, but warn in the latter case
						if (err == GSL_ELOSS)
						{
							if (warn)
							{
								fprintf(stderr, "\t\t\tG[%d](%g,%g) inaccurate.\n", ln, kn, r);
								warn = false;
							}
							
							abort();
						}
						
						// scale G if needed
						if (i > n1)
							yg[i][0] = const_re * G.val;
					}
					else
					{
						// some other more peculiar error
						if (warn)
						{
							fprintf (
								stderr,
								"[distort_irregular_allowed] Error %d (\"%s\") while evaluating G[%d](%g,%g).\n",
								err, gsl_strerror(err), ln, kn, r
							);
							warn = false;
						}
						
						abort();
						
						yg[i][0] = 0.;
					}
				}
			}
		}
		
		// save
		array = rArray(samples);
		array[0] = 0.;
		for (int i = 1; i < samples; i++)
			array[i] = yg[samples-i-1][0];
		
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
	dydx[1] = (x == 0.) ? 0. : (kn*kn + 2.*U(x) + ln*(ln+1)/(x*x)) * y[0];
	return 0;
}

double HyperbolicWave::operator()(double x) const
{
	if (x > grid.back())
		
		// extrapolate
		return 2. / M_PI * ric_k(ln,kn*x);
	
	else
		
		// interpolate
		return interpolator.interp(x);
}

Complex HyperbolicWave::getPhasef() const
{
	return pow(Complex(0.,1.), -ln);
}

double HyperbolicWave::getPhase() const
{
	return M_PI * (ln % 4) / 2.;
}

void HyperbolicWave::toFile(const char* filename) const
{
	write_array(grid, array, filename);
}

