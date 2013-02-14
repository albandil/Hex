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

#ifndef HEX_HYPERBOLIC_WAVE
#define HEX_HYPERBOLIC_WAVE

#include "complex.h"
#include "potential.h"
#include "ode.h"
#include "specf.h"

/**
 * Irregular forbidden wave information.
 */
class HyperbolicWave : public RadialFunction<double>
{
public:
	
	HyperbolicWave(double _kn, int _ln, DistortingPotential const & _U);
	HyperbolicWave(HyperbolicWave const& W) { *this = W; }
	
	HyperbolicWave operator= (HyperbolicWave const& W);
	
	/**
	 * Evaluate distorted wave.
	 * \param x Coordinate where to evaluate.
	 * \return Real distorted wave.
	 */
	double operator()(double x) const;
	
	/**
	 * Return the phase factor.
	 */
	Complex getPhasef() const;
	
	/**
	 * Return the phase.
	 */
	double getPhase() const;
	
	/**
	 * Return derivatives from the distorted wave equation.
	 */
	int derivs(double x, size_t nv, const o2scl::ovector_base& y, o2scl::ovector_base& dydx);
	
	void toFile(const char * filename) const;
	
	mutable unsigned Evaluations;
	
private:
	
	// distorting potential
	DistortingPotential U;
	
	// interpolator
	o2scl::interp_cspline<rArray> interpolator;
	
	// distorted wave input parameters
	double kn;		// wavenumber of the distorted wave
	int ln;			// angular momentum of the distorted wave
	
	// distorted wave computed attributes
	int samples;	// sample count
	double h;		// discretization step
	
	// distorted wave data
	rArray grid;	// grid
	rArray array;	// samples
};

#endif
