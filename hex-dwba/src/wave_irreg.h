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

#ifndef HEX_IRREGULAR_WAVE
#define HEX_IRREGULAR_WAVE

#include "complex.h"
#include "potential.h"
#include "ode.h"
#include "specf.h"

/**
 * Irregular wave information.
 */
class IrregularWave : public RadialFunction<Complex>
{
public:
	
	IrregularWave(double _kn, int _ln, DistortingPotential const & _U);
	IrregularWave(IrregularWave const& W) { *this = W; }
	
	IrregularWave operator= (IrregularWave const& W);
	
	/**
	 * Evaluate distorted wave.
	 * \param x Coordinate where to evaluate.
	 * \return Real distorted wave.
	 */
	Complex operator()(double x) const;
	
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
	o2scl::interp_cspline<rArray> interpolator_re, interpolator_im;
	
	// distorted wave input parameters
	double kn;		// wavenumber of the distorted wave
	int ln;			// angular momentum of the distorted wave
	
	// distorted wave computed attributes
	int samples;	// sample count
	double h;		// discretization step
	
	// distorted wave data
	rArray grid;	// grid
	rArray array_re, array_im;	// samples
};

#endif
