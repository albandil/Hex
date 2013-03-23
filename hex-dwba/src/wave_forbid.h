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

#ifndef HEX_FORBIDDEN_WAVE
#define HEX_FORBIDDEN_WAVE

#include "complex.h"
#include "potential.h"
#include "ode.h"
#include "specf.h"

/**
 * \brief Irregular forbidden wave information.
 * 
 */
class ForbiddenWave : public RadialFunction<double>
{
public:
	
	// constructors
	// @{
	ForbiddenWave(double _kn, int _ln, DistortingPotential const & _U);
	ForbiddenWave(ForbiddenWave const& W) { *this = W; }
	// @}
	
	// assignment
	ForbiddenWave operator= (ForbiddenWave const& W);
	
	/**
	 * Evaluate distorted wave.
	 */
	double operator()(double x) const;
	
	/**
	 * \brief Return derivatives from the distorted wave equation.
	 * 
	 */
	int derivs(double x, size_t nv, const o2scl::ovector_base& y, o2scl::ovector_base& dydx);
	
	/**
	 * Export data to file using \ref write_array.
	 */
	void toFile(const char * filename) const;
	
	void scale(bool s);
	
	/// (debuging parameter) number of evaluations
	mutable unsigned Evaluations;
	
private:
	
	/// distorting potential
	DistortingPotential U;
	
	/// interpolator
	o2scl::interp_cspline<rArray> interpolator;
	
	/// wavenumber of the distorted wave
	double kn;
	
	/// angular momentum of the distorted wave
	int ln;
	
	/// sample count
	int samples;
	
	/// discretization step
	double h;
	
	/// grid
	rArray grid;
	
	/// samples
	rArray array;
	
	bool Scaled;
};

#endif
