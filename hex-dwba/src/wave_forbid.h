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
 * \brief Regular forbidden wave information.
 * 
 * The forbidden wave \f$ \theta_{l_n}(k_n, r) \f$ is a solution of the equation
 * \f[
 *     \left(-\frac{\mathrm{d}^2}{\mathrm{d}r^2} + \frac{l_n(l_n+1)}{r^2} + 2U_n(r) \right)
 *     \theta_{l_n} (k_n, r) = -k_n^2 \theta_{l_n} (k_n, r) 
 * \f]
 * which is regular at origin and which satisfies the boundary condition
 * \f[
 *     \theta_{l_n} (k_n, r) \propto \mathrm{i}^{l_n + 1} \hat{i}_{l_n}(k_n r)
 *          + t_{l_n}(k_n) \mathrm{i}^{-l_n} \frac{2}{pi} \hat{k}_{l_n}(k_n r) \ .
 * \f]
 * Most of the work is being done in the constructor. The radial interval between 
 * the origin and the far radius of DistortingPotential is divided into an equidistant
 * grid. Spacing doesn't depend on energy, as this dependence is unclear.
 * 
 * \todo Forbidden wave discretization.
 * 
 * Being regular function, the initial condition is stated at origin as
 * \f[
 *     \theta_{l_n} (k_n, h) = \theta_{l_n}' (k_n, h) = h
 * \f]
 * and normalized on the fly to not exceed specified threshold. When the solution
 * is found, it is renormalized using the correct asymptotics:
 */
class ForbiddenWave : public RadialFunction<Complex>
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
	Complex operator()(double x) const;
	
	/**
	 * \brief Return derivatives from the distorted wave equation.
	 * 
	 * The defining differential equation can be recast into a set of two linear
	 * differential euqations of first order:
	 * \f[
	 *     \frac{\mathrm{d}}{\mathrm{d}r} \theta_{l_n}(k_n,r) = \theta_{l_n}'(k_n,r) \ ,
	 * \f]
	 * \f[
	 *     \frac{\mathrm{d}}{\mathrm{d}r} \theta_{l_n}'(k_n,r) = \left( \frac{l_n(l_n+1)}{r^2}
	 *     + 2U(r) + k^2 \right) \theta_{l_n}(k_n,r) \ .
	 * \f]
	 * This function computes vector \f$ ( \theta_{l_n}'(k_n, r), \theta_{l_n}''(k_n, r) ) \f$
	 * from the vector \f$ ( \theta_{l_n}(k_n, r), \theta_{l_n}'(k_n, r) ) \f$ and the independent
	 * coordinate \f$ x \f$.
	 */
	int derivs(double x, size_t nv, const o2scl::ovector_base& y, o2scl::ovector_base& dydx);
	
	/**
	 * Export data to file using \ref write_array.
	 */
	void toFile(const char * filename) const;
	
	/// (debuging parameter) number of evaluations
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
