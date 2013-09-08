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
 * Regular forbidden wave information.
 * 
 * The hyperbolic wave \f$ \theta_{l_n}(k_n, r) \f$ is a solution of the equation
 * \f[
 *     \left(-\frac{\mathrm{d}^2}{\mathrm{d}r^2} + \frac{l_n(l_n+1)}{r^2} + 2U_n(r) \right)
 *     \theta_{l_n} (k_n, r) = -k_n^2 \theta_{l_n} (k_n, r) 
 * \f]
 * which is regular at origin and which satisfies the boundary condition
 * \f[
 *     \theta_{l_n} (k_n, r) \propto \hat{i}_{l_n}(k_n r) \ .
 * \f]
 * Most of the work is being done in the constructor. The radial interval between 
 * the origin and the far radius of DistortingPotential is divided into an equidistant
 * grid. Spacing doesn't depend on energy, as this dependence is unclear.
 * 
 * \todo Forbidden wave discretization.
 * 
 */
class HyperbolicWave : public RadialFunction<double>
{
public:
	
	HyperbolicWave(double _kn, int _ln, DistortingPotential const & _U);
	HyperbolicWave(HyperbolicWave const& W) { *this = W; }
	
	~HyperbolicWave();
	
	HyperbolicWave operator= (HyperbolicWave const& W);
	
	/**
	 * Evaluate distorted wave.
	 * \param x Coordinate where to evaluate.
	 * \return Real distorted wave.
	 */
	double operator()(double x) const;
	
    /// Classical turning point.
    double getTurningPoint () const { return r0; }
    
    /// Near-zero asymptotic behaviour.
    std::pair<double,int> getZeroAsymptotic (double x) const;
    
	/**
	 * \brief Return derivatives from the distorted wave equation.
	 * 
	 * The defining differential equation can be recast into a set of two linear
	 * differential euqations of first order:
	 * \f[
	 *     \frac{\mathrm{d}}{\mathrm{d}r} \tilde{theta}_{l_n}(k_n,r) = \tilde{\theta}_{l_n}'(k_n,r) \ ,
	 * \f]
	 * \f[
	 *     \frac{\mathrm{d}}{\mathrm{d}r} \tilde{\theta}_{l_n}'(k_n,r) = \left( \frac{l_n(l_n+1)}{r^2}
	 *     + 2U(r)\right) \tilde{theta}_{l_n}(k_n,r) - 2k\tilde{\theta}_{l_n}'(k_n,r) \ .
	 * \f]
	 * This function computes vector \f$ ( \theta_{l_n}'(k_n, r), \theta_{l_n}''(k_n, r) ) \f$
	 * from the vector \f$ ( \theta_{l_n}(k_n, r), \theta_{l_n}'(k_n, r) ) \f$ and the independent
	 * coordinate \f$ x \f$.
	 */
	int derivs(double x, size_t nv, const o2scl::ovector_base& y, o2scl::ovector_base& dydx);
    int derivs0(double x, size_t nv, const o2scl::ovector_base& y, o2scl::ovector_base& dydx);
	
	void toFile(const char * filename) const;
	
	void scale(bool s);
	
	mutable unsigned Evaluations;
	
private:
	
	/// distorting potential
	DistortingPotential U;
	
	/// interpolator
	gsl_interp *interpolator, *interpolator0;
	
	/// wavenumber of the distorted wave
	double kn;
	
	/// angular momentum of the distorted wave
	int ln;
	
	/// sample count
	int samples, samples0;
	
	/// discretization step
	double h, h0;
	
	/// grid
	rArray grid, grid0;
	
	/// samples
	rArray array, array0;
	
	bool Scaled;
    
    /// classical turning point, far radius
    double r0, rf;
};

#endif
