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

#ifndef HEX_DISTORTED_WAVE
#define HEX_DISTORTED_WAVE

#include "complex.h"
#include "potential.h"
#include "specf.h"

/**
 * \brief Distorted wave information.
 * 
 * The distorted wave \f$ \phi_{l_n}(k_n, r) \f$ is a solution of the equation
 * \f[
 *     \left(-\frac{\mathrm{d}^2}{\mathrm{d}r^2} + \frac{l_n(l_n+1)}{r^2} + 2U_n(r) \right)
 *     \phi_{l_n} (k_n, r) = k_n^2 \phi_{l_n} (k_n, r) 
 * \f]
 * which is regular at origin and which satisfies the boundary condition
 * \f[
 *     \phi_{l_n} (k_n, r) \propto \hat{j}_{l_n}(k_n r) + t_{l_n}(k_n) \mathrm{i}\hat{h}_{l_n}^{(+)}(k_n r) \ .
 * \f]
 * Most of the work is being done in the constructor. The radial interval between 
 * the origin and the far radius of DistortingPotential is divided into an equidistant
 * grid (with spacing adjusted to the energy of projectile),
 * on which the solution is sought. Integration starts from \f$ r = h \f$, where \f$ h \f$
 * is a small number compared to grid size. The initial conditions here are chosen as
 * \f[
 *     \phi_{l_n}(k_n, h) \simeq h \cdot 10^{-5} \ ,
 * \f]
 * \f[
 *     \phi_{l_n}'(k_n, h) \simeq (l_n + 1) \cdot 10^{-5} \ ,
 * \f]
 * which is in accord with the behaviour of Ricatti-Bessel functions near origin. Resulting
 * solution is real and using the standard formula
 * \f[
 *     \tan \delta_{l_n} = \frac{k_n \cos(k_n R - l_n\pi/2) - D \sin(k_n R - l_n\pi/2)}
                                {k_n \sin(k_n R - l_n\pi/2) + D \cos(k_n R - l_n\pi/2)}
 * \f]
 * the phase shift is determined. This number sets the phase factor \f$ \mathrm{e}^{\mathrm{i}\delta_{l_n}} \f$
 * of the whole wave.
 */
class DistortedWave : public RadialFunction<double>
{
public:
	
	// constructors
	// @{
	DistortedWave(double _kn, int _ln, DistortingPotential const & _U);
	DistortedWave(DistortedWave const& W) { *this = W; }
	// @}
	
	// assignment
	DistortedWave operator= (DistortedWave const& W);
	
	/**
	 * Evaluate distorted wave.
	 */
	double operator() (double x) const;
	
	/**
	 * Return the phase factor \f$ \mathrm{e}^{\mathrm{i}\delta_{l_n}} \f$.
	 */
	Complex getPhasef() const;
	
	/**
	 * Return the phase \f$ \delta_{l_n} \f$.
	 */
	double getPhase() const;
	
	/// Wavenumber.
	double k() const { return kn; }
	
	/// Angular momentum.
	int l() const { return ln; }
	
	/**
	 * \brief Return derivatives from the distorted wave equation.
	 * 
	 * The defining differential equation can be recast into a set of two linear
	 * differential euqations of first order:
	 * \f[
	 *     \frac{\mathrm{d}}{\mathrm{d}r} \phi_{l_n}(k_n,r) = \phi_{l_n}'(k_n,r) \ ,
	 * \f]
	 * \f[
	 *     \frac{\mathrm{d}}{\mathrm{d}r} \phi_{l_n}'(k_n,r) = \left( \frac{l_n(l_n+1)}{r^2}
	 *     + 2U(r) - k^2 \right) \phi_{l_n}(k_n,r) \ .
	 * \f]
	 * This function computes vector \f$ ( \phi_{l_n}'(k_n, r), \phi_{l_n}''(k_n, r) ) \f$
	 * from the vector \f$ ( \phi_{l_n}(k_n, r), \phi_{l_n}'(k_n, r) ) \f$ and the independent
	 * coordinate \f$ x \f$.
	 */
	int derivs(double x, size_t nv, const o2scl::ovector_base& y, o2scl::ovector_base& dydx);
	
	/**
	 * Export data to file using \ref write_array.
	 */
	void toFile(const char * filename) const;
	
	/**
	 * Return radius from which the asymptotic form is used.
	 */
    double farRadius() const;
	
	size_t sampleCount() const;
	
	/// (debuging parameter) number of evaluations
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
	double phase;	// phase shift
	rArray grid;	// grid
	rArray array;	// samples
};

#endif
