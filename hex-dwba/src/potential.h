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

#ifndef HEX_DISTORTING_POTENTIAL
#define HEX_DISTORTING_POTENTIAL

#include <o2scl/interp.h>
#include <o2scl/ode_funct.h>

#include "arrays.h"
#include "hydrogen.h"
#include "specf.h"

class DistortedWave;
class IrregularWave;
class ForbiddenWave;
class HyperbolicWave;

/**
 * \brief Distorting potential information.
 * 
 * Distorting potential is a spherically averaged potential of atomic particles.
 * Here, it is generated by the hydrogen proton and electron. The formula for
 * the distorting potential says
 * \f[
 *     U_n(r) = \int_{r}^{\infty} P_{n0}(r) \left( \frac{1}{r'} - \frac{1}{r} \right) \mathrm{d}r' \ ,
 * \f]
 * where \f$ P_{nl}(r) = rR_{nl}(r) \f$ and \f$ R_{nl}(r) \f$ is the common hydrogen
 * radial function (see Hydrogen).
 * 
 * This class is used to evaluate the distorting potential
 * \f$ U_n(r) \f$ and also to compute distorted waves \f$ \phi_{l_n}(k_n,r) \f$,
 * \f$ \eta_{l_n}(k_n, r) \f$, \f$ \theta_{l_n}(k_n, r) \f$ and
 * \f$ \zeta_{l_n}(k_n, r) \f$ from their respective definition differential
 * equations, i.e.
 * \f[
 *     \left(-\frac{\mathrm{d}^2}{\mathrm{d}r^2} + \frac{l_n(l_n+1)}{r^2} + 2U_n(r) \right)
 *     F_{l_n} (k_n, r) = k_n^2 F_{l_n} (k_n, r)
 * \f]
 * for allowed waves ( DistortedWave \f$ \phi_{l_n}(k_n,r) \f$ and IrregularWave
 * \f$ \eta_{l_n}(k_n, r) \f$ ) and
 * \f[
 *     \left(-\frac{\mathrm{d}^2}{\mathrm{d}r^2} + \frac{l_n(l_n+1)}{r^2} + 2U_n(r) \right)
 *     G_{l_n} (k_n, r) = -k_n^2 G_{l_n} (k_n, r)
 * \f]
 * for forbidden waves ( ForbiddenWave \f$ \theta_{l_n}(k_n, r) \f$ and HyperbolicWave
 * \f$ \zeta_{l_n}(k_n, r) \f$ ).
 */
class DistortingPotential : public RadialFunction<double>
{
public:
	// constructors
	// @{
	DistortingPotential() : n(0), k(0.) {}
	DistortingPotential(int _n) : n(_n), k(0.) {}
	DistortingPotential(double _k) : n(0), k(_k) {}
	DistortingPotential(HydrogenFunction const& psi) : n(psi.getN()), k(psi.getK()) {}
	DistortingPotential(DistortingPotential const& U) : n(U.n), k(U.k) {}
	// @}
	
	// socialize
	friend class DistortedWave;
	friend class IrregularWave;
	friend class ForbiddenWave;
	friend class HyperbolicWave;
	
	/**
	 * \brief Assignment
	 */
	DistortingPotential operator= (DistortingPotential const& V);
	
	/**
	 * \brief Comparison
	 */
	bool operator== (DistortingPotential const & V) const;
	
	/**
	 * \brief Evaluate the distorting potential.
	 * 
	 * Evaluate the distorting potential. At the moment, it is hard-coded,
	 * and only for \f$ n = 1 \f$ state, which is
	 * \f[
	 *      U(r) = -\left(1 + \frac{1}{r}\right) \mathrm{e}^{-2r} \ .
	 * \f]
	 * \todo Hard-code more distorting potentials.
	 * \todo Implement runtime generation of distorting potentials.
	 * \param x Coordinate where to evaluate.
	 */
	double operator() (double x) const;
	
	/**
	 * \brief Add multipole field potential to the distorting potential.
	 * 
	 * Returns
	 * \f[
	 *     U'(r) = U(r) + \frac{1}{r} \ .
	 * \f]
	 * The function handles correctly the input \f$ r = 0 \f$.
	 */
	double plusMonopole(double x) const;
	
	/**
	 * \brief Return the zero limit.
	 * 
	 * Return the asymptotic constant around zero,
	 * \f[
	 *     a = \lim_{r \rightarrow 0+} \left( \frac{1}{r} + U(r) \right)
	 * \f]
	 * At the moment, it is hard-coded, and only for \f$ n = 1 \f$ state
	 * (for which \f$ a = 1 \f$.
	 * \todo Hard-code more distorting potential limits.
	 * \todo Implement runtime generation of distorting potential limits.
	 */
	double getConstant() const;
	
	/**
	 * \brief Return largest evaluated coordinate.
	 * 
	 * Return a radius sufficiently far from the atom. The radial
	 * orbital ought to be small here. At the moment it is hard-coded
	 * with a value \f$ R_{\mathrm{far}} = 200 \f$, which is sufficient
	 * for \f$ n = 1 \f$ distorting state, but less for higher states.
	 * \todo Implement runtime computation of far radius based on
	 *       Laguerre polynomials behaviour.
	 */
	double getFarRadius() const;
	
	/**
	 * \brief Compute \f$ \phi_{l_n}(k_n, r) \f$.
	 * 
	 * Return the regular distorted wave (kn,ln). The DistortedWave
	 * \f$ \phi_{l_n}(k_n, r) \f$ has the following properties:
	 * \f[
	 *     \phi_{l_n}(k_n, 0) = 0 \ ,
	 * \f]
	 * \f[
	 *     \phi_{l_n}(k_n, r) \propto \hat{j}_{l_n}(k_n r) + t_{l_n}(k_n) \mathrm{i}\hat{h}_{l_n}^{(+)}(k_n r) \ .
	 * \f]
	 * \param kn Wave number of the wave.
	 * \param ln Angular momentum (partial wave).
	 */
	DistortedWave getDistortedWave(double kn, int ln) const;
	
	/**
	 * \brief Compute \f$ \eta_{l_n}(k_n, r) \f$
	 * 
	 * Return the irregular distorted wave (kn,ln). The IrregularWave
	 * \f$ \eta_{l_n}(k_n, r) \f$ has the following asymptotics:
	 * \f[
	 *     \eta_{l_n}(k_n, r) \propto \mathrm{i}\hat{h}_{l_n}^{(+)}(k_n r)
	 * \f]
	 * and is nonzero at origin.
	 * \param kn Wave number of the wave.
	 * \param ln Angular momentum (partial wave).
	 */
	IrregularWave getIrregularWave(double kn, int ln) const;
	
	/**
	 * \brief Compute \f$ \theta_{l_n}(k_n, r) \f$.
	 * 
	 * Return the regular forbidden wave (i·kn,ln). The properties of the ForbiddenWave are
	 * \f[
	 *     \theta_{l_n}(k_n, 0) = 0 \ ,
	 * \f]
	 * \f[
	 *     \theta_{l_n}(\mathrm{i}k_n, r) \propto \mathrm{i}^{l_n + 1} \hat{i}_{l_n}(k_n r)
	 *                   + t_{l_n}(\mathrm{i}k_n) \mathrm{i}^{-l_n} \frac{2}{\pi} \hat{k}_{l_n}(k_n r) \ .
	 * \f]
	 * \param kn Absolute value of the imaginary wave number of the wave.
	 * \param ln Angular momentum (partial wave).
	 */
	ForbiddenWave getForbiddenWave(double kn, int ln) const;
	
	/**
	 * \brief Compute \f$ \zeta_{l_n}(k_n, r) \f$.
	 * 
	 * Return the irregular forbidden wave (i·kn,ln). The asymptotics of
	 * the HyperbolicWave is
	 * \f[
	 *     \zeta_{l_n}(\mathrm{i}k_n, r) \propto \mathrm{i}^{-l_n} \frac{2}{\pi} \hat{k}_{l_n}(k_n r) \ .
	 * \f]
	 * \param kn Absolute value of the imaginary wave number of the wave.
	 * \param ln Angular momentum (partial wave).
	 */
	HyperbolicWave getHyperbolicWave(double kn, int ln) const;
	
	void toFile(const char * filename) const;
	
private:
	int n;		// principal quantum number of distorting state
	double k;	// wavenumber of distorting state
};

#include "wave_distort.h"
#include "wave_irreg.h"
#include "wave_forbid.h"
#include "wave_hyperb.h"

#endif
