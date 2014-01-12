/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2014                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HEX_DISTORTING_POTENTIAL
#define HEX_DISTORTING_POTENTIAL

#include "arrays.h"
#include "hydrogen.h"
#include "specf.h"

/**
 * @brief U polynomial coefficients.
 * 
 * Can be generated in Maxima CAS.
 * @code
 *     R(n,l,r) := sqrt((2/n)^3 * (n-l-1)!/(2*n*(n+l)!)) * gen_laguerre(n-l-1,2*l+1,2*r/n) * (2*r/n)^l * exp(-r/n);
 *     Ucoeffs(n) := block ( [ tmp ],
 *         tmp : expand(integrate(R(n,0,r)^2*r^2*(1/r-1/x), r, x, inf)*exp(2*x/n)+1/x),
 *         reverse(makelist(float(coeff(tmp, x, i)), i, 0, length(tmp)-1))
 *     );
 *     Ucoeffs(1);
 *     Ucoeffs(2);
 *       ...
 * @endcode
 */
extern const rArrays Ucoeffs;

/**
 * @brief Distorting potential information.
 * 
 * Distorting potential is a spherically averaged potential of atomic particles.
 * Here, it is generated by the hydrogen proton and electron. The formula for
 * the distorting potential says
 * @f[
 *     U_n(r) = \int_{r}^{\infty} P_{n0}(r) \left( \frac{1}{r'} - \frac{1}{r} \right) \mathrm{d}r' \ ,
 * @f]
 * where @f$ P_{nl}(r) = rR_{nl}(r) @f$ and @f$ R_{nl}(r) @f$ is the common hydrogen
 * radial function (see Hydrogen).
 * 
 * This class is used to evaluate the distorting potential
 * @f$ U_n(r) @f$ and also to compute distorted waves @f$ \phi_{l_n}(k_n,r) @f$,
 * @f$ \eta_{l_n}(k_n, r) @f$, @f$ \theta_{l_n}(k_n, r) @f$ and
 * @f$ \zeta_{l_n}(k_n, r) @f$ from their respective definition differential
 * equations, i.e.
 * @f[
 *     \left(-\frac{\mathrm{d}^2}{\mathrm{d}r^2} + \frac{l_n(l_n+1)}{r^2} + 2U_n(r) \right)
 *     F_{l_n} (k_n, r) = k_n^2 F_{l_n} (k_n, r)
 * @f]
 * for allowed waves ( DistortedWave @f$ \phi_{l_n}(k_n,r) @f$ and IrregularWave
 * @f$ \eta_{l_n}(k_n, r) @f$ ) and
 * @f[
 *     \left(-\frac{\mathrm{d}^2}{\mathrm{d}r^2} + \frac{l_n(l_n+1)}{r^2} + 2U_n(r) \right)
 *     G_{l_n} (k_n, r) = -k_n^2 G_{l_n} (k_n, r)
 * @f]
 * for forbidden waves ( ForbiddenWave @f$ \theta_{l_n}(k_n, r) @f$ and HyperbolicWave
 * @f$ \zeta_{l_n}(k_n, r) @f$ ).
 */
class DistortingPotential : public RadialFunction<double>
{
public:
    // constructors
    // @{
    DistortingPotential() : n_(0), k_(0.), rmax_(0) {}
    DistortingPotential(int n, double rmax = 0.) : n_(n), k_(0.), rmax_(rmax) {}
    DistortingPotential(double k, double rmax = 0.) : n_(0), k_(k), rmax_(rmax) {}
    DistortingPotential(DistortingPotential const & U) : n_(U.n_), k_(U.k_), rmax_(U.rmax_) {}
    // @}
    
    /**
     * @brief Assignment
     */
    DistortingPotential operator= (DistortingPotential const& V);
    
    /**
     * @brief Comparison
     */
    bool operator== (DistortingPotential const & V) const;
    
    /**
     * @brief Evaluate the distorting potential.
     * 
     * Evaluate the distorting potential. At the moment, it is hard-coded,
     * and only for @f$ n = 1 @f$ state, which is
     * @f[
     *      U(r) = -\left(1 + \frac{1}{r}\right) \mathrm{e}^{-2r} \ .
     * @f]
     * @todo Implement runtime generation of distorting potentials.
     * @param x Coordinate where to evaluate.
     */
    double operator() (double x) const;
    
    /**
     * @brief Classical turning point.
     * 
     * This is only a compulsory pure-virtual function of the base class;
     * it has no physical meaning in DistortingPotential.
     */
    double getTurningPoint () const { return 0.; }
    
    /**
     * @brief Near-zero asymptotic behaviour.
     * @param x Evaluation radius.
     * @return Pair (y,k) that can be used to reconstruct the value of the potential,
     * @f[
     *     U(x) = y \cdot r^k\ .
     * @f]
     */
    std::pair<double,int> getZeroAsymptotic (double x) const { return std::make_pair(-1,-1); }
    
    /**
     * @brief Add multipole field potential to the distorting potential.
     * 
     * Returns
     * @f[
     *     U'(r) = U(r) + \frac{1}{r} \ .
     * @f]
     * The function handles correctly the input @f$ r = 0 @f$.
     */
    double plusMonopole(double x) const;
    
    /**
     * @brief Return the zero limit.
     * 
     * Return the asymptotic constant around zero,
     * @f[
     *     a = \lim_{r \rightarrow 0+} \left( \frac{1}{r} + U(r) \right) \ .
     * @f]
     */
    double getConstant() const;
    
    /**
     * @brief Return largest evaluated coordinate.
     * 
     * Return a radius sufficiently far from the atom. The radial
     * orbital ought to be small here.
     */
    double getFarRadius() const;
    
    void toFile(const char * filename) const;
    int n() const { return n_; }
    double k() const { return k_; }
    
private:
    int n_;        // principal quantum number of distorting state
    double k_;     // wavenumber of distorting state
    double rmax_;  // far radius
};

#endif
