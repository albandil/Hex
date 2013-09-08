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
 * \* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HEX_DWBA_MULTI
#define HEX_DWBA_MULTI

#include "arrays.h"
#include "hydrogen.h"
#include "potential.h"
#include "chebyshev.h"
#include "specf.h"

/**
 * @brief Base class for all multipole functions (@f$ \varphi_{\alpha\beta}^\lambda @f$).
 */
class PhiFunction : public RadialFunction<double>
{
public:
    
    virtual double operator() (double x) const = 0;
    virtual double getTurningPoint() const = 0;
    virtual std::pair<double,int> getZeroAsymptotic (double x) const = 0;
    virtual bool isZero() const = 0;
    virtual double k() const = 0;
};

/**
 * \brief Multipole function (@f$ \varphi_{\alpha\beta}^\lambda @f$) of the "direct" type.
 * 
 * Auxiliary class holding information about intermediate integration.
 * When evaluated, it shall return
 * \f[
 *     \int 
 *         \psi_n(r_1) \left(
 *             \frac{r_<^\lambda}{r_>^{\lambda+1}} 
 *           - \frac{\delta_{\lambda 0}}{r_2} 
 *           - \delta_{\lambda 0} U_\alpha(r_2)
 *         \right) \psi_\alpha(r_1) \mathrm{d}r_1  \ .
 * \f]
 * 
 * The integral is split into two integrals with a fixed integrand (no \f$ r_< \f$,
 * \f$ r_> \f$ anymore). The constructor tries to approximate the integrands
 * using a standard Chebyshev approximation. If successfull, the expansions are
 * trivially integrated and stored for use in evaluation by the operator().
 * 
 * If the maximal given number of allowed %Chebyshev nodes to use is hit, then
 * a different method is used. Again, a %Chebyshev approximation is sought, but
 * now of the whole integral as a function of \f$ r_2 \f$. The integral is done
 * either by ordinary integration along the real radius or, for large momenta,
 * by integration in complex plane along a contour where the oscillating free
 * intermediate state quickly damps.
 * 
 * \param psin Intermediate hydrogen function (bound or free).
 * \param lam Multipole moment \f$ \lambda \ge 0 \f$.
 * \param U Distorting potential.
 * \param psi Initial (or final) atomic bound state.
 * \param cblimit Maximal number of %Chebyshev evaluation points used. If the
 *     routine doesn't find a satisfactory approximation even after evaluating
 *     the function at 'cblimit' %Chebyshev roots, the integrand is assumed
 *     to be highly oscillatory and an alternative method is used. Minus one disables
 *     the limit.
 */
class PhiFunctionDir : public PhiFunction
{
public:
    
    PhiFunctionDir (
        HydrogenFunction const & psin, 
        int lam, 
        DistortingPotential const & U, 
        HydrogenFunction const & psi,
        int cblimit = 1024
    );
    
    /**
     * @brief Evaluate the function.
     * 
     * This function will interpolate the stored arrays to return
     * as precise result as possible.
     */
    double operator() (double x) const;
    
    /**
     * @brief Classical turning point.
     * 
     * This routine returns the effective classical turning point for the
     * multipole function, though it is not a radial partial wave for which
     * a turning point were defined.
     */
    double getTurningPoint() const;
    
    /**
     * @brief Near-zero asymptotic behaviour.
     * 
     * In the vicinity of @f$ r = 0 @f$ is
     * @f[
     *     \varphi_{\alpha\beta}^\lambda(r_2)
     *     =
     *     r_2^\lambda \Phi_{\alpha\beta}^\lambda
     *     -
     *     \frac{(2\lambda+1)\overline{\psi}_\alpha\overline{\psi}_\beta}{\left(L_\alpha+L_\beta+\frac{5}{2}\right)^2 - \left(\lambda+\frac{1}{2}\right)^2}
     *     r^{L_\alpha+L_\beta+2} \ ,
     * @f]
     * where
     * @f[
     *     \Phi_{\alpha\beta}^\lambda = \int_0^\infty \psi_\alpha(r) \psi_\beta(r) r^{-\lambda-1} \mathrm{d}r \ ,
     * @f]
     * @f[
     *     \overline{\psi}_\alpha = \lim_{r \rightarrow 0} r^{L_\alpha + 1} \psi_\alpha(r) \ .
     * @f]
     * Thus, the function is asymptotically equal to @f$ r^\lambda \Psi_{\alpha\beta}^\lambda @f$,
     * as the multipole moment @f$ \lambda @f$ needs to be smaller than @f$ L_\alpha + L_\beta + 2 @f$
     * in order to make @f$ \varphi_{\alpha\beta}^\lambda @f$ integrable.
     */
    std::pair<double,int> getZeroAsymptotic (double x) const;
    
    /**
     * @brief Whether the function is identical zero.
     * 
     * If the @f$ \alpha @f$ and @f$ \beta @f$ states correspond to the states
     * used to compute the distorting potential @f$ U(r) @f$, then the function
     * is identically zero due to exact subtraction of the potentials.
     */
    bool isZero() const { return Zero; }
    
    /**
     * @brief Return effective wavenumber.
     * 
     * Though not being a proper radial part of a partial wave, still the oscillations
     * of the multipole function have an effective wave number that can be computed
     * from the wave numbers of the function in the integrand.
     */
    double k() const { return Wavenum; }
    
private:
    
    /// Compute Chebyshev expansion coefficients by evaluating the integrands.
    void tryFullRealChebyshev (
        int cblimit,
        bool & Cheb_L_conv,
        bool & Cheb_mLm1_conv
    );
    
    void tryFrontRealChebyshev (
        int cblimit,
        bool & Cheb_L_conv,
        bool & Cheb_mLm1_conv
    );
    
    /// Compute Chebyshev expansion coefficients by computing the complex integral.
    void tryFullComplexChebyshev (
        int cblimit,
        bool & Cheb_L_conv,
        bool & Cheb_mLm1_conv
    );
    
    HydrogenFunction psin;
    HydrogenFunction psi;
    
    int Lam;
    DistortingPotential U;
    
    /// whether psin == psi
    bool Diag;
    
    /// whether the integral is identical zero for all "x2"
    bool Zero;
    
    /// Wavenumber of the intermediate state.
    double Wavenum;
    
    /// Use frontal approximation.
    bool UseFront;
    
    /// Chebyshev approximations of the integrand \f$ \psi_n(r) r^\lambda \psi_i(r) \f$.
    Chebyshev<double,double> Cheb_L, Cheb0_L;
    
    /// Chebyshev approximations of the integrand \f$ \psi_n(r) r^{-\lambda-1} \psi_i(r) \f$.
    Chebyshev<double,double> Cheb_mLm1, Cheb0_mLm1;
    
    Chebyshev<double,double> ACoeff, BCoeff;
    int ATail, BTail;
    
    int Cheb_mLm1_tail, Cheb0_mLm1_tail;
    int Cheb_L_tail, Cheb0_L_tail;
    
    double r0;
};

#endif
