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

#ifndef HEX_PWBA2_MULTIPOT
#define HEX_PWBA2_MULTIPOT

#include <gsl/gsl_interp.h>

#include "arrays.h"

/**
 * @brief Multipole potential.
 * 
 * Defined as
 * @f[
 *     V_{\alpha \beta}^\lambda (r) = \int_0^\infty
 *     \left<\psi_\alpha(r')\right|
 *     \left(\frac{r_<^\lambda}{r_>^{\lambda+1}} - \frac{1}{r}\right)
 *     \left|\psi_\beta(r')\right> \,
 *     \mathrm{d}r' \ ,
 * @f]
 * where @f$ \psi_\alpha(r') @f$ is either @f$ P_\alpha(r') @f$ (the hydrogenic
 * orbital for bound states) or @f$ F_\alpha(r') @f$ (the Coulomb partial wave
 * for free states).
 */
class MultipolePotential
{
    public:
        
        /**
         * @brief Bound-bound constructor
         * 
         * This constructor will setup the bound-bound matrix element of the multipole
         * potential. The integral to be computed is
         * @f[
         *     V_{ab}^\lambda(r') = \int_0^\infty R_{N_a L_a}(r) \left(
         *         \frac{r_<^\lambda}{r_>^{\lambda+1}} - \frac{1}{r'}
         *     \right) R_{N_b L_b}(r) r^2 \, \mathrm{d}r \ ,
         * @f]
         * where @f$ R_{nl}(r) @f$ are the radial hydrogen functions (see @ref Hydrogen::P
         * for definition of these functions).
         * The integrand is integrated term by term using the identity
         * @f[
         *     \int_0^x t^{a-1} \mathrm{e}^{-t} \mathrm{d}t = \Gamma(a) P(a,x) \ ,
         * @f]
         * @f[
         *     \int_x^\infty t^{a-1} \mathrm{e}^{-t} \mathrm{d}t = \Gamma(a) Q(a,x) \ ,
         * @f]
         * where @f$ Q(a,x) \equiv \Gamma(a,x)/\Gamma(a) @f$ is the normalized incomplete
         * Gamma function and @f$ P(a,x) = 1 - Q(a,x) @f$ is the complementary scaled
         * incomplete Gamma function.
         */
        MultipolePotential (int lambda, int Na,    int La, int Nb,    int Lb);
        
        /**
         * @brief Bound-bound constructor
         * 
         * This constructor will setup the bound-bound matrix element of the multipole
         * potential. The integral to be computed is
         * @f[
         *     V_{ab}^\lambda(r') = \int_0^\infty R_{N_a L_a}(r) \left(
         *         \frac{r_<^\lambda}{r_>^{\lambda+1}} - \frac{1}{r'}
         *     \right) R_{N_b L_b}(r) r^2 \, \mathrm{d}r \ ,
         * @f]
         * where @f$ R_{nl}(r) @f$ are the radial hydrogen functions and @f$ F_l(k,r) @f$
         * the Coulomb partial waves. See @ref Hydrogen::P and @ref Hydrogen::F for definitions
         * of these functions.
         */
        MultipolePotential (int lambda, double Ka, int La, int Nb,    int Lb);
        
        MultipolePotential (int lambda, int Na,    int La, double Kb, int Lb);
        
        /**
         * @brief Free-free constructor.
         * 
         * This function does nothing at the moment, because free-free elements are
         * used only in exchange integrals, which are not yet supported in Hex-PWBA2.
         */
        MultipolePotential (int lambda, double Ka, int La, double Kb, int Lb);
        
        /**
         * @brief Evaluates the potential.
         * 
         * This is the evaluation function (overloaded ()-operator). Depending on the
         * internal flag (bound-bound, etc.) it will use the attributes to evaluate
         * the matrix element. The usage is straightforward:
           @code
                // create a new class instance for bound-bound potential
                MultipolePotential V (lambda, na, la, nb, lb);
                
                // evaluate the potential for r = 1.0
                std::cout << "Vab(1.0) = " << V(1.0) << std::endl;
           @endcode
         */
        double operator() (double x) const;
        
        /**
         * @brief Different types of the potential.
         * 
         * These flags determine which of the internal constants are being used
         * (which potential is stored in the specific class instance).
         */
        typedef enum
        {
            /// @f$ V_{\alpha\beta}^\lambda(r') = \left<P_\alpha(r)\right|V^\lambda(r,r')\left|P_\beta(r)\right> @f$
            bound_bound,
            
            /// @f$ V_{\alpha\beta}^\lambda(r') = \left<F_\alpha(r)\right|V^\lambda(r,r')\left|P_\beta(r)\right> @f$
            free_bound,
            
            /// @f$ V_{\alpha\beta}^\lambda(r') = \left<P_\alpha(r)\right|V^\lambda(r,r')\left|F_\beta(r)\right> @f$
            bound_free,
            
            /// @f$ V_{\alpha\beta}^\lambda(r') = \left<F_\alpha(r)\right|V^\lambda(r,r')\left|F_\beta(r)\right> @f$
            free_free
        }
        Type;
        
        /**
         * @brief Check potential type.
         */
        Type getType() const
        {
            return type_;
        }
        
    private:
        
        // which type
        Type type_;
        
        // multipole moment (transferred angular momentum)
        int lambda_;
        
        // final or initial bound state quantum numbers
        int Na_, La_, Nb_, Lb_;
        
        // final or initial free state wavenumbers
        double Ka_, Kb_;
        
        // bound-bound specific variables
        double c_;
        rArray coefs_;
        double norm_;
        
        // free-bound and bound-free specific variables
        rArray r_, V_;
        gsl_interp * interp_;
        gsl_interp_accel * accel_;
};

#endif
