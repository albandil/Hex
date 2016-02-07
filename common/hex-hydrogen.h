//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2016, Jakub Benda, Charles University in Prague                    //
//                                                                                   //
// MIT License:                                                                      //
//                                                                                   //
//  Permission is hereby granted, free of charge, to any person obtaining a          //
// copy of this software and associated documentation files (the "Software"),        //
// to deal in the Software without restriction, including without limitation         //
// the rights to use, copy, modify, merge, publish, distribute, sublicense,          //
// and/or sell copies of the Software, and to permit persons to whom the             //
// Software is furnished to do so, subject to the following conditions:              //
//                                                                                   //
//  The above copyright notice and this permission notice shall be included          //
// in all copies or substantial portions of the Software.                            //
//                                                                                   //
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS          //
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       //
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE       //
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, //
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF         //
// OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.  //
//                                                                                   //
//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //

#ifndef HEX_HYDROGEN
#define HEX_HYDROGEN

#include <vector>

#include "hex-special.h"

/// Default range parameter in the Sturmian functions.
#define DEFAULT_LAMBDA      1

// Default maximal iteration count when searching for nodes of the wave functions.
#define DEFAULT_MAXSTEPS    1000

/**
 * @brief Hydrogen atom.
 * 
 * Namespace holding routines concerned with the hydrogen atom. Only a narrow
 * subset of the members is used.
 */
namespace Hydrogen
{

/**
 * @brief Hydrogen bound wave function in Carthesian coordinates.
 */
class CarthesianBoundWaveFunction
{
    public:
        
        CarthesianBoundWaveFunction (int N, int L, int M);
        
        ~CarthesianBoundWaveFunction ();
        
        /**
         * @brief Term of the multidimensional hydrogen bound function.
         * 
         * Term in the form
         * @f[
         *     c x^u y^v z^w r^n
         * @f]
         */
        typedef struct
        {
            Complex c;
            unsigned u, v, w, n;
        }
        Term;
        
        /**
         * @brief Evaluate the function.
         * 
         * Evaluates the hydrogen function at @f$ \mathbf{r} = (x,y,z) @f$.
         */
        Complex operator() (double x, double y, double z) const;
        
        /**
         * @brief Access the overall normalization factor.
         */
        double norm () const { return norm_; }
        
        /**
         * @brief Access the list of terms.
         */
        std::vector<Term> const & terms () const { return terms_; }
        
    private:
        
        int N_, L_, M_;
        double norm_;
        std::vector<Term> terms_;
};

/**
 * @brief Compose name of the hydrogen state.
 * 
 * Return the spectroscopic notation for given quantum numbers.
 * Example are: 1s(0), 2s(0), 2p(0), 2p(1), 2p(-1).
 */
std::string stateName (int n, int l, int m);

/**
 * @brief Compose name of the semi-relativistic hydrogen state.
 * 
 * Return the spectroscopic notation for given quantum numbers.
 * Example are: 1s1/2(0), 2s1/2(0), 2p1/2(0), 2p3/2(0), 2p1/2(-1/2).
 */
std::string stateRName (int n, int l, int two_j, int two_m);

/**
 * @brief Hydrogen bound radial orbital.
 * 
 * The hydrogen radial functions @f$ R_{nl} @f$ are defined as
 * @f[
 *     R_{nl}(r) = \sqrt{\left(\frac{2}{n}\right)^3 \frac{(n-l-1)!}{2n(n+l)!}}
 *     \left(\frac{2r}{n}\right)^l L_{n-l-1}^{(2l+1)}\left(\frac{2r}{n}\right)
 *     \mathrm{e}^{-r/n} \ ,
 * @f]
 * where @f$ L_n^{(\alpha)}(r) @f$ are the generalized Laguerre polynomials.
 * This function, however, returns the "radial orbital", which is understood
 * here to be the radial function multiplied by the radius,
 * @f[
 *     P_{nl}(r) = r R_{nl}(r) \ .
 * @f]
 */
double P (unsigned n, unsigned l, double r, double Z = 1);

/**
 * @brief Get last node of the bound radial wave function.
 */
double lastZeroBound (int n, int l);

/**
 * @brief Return all constant factors of the bound state.
 * 
 * This function doesn't compute only the normalization
 * factor @f$ N_{nl} @f$ but also multiplied the angular-momentum
 * factor @f$ (2/n)^l @f$.
 */
double getBoundN (int n, int l);

/**
 * @brief Sturmian wave function
 * 
 * Evaluate the Sturmian wave function from the formula
 * \f[
 * S_{n\ell}(r) = \left(\frac{\lambda_\ell (k-1)!}{(2\ell+1+k)!}\right)^{1/2} 
 * (\lambda_\ell r)^{\ell+1} \exp(-\lambda_\ell r/2) L_{k-1}^{2\ell+2}(\lambda_\ell r)
 * \f]
 */
double S (int n, int l, double r, double lambda = DEFAULT_LAMBDA);

/**
 * Return radial distance in the exponential decreasing regions,
 * for which the radial function is equal to "eps".
 * 
 * \warning A naive hunt & bisection algorithm is used, which will collapse
 * if any of the roots lies in the vicinity of \f$ r_k = 2^k \f$.
 */
double getBoundFar (int n, int l, double eps, int max_steps = DEFAULT_MAXSTEPS);

/**
 * Return radial distance in the exponential decreasing regions,
 * for which the radial function is equal to "eps".
 * 
 * \warning A naive hunt & bisection algorithm is used, which will collapse
 * if any of the roots lies in the vicinity of \f$ r_k = 2^k \f$.
 */
double getSturmFar (int n, int l, double lambda, double eps, int max_steps = DEFAULT_MAXSTEPS);

/**
 * @brief Free hydrogen state wave function.
 * 
 * Evaluate free state
 * \f[
 *     \psi_{\mathbf{k}lm}(\mathbf{r}) = \frac{1}{k} \sqrt{\frac{2}{\pi}} F_l(k,r) 
 * 			Y_{lm} (\mathbf{\hat{r}})
 * 			Y_{lm}^\ast(\mathbf{\hat{k}}) \,,
 * \f]
 * where the Coulomb wave function is computed by GSL, which uses the formula
 * \f[
 *     F_l(k,r) = \mathrm{e}^{\pi/2k} \frac{|\Gamma(l+1-\mathrm{i}/k)|}{(2l+1)!}
 *     \mathrm{e}^{-\mathrm{i}kr} \,{}_1\!F_1\left(
 *         \left.\matrix{l+1-\mathrm{i}/k \cr 2l+2}\right|
 *         2\mathrm{i}kr
 *     \right) .
 * \f]
 * Use precomputed value of the Coulomb phase shift "sigma" if available.
 * Complete wave function is further modified by a complex unit factor,
 * \f[
 *     \Psi_{\mathbf{k}lm}(\mathbf{r}) = \mathrm{i}^l 
 * 			\mathrm{e}^{\mathrm{i}\sigma_l(k)}
 * 			\psi_{\mathbf{k}lm}(\mathbf{r}) \ .
 * \f]
 * The missing phase (as a real number) can be retrieved by \ref evalFreeStatePhase .
 * 
 * \note Note that the factor \f$ \mathrm{i}^l \f$ is used rather inconsistently
 * at the moment: sometimes as \f$ (-\mathrm{i})^l \f$ instead. This just flips
 * the direction of the momentum vector, which is always integrated over, anyway,
 * so it doesn't matter.
 */
double F (double k, int l, double r, double sigma = special::constant::Nan);

/**
 * Evaluate phase of the free state function.
 * Use precomputed value of the Coulomb phase shift "sigma" if available.
 */
double evalFreeStatePhase (double k, int l, double sigma = special::constant::Nan);

/**
 * Evaluate free state asymptotics \f$ \sin (kr - \pi l / 2 + \sigma_l) \f$.
 */
double evalFreeState_asy (double k, int l, double r, double sigma);

/**
 * Find zeros of the free function asymptotics.
 */
double getFreeAsyZero (double k, int l, double Sigma, double eps, int max_steps, int nzero);

/**
 * Find local maxima of the free function asymptotics.
 */
double getFreeAsyTop (double k, int l, double Sigma, double eps, int max_steps, int ntop);

/**
 * \brief Return sufficiently far radius for using the asymptotic form of the free state.
 * 
 * Return radial distance in the oscillating region,
 * for which the radial function less than "eps" in some zero-node of the
 * asymptotical form. The asymptotic form is
 * \f[
 *      F_\ell(k,r) \propto \sin \left(kr - \frac{\ell\pi}{2} + \frac{1}{k}\log 2k + \sigma_\ell(k)\right) \ ,
 * \f]
 * so the free state will be evaluated in such radii that the following
 * is fulfilled:
 * \f[
 *      n\pi = kr - \frac{\ell\pi}{2} + \frac{1}{k}\log 2k + \sigma_\ell(k) \ .
 * \f]
 */
double getFreeFar
(
    double k, int l,
    double Sigma = special::constant::Nan,
    double eps = 1e-10,
    int max_steps = DEFAULT_MAXSTEPS
);

} // end of namespace hydrogen

/**
 * Hydrogen radial function.
 */
class HydrogenFunction : public special::RadialFunction<double>
{
public:
    
    /// Constructor for bound state
    HydrogenFunction (int n, int l)
        : n_(n), l_(l) {}
    
    /**
     * \brief Get far radius.
     * 
     * Compute far radius \f$ R \f$. For \f$ r > R \f$ the hydrogen radial function
     * will always be less than 'eps'. Or, if the function is a free state
     * wave function, return the smallest radius such that the value of the
     * precise free state is less than "eps" and the value of the asymptotic form
     * is zero. (I.e. the free state will be evaluated at the zeros of the asymptotic
     * form.)
     */
    inline double far (double eps = 1e-10, int max_steps = 1000) const
    {
        return Hydrogen::getBoundFar(n_,l_,eps,max_steps);
    };
    
    /// Get principal quantum number.
    inline int n () const { return n_; }
    
    /// Get orbital quantum number.
    inline int l () const { return l_; }
    
    /// Evaluate the function.
    double operator() (double r) const;
    
    /// Classical turning point.
    double getTurningPoint () const;
    
    /// Comparison
    inline bool operator== (HydrogenFunction const & psi) const
    {
        return n_ == psi.n_ and l_ == psi.l_;
    }
    
private:
    
    /// Principal quantum number of bound state.
    int n_;
    
    /// Angular momentum.
    int l_;
};

#endif
