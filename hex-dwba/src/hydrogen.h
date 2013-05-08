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

#ifndef HEX_HYDROGEN
#define HEX_HYDROGEN

#include <vector>

#include "arrays.h"
#include "coulomb.h"
#include "specf.h"

#define DEFAULT_LAMBDA		1
#define DEFAULT_MAXSTEPS	1000

/**
 * Namespace holding routines concerned with the hydrogen atom.
 */
namespace Hydrogen
{
	/**
	 * Bound state \f$ P_{nl}(r) = rR_{nl}(r) \f$.
	 */
	double evalBoundState(int n, int l, double r);
	double lastZeroBound(int n, int l);
	double getBoundN(int n, int l);
	
	/**
	 * Sturmian wave function
	 * \f[
	 * S_{n\ell}(r) = \left(\frac{\lambda_\ell (k-1)!}{(2\ell+1+k)!}\right)^{1/2} 
	 * (\lambda_\ell r)^{\ell+1} \exp(-\lambda_\ell r/2) L_{k-1}^{2\ell+2}(\lambda_\ell r)
	 * \f]
	 */
	double evalSturmian(int n, int l, double r, double lambda = DEFAULT_LAMBDA);
	
	/**
	 * Return radial distance in the exponential decreasing regions,
	 * for which the radial function is equal to "eps".
	 * 
	 * \warning A naive hunt & bisection algorithm is used, which will collapse
	 * if any of the roots lies in the vicinity of \f$ r_k = 2^k \f$.
	 */
	double getBoundFar(int n, int l, double eps, int max_steps = DEFAULT_MAXSTEPS);
	
	/**
	 * Return radial distance in the exponential decreasing regions,
	 * for which the radial function is equal to "eps".
	 * 
	 * \warning A naive hunt & bisection algorithm is used, which will collapse
	 * if any of the roots lies in the vicinity of \f$ r_k = 2^k \f$.
	 */
	double getSturmFar(int n, int l, double lambda, double eps, int max_steps = DEFAULT_MAXSTEPS);
	
	/**
	 * Evaluate free state
	 * \f[
	 *     \psi_{\mathbf{k}lm}(\mathbf{r}) = \frac{1}{k} \sqrt{\frac{2}{\pi}} F_l(k,r) 
	 * 			Y_{lm} (\mathbf{\hat{r}})
	 * 			Y_{lm}^\ast(\mathbf{\hat{k}})
	 * \f]
	 * Use precomputed value of the Coulomb phase shift "sigma" if available.
	 * Complete wave function is further modified by a complex unit factor,
	 * \f[
	 *     \Psi_{\mathbf{k}lm}(\mathbf{r}) = \mathrm{i}^l 
	 * 			\mathrm{e}^{\mathrm{i}\sigma_l(k)}
	 * 			\psi_{\mathbf{k}lm}(\mathbf{r}) \ .
	 * \f]
	 * The missing phase (as a real number) can be retrieved by \ref evalFreeStatePhase .
	 */
	double evalFreeState(double k, int l, double r, double sigma = Nan);
	
	/**
	 * Evaluate phase of the free state function.
	 * Use precomputed value of the Coulomb phase shift "sigma" if available.
	 */
	double evalFreeStatePhase(double k, int l, double sigma = Nan);
	
	/**
	 * Evaluate free state asymptotics \f$ \sin (kr - \pi l / 2 + \sigma_l) \f$.
	 */
	double evalFreeState_asy(double k, int l, double r, double sigma);
	
	/**
	 * Find zeros of the free function asymptotics.
	 */
	double getFreeAsyZero(double k, int l, double Sigma, double eps, int max_steps, int nzero);
	
	/**
	 * Find local maxima of the free function asymptotics.
	 */
	double getFreeAsyTop(double k, int l, double Sigma, double eps, int max_steps, int ntop);
	
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
	double getFreeFar(double k, int l, double Sigma = Nan, double eps = 1e-10, int max_steps = DEFAULT_MAXSTEPS);
};

/**
 * Hydrogen radial function.
 */
class HydrogenFunction : public RadialFunction<double>, public Coulomb_wave_functions
{
public:
	
	/// Constructor for bound state
	HydrogenFunction(int n, int l)
		: Coulomb_wave_functions(), n(n), k(0), l(l), Sigma(0), Verbose(false), Far(far()) {}
	
	/// Constructor for free state
	HydrogenFunction(double k, int l)
		: Coulomb_wave_functions(true,l,-1./k), n(0), k(k), l(l), Sigma(F_sigma(l,k)), Verbose(false), Far(far()) {}
	
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
		if (n != 0)
			return Hydrogen::getBoundFar(n,l,eps,max_steps);
		else
			return Hydrogen::getFreeFar(k,l,Sigma,eps,max_steps);
	};
	
	/// Get principal quantum number.
	inline int getN () const { return n; }
	
	/// Get orbital quantum number.
	inline int getL () const { return l; }
	
	/// Get Coulomb wave momentum.
	inline double getK () const { return k; }
	
	/// Get precomputed far radius.
	inline double getFar() const { return Far; }
	
	/// Verbosity control
	bool verbose() const { return Verbose; }
	void verbose(bool b) { Verbose = b; }
	
	/// Evaluate the function.
	double operator() (double r) const;
	
	/// Evaluate the function for complex radius.
	/// \note Only purely imaginary radius \f$ r \f$ expected for the free function.
	Complex operator() (Complex r) const;
	
	/// Comparison
	inline bool operator== (HydrogenFunction const & psi) const
	{
		return n == psi.n and k == psi.k and l == psi.l;
	}

	/// Complex regular Coulomb function for the parameters "l" and "k".
	inline Complex F(Complex r) const { Complex f, df; F_dF(k*r,f,df); return f; }
	
	/// Complex outgoing Coulomb function for the parameters "l" and "k".
	inline Complex Hplus(Complex r) const { Complex h, dh; H_dH(1,k*r,h,dh); return h; }
	
private:
	
	/// Principal quantum number of bound state.
	int n;
	
	/// Wavenumber of free state.
	double k;
	
	/// Angular momentum.
	int l;
	
	/// Coulomb phase shift.
	double Sigma;
	
	/// Verbosity
	bool Verbose;
	
	/// Far value.
	double Far;
};

/**
 * Hydrogen radial function.
 */
class SturmianFunction : public RadialFunction<double>
{
public:
	
	SturmianFunction() : n(0), l(0), lambda(0) {}
	SturmianFunction(int n, int l, double lambda = DEFAULT_LAMBDA)
	    : n(n), l(l), lambda(lambda) {}
	
	/**
	 * Compute far radius \f$ R \f$. For \f$ r > R \f$ the hydrogen Sturmian function
	 * will always be less than 'eps'.
	 */
	inline double far (double eps, int max_steps = 1000) const
	{
		return Hydrogen::getSturmFar(n,l,lambda,eps,max_steps);
	};
	
	/**
	 * Get principal quantum number.
	 */
	inline int getN () const { return n; }
	
	/**
	 * Get orbital quantum number.
	 */
	inline int getL () const { return l; }
	
	/**
	 * Get shielding constant.
	 */
	inline double getLambda () const { return lambda; }
	
	/**
	 * Evaluate the function.
	 */
	inline double operator() (double r) const
	{
		return Hydrogen::evalSturmian(n, l, r, lambda);
	}
	
private:
	
	/// \f$ T_3 \f$ operator eigenvalue.
	int n;
	
	/// Angular momentum.
	int l;
	
	/// Scaling (shielding) factor.
	double lambda;
};

#endif
