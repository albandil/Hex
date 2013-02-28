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
#include "specf.h"

#define DEFAULT_LAMBDA		1
#define DEFAULT_MAXSTEPS	1000

/**
 * Utility class holding information about the hydrogen atom.
 */
class Hydrogen
{
public:
	
	/**
	 * Bound state \f$ P_{nl}(r) = rR_{nl}(r) \f$.
	 */
	static double evalBoundState(int n, int l, double r);
	static double lastZeroBound(int n, int l);
	
	/**
	 * Sturmian wave function
	 * \f[
	 * S_{n\ell}(r) = \left(\frac{\lambda_\ell (k-1)!}{(2\ell+1+k)!}\right)^{1/2} 
	 * (\lambda_\ell r)^{\ell+1} \exp(-\lambda_\ell r/2) L_{k-1}^{2\ell+2}(\lambda_\ell r)
	 * \f]
	 */
	static double evalSturmian(int n, int l, double r, double lambda = DEFAULT_LAMBDA);
	
	/**
	 * Return radial distance in the exponential decreasing regions,
	 * for which the radial function is equal to "eps".
	 * 
	 * \warning A naive hunt & bisection algorithm is used, which will collapse
	 * if any of the roots lies in the vicinity of \f$ r_k = 2^k \f$.
	 */
	//@{
	static double getFarRadius(int n, int l, double eps, int max_steps = DEFAULT_MAXSTEPS);
	static double getFarRadiusS(int n, int l, double lambda, double eps, int max_steps = DEFAULT_MAXSTEPS);
	//@}
	
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
	static double evalFreeState(double k, int l, double r, double sigma = Nan);
	
	/**
	 * Evaluate phase of the free state function.
	 * Use precomputed value of the Coulomb phase shift "sigma" if available.
	 */
	static double evalFreeStatePhase(double k, int l, double sigma = Nan);
	
	/**
	 * Evaluate free state asymptotics \f$ \sin (kr - \pi l / 2 + \sigma_l) \f$.
	 */
	static double evalFreeState_asy(double k, int l, double r);
};

/**
 * Hydrogen radial function.
 */
class HydrogenFunction : public RadialFunction<double>
{
public:
	
	HydrogenFunction() : n(0), k(0), l(0), Sigma(0) {}
	HydrogenFunction(int n, int l) : n(n), k(0), l(l), Sigma(0) {}
	HydrogenFunction(double k, int l) : n(0), k(k), l(l), Sigma(F_sigma(l,k)) {}
	
	/**
	 * \brief Get far radius.
	 * 
	 * Compute far radius \f$ R \f$. For \f$ r > R \f$ the hydrogen radial function
	 * will always be less than 'eps'.
	 */
	inline double far (double eps, int max_steps = 1000) const
	{
		if (n != 0)
		{
			return Hydrogen::getFarRadius(n,l,eps,max_steps);
		}
		else
		{
			throw "Coulomb function has no \"far\"!";
		}
	};
	
	/// Get principal quantum number.
	inline int getN () const { return n; }
	
	/// Get orbital quantum number.
	inline int getL () const { return l; }
	
	/// Get Coulomb wave momentum.
	inline double getK () const { return k; }
	
	/// Evaluate the function.
	inline double operator() (double r) const
	{
		if (n != 0)
			return Hydrogen::evalBoundState(n, l, r);
		else
			return Hydrogen::evalFreeState(k, l, r, Sigma == 0 ? std::numeric_limits<double>::quiet_NaN() : Sigma);
	}
	
	/// Comparison
	inline bool operator== (HydrogenFunction const & psi) const
	{
		return n == psi.n and k == psi.k and l == psi.l;
	}

private:
	
	/// Principal quantum number of bound state.
	int n;
	
	/// Wavenumber of free state.
	double k;
	
	/// Angular momentum.
	int l;
	
	/// Coulomb phase shift.
	double Sigma;
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
	inline double far (double eps, int max_steps = 1000) const { return Hydrogen::getFarRadiusS(n,l,lambda,eps,max_steps); };
	
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
	inline double operator() (double r) const { return Hydrogen::evalSturmian(n, l, r, lambda); }
	
private:
	
	/// \f$ T_3 \f$ operator eigenvalue.
	int n;
	
	/// Angular momentum.
	int l;
	
	/// Scaling (shielding) factor.
	double lambda;
};

#endif
