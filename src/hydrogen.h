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
	 * for which the radial function is equal to \ref eps.
	 * 
	 * \warning A naive hunt & bisection algorithm is used, which will collapse
	 * if any of the roots lies in the vicinity of \f$ r_k = 2^k \f$.
	 */
	//@{
	static double getFarRadius(int n, int l, double eps, int max_steps = DEFAULT_MAXSTEPS);
	static double getFarRadiusS(int n, int l, double lambda, double eps, int max_steps = DEFAULT_MAXSTEPS);
	//@}
	
	/**
	 * Evaluate free state (Coulomb wave function) \f$ F_{kl}(r) \f$.
	 */
	static double evalFreeState(double k, int l, double r);
	
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
	
	HydrogenFunction() : n(0), l(0) {}
	HydrogenFunction(int n, int l) : n(n), l(l) {}
	
	/**
	 * Compute far radius \f$ R \f$. For \f$ r > R \f$ the hydrogen radial function
	 * will always be less than 'eps'.
	 */
	inline double far (double eps, int max_steps = 1000) const { return Hydrogen::getFarRadius(n,l,eps,max_steps); };
	
	/**
	 * Get principal quantum number.
	 */
	inline int getN () const { return n; }
	
	/**
	 * Get orbital quantum number.
	 */
	inline int getL () const { return l; }
	
	/**
	 * Evaluate the function.
	 */
	inline double operator() (double r) const { return Hydrogen::evalBoundState(n,l,r); }

private:
	
	int n, l;
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
	
	int n, l;
	double lambda;
};

#endif
