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

#ifndef HEX_SPECF
#define HEX_SPECF

#include <gsl/gsl_sf.h>
	/// second power
	#define sqr(x) (gsl_sf_pow_int((x),2))

#include "complex.h"
#include "arrays.h"

/// Shifting coefficients for Sturmian T-operators
#define ALPHA_PLUS(n,l)  (sqrt(((n)-(l))*((n)+(l)+1)))
#define ALPHA_MINUS(n,l) (sqrt(((n)-(l)-1)*((n)+(l))))

/// Kronecker delta
#define DELTA(a,b)       ((a) == (b) ? 1 : 0)

/**
 * \brief Ricatti-Bessel function of the first kind, \f$ \hat{j}_n(x) \f$.
 */
#define ric_j(n,x)			((x)*gsl_sf_bessel_jl(n,x))
/**
 * \brief Ricatti-Bessel function of the second kind, \f$ \hat{y}_n(x) \f$.
 */
#define ric_n(n,x)			((x)*gsl_sf_bessel_yl(n,x))
/**
 * \brief Scaled modified Ricatti-Bessel function of the first kind, \f$ \mathrm{e}^{-x} \hat{i}_n(x) \f$.
 */
#define ric_i_scaled(n,x)	((x)*gsl_sf_bessel_il_scaled(n,x))
/**
 * \brief Scaled modified Ricatti-Bessel function of the second kind, \f$ \mathrm{e}^x \hat{k}_n(x) \f$.
 */
#define ric_k_scaled(n,x)	((x)*gsl_sf_bessel_kl_scaled(n,x))
/**
 * \brief Modified Ricatti-Bessel function of the first kind, \f$ \hat{i}_n(x) \f$.
 */
#define ric_i(n,x)			(exp(fabs(x))* ric_i_scaled(n,x))
/**
 * \brief Modified Ricatti-Bessel function of the second kind, \f$ \hat{k}_n(x) \f$.
 */
#define ric_k(n,x)			(exp(-x)     * ric_k_scaled(n,x))
/**
 * \brief Ricatti-Hankel function of the first kind, \f$ \hat{h}_n^{(+)}(x) \f$.
 */
#define ric_h_plus(n,x) 	Complex( ric_n(n,x), ric_j(n,x))
/**
 * \brief Derivative of Ricatti-Hankel function of the first kind, \f$ {\hat{h}_n^{(+)}}'(x) \f$.
 */
#define dric_h_plus(n,x)	Complex(dric_n(n,x),dric_j(n,x))

/**
 * \brief Spherical harmonic function.
 */
Complex sphY(int l, int m, double theta, double phi);

/**
 * \brief Base class for radial distorted functions.
 * 
 * This class serves as a heritage base for
 * \ref DistortingPotential,
 * \ref DistortedWave,
 * \ref HydrogenFunction
 * and others.
 */
template <typename T> class RadialFunction
{
public:
	
	/// Evaluate the function.
	virtual T operator() (double x) const = 0;
};

/**
 * \brief Derivative of Ricatti-Bessel function of the first kind, \f$ \hat{j}_n'(x) \f$.
 */
inline double dric_j(int n, double x)
{
	if (n == 0)
		return cos(x);
	
	return gsl_sf_bessel_jl(n,x) + (n * ric_j(n-1,x) - (n+1) * ric_j(n+1,x)) / (2*n+1);
}

/**
 * \brief Derivative of Ricatti-Bessel function of the second kind, \f$ \hat{y}_n'(x) \f$.
 */
inline double dric_n(int n, double x)
{
	if (n == 0)
		return sin(x);
	
	return gsl_sf_bessel_yl(n,x) + (n * ric_n(n-1,x) - (n+1) * ric_n(n+1,x)) / (2*n+1);
}

/**
 * \brief Scaled derivative of modified Ricatti-Bessel function of the second kind, \f$ \mathrm{e}^{-x} \hat{i}_n'(x) \f$.
 */
inline double dric_i_scaled(int n, double x)
{
	if (n == 0)
		return 0.5 * (1 + exp(-2*x)); // = exp(-x) cosh(x)
	
	return gsl_sf_bessel_il_scaled(n,x) + (n * ric_i_scaled(n-1,x) - (n+1) * ric_i_scaled(n+1,x)) / (2*n+1);
}

/**
 * \brief Scaled derivative of modified Ricatti-Bessel function of the second kind, \f$ \hat{i}_n'(x) \f$.
 */
inline double dric_i(int n, double x)
{
	if (n == 0)
		return cosh(x);
	
	return exp(fabs(x)) * (gsl_sf_bessel_il_scaled(n,x) + (n * ric_i_scaled(n-1,x) - (n+1) * ric_i_scaled(n+1,x)) / (2*n+1));
}

/**
 * \brief Scaled derivative of modified Ricatti-Bessel function of the second kind, \f$ \mathrm{e}^x \hat{k}_n'(x) \f$.
 */
inline double dric_k_scaled(int n, double x)
{
	if (n == 0)
		return 1/x;
	
	return gsl_sf_bessel_kl_scaled(n,x) + (n * ric_k_scaled(n-1,x) - (n+1) * ric_k_scaled(n+1,x)) / (2*n+1);
}

/**
 * \brief Derivative of modified Ricatti-Bessel function of the second kind, \f$ \hat{k}_n'(x) \f$.
 */
inline double dric_k(int n, double x)
{
	if (n == 0)
		return exp(-x)/x;
	
	return exp(-x) * (gsl_sf_bessel_kl_scaled(n,x) + (n * ric_k_scaled(n-1,x) - (n+1) * ric_k_scaled(n+1,x)) / (2*n+1));
}

/**
 * \brief Asymptotic form of the regular Coulomb wave.
 * \param l Angular momentum.
 * \param k Wavenumber.
 * \param r Radial coordinate.
 * \param sigma Optionally, the precomputed Coulomb phase shift.
 */
double F_asy(int l, double k, double r, double sigma = Nan);

/**
 * \brief Coulomb phase shift.
 * \param l Angular momentum.
 * \param k Wavenumber.
 */
double F_sigma(int l, double k);

/**
 * \brief Complex upper incomplete Gamma function \f$ \Gamma(a,z) \f$.
 * 
 * Computes the analytically extended upper incomplete Gamma function
 * \f[
 *     \Gamma(a,x) = \int_z^\infty t^{a-1} \mathrm{e}^{-t} \mathrm{d}t \ ,
 * \f]
 * \f[
 *     \Gamma(a,x) \ [x \in \mathbb{R}] \rightarrow \Gamma(a,z) \ [z \in \mathbb{C}] \ .
 * \f]
 * 
 */
Complex cgamma (double a, Complex z, int max_iter = 1000, int* iter = 0, double eps = 1e-10);

/**
 * \brief Complex upper incomplete Gamma function (continued fraction).
 * 
 * Uses the method of continued fraction to compute the function value.
 * The formula is
 * \f[
 *     \Gamma(a,x) = \mathrm{e}^{-x} x^a \left( \frac{1}{x+1-a-{}}
 *         \frac{1\cdot(1-a)}{x+3-a-{}} \frac{2\cdot(2-a)}{x+5-a-{}} \dots \right) \ ,
 * \f]
 * which converges well for \f$ x > a + 1 \f$ .
 */
Complex cgamma_cfrac (double a, Complex z, int max_iter = 1000, int* iter = 0, double eps = 1e-10);

/**
 * \brief Complex upper incomplete Gamma function (infinite series).
 * 
 * Uses the method of infinite sries to compute the function value.
 * The formula is
 * \f[
 *     \gamma(a,x) = \Gamma(a) - \Gamma(a,x) = \mathrm{e}^{-x} x^a
 *         \sum_{n=0}^\infty \frac{\Gamma(a)}{\Gamma(a+1+n)} x^n \ ,
 * \f]
 * which converges well for \f$ x < a + 1 \f$.
 */
Complex cgamma_series (double a, Complex z, int max_iter = 1000, int* iter = 0, double eps = 1e-10);

//@}

#endif
