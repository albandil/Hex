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

#include <cmath>
#include <gsl/gsl_sf.h>

/// factorial
#define fac(n) gsl_sf_fact(n)

/// second power
#define sqr(x) (gsl_sf_pow_int((x),2))

#include "complex.h"
#include "misc.h"

/// Shifting coefficients for Sturmian T-operators
#define ALPHA_PLUS(n,l)  (sqrt(((n)-(l))*((n)+(l)+1)))
#define ALPHA_MINUS(n,l) (sqrt(((n)-(l)-1)*((n)+(l))))

/// Kronecker delta
#define DELTA(a,b)       ((a) == (b) ? 1 : 0)

//
// Hydrogen radial orbital
//

inline double hydro_P(unsigned n, unsigned l, double z)
{
	return gsl_sf_hydrogenicR(n, l, 1, z);
}

/** Hydrogen radial function (radius-multiplied)
 * Evaluate hydrogen radial function (radius-multiplied) for complex argument
 * \param n Principal quantum number.
 * \param l Orbital quantum number.
 * \param z Radius: real or complex argument.
 */
Complex hydro_P(unsigned n, unsigned l, Complex z);

/** Derivative of Pnl.
 * \param n
 * \param l
 * \param z
 */
Complex dhydro_P(unsigned n, unsigned l, Complex z);

//
// Ricatti-Bessel functions.
//

/** Riccati-Bessel function
 * 
 * Evaluate Riccati-Bessel function for complex argument. Function is not suitable for
 * large degrees, it uses the most naïve (and least stable) evaluation method.
 * Starting from the expressions for zeroth and first Riccati-Bessel function
 * \f[
 *      j_0(z) = \sin z, \qquad j_1(z) = \frac{\sin z}{z} - \cos z
 * \f]
 * the function employs the forward(!) recurrence relation
 * \f[
 *      j_{n+1}(z) = \frac{2n+1}{z} j_n(z) - j_{n-1}(z) .
 * \f]
 * 
 * \param n Degree of the Riccati-Bessel function.
 * \param z Complex argument.
 */
LComplex ric_j(int n, LComplex z);

/** Derivative of Riccati-Bessel function
 * 
 * \param n Degree of the function.
 * \param z Complex argument.
 */
LComplex dric_j(int n, LComplex z);

/**
 * \brief Ricatti-Bessel function of the first kind, \f$ \hat{j}_n(x) \f$.
 * 
 * \f$ \hat{j}_0(x) = \sin x \f$
 */
inline double ric_j (int n, double x)
{
	return x * gsl_sf_bessel_jl(n,x);
}

/**
 * \brief Ricatti-Bessel function of the second kind, \f$ \hat{y}_n(x) \f$.
 * 
 * \f$ \hat{n}_0(x) = \cos x \f$
 */
inline double ric_n (int n, double x)
{
	return x * gsl_sf_bessel_yl(n,x);
}

//
// Modified Bessel function of the first kind.
//

/**
 * \brief Scaled modified spherical Bessel function of the first kind, \f$ \mathrm{e}^{-x} i_n(x) \f$.
 * 
 * \f$ \mathrm{e}^{-x} i_0(x) = \mathrm{e}^{-x} \frac{1}{x} \sinh x \f$
 */
inline double sph_i_scaled (int n, double x)
{
	return gsl_sf_bessel_il_scaled(n,x);
}

/**
 * \brief Scaled modified Ricatti-Bessel function of the first kind, \f$ \mathrm{e}^{-x} \hat{i}_n(x) \f$.
 * 
 * \f$ \mathrm{e}^{-x} \hat{i}_0(x) = \mathrm{e}^{-x} \sinh x \f$
 */
inline double ric_i_scaled (int n, double x)
{
	return x * gsl_sf_bessel_il_scaled(n,x);
}

/**
 * \brief Modified spherical Bessel function of the first kind, \f$ i_n(x) \f$.
 * 
 * \f$ i_0(x) = \frac{1}{x} \sinh x \f$
 */
inline double sph_i (int n, double x)
{
	return exp(x) * sph_i_scaled(n,x);
}

/**
 * \brief Modified Ricatti-Bessel function of the first kind, \f$ \hat{i}_n(x) \f$.
 * 
 * \f$ \hat{i}_0(x) = \sinh x \f$
 */
inline double ric_i (int n, double x)
{
	return exp(x)*ric_i_scaled(n,x);
}


//
// Modified Bessel function of the second kind.
//

/**
 * \brief Scaled modified spherical Bessel function of the second kind, \f$ \mathrm{e}^x k_n(x) \f$.
 * 
 * \f$ \mathrm{e}^{x} k_0(x) = \frac{1}{x} \f$.
 */
inline double sph_k_scaled (int n, double x)
{
	return gsl_sf_bessel_kl_scaled(n,x) * M_2_PI;
}

/**
 * \brief Scaled modified Ricatti-Bessel function of the second kind, \f$ \mathrm{e}^x \hat{k}_n(x) \f$.
 * 
 * \f$ \mathrm{e}^{x} \hat{k}_0(x) = 1 \f$.
 */
inline double ric_k_scaled (int n, double x)
{
	return x * sph_k_scaled(n,x);
}

/**
 * \brief Modified spherical Bessel function of the second kind, \f$ k_n(x) \f$.
 * 
 * \f$ k_0(x) = \frac{1}{x} \mathrm{e}^{-x} \f$.
 */
inline double sph_k (int n, double x)
{
	return exp(-x) * sph_k_scaled(n,x);
}

/**
 * \brief Modified Ricatti-Bessel function of the second kind, \f$ \hat{k}_n(x) \f$.
 * 
 * \f$ \hat{k}_0(x) = \mathrm{e}^{-x} \f$.
 */
inline double ric_k (int n, double x)
{
	return exp(-x) * ric_k_scaled(n,x);
}

//
// Ricatti-Hankel functions.
//

/**
 * \brief Ricatti-Hankel function of the first kind, \f$ \hat{h}_n^{(+)}(x) \f$.
 * 
 * \f$ \hat{h}_0^{(+)}(x) = -\mathrm{i} \exp (\mathrm{i}x) \f$.
 */
inline Complex ric_h_plus (int n, double x)
{
	return Complex(ric_j(n,x), ric_n(n,x));
}

//
// Other special functions.
//

/**
 * \brief Spherical harmonic function.
 */
Complex sphY(int l, int m, double theta, double phi);

/**
 * \brief Base class for real functions.
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
 * \brief Derivative of the modified Ricatti-Bessel function of the second kind, \f$ \hat{i}_n'(x) \f$.
 */
inline double dric_i(int n, double x)
{
	if (n == 0)
		return cosh(x);
	
	return exp(x) * (sph_i_scaled(n,x) + (n * ric_i_scaled(n-1,x) + (n+1) * ric_i_scaled(n+1,x)) / (2*n+1));
}

/**
 * \brief Derivative of the modified Ricatti-Bessel function of the second kind, \f$ \hat{k}_n'(x) \f$.
 */
inline double dric_k(int n, double x)
{
	if (n == 0)
		return -exp(-x);
	
	return exp(-x) * (sph_k_scaled(n,x) - (n * ric_k_scaled(n-1,x) + (n+1) * ric_k_scaled(n+1,x)) / (2*n+1));
}

/**
 * \brief Derivative of scaled modified Ricatti-Bessel function of the second kind, \f$ (\mathrm{e}^{-x} \hat{i}_n(x))' \f$.
 */
inline double dric_i_scaled(int n, double x)
{
	if (n == 0)
		return exp(-2*x);
	
	return sph_i_scaled(n,x) - ric_i_scaled(n,x) + (n * ric_i_scaled(n-1,x) + (n+1) * ric_i_scaled(n+1,x)) / (2*n+1);
}

/**
 * \brief Scaled derivative of modified Ricatti-Bessel function of the second kind, \f$ (\mathrm{e}^x \hat{k}_n(x)) \f$.
 */
inline double dric_k_scaled(int n, double x)
{
	if (n == 0)
		return 0;
	
	return sph_k_scaled(n,x) + ric_k_scaled(n,x) - (n * ric_k_scaled(n-1,x) + (n+1) * ric_k_scaled(n+1,x)) / (2*n+1);
}

/**
 * \brief Derivative of Ricatti-Hankel function of the first kind, \f$ {\hat{h}_n^{(+)}}'(x) \f$.
 */
inline Complex dric_h_plus (int n, double x)
{
	return Complex(dric_j(n,x),dric_n(n,x));
}

/**
 * \brief Uniform approximation to the Coulomb wave function.
 * 
 * This routine uses algorithm from the following article:
 *    Michel N, Uniform WKB approximation of Coulomb wave functions for arbitrary partial wave, EPL, 83 (2008) 10002.
 * 
 * The method is asymptotically valid for high energies and partial waves.
 * 
 * \param l Angular momentum.
 * \param k Wavenumber.
 * \param r Radial coordinate.
 * \param F Output reference for resulting value.
 * \param Fp Output reference for resulting derivative.
 */
int coul_F_michel(int l, double k, double r, double& F, double& Fp);

/**
 * \brief Evaluate Coulomb wave function (and its derivative).
 * \param l Angular momentum.
 * \param k Wavenumber.
 * \param r Radial coordinate.
 * \param F Output reference for resulting value.
 * \param Fp Output reference for resulting derivative.
 */
int coul_F (int l, double k, double r, double & F, double & Fp);

/**
 * \brief Asymptotic form of the regular Coulomb wave.
 * \param l Angular momentum.
 * \param k Wavenumber.
 * \param r Radial coordinate.
 * \param sigma Optionally, the precomputed Coulomb phase shift.
 */
double coul_F_asy(int l, double k, double r, double sigma = Nan);

/**
 * \brief Coulomb phase shift.
 * \param l Angular momentum.
 * \param k Wavenumber.
 */
double coul_F_sigma(int l, double k);

#endif
