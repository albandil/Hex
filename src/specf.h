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

#ifndef HEX_SPECF
#define HEX_SPECF

#include <cmath>
#include <gsl/gsl_sf.h>

/// factorial
#define fac(n) gsl_sf_fact(n)

/// second power
#define sqr(x) (gsl_sf_pow_int((x),2))

#include "arrays.h"
#include "complex.h"
#include "misc.h"

/// Shifting coefficients for Sturmian T-operators
#define ALPHA_PLUS(n,l)  (sqrt(((n)-(l))*((n)+(l)+1)))
#define ALPHA_MINUS(n,l) (sqrt(((n)-(l)-1)*((n)+(l))))

/// Kronecker delta
#define DELTA(a,b)       ((a) == (b) ? 1 : 0)

namespace special
{

namespace constant
{
    
const double pi_half = 1.57079632679489661923;
const double sqrt_half = 0.70710678118654752440;
    
}; // end of namespace "special::constant"

namespace integral
{

/**
 * @brief Compute integral of the confluent hypergeometric function.
 * 
 * The function template will return the scaled value of the indefinite integral
 * @f[
 *     I(a,b,c,u,v;z) = \frac{1}{x^{a+1}} \int x^a \mathrm{e}^{bx}
 *     {}_1\!F_1\left(\matrix{u \cr v}\Big|\,cx\right)
 *     \mathrm{d}x \ .
 * @f]
 * The series representation, which is used here, is obtained by
 * a term-by-term integration. The result is
 * @f[
 *     I(a,b,c,u,v;z) = \sum_{k=0}^\infty \frac{1}{a+1+k}
 *     \sum_{n=0}^k \frac{(u)_n}{(v)_n} \frac{(bx)^{k-n}}{(k-n)!}
 *     \frac{(cx)^n}{n!} \ .
 * @f]
 * This is an analytic function in all arguments except for @f$ a \in \mathbb{Z}^- @f$.
 * 
 * The function also accepts optional arguments that govern the precision
 * and the maximal number of iterations (terms of the outer sum).
 * Due to the template character the arguments "a", "b", "c", "u", "v" and "z"
 * are expected to be either all real or all complex.
 * The result has the corresponding type as well.
 */
template <class T> T pow_exp_hyperg1F1 (T a, T b, T c, T u, T v, T x, double epsrel = 1e-10, unsigned maxiter = 1000)
{
    // last added term and up-to-now sum of terms
    T term = 1. / (a + 1.);
    T bfactor = 1.;
    T suma = term;
    
    // auxiliary variables
    T c_over_b = c / b;
    
    // the outer sum
    for (unsigned n = 1; std::abs(term) > epsrel * std::abs(suma); n++)
    {
        // compute the inner sum
        T iterm = 1, isum = 0;
        for (unsigned k = 0; k <= n; k++)
        {
            // update the inner sum
            isum += iterm;
            
            // update term for the next k
            iterm *= (u + T(k)) / (v + T(k)) * ((n - k) / (k + 1.)) * c_over_b;
        }
        
        // update the term and the outer sum
        bfactor *= b * x / T(n);
        term = bfactor * isum / (a + T(n + 1));
        suma += term;
        
        std::cout << "n = " << n << std::endl;
        std::cout << "\tisum = " << isum << std::endl;
        std::cout << "\tbfac = " << bfactor << std::endl;
        std::cout << "\tterm = " << term << std::endl;
        std::cout << "\tsuma = " << suma << std::endl;
        
        // check if we run out of allowed iterations
        if (n == maxiter)
        {
            throw exception
            (
                "Maximal number of iterations (%d) reached in pow_exp_hyperg1F1.",
                maxiter
            );
        }
    }
    
    // return the result
    return suma;
}
// template <class T> T pow_exp_hyperg1F1 (T a, T b, T c, T u, T v, T x, double epsrel = 1e-10, unsigned maxiter = 1000)
// {
//     // last added term and up-to-now sum of terms
//     T term = 1. / (a + 1.);
//     T suma = term;
//     
//     // arrays containing (bx)^k/k! and (cx)^k/k!
//     NumberArray<T> bterms = { T(1) };
//     NumberArray<T> cterms = { T(1) };
//     
//     // the outer sum
//     for (unsigned k = 1; std::abs(term) > epsrel * std::abs(suma); k++)
//     {
//         // add new term of Taylor expansion of exponential
//         T x_over_k = x / T(k);
//         bterms.push_back (bterms.back() * b * x_over_k);
//         cterms.push_back (cterms.back() * c * x_over_k);
//         
//         // compute the inner sum
//         T poch = 1, isum = 0;
//         for (unsigned n = 0; n <= k; n++)
//         {
//             // update the inner sum
//             isum += poch * bterms[k-n] * cterms[n];
//             
//             // update Pochhammer symbols for next n
//             poch *= (u + T(n)) / (v + T(n));
//         }
//         
//         // update the term and the outer sum
//         term = isum / (a + T(k + 1));
//         suma += term;
//         
//         // check if we run out of allowed iterations
//         if (k == maxiter)
//         {
//             throw exception
//             (
//                 "Maximal number of iterations (%d) reached in pow_exp_hyperg1F1.",
//                 maxiter
//             );
//         }
//     }
//     
//     // return the result
//     return suma;
// }

};

}; // end of namespace "special"

//
// Hydrogen radial orbital
//

inline double hydro_P (unsigned n, unsigned l, double z)
{
    return z * gsl_sf_hydrogenicR(n, l, 1, z);
}

inline double hydro_F (double k, unsigned l, double z)
{
    gsl_sf_result F, G, Fp, Gp;
    double expF, expG;
    
    gsl_sf_coulomb_wave_FG_e (-1./k, k*z, l, 0, &F, &Fp, &G, &Gp, &expF, &expG);
    
    return special::constant::pi_half * special::constant::sqrt_half * F.val;
}

/** Hydrogen radial function (radius-multiplied)
 * Evaluate hydrogen radial function (radius-multiplied) for complex argument
 * @param n Principal quantum number.
 * @param l Orbital quantum number.
 * @param z Radius: real or complex argument.
 */
Complex hydro_P(unsigned n, unsigned l, Complex z);

/** Derivative of Pnl.
 * @param n
 * @param l
 * @param z
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
 * @f[
 *      j_0(z) = \sin z, \qquad j_1(z) = \frac{\sin z}{z} - \cos z
 * @f]
 * the function employs the forward(!) recurrence relation
 * @f[
 *      j_{n+1}(z) = \frac{2n+1}{z} j_n(z) - j_{n-1}(z) .
 * @f]
 * 
 * \note Forward recurrence is unstable for small arguments. In the present
 *       implementation those are expected to be real. For real arguments the
 *       GSL routine is called, which doesn't suffer from this instability.
 * 
 * @param n Degree of the Riccati-Bessel function.
 * @param z Complex argument.
 */
Complex ric_j(int n, Complex z);

/**
 * Vectorized interface for \ref ric_j.
 * 
 * Returns values of all Riccati-Bessel functions of order less than or equal to lmax.
 * 
 * @param lmax Angular momentum limit.
 * @param z Complex argument.
 */
cArray ric_jv(int lmax, Complex z);

/** Derivative of Riccati-Bessel function
 * 
 * @param n Degree of the function.
 * @param z Complex argument.
 */
Complex dric_j(int n, Complex z);

/**
 * Vectorized interface for \ref dric_j.
 * 
 * Returns values of all derivatives of Riccati-Bessel functions of order less than or equal to lmax.
 * 
 * @param lmax Angular momentum limit.
 * @param z Complex argument.
 */
cArray dric_jv(int lmax, Complex z);

/**
 * @brief Ricatti-Bessel function of the first kind, @f$ \hat{j}_n(x) @f$.
 * 
 * @f$ \hat{j}_0(x) = \sin x @f$
 */
inline double ric_j (int n, double x)
{
    return x * gsl_sf_bessel_jl(n,x);
}

/**
 * @brief Ricatti-Bessel function of the second kind, @f$ \hat{y}_n(x) @f$.
 * 
 * @f$ \hat{n}_0(x) = \cos x @f$
 */
inline double ric_n (int n, double x)
{
    return x * gsl_sf_bessel_yl(n,x);
}

//
// Modified Bessel function of the first kind.
//

/**
 * @brief Scaled modified spherical Bessel function of the first kind, @f$ \mathrm{e}^{-x} i_n(x) @f$.
 * 
 * @f$ \mathrm{e}^{-x} i_0(x) = \mathrm{e}^{-x} \frac{1}{x} \sinh x @f$
 */
inline double sph_i_scaled (int n, double x)
{
    return gsl_sf_bessel_il_scaled(n,x);
}

/**
 * @brief Scaled modified Ricatti-Bessel function of the first kind, @f$ \mathrm{e}^{-x} \hat{i}_n(x) @f$.
 * 
 * @f$ \mathrm{e}^{-x} \hat{i}_0(x) = \mathrm{e}^{-x} \sinh x @f$
 */
inline double ric_i_scaled (int n, double x)
{
    return x * gsl_sf_bessel_il_scaled(n,x);
}

/**
 * @brief Modified spherical Bessel function of the first kind, @f$ i_n(x) @f$.
 * 
 * @f$ i_0(x) = \frac{1}{x} \sinh x @f$
 */
inline double sph_i (int n, double x)
{
    return exp(x) * sph_i_scaled(n,x);
}

/**
 * @brief Modified Ricatti-Bessel function of the first kind, @f$ \hat{i}_n(x) @f$.
 * 
 * @f$ \hat{i}_0(x) = \sinh x @f$
 */
inline double ric_i (int n, double x)
{
    return exp(x)*ric_i_scaled(n,x);
}


//
// Modified Bessel function of the second kind.
//

/**
 * @brief Scaled modified spherical Bessel function of the second kind, @f$ \mathrm{e}^x k_n(x) @f$.
 * 
 * @f$ \mathrm{e}^{x} k_0(x) = \frac{1}{x} @f$.
 */
inline double sph_k_scaled (int n, double x)
{
    return gsl_sf_bessel_kl_scaled(n,x) * M_2_PI;
}

/**
 * @brief Scaled modified Ricatti-Bessel function of the second kind, @f$ \mathrm{e}^x \hat{k}_n(x) @f$.
 * 
 * @f$ \mathrm{e}^{x} \hat{k}_0(x) = 1 @f$.
 */
inline double ric_k_scaled (int n, double x)
{
    return x * sph_k_scaled(n,x);
}

/**
 * @brief Modified spherical Bessel function of the second kind, @f$ k_n(x) @f$.
 * 
 * @f$ k_0(x) = \frac{1}{x} \mathrm{e}^{-x} @f$.
 */
inline double sph_k (int n, double x)
{
    return exp(-x) * sph_k_scaled(n,x);
}

/**
 * @brief Modified Ricatti-Bessel function of the second kind, @f$ \hat{k}_n(x) @f$.
 * 
 * @f$ \hat{k}_0(x) = \mathrm{e}^{-x} @f$.
 */
inline double ric_k (int n, double x)
{
    return exp(-x) * ric_k_scaled(n,x);
}

//
// Ricatti-Hankel functions.
//

/**
 * @brief Ricatti-Hankel function of the first kind, @f$ \hat{h}_n^{(+)}(x) @f$.
 * 
 * @f$ \hat{h}_0^{(+)}(x) = -\mathrm{i} \exp (\mathrm{i}x) @f$.
 */
inline Complex ric_h_plus (int n, double x)
{
    return Complex(ric_j(n,x), ric_n(n,x));
}

//
// Other special functions.
//

/**
 * @brief Spherical harmonic function.
 */
Complex sphY(int l, int m, double theta, double phi);

/**
 * @brief Bi-polar spherical harmonic function.
 */
Complex sphBiY(int l1, int l2, int L, int M, double theta1, double phi1, double theta2, double phi2);

/**
 * @brief Base class for radial functions.
 */
template <typename T> class RadialFunction
{
public:
    
    /// Evaluate the function.
    virtual T operator() (double x) const = 0;
};

/**
 * @brief Derivative of Ricatti-Bessel function of the first kind, @f$ \hat{j}_n'(x) @f$.
 */
inline double dric_j(int n, double x)
{
    if (n == 0)
        return cos(x);
    
    return gsl_sf_bessel_jl(n,x) + (n * ric_j(n-1,x) - (n+1) * ric_j(n+1,x)) / (2*n+1);
}

/**
 * @brief Derivative of Ricatti-Bessel function of the second kind, @f$ \hat{y}_n'(x) @f$.
 */
inline double dric_n(int n, double x)
{
    if (n == 0)
        return sin(x);
    
    return gsl_sf_bessel_yl(n,x) + (n * ric_n(n-1,x) - (n+1) * ric_n(n+1,x)) / (2*n+1);
}

/**
 * @brief Derivative of the modified Ricatti-Bessel function of the second kind, @f$ \hat{i}_n'(x) @f$.
 */
inline double dric_i(int n, double x)
{
    if (n == 0)
        return cosh(x);
    
    return exp(x) * (sph_i_scaled(n,x) + (n * ric_i_scaled(n-1,x) + (n+1) * ric_i_scaled(n+1,x)) / (2*n+1));
}

/**
 * @brief Derivative of the modified Ricatti-Bessel function of the second kind, @f$ \hat{k}_n'(x) @f$.
 */
inline double dric_k(int n, double x)
{
    if (n == 0)
        return -exp(-x);
    
    return exp(-x) * (sph_k_scaled(n,x) - (n * ric_k_scaled(n-1,x) + (n+1) * ric_k_scaled(n+1,x)) / (2*n+1));
}

/**
 * @brief Derivative of scaled modified Ricatti-Bessel function of the second kind, @f$ (\mathrm{e}^{-x} \hat{i}_n(x))' @f$.
 */
inline double dric_i_scaled(int n, double x)
{
    if (n == 0)
        return exp(-2*x);
    
    return sph_i_scaled(n,x) - ric_i_scaled(n,x) + (n * ric_i_scaled(n-1,x) + (n+1) * ric_i_scaled(n+1,x)) / (2*n+1);
}

/**
 * @brief Scaled derivative of modified Ricatti-Bessel function of the second kind, @f$ (\mathrm{e}^x \hat{k}_n(x)) @f$.
 */
inline double dric_k_scaled(int n, double x)
{
    if (n == 0)
        return 0;
    
    return sph_k_scaled(n,x) + ric_k_scaled(n,x) - (n * ric_k_scaled(n-1,x) + (n+1) * ric_k_scaled(n+1,x)) / (2*n+1);
}

/**
 * @brief Derivative of Ricatti-Hankel function of the first kind, @f$ {\hat{h}_n^{(+)}}'(x) @f$.
 */
inline Complex dric_h_plus (int n, double x)
{
    return Complex(dric_j(n,x),dric_n(n,x));
}

/**
 * @brief Uniform approximation to the Coulomb wave function.
 * 
 * This routine uses algorithm from the following article:
 *    Michel N, Uniform WKB approximation of Coulomb wave functions for arbitrary partial wave, EPL, 83 (2008) 10002.
 * 
 * The method is asymptotically valid for high energies and partial waves.
 * 
 * @param l Angular momentum.
 * @param k Wavenumber.
 * @param r Radial coordinate.
 * @param F Output reference for resulting value.
 * @param Fp Output reference for resulting derivative.
 */
int coul_F_michel(int l, double k, double r, double& F, double& Fp);

/**
 * @brief Evaluate Coulomb wave function (and its derivative).
 * @param l Angular momentum.
 * @param k Wavenumber.
 * @param r Radial coordinate.
 * @param F Output reference for resulting value.
 * @param Fp Output reference for resulting derivative.
 */
int coul_F (int l, double k, double r, double & F, double & Fp);

/**
 * @brief Asymptotic form of the regular Coulomb wave.
 * @param l Angular momentum.
 * @param k Wavenumber.
 * @param r Radial coordinate.
 * @param sigma Optionally, the precomputed Coulomb phase shift.
 */
double coul_F_asy(int l, double k, double r, double sigma = Nan);

/**
 * @brief Coulomb phase shift.
 * @param l Angular momentum.
 * @param k Wavenumber.
 */
double coul_F_sigma(int l, double k);

/**
 * @brief Check triangle inequality.
 * 
 * Verify that the three angular momenta satisfy triangle inequalities, i.e.
 * @f[
 *       |j_1 - j_2 | \le j_3 \le j_1 + j_2
 * @f]
 * and analogously for all (cyclical) permutations of indices.
 */
bool makes_triangle (int two_j1, int two_j2, int two_j3);

/**
 * @brief Auxiliary function used in coupling coefficient computations.
 * 
 * The triangle function, defined as
 * @f[
 *     \Delta(j_1,j_2,j_3) = \sqrt{
 *         \frac{(j_1+j_2-j_3)!(j_1-j_2+j_3)!(-j_1+j_2+j_3)!}{(j_1+j_2+j_3+1)!}
 *     } \ .
 * @f]
 */
double logdelta (int two_j1, int two_j2, int two_j3);

/**
 * @brief Wigner 3j coefficient.
 * 
 * Compute the Wigner 3j coefficient. Even though there is a routine "gsl_sf_coupling_3j"
 * in GSL library, it has to be implemented anew, because of factorial overflows. This
 * routine computes all factorials only in logarithms using the standard function "lgamma".
 * The formula is taken from Edmonds, A. R.: Angular momentum in quantum mechanics, Princeton
 * 1968.
 * @f[
 *     \left( \matrix{j_1 7 j_2 & j_3 \cr m_1 & m_2 & m_3} \right)
 *     = \epsilon(j_1,j_2,j_3)\Delta(j_1,j_2,j_3) \delta_{m_1+m_2+m_3}^0 (-1)^{j_1-j_2-m_3}
 *     \sqrt{(j_1+m_1)! (j_1-m_1)! (j_2+m_2)! (j_2-m_2)! (j_3+m_3)! (j_3-m_3)!}
 *     \sum_k \frac{(-1)^2}{k! (j_1+j_2-j_3-k)! (j_1-m_1-k)! (j_2+m_2-k)! (j_3-j_2+m_1+k)!
 *     (j_3-j_1-m_2+k)!} \ ,
 * @f]
 * @f[
 *     \epsilon(j_1,j_2,j_3) = \cases{1 & triangle inequality satisfied \cr 0 & otherwise} \ .
 * @f]
 * See @ref logdelta for definition of the triangle function @f$ \Delta(a,b,c) @f$.
 * Note that the arguments @f$ j_1, j_2, j_3, m_1, m_2, m_3 @f$ need to be supplied
 * doubled, as @f$ 2j_1, 2j_2, 2j_3, 2m_1, 2m_2, 2m_3 @f$ so that the parameters can
 * be considered integral even though the half-integral angular momentum is allowed.
 */
//@{
double Wigner3j_2 (int two_j1, int two_j2, int two_j3, int two_m1, int two_m2, int two_m3);
#define Wigner3j(a,b,c,d,e,f) Wigner3j_2(2*(a),2*(b),2*(c),2*(d),2*(e),2*(f))
//@}

/**
 * @brief Wigner 6j coefficient.
 * 
 * Compute the Wigner 6j coefficient. Even though there is a routine "gsl_sf_coupling_6j"
 * in GSL library, it has to be implemented anew, because of factorial overflows. This
 * routine computes all factorials only in logarithms using the standard function "lgamma".
 * The formula is taken from Edmonds, A. R.: Angular momentum in quantum mechanics, Princeton
 * 1968.
 * @f[
 *     \left\{ \matrix{j_1 & j_2 & j_3 \cr j_4 & j_5 & j_6} \right\}
 *     =
 *     \Delta(j_1,j_2,j_3) \Delta(j_3,j_4,j_5) \Delta(j_1,j_5,j_6) \Delta(j_2,j_4,j_6)
 *     \sum_k \frac{(-1)^k (k+1)!}{(k-j_1-j_2-j_3)! (k-j_1-j_5-j_6)! (k-j_2-j_4-j_6)!
 *     (k-j_3-j_4-j_5)! (j_1+j_2+j_4+j_5-k)! (j_1+j_3+j_4+j_6-k)! (j_2+j_3+j_5+j_6-k)!} \ .
 * @f]
 * See @ref logdelta for definition of the triangle function @f$ \Delta(a,b,c) @f$.
 * Note that the arguments @f$ j_1, j_2, j_3, j_4, j_5, j_6 @f$ need to be supplied
 * doubled, as @f$ 2j_1, 2j_2, 2j_3, 2j_4, 2j_5, 2j_6 @f$ so that the parameters can
 * be considered integral even though the half-integral angular momentum is allowed.
 */
//@{
double Wigner6j_2 (int two_j1, int two_j2, int two_j3, int two_j4, int two_j5, int two_j6);
#define Wigner6j(a,b,c,d,e,f) Wigner6j_2(2*(a),2*(b),2*(c),2*(d),2*(e),2*(f))
//@}

/**
 * @return Value of @f$ f(\lambda,l_1,l_2,l_1',l_2',L) @f$.
 */
double computef (int lambda, int l1, int l2, int l1p, int l2p, int L);

/**
 * @brief Clebsch-Gordan coefficient.
 * 
 * @note In present implementation valid only for
 * integer (not half-integer) angular momenta.
 * [Otherwise one needs to correct the signs.]
 */
double ClebschGordan(int l1, int m1, int l2, int m2, int L, int M);

/**
 * @brief Gaunt's integral.
 * 
 * Computes the integral of three spherical harmonic functions.
 * @f[
 * \int_{4\pi} Y_{l_1m_1} Y_{l_2m_2} Y^{\ast}_{lm} \mathrm{d}\Omega \ .
 * @f]
 */
double Gaunt(int l1, int m1, int l2, int m2, int l, int m);

/**
 * Compute number of angular momenta pairs @f$ \ell_1 @f$ and @f$ \ell_2 @f$ that
 * are less than or equal to "maxell" and compose the total angular momentum @f$ L @f$.
 */
int triangle_count(int L, int maxell);

/**
 * @brief The hypergeometric function @f$ {}_2F_1 @f$.
 * 
 * Evaluates the Gauss hypergeometric function
 * @f[
 *     {}_2F_1(a,b;c;x)
 * @f]
 * for arbitrary real argument @f$ x @f$. It uses the routine gsl_sf_hyperg_2F1 from
 * GSL library, which is defined for |x| < 1. Using the transformations Abramowitz & Stegun
 * 15.3.3-9 the hypergeometric function is transformed so that is can be evaluated by
 * the library function.
 * 
 * @todo Document transformations.
 */
double Hyper2F1 (double a, double b, double c, double x);

#endif
