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

#ifndef HEX_SYMBOLIC
#define HEX_SYMBOLIC

#include <iostream>
#include <vector>

#include <cln/cln.h>

#include "misc.h"

#define GF_NONE   0
#define GF_SIN    1
#define GF_COS    2

namespace symbolic
{


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                             * 
 *                              Rational numbers                               *
 *                                                                             * 
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/// Rational number (CLN)
typedef cln::cl_RA rational;

/**
 * @brief Construct a rational number.
 * 
 * Return a @ref rational object, representing the fraction "num"/"den".
 */
inline rational make_rational (int num, int den) { return rational(num)/cln::cl_I(den); }

/// Convert rational to double.
inline double double_approx (rational const & r) { return cln::double_approx(r); }

/// Rational constant 1/2
extern rational onehalf;


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                             * 
 *                              Symbolic polynomials                           *
 *                                                                             * 
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/**
 * \brief Symbolic term structure
 * 
 * Contains information about one term of the form
 * \f[
 *       k x^a \mathrm{gon}\,(bx) \exp(-cx)
 * \f]
 * The function ("gon") is identified by a unique number and presently
 * can be one of the following:
 * - none (0)
 * - sine (1)
 * - cosine (2)
 */
class term
{
public:
    
    /**
     * @brief Default constructor.
     * 
     * Create a new symbolic term that will evaluate to zero.
     */
    term () : ki(1), kr(0), a(0), gf(GF_NONE), b(0), c(0) {}
    
    /**
     * @brief Construct from fraction.
     * 
     * Create a new symbolic term that will evaluate to a rational number.
     */
    term (rational x) : ki(1), kr(x), a(0), gf(GF_NONE), b(0), c(0) {}
    
    /**
     * @brief Conversion to plain number.
     * 
     * Returns the double representation of the term if all symbolic
     * factors are equal to 1. Otherwise throws an exception.
     */
    double todouble() const;
    
    /**
     * @brief Ordering function.
     * 
     * This function is ised in std::sort for ordering of terms. The ordering
     * is done according to multiple fields. In decreasing order of importance
     * those fields are: power, frequency, scale, goniometric function and
     * finally the irrational factor. The ordering is descending, so largest numbers
     * come first.
     */
    static bool ordering (term const & a, term const & b);
    
    //
    // data fields
    //
    
    /// positive (!) irrational multiplication factor
    double ki;
    
    /// rational multiplication factor
    rational kr;
    
    /// exponent
    int a;
    
    /// goniometric function (none, sin or cos)
    int gf;
    
    /// goniometric function wave number
    double b;
    
    /// exponential factor
    rational c;
    
};

/**
 * @brief Symbolic addition.
 * 
 * This overloaded operator will sum two simillar symbolic terms. Two symbolic
 * terms are considered simillar here if they share the power, exponential argument
 * and the whole goniometric function. In that case only the multiplicative constants
 * need to be added.
 */
term operator + (term const & A, term const & B);

/**
 * @brief Symbolic division.
 * 
 * This will try to compute a factor of two symbolic terms. Because the result has to
 * be a symbolic term defined above, only a subset of possible divisions is allowed:
 * The goniometric functions can either be identical in both terms or the second term
 * has to contain none goniometric function. If this is violated, the function will throw.
 */
term operator / (term const & A, term const & B) throw (exception);

/**
 * @brief Symbolic exponential.
 * 
 * Construct a symbolic term containing just the exponential function.
 * @f[
 *     \exp(-cx)
 * @f]
 */
term expm (rational c);

/**
 * \brief Sum of symbolic terms.
 * 
 * This class holds a representation of a polynomial as an array of symbolic
 * terms (@ref term). It also posesses a STL-like interface for efficient use
 * in range-based loops.
 */
class poly
{
public:
    
    //
    // constructors
    //
    
    /**
     * @brief Default constructor.
     * 
     * Initialize the object by an empty polynomial. It contains
     * no terms and evaluates to zero.
     */
    poly() {}
    
    /**
     * @brief Size constructor.
     * 
     * Initialize the array of terms to a given size. The terms
     * will evaluate to zero unless modified.
     */
    poly(int n) : terms(n) {}
    
    /**
     * @brief Single-term constructor.
     * 
     * This contructor will create a new polynomial containing
     * just a single term, which is copied from the argument.
     */
    poly(term const & T) : terms({T}) {}
    
    /**
     * @brief Copy constructor.
     * 
     * This constructor will initialize the array with existing array
     * of symbolic terms. The resulting polynomial will be identical to
     * the argument.
     */
    poly(poly const & P) : terms(P.terms) {}
    
    //
    // data storage
    //
    
    /**
     * @brief Terms array.
     * 
     * This array contains all symbolic terms that make the polynomial.
     * Every term is represented by a single item in the array.
     */
    std::vector<term> terms;
    
    /**
     * @brief Optimize the polynomial.
     * 
     * Sums similar terms. The terms are simillar if they agree in the power,
     * exponential argument and the whole goniometric function. Only the
     * multiplicative constant is summed. If the result is zero, the term
     * is utterly omitted.
     */
    void optimize();
    
    
    //
    // std::vector-like interface
    //
    
    typedef std::vector<term>::iterator       iterator;
    typedef std::vector<term>::const_iterator const_iterator;
    inline iterator begin() { return terms.begin(); }
    inline iterator end()   { return terms.end();   }
    inline const_iterator begin() const { return terms.begin(); }
    inline const_iterator end()   const { return terms.end();   }
    inline term & front() { return terms.front(); }
    inline term & back()  { return terms.back();  }
    inline term const & front() const { return terms.front(); }
    inline term const & back()  const { return terms.back();  }
    inline void push_back(term const & T) { terms.push_back(T); }
    inline void insert(iterator a, const_iterator u, const_iterator v) { terms.insert(a,u,v); }
    inline size_t size() const { return terms.size(); }
    inline term & operator [] (size_t i) { return terms[i]; }
    inline term const & operator [] (size_t i) const { return terms[i]; }
    
    void clear () { terms.clear(); }
    
    //
    // stream output
    //
    
    /**
     * @brief Formatted output.
     * 
     * Write polynomial to standard output in the natural format
     * <pre>
     * k x^a exp(bx)
     * k x^a exp(bx) sin(cx)
     * k x^a exp(bx) cos(cx)
     * </pre>
     */
    friend std::ostream & operator << (std::ostream & os, poly const & P);
};

/**
 * @brief Evaluate the polynomial.
 * 
 * Evaluate the polynomial for a given value of the variable. The evaluation is
 * done term after term and the results are summed afterwards.
 */
double eval (poly const & P, double r);

/**
 * @brief Collect common term.
 * 
 * Split polynomial into two polynomials that do or do not contain the given factor in the terms.
 * One of the following terms is expected: exp(-cx), sin(bx), cos(bx), sin(bx)exp(-cx),
 * cos(bx)exp(-cx). Other terms than those listed cannot be collected.
 * 
 * The resulting splitting has the form
 * @f[
 *     P(x) = p(x)Q(x) + R(x)
 * @f]
 * @param P The original polynomial.
 * @param p Factor to be collected.
 * @param Q Collected terms with @f$ p(x) @f$ factored out.
 * @param R The uncollected terms.
 */
void collect (poly const & P, term const & p, poly & Q, poly & R);

/**
 * @brief Write to string.
 */
std::string tostring (poly const & P);

/**
 * @brief Write to stream.
 */
std::ostream & operator << (std::ostream & os, poly const & P);

/**
 * @brief Symbolic addition.
 * 
 * Overloaded summation operator that will add two symbolic polynomials.
 * The addition is simple concatenation of the terms followed by optimization
 * of the resulting polynomial (function @ref optimize).
 */
poly operator + (poly const & P, poly const & Q);

/**
 * @brief Symbolic subtraction.
 * 
 * Overloaded subtraction operator that will subtract two symbolic polynomials.
 * The subtraction is simple concatenation of the terms (some will invert sign)
 * followed by optimization of the resulting polynomial (function @ref optimize).
 */
poly operator - (poly const & P, poly const & Q);

/**
 * @brief Symbolic multiplication.
 * 
 * This overloaded operator uses the fact that operation on symbolic terms
 * that consist only of power, exponential and goniometric function are closed,
 * i.e., result of multiplication is again a symbolic polynomial of this structure.
 * This holds due to product identities for goniometric functions.
 * @f[
 *     \sin(a+b) = \sin a \cos b + \cos a \sin b \ ,
 * @f] 
 * @f[
 *     \cos(a+b) = \cos a \cos b - \sin a \sin b \ .
 * @f]
 */
poly operator * (poly const & P, poly const & Q);

/**
 * @brief Division by a rational number.
 */
poly operator / (poly const & P, rational const & q);


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                             * 
 *                              Special functions                              *
 *                                                                             * 
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/**
 * @brief Associated Laguerre polynomial.
 * 
 * Return the symbolic representation of the associated Laguerre polynomial.
 * @f[
 *      L_k^{s}(x) = \sum_{j = s}^k (-1)^{j+1} \frac{(k!)^2 x^{j-s}}{(k-j)! j! (j-s)! } \ .
 * @f]
 * The associated Laguerre polynomials can be used to compute the hydrogen
 * wave function
 * @f[
 *      \psi_{nlm}(\mathbf{r})
 *      =
 *      \sqrt{\left(\frac{2Z}{a_0n}\right)^3 \frac{(n-l-1)!}{2n(n+l)!^3}}
 *      \left(\frac{2Zr}{na_0}\right)^l L_{n+l}^{2l+1}(2Zr/na_0) \mathrm{e}^{-Zr/na_0}
 *      Y_{lm}(\hat{\mathbf{r}}) \ .
 * @f]
 * Note the different normalization factor than in the expression with
 * the generalized Laguerre polynomials (see symbolic::GeneralizedLaguerre).
 */
poly AssociatedLaguerre (int k, int s);

/**
 * @brief Generalized Laguerre polynomial.
 * 
 * Return the symbolic representation of the generalized Laguerre polynomial.
 * @f[
 *      L_n^{(\alpha)}(x) = \sum_{j = 0}^n (-1)^j {n+\alpha \choose n-j} \frac{x^j}{j!} \ .
 * @f]
 * The generalized Laguerre polynomials can be used to compute the hydrogen
 * wave function
 * @f[
 *      \psi_{nlm}(\mathbf{r})
 *      =
 *      \sqrt{\left(\frac{2Z}{a_0n}\right)^3 \frac{(n-l-1)!}{2n(n+l)!}}
 *      \left(\frac{2Zr}{na_0}\right)^l L_{n-l-1}^{(2l+1)}(2Zr/na_0) \mathrm{e}^{-Zr/na_0}
 *      Y_{lm}(\hat{\mathbf{r}}) \ .
 * @f]
 * Note the different normalization factor than in the expression with
 * the associated Laguerre polynomials (see symbolic::AssociatedLaguerre).
 */
poly GeneralizedLaguerre (int n, int a);

/**
 * @brief Hydrogen radial function normalization factor.
 * 
 * Return the symbolic representation of the hydrogen radial function normalization
 * factor compatible with symbolic::AssociatedLaguerre,
 * @f[
 *     N_{nl} = \sqrt{\left(\frac{2}{n}\right)^3\frac{(n-l-1)!}{2n (n+l)!^3}} \ .
 * @f]
 */
term HydrogenN (int n, int l);

/**
 * @brief Un-normalized hydrogen radial function.
 * 
 * Return the symbolic representation of the hydrogen radial function that uses
 * symbolic::AssociatedLaguerre,
 * @f[
 *     \pi_{nl}(r) = \frac{ P_{nl}(r) }{ N_{nl} } = r \left(\frac{2r}{n}\right)^l L_{n+l}^{2l+1}(2r/n) \exp(-r/n) \ ,
 * @f]
 */
poly HydrogenNP (int n, int l);

/**
 * @brief Normalized hydrogen radial function.
 * 
 * Return the symbolic representation of the hydrogen radial function,
 * @f[
 *     P_{nl}(r) = r N_{nl} \left(\frac{2r}{n}\right)^l L_{n+l}^{2l+1}(2r/n) \exp(-r/n) \ ,
 * @f]
 * where the constant @f$ N_{nl} @f$ is computed by @ref HydrogenN ans the function
 * @f$ L_{n+l}^{2l+1} @f$ stands for the associated Laguerre polynomial (see definition
 * in symbolic::AssociatedLaguerre). Alternative formula uses the generalized Laguerre
 * polynomial @f$ L_{n-l-1}^{(2l+1)}(x) @f$; the resulting expressions differ only in the
 * normalization factor.
 */
poly HydrogenP (int n, int l);

/**
 * @brief Hydrogen sturmian function.
 * 
 * Return the symbolic representation of the hydrogen Sturmian function.
 * @f[
 *     S_n^{l}(r) = n^2 N_{nl} (2\lambda r)^l exp(-\lambda r) L_{n+l}^{2l+1}(2\lambda r) \ ,
 * @f]
 * where the constant @f$ N_{nl} @f$ is computed by @ref HydrogenN. As in @ref HydrogenP
 * do not interchange the associated Laguerre polynomial @f$ L_{n+l}^{2l+1} @f$ (used here)
 * and the generalized Laguerre polynomial @f$ L_{n-l-1}^{(2l+1)} @f$ used by other authors.
 * The formulas are identical when an appropriate normalization constant is used (see
 * symbolic::AssociatedLaguerre and symbolic::GeneralizedLaguerre).
 */
poly HydrogenS (int n, int l, cln::cl_RA lambda);

/**
 * @brief Return basis Laguerre orbital.
 * 
 * Return the symbolic reprezentation of the requested Laguerre basis orbital.
 * The orbital is defined by the expression
 * @f[
 *     \xi_N^{(L)}(x) = \sqrt{\frac{\lambda_L (N-1)!}{(2L + N - 1)!}}
 *     (\lambda_L x)^{L + 1} \exp(-\lambda_L x/2) L_{N-1}^{(2L+2)}(\lambda_L x)
 * @f]
 * The screening parameter @f$ \lambda_L @f$ needs to be given as a CLN's rational
 * number. Here @f$ L_{N-1}^{2L+2} @f$ is the generalized Laguerre polynomial.
 */
poly LaguerreBasisFunction (int N, int L, rational lambda);

/**
 * @brief Squared normalization factor.
 */
rational LaguerreBasisFunctionNsqr (int N, int L, rational lambda);

/**
 * @brief Riccati-Bessel function.
 * 
 * Constructs poly object that will represent the Riccat-Bessel function.
 * According to Abramovitz & Stegun §10.1.8 it is
 * @f[
 *     j_l (z) = z^{-1} \left[
 *          P(n+1/2,z) \sin(z - n\pi/2) + Q(n+1/2,z) \cos(z - n\pi/2)
 *     \right] \ ,
 * @f]
 * @f[
 *     P(n+1/2,z) = \sum_{k = 0}^{[n/2]} (-1)^k
 *     \frac{(n+2k)!}{(2k)!(n-2k+1)!}
 *     (2z)^{-2k} \ ,
 * @f]
 * @f[
 *     Q(n+1/2,z) = \sum_{k = 0}^{[(n-1)/2]} (-1)^k
 *     \frac{(n+2k+1)!}{(2k+1)!(n-2k)!}
 *     (2z)^{-2k-1} \ .
 * @f]
 * The arguments of the goniometric functions can be further simplified by the
 * summation formulas
 * @f[
 *     \sin(z-n\pi/2) = \sin(z) \cos(n\pi/2) - \cos(z)\sin(n\pi/2)
 *     = \cases{(-1)^{n/2}\sin(z) & $n$ even \cr (-1)^{(n+1)/2}\cos(z) & $n$ odd}
 * @f]
 * @f[
 *     \cos(z-n\pi/2) = \cos(z) \cos(n\pi/2) + \sin(z)\sin(n\pi/2)
 *     = \cases{(-1)^{n/2}\cos(z) & $n$ even \cr (-1)^{(n-1)/2}\sin(z) & $n$ odd}
 * @f]
 */
poly RiccatiBessel (int l, double k);

/**
 * @brief Binomial coefficient.
 * 
 * Return the exact rational representation of the binomial coefficient.
 * @f[
 *     \left(
 *         \matrix{n \cr k}
 *     \right) = \frac{n!}{k! (n-k)!} \ .
 * @f]
 */
rational combination (rational const & alpha, int k);


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                             * 
 *                            Symbolic integration                             *
 *                                                                             * 
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/**
 * @brief Integrate to infinity.
 * 
 * Returns a symbolic reprezentation of the semi-definite integral from "x" to infinity
 * (which is actually the indefinite integral subtracted from the definite integral
 * from zero to infinity),
 * @f[
 *     \int_x^\infty f(t) \mathrm{d}t \ .
 * @f]
 * The following identities are used:
 * <center>
 * <table border = "0">
 * <tr>
 *     <td width = "0" align = "right">@f$ \displaystyle \int_x^\infty t^a \mathrm{e}^{-ct} \mathrm{d}t @f$</td>
 *     <td width = "0">@f$ \displaystyle \ =\ @f$</td>
 *     <td width = "0" align = "left">@f$ \displaystyle \frac{a!}{c^{a+1}} \mathrm{e}^{-cx} \sum_{j = 0}^a \frac{c^j}{j!} @f$</td>
 *     <td width = "0" slign = "left">@f$ \displaystyle (a \in \mathbb{Z}_0^+, c \in \mathbb{R}^+) @f$</td>
 * </tr>
 * <tr>
 *     <td width = "0" align = "right">@f$ \displaystyle \int_x^\infty t^a \cos bt \,\mathrm{e}^{-ct} \mathrm{d}t @f$</td>
 *     <td width = "0">@f$ \displaystyle \ =\ @f$</td>
 *     <td width = "0" align = "left">@f$ \displaystyle \mathrm{Re}\;\left(\frac{a!}{(c-\mathrm{i}|b|)^{a+1}} \mathrm{e}^{-cx} \sum_{j = 0}^a \frac{(c-\mathrm{i}|b|)^j}{j!} \cos bx\right) - \mathrm{Im}\;\left(\frac{a!}{(c-\mathrm{i}|b|)^{a+1}} \mathrm{e}^{-cx} \sum_{j = 0}^a \frac{(c-\mathrm{i}|b|)^j}{j!} \sin bx\right) @f$</td>
 *     <td width = "0" slign = "left">@f$ \displaystyle (a \in \mathbb{Z}_0^+, b \in \mathbb{R}, c \in \mathbb{R}^+) @f$</td>
 * </tr>
 * <tr>
 *     <td width = "0" align = "right">@f$ \displaystyle \int_x^\infty t^a \sin bt \,\mathrm{e}^{-ct} \mathrm{d}t @f$</td>
 *     <td width = "0">@f$ \displaystyle \ =\ @f$</td>
 *     <td width = "0" align = "left">@f$ \displaystyle \left[\mathrm{Im}\;\left(\frac{a!}{(c-\mathrm{i}|b|)^{a+1}} \mathrm{e}^{-cx} \sum_{j = 0}^a \frac{(c-\mathrm{i}|b|)^j}{j!} \cos bx\right) + \mathrm{Re}\;\left(\frac{a!}{(c-\mathrm{i}|b|)^{a+1}} \mathrm{e}^{-cx} \sum_{j = 0}^a \frac{(c-\mathrm{i}|b|)^j}{j!} \sin bx \right)\right]\; \mathrm{sign}\;b\qquad @f$</td>
 *     <td width = "0" slign = "left">@f$ \displaystyle (a \in \mathbb{Z}_0^+, b \in \mathbb{R} \setminus \{0\}, c \in \mathbb{R}^+) @f$</td>
 * </tr>
 * </table>
 * </center>
 */
poly integrate_inf (poly const & P);

/**
 * @brief Integrate from zero.
 * 
 * Returns a symbolic reprezentation of the semi-definite integral from zero to "x"
 * (which is actually the indefinite integral),
 * @f[
 *     \int_0^x f(t) \mathrm{d}t \ .
 * @f]
 * This integral is computes simply as a difference of the results of @ref integrate_full and @ref integrate_inf.
 */
poly integrate_low (poly const & P);

/**
 * @brief Integrate over positive real numbers.
 * 
 * Integrate the function represented as a symbolic polynomial from 0 to inifinity.
 * The result is returned as a single term with the multiplicative constant set to
 * the value of the integral. The following identities are used:
 * <center>
 * <table width = "0" border = "0">
 * <tr>
 *     <td width = "0" align = "right">@f$ \displaystyle \int_0^\infty x^a \mathrm{e}^{-cx} \mathrm{d}x @f$</td>
 *     <td width = "0">@f$ \displaystyle \ =\  @f$</td>
 *     <td>@f$ \displaystyle \frac{a!}{c^{a+1}} @f$</td>
 *     <td align = "left">@f$ \displaystyle (a \in \mathbb{Z}_0^+, c \in \mathbb{R}^+) @f$</td>
 * </tr>
 * <tr>
 *     <td align = "right">@f$ \displaystyle \int_0^\infty x^a \sin bx\; \mathrm{e}^{-cx} \mathrm{d}x @f$</td>
 *     <td>@f$ \displaystyle \ =\  @f$</td>
 *     <td>@f$ \displaystyle \mathrm{Im}\,\frac{a!}{(c-\mathrm{i}|b|)^{a+1}} \mathrm{sign}\; b \qquad @f$</td>
 *     <td align = "left">@f$ \displaystyle (a \in \mathbb{Z}_0^+, b \in \mathbb{R} \setminus \{0\}, c \in \mathbb{R}^+) @f$</td>
 * </tr>
 * <tr>
 *     <td align = "right">@f$ \displaystyle \int_0^\infty x^a \cos bx\; \mathrm{e}^{-cx} \mathrm{d}x @f$</td>
 *     <td>@f$ \displaystyle \ =\  @f$</td>
 *     <td>@f$ \displaystyle \mathrm{Re}\,\frac{a!}{(c-\mathrm{i}|b|)^{a+1}} @f$</td>
 *     <td align = "left">@f$ \displaystyle (a \in \mathbb{Z}_0^+, b \in \mathbb{R}, c \in \mathbb{R}^+) @f$</td>
 * </tr>
 * <tr>
 *     <td align = "right">@f$ \displaystyle \int_0^\infty \frac{1}{x} \sin bx\; \mathrm{e}^{-cx} \mathrm{d}x @f$</td>
 *     <td>@f$ \displaystyle \ = \ @f$</td>
 *     <td>@f$ \displaystyle \left(\frac{\pi}{2} - \arctan \frac{c}{b}\right) \mathrm{sign}\;b @f$</td>
 *     <td align = "left">@f$ \displaystyle (b \in \mathbb{R} \setminus \{0\}, c \in \mathbb{R}_0^+) @f$</td>
 * </tr>
 * <tr>
 *     <td align = "right">@f$ \displaystyle \int_0^\infty \frac{1}{x} (1 - \cos bx)\; \mathrm{e}^{-cx} \mathrm{d}x @f$</td>
 *     <td>@f$ \displaystyle \ =\  @f$</td>
 *     <td>@f$ \displaystyle \frac{1}{2} \ln \left(1 + \frac{b^2}{c^2}\right) @f$</td>
 *     <td align = "left">@f$ \displaystyle (b \in \mathbb{R}, c \in \mathbb{R}^+) @f$</td>
 * </tr>
 * </table>
 * </center>
 * The value @f$ b = 0 @f$ results in zero integral for integrands containing sine.
 */
term integrate_full (poly const & P);

} // endof namespace symbolic

#endif
