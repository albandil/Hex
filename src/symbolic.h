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

#define GF_NONE   0
#define GF_SIN    1
#define GF_COS    2

namespace symbolic
{

/// Rational number (CLN)
typedef cln::cl_RA rational;

/// Rational constant 1/2
extern rational onehalf;

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
    
    //
    // constructors
    //
    
    /**
     * @brief Default constructor.
     * 
     * Create a new symbolic term that will evaluate to zero.
     */
    term() : ki(1), kr(0), a(0), gf(GF_NONE), b(0), c(0) {}
    
    /**
     * @brief Constructo from fraction.
     * 
     * Create a new symbolic term that will evaluate to a rational number.
     */
    term(rational x) : ki(1), kr(x), a(0), gf(GF_NONE), b(0), c(0) {}
    
    //
    // data fields
    //
    
    /// constant irrational multiplication factor
    double ki;
    
    /// constant rational multiplication factor
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
    std::ostream & operator << (std::ostream & os) const;
};

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
 * @brief Associated laguerre polynomial.
 * 
 * Return the symbolic representation of the associated Laguerre polynomial.
 * @f[
 *      L_k^{(s)}(x) = \sum_{j = s}^k (-1)^j \frac{(k!)^2 x^{j-s}}{(k-j)! j! (j-s)! } \ .
 * @f]
 */
poly Laguerre (int k, int s);

/**
 * @brief Hydrogen radial function normalization factor.
 * 
 * Return the symbolic representation of the hydrogen radial function normalization
 * factor,
 * @f[
 *     N_{nl} = \sqrt{\left(\frac{2}{n}\right)^3\frac{(n-l-1)!}{2n (n+l)!^3}} \ .
 * @f]
 */
term HydrogenN (int n, int l);

/**
 * @brief Hydrogen radial function.
 * 
 * Return the symbolic representation of the hydrogen radial function,
 * @f[
 *     P_{nl}(r) = r N_{nl} \left(\frac{2r}{n}\right)^l L_{n+l}^{2l+1}(2r/n) \exp(-r/n) \ ,
 * @f]
 * where the constant @f$ N_{nl} @f$ is computed by @ref HydrogenN. Alternative formula
 * uses the <i>generalized</i> Laguerre polynomial @f$ L_{n-l-1}^{(2l+1)}(x) @f$;
 * however, the resulting expressions are identical.
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
 * The formulas are identical, nevertheless.
 */
poly HydrogenS (int n, int l, cln::cl_RA lambda);

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
poly LaguerreBasisFunction (int N, int L, cln::cl_RA lambda);

//
// integrals
//

/**
 * @brief Integrate to infinity.
 * 
 * Returns a symbolic reprezentation of the semi-definite integral from "x" to infinity
 * (which is actually the indefinite integral subtracted from the definite integral
 * from zero to infinity),
 * @f[
 *     \int_x^\infty f(x') \mathrm{d}x' \ .
 * @f]
 */
poly integrate_inf (poly const & P);

/**
 * @brief Integrate from zero.
 * 
 * Returns a symbolic reprezentation of the semi-definite integral from zero to "x"
 * (which is actually the indefinite integral),
 * @f[
 *     \int_0^x f(x') \mathrm{d}x' \ .
 * @f]
 */
poly integrate_low (poly const & P);

/**
 * @brief Integrate over positive real numbers.
 * 
 * Integrate the function represented as a symbolic polynomial from 0 to inifinity.
 * The result is returned as a single term with the multiplicative constant set to
 * the value of the integral.
 */
term integrate_full (poly const & P);

/**
 * @brief Evaluate the polynomial.
 * 
 * Evaluate the polynomial for a given value of the variable. The evaluation is
 * done term after term and the results are summed afterwards.
 */
double eval (poly const & P, double r);

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

} // endof namespace symbolic

#endif
