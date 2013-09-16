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

#ifndef HEX_CHEBYSHEV
#define HEX_CHEBYSHEV

#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

#include <fftw3.h>

#include "arrays.h"
#include "misc.h"

/**
 * @brief Chebyshev approximation.
 * 
 * This class manages a Chebyshev approximation
 * @f[
 *     F(x) = \frac{c_0}{2} + \sum_{k=1}^N c_k T_k(x)
 * @f]
 * of a function F with signature
 * @code
 * Tout F(double Tin)
 * @endcode
 */
template <typename Tin, typename Tout> class Chebyshev
{
public:
    
    Chebyshev () : N(0), C(), xt(0.), m(0.) {}
    Chebyshev (Chebyshev const & cb) : N(cb.N), C(cb.C), xt(cb.xt), m(cb.m) {}
    
    /**
     * @brief Constructor.
     * 
     * @param f Function to approximate.
     * @param n How many Chebyshev nodes to use.
     * @param a Left boundary of the approximation interval.
     * @param b Right boundary of the approximation interval.
     */
    template <class Functor> Chebyshev (Functor const & f, int n, Tin a, Tin b)
    {
        generate(f, n, a, b);
    }
    
    /**
     * @brief Constructor from reference to array.
     * 
     * @param array Array of precomputed Chebyshev coefficients.
     * @param a Left boundary of the approximation interval.
     * @param b Right boundary of the approximation interval.
     */
    Chebyshev (ArrayView<Tout> const & array, Tin a, Tin b)
    {
        N  = array.size();
        xt = 0.5 * (b + a);
        m  = 0.5 * (b - a);
        
        // copy coefficients
        C = array;
    }
    
    /**
     * @brief Compute the transformation.
     * 
     * Generate and store the Chebyshev coefficients.
     * @param f Function of the signature Tout(*)(Tin) to be approximated.
     * @param n Number of the coefficients to compute.
     * @param a Left boundary of the approximation inerval.
     * @param b Right boundary of the approximation inerval.
     */
    template <class Functor>
    void generate (Functor const & f, int n, Tin a, Tin b);
    
    /**
     * Return full approximation value.
     */
    Tout operator() (Tin const & x) const
    {
        Tout ret = 0.5 * C[0];
        const double xp = scale (x);
        
        for (int k = 1; k < N; k++)
        {
            double Tk_x = cos(k * acos(xp));
            ret += C[k] * Tk_x;
        }
        return ret;
    }
    
    /**
     * Use Clenshaw recurrence formula for evaluation of 'm' terms.
     * The formula has the advantage of not evaluating goniometric funtions.
     */
    Tout clenshaw (Tin const & x, int const & m) const
    {
        Tout d_j = 0, d_jp1 = 0, d_jp2 = 0;
        const double one_x = scale(x);
        const double two_x = 2 * one_x; // due to linearity of 'scale'
        
        for (int j = m - 1; j >= 1; j--)
        {
            d_j   = two_x * d_jp1 - d_jp2 + C[j];
            
            d_jp2 = d_jp1;
            d_jp1 = d_j;
        }
        
        d_j = one_x * d_jp1 - d_jp2 + 0.5 * C[0];
        
        return d_j;
    }
    
    /**
     * Get index of the first Chebyshev approximation coefficient
     * that is smaller than 'eps' times the sum of fabs(C[k])
     * truncated after this term. If no such term exists, the total
     * count of terms is returned.
     *
     * @note The Chebyshev polynomial corresponding to the last considered
     * term can be negligible near some evaluation point 'x'. In that case,
     * its contribution might be shadowed by the contribution of the following
     * polynomial. So, the evaluated result can have worse precision than the
     * requested 'eps'. Nevertheless, the more polynomials get involved, the
     * less is this fact problematic.
     */
    int tail (double eps) const
    {
        double sum = std::abs(0.5 * C[0]);
        double abs_Ck;
        
        for (int k = 1; k < N; k++)
        {
            abs_Ck = std::abs(C[k]);
            sum += abs_Ck;
            
            if (abs_Ck <= eps * sum)
                return k + 1;
        }
        
        return N;
    }
    
    /**
     * Return approximated value where the terms of magnitude less than 'eps'
     * are discarded.
     */
    Tout approx (double x, double eps, int * n = 0) const
    {
        Tout ret = 0.5 * C[0];
        double xp = scale(x);
        
        int k;
        for (k = 1; k < N; k++)
        {
            double Tk_x = cos(k * acos(xp));
            Tout delta = C[k] * Tk_x;
            ret += delta;
            
            if (std::abs(delta) <= eps * std::abs(ret))
            {
                k++;
                break;
            }
        }
        
        if (n != 0)
            *n = k;
        
        return ret;
    }
    
    /**
     * Integration types.
     */
    typedef enum {
        Integ_Indef,
        Integ_Low,
        Integ_High
    } Integration;
    
    /**
     * Return Chebyshev aproximation of the function primitive to the
     * stored Chebyshev approximation.
     * 
     * @param itype Whether to return general indefinite integral expansion
     * @f[
     *     \int_a^x f(x) \mathrm{d}x \ ,
     * @f]
     * (corresponds to "indef") or definite low
     * @f[
     *     \int_0^x f(x) \mathrm{d}x \ ,
     * @f]
     * (corresponds to "def_low") or definite high
     * @f[
     *     \int_x^\infty f(x) \mathrm{d}x \ .
     * @f]
     */
    Chebyshev integrate (Integration itype = Integ_Indef) const
    {
        Chebyshev ret;
        ret.xt = xt;
        ret.m  = m;
        ret.N  = N;
        
        ret.C.resize(ret.N);
        ret.C[0] = 0;
        ret.C[N-1] = ret.m * C[N-2] / (2.*(N-2.));
        
        for (int i = 1; i < N - 1; i++)
            ret.C[i] = ret.m * (C[i-1] - C[i+1]) / (2.*i);
        
        // transform coefficients to get a definite integral, if requested
        switch (itype)
        {
            case Integ_Indef:
            {
                // already done
                break;
            }
            case Integ_Low:
            {
                // compute C[0] as an alternating sum of other coefficients
                for (int i = 1; i < N; i += 2)
                    ret.C[0] += ret.C[i];
                for (int i = 2; i < N; i += 2)
                    ret.C[0] -= ret.C[i];
                break;
            }
            case Integ_High:
            {
                // compute C[0] as a sum of other coefficients and negate other coefficients
                for (int i = 1; i < N; i++)
                {
                    ret.C[0] += ret.C[i];
                    ret.C[i] = -ret.C[i];
                }
                break;
            }
        }
        
        // take into account the normalization (F = c_0/2 + ...)
        ret.C[0] *= 2.;
        
        return ret;
    }
    
    /**
     * Get Chebyshev root in the interval (x1,x2).
     * @param N Order of the polynomial.
     * @param k index of the root.
     * @param x1 Left bound of the interval.
     * @param x2 Right bound of the interval.
     */
    static Tin root (int N, int k, Tin x1 = 0., Tin x2 = 1.)
    {
        return x1 + 0.5 * (1. + cos(M_PI * (k + 0.5) / N)) * (x2 - x1);
    }
    
    /**
     * Write out the coefficients.
     */
    std::string str() const
    {
        std::ostringstream out;
        
        out << "[";
        if (N > 0)
        {
            for (int i = 0; i < N - 1; i++)
                out << C[i] << ", ";
            out << C[N-1];
        }
        out << "]" << std::endl;
        
        return out.str();
    }
    
    NumberArray<Tout> const & coeffs() const
    {
        return C;
    }
    
    /**
     * Get k-th root of the N-order Chebyshev polynomial.
     */
    inline static double node (int k, int N)
    {
        return cos(M_PI * (k + 0.5) / N);
    }
    
private:
    
    /// map interval (xt-m,xt+m) to (-1,1)
    inline double scale(Tin x) const {
        return (x - xt) / m;
    }
    
    /// map interval (-1,1) to (xt-m,xt+m)
    inline Tin unscale(double x) const {
        return (xt + m*x);
    }
    
    /// coefficient number
    int N;
    
    /// Chebyshev coefficients
    NumberArray<Tout> C;
    
    /// approximation interval center
    Tin xt;
    
    /// approximation interval half-width
    Tin m;
};

//
// template member specializations
//

/**
 * @brief Chebyshev approximation of a given real function.
 */
template<> template <class Functor> 
void Chebyshev<double,double>::generate (Functor const & f, int n, double a, double b)
{
    N  = n;
    xt = 0.5 * (b + a);
    m  = 0.5 * (b - a);
    C.resize(N);
    
    // input array
    rArray fvals(N);
    
    // create the FFTW plan
    fftw_plan plan = fftw_plan_r2r_1d (N, &fvals[0], &C[0], FFTW_REDFT10, 0);
    
    // evaluate nodes and function
    double pi_over_N = M_PI / N;
    for (int k = 0; k < N; k++)
    {
        double xk = cos(pi_over_N * (k + 0.5));
        fvals[k] = f(unscale(xk));
    }
    
    // execute the transform
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    // normalize
    double scal = 1./N;
    for (double & c : C)
        c *= scal;
}

/**
 * @brief Chebyshev approximation of a given complex function.
 */
template<> template <class Functor>
void Chebyshev<double,Complex>::generate (Functor const & f, int n, double a, double b)
{
    N  = n;
    xt = 0.5 * (b + a);
    m  = 0.5 * (b - a);
    C.resize(N);
    
    // input array
    cArray fvals(4*N), ftraf(4*N);
    
    // create the FFTW plan
    fftw_plan plan = fftw_plan_dft_1d (
        4*N,
        reinterpret_cast<fftw_complex*>(&fvals[0]),
        reinterpret_cast<fftw_complex*>(&ftraf[0]),
        FFTW_FORWARD,
        0
    );
    
    // evaluate nodes and function
    double pi_over_N = M_PI / N;
    for (int k = 0; k < N; k++)
    {
        double xk = cos(pi_over_N * (k + 0.5));
        fvals[2*k+1] = fvals[4*N-(2*k+1)] = f(unscale(xk));
    }
    
    // execute the transform
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    // copy normalized coefficients
    double scal = 1./N;
    for (int i = 0; i < N; i++)
        C[i] = ftraf[i] * scal;
}

#endif
