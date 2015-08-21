//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2015, Jakub Benda, Charles University in Prague                    //
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

#ifndef HEX_CHEBYSHEV
#define HEX_CHEBYSHEV

#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_complex.h>

#include "arrays.h"
#include "misc.h"
#include "special.h"

/**
 * @brief Chebyshev approximation.
 * 
 * This class manages a Chebyshev approximation
 * @f[
 *     F(x) = \frac{c_0}{2} + \sum_{k=1}^N c_k T_k(x)
 * @f]
 * of a function F with signature
 * @code
 * Tout F (Tin x)
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
     * Generate and store the Chebyshev coefficients. After evaluation
     * of the approximated function in the Chebyshev nodes, the computation
     * of Chebyshev expansion coefficients is done using the fast Fourier
     * transform provided by FFTW.
     * 
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
        const double xp = scale(x);
        
        for (int k = 1; k < N; k++)
        {
            double Tk_x = cos(k * acos(xp));
            ret += C[k] * Tk_x;
        }
        return ret;
    }
    
    /**
     * @brief Evaluate expansion.
     * 
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
     * @brief Get significant number of coefficients.
     * 
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
     * @brief Evaluate expansion.
     * 
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
            double Tk_x = std::cos(k * std::acos(xp));
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
    typedef enum
    {
        Integ_Indef,
        Integ_Low,
        Integ_High
    }
    Integration;
    
    /**
     * @brief Get expansion of integral of the approximated function.
     * 
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
     * @brief Chebyshev node in interval.
     * 
     * Get Chebyshev root in the interval (x1,x2).
     * @param N Order of the polynomial.
     * @param k index of the root.
     * @param x1 Left bound of the interval.
     * @param x2 Right bound of the interval.
     */
    static Tin root (int N, int k, Tin x1 = 0., Tin x2 = 1.)
    {
        return x1 + 0.5 * (1. + std::cos(special::constant::pi * (k + 0.5) / N)) * (x2 - x1);
    }
    
    /**
     * @brief Convert to string.
     * 
     * Write out the coefficients to std::string.
     */
    std::string str () const
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
    
    /// Return reference to the coefficient array.
    NumberArray<Tout> const & coeffs () const
    {
        return C;
    }
    
    /**
     * @brief Chebyshev node from (0,1).
     * 
     * Get k-th root of the N-order Chebyshev polynomial.
     */
    inline static double node (int k, int N)
    {
        return cos(special::constant::pi * (k + 0.5) / N);
    }
    
private:
    
    /// map interval (xt-m,xt+m) to (-1,1)
    inline double scale (Tin x) const
    {
        return (x - xt) / m;
    }
    
    /// map interval (-1,1) to (xt-m,xt+m)
    inline Tin unscale (double x) const
    {
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
 * 
 * Uses the function "fftw_plan_r2r_1d" for real data.
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
    
    // evaluate nodes and function
    double pi_over_N = special::constant::pi / N;
    for (int k = 0; k < N; k++)
    {
        double xk = std::cos(pi_over_N * (k + 0.5));
        fvals[k] = f(unscale(xk));
    }
    
    // execute the transform
    gsl_fft_real_radix2_transform(fvals.data(), 1, N);
    
    // normalize
    C *= 1. / N;
}

/**
 * @brief Chebyshev approximation of a given complex function.
 * 
 * Uses the function "fftw_plan_dft_1d" for complex data. Slower than the
 * real "generate".
 */
template<> template <class Functor>
void Chebyshev<double,Complex>::generate (Functor const & f, int n, double a, double b)
{
    N  = n;
    xt = 0.5 * (b + a);
    m  = 0.5 * (b - a);
    C.resize(N);
    
    // input array
    cArray fvals(4*N);
    
    // evaluate nodes and function
    double pi_over_N = special::constant::pi / N;
    for (int k = 0; k < N; k++)
    {
        double xk = std::cos(pi_over_N * (k + 0.5));
        fvals[2*k+1] = fvals[4*N-(2*k+1)] = f(unscale(xk));
    }
    
    // execute the transform
    gsl_fft_complex_radix2_forward(reinterpret_cast<double*>(fvals.data()), 1, 4 * N);
    
    // copy normalized coefficients
    C = fvals / double(N);
}

#endif /* HEX_CHEBYSHEV */
