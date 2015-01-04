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

#ifndef HEX_CFRACT
#define HEX_CFRACT

#include <cmath>
#include <complex>
#include <vector>

#include <gsl/gsl_sf.h>

/**
 * @brief Evaluate polynomial.
 * 
 * Evaluates polynomial specified as array of coefficients. The order of the
 * oefficients is expected from highest order (a[N]) to lowest order (a[0]).
 * 
 * This is a function template, which accepts arbitrary data types (supporting basic
 * arithmetic).
 * 
 * @param poly Array of polynomials.
 * @param z Evaluation point.
 */
template <class T, class U> auto eval_poly (std::vector<T> poly, U z) -> decltype(T(0)*U(0))
{
    // polynomial length (highest order)
    unsigned N = poly.size();
    
    // skip empty polynomials
    if (N == 0)
        return 0.;
    
    // initialize by highest-order coefficient
    decltype(T(0)*U(0)) result = poly[0];
    
    // if N > 1, iterate over the remaining coefficients
    for (unsigned i = 1; i < N; i++)
        result = result * z + poly[i];
    
    return result;
}

/**
 * @brief Evaluate power series using the method of continues fractions.
 * 
 * This function template implements the CFRACT algorithm as presented in
 * Nesbet R. K., Analytical evaluation of integrals over Coulomb wave functions, Comput. Phys. Commun. 52 (1988) 29--33.
 * 
 * @param GetCoeff Function of functor (signature "unsigned int -> T") that will return i-th coefficient of the
 *                 power series to evaluate. The calls are consecutive in the order (0, 1, 2, ...).
 * @param z Evaluation point.
 * @param tolelance Absolute convergence criterium.
 * @param maxiter Maximal number of iterations.
 * @param pn Optional pointer to an unsigned variable. If available, the actual number of iterations will be stored there.
 */
template <class T, class TGetCoeff> T cfract (TGetCoeff GetCoeff, T z, double tolerance = 1e-5, unsigned maxiter = 1000, unsigned * pn = nullptr)
{
    // polynomial coefficients
    std::vector<T> G = { GetCoeff(0) };
    
    // continued fraction coefficients
    std::vector<T> A = { G[0] };
    
    // collected numerator and denominator of the continued fraction
    std::vector<T> C_prev = { 0. }, C_curr = { G[0] };
    std::vector<T> D_prev = { 1. }, D_curr = { 1. };
    
    // remainders
    T r_curr, r_prev = G[0];
    
    // convergence loop
    for (unsigned n = 1; ; n++)
    {
        // get next polynomial coefficient
        G.push_back(GetCoeff(n));
        
        // calculate remainder
        r_curr = 0.;
        for (unsigned i = 0; i <= n / 2; i++)
            r_curr += G[n-i] * (n-i < n ? D_curr[D_curr.size()-1-i] : T(1));
        
        // add new continued fraction coefficient
        T a = -r_curr / r_prev;
        
        // multiply C_prev and D_prev by 'z'
        C_prev.push_back(0);
        D_prev.push_back(0);
        
        // multiply C_prev and D_prev by 'a'
        for (unsigned i = 0; i < C_prev.size(); i++)
            C_prev[i] *= a;
        for (unsigned i = 0; i < D_prev.size(); i++)
            D_prev[i] *= a;
        
        // new polynomials C and D
        std::vector<T> C_new(std::max(C_curr.size(),C_prev.size()));
        std::vector<T> D_new(std::max(D_curr.size(),D_prev.size()));
        
        // C_new = C_curr + (a * z * C_prev), end-aligned sum
        for (unsigned i = 0; i < C_new.size(); i++)
        {
            C_new[C_new.size() - 1 - i] = (i < C_curr.size() ? C_curr[C_curr.size() - 1 - i] : T(0))
                                        + (i < C_prev.size() ? C_prev[C_prev.size() - 1 - i] : T(0));
        }
        
        // D_new = D_curr + (a * z * D_prev), end-aligned sum
        for (unsigned i = 0; i < D_new.size(); i++)
        {
            D_new[C_new.size() - 1 - i] = (i < D_curr.size() ? D_curr[D_curr.size() - 1 - i] : T(0))
                                        + (i < D_prev.size() ? D_prev[D_prev.size() - 1 - i] : T(0));
        }
        
        // store previous polynomials
        C_prev = C_curr; C_curr = C_new;
        D_prev = D_curr; D_curr = D_new;
        
        // evaluate polynomials
        T c_curr = eval_poly(C_curr,z), c_prev = eval_poly(C_prev,z);
        T d_curr = eval_poly(D_curr,z), d_prev = eval_poly(D_prev,z);
        
        // check numerics
        if (not std::isfinite(std::abs(c_curr)) or not std::isfinite(std::abs(d_curr)))
        {
            if (pn != nullptr)
                *pn = maxiter;
            return c_prev / d_prev;
        }
        
        // check convergence and iteration limit
        double diff = std::abs(c_curr * d_prev - c_prev * d_curr) / std::abs(d_curr * d_prev);
        if (diff < tolerance or n >= maxiter)
        {
            if (pn != nullptr)
                *pn = n;
            return c_curr / d_curr;
        }
        
        // save remainder
        r_prev = r_curr;
    }
}

/**
 * @brief Evaluation of power series.
 * 
 * This function is just for confronting with the @ref cfract function. It is not intended for usage,
 * because the other is superior to this one. The interface is identical as for @ref cfract.
 */
template <class TGetCoeff, class T> T series (TGetCoeff GetCoeff, T z, double tolerance = 1e-5, unsigned maxiter = 1000, unsigned * pn = nullptr)
{
    std::vector<T> coeffs;
    T prev_eval = 0;
    
    for (unsigned n = 0; ; n++)
    {
        coeffs.resize(n + 1);
        for (int i = n; i >= 1; i--)
            coeffs[i] = coeffs[i-1];
        coeffs[0] = GetCoeff(n);
        
        T eval = eval_poly(coeffs,z);
        double diff = std::abs(prev_eval - eval);
        prev_eval = eval;
        
        if (diff < tolerance or n >= maxiter)
        {
            if (pn != nullptr)
                *pn = n;
            return eval;
        }
    }
}

#endif
