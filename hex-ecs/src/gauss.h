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

#ifndef HEX_GAUSS
#define HEX_GAUSS

#include <map>
#include <vector>

#include "bspline.h"

class GaussLegendreData
{
    public:
        
        /**
         * @brief Retrieve Gauss-Legendre data.
         * 
         * Precompute Gauss-Legendre quadrature data (if not already done
         * is some previous call) and return pointers to the chache table.
         * 
         * @param points Gauss-Legendre points half-count. If too low/high, the return value
         *               will contain the (used) nearest implemented value.
         * @param vx     On return, the Gauss-Legendre nodes (nonnegative half of them).
         * @param vw     On return, the corresponding Gauss-Legendre weights.
         */
        int gauss_nodes_and_weights (int points, const double* &vx, const double* &vw) const;
    
    private:
        
        // precomputed nodes and weights, common to all instances
        static std::vector<std::pair<double*,double*>> data_;
};

/**
 * @brief Gauss-Legendre quadrature.
 * 
 * This class's purpose is aid in fixed-order Gauss-Legendre quadrature.
 * It is constructed above a B-spline environment, which is passed to the
 * constructor. The functions "p_points" and "p_weights" return evaluation nodes
 * and weights, respectively, and the functions "quad" do the fixed-order
 * quadrature itself. Every call to p_points or p_weights resuts in call
 * to gauss_nodes_and_weights, which uses GSL to get the requested data.
 * The computed nodes and weights are stored in a cache table to allow
 * faster subsequent computations.
 */
class GaussLegendre : public GaussLegendreData
{
    public:
        
        // constructor
        GaussLegendre(Bspline const & bspline)
            : bspline_(bspline) {}
            
        /**
         * @brief Get Gauss-Legendre quadrature points in interval.
         * 
         * Map Gauss-Legendre points provided by @ref gauss_nodes_and_weights
         * to a complex interval @f$ (x_1,x_2) @f$.
         * 
         * @param points Number of Gauss-Legendre points.
         * @param x1 Left boundary of the interval.
         * @param x2 Right boundary of the interval.
         */
        cArray p_points (int points, Complex x1, Complex x2) const;

        /**
         * @brief Get Gauss-Legendre quadrature weights in interval.
         * 
         * Map Gauss-Legendre weights provided by @ref gauss_nodes_and_weights
         * to a complex interval @f$ (x_1,x_2) @f$.
         * 
         * @param points Number of Gauss-Legendtre points.
         * @param x1 Left boundary of the interval.
         * @param x2 Right boundary of the interval.
         */
        cArray p_weights (int points, Complex x1, Complex x2) const;

        /**
         * @brief Gauss-Legendre integrator
         * 
         * Integrate the given function.
         * 
         * @param f Function of type (void (*) (unsigned, complex*, complex*, data...)),
         *          where the first parameter is number of points to evaluate, the second parameter
         *          is pointer to an array of double complex values at which to evaluate, the third
         *          is pointer to an array into which the data will be written (responsibility for
         *          memory management is on callers side) and finally the fourth (and further) parameter is 
         *          any other data to be suplied to the function.
         * @param data Data to pass to the function.
         * @param points Gauss-Legendre points count.
         * @param iknot Knot index.
         * @param x1 Left integration boundary.
         * @param x2 Right integration boundary.
         * 
         * It must be
         * @code
         *    t[iknot] <= x1 <= t[iknot+1]
         *    t[iknot] <= x2 <= t[iknot+1]
         * @endcode
         */
        template <class Functor, class... Data> Complex quad
        (
            Functor f,
            int points, int iknot, Complex x1, Complex x2,
            Data... data
        ) const
        {
            // check boundaries
            if (x1.real() < bspline_.t(iknot).real() or bspline_.t(iknot+1).real() < x1.real() or
                x2.real() < bspline_.t(iknot).real() or bspline_.t(iknot+1).real() < x2.real())
            {
                HexException("[quad] Error: boundaries not for this iknot!");
            }
            
            // get evaluation points and weights
            cArray xs = p_points(points, x1, x2);
            cArray ws = p_weights(points, x1, x2);
            cArray values (points);
            
            // evaluate the function
            f(points, xs.data(), values.data(), data...);
            
            // pointers for fast restricted access
            Complex const * const restrict pw = ws.data();
            Complex const * const restrict py = values.data();
            
            // sum the results
            Complex result = 0.;
            for (int ipt = 0; ipt < points; ipt++)
                result += pw[ipt] * py[ipt];
            
            return result;
        }
        
        /**
         * @brief Gauss-Legendre integrator.
         * 
         * This is a variant of the other quad function that can be used together with
         * the pointer-to-member-function arguments. A typical usage woud be
         * @code
         * Class A { ... };
         * A a;
         * result = quad (&a, &a::integrand, ...);
         * @endcode
         */
        template <class ClassPtr, class Functor, class... Data> Complex quadMFP
        (
            ClassPtr ptr, Functor f,
            int points, int iknot, Complex x1, Complex x2,
            Data... data
        ) const
        {
            // check boundaries
            if (x1.real() < bspline_.t(iknot).real() or bspline_.t(iknot+1).real() < x1.real() or
                x2.real() < bspline_.t(iknot).real() or bspline_.t(iknot+1).real() < x2.real())
            {
                HexException("[quad] Error: boundaries not for this iknot!");
            }
            
            // get evaluation points and weights
            cArray xs = p_points(points, x1, x2);
            cArray ws = p_weights(points, x1, x2);
            cArray values (points);
            
            // evaluate the function
            (ptr->*f)(points, xs.data(), values.data(), data...);
            
            // pointers for fast restricted access
            Complex const * const restrict pw = ws.data();
            Complex const * const restrict py = values.data();
            
            // sum the results
            Complex result = 0.;
            for (int ipt = 0; ipt < points; ipt++)
                result += pw[ipt] * py[ipt];
            
            return result;
        }
        
    private:

        // B-spline environment
        Bspline const & bspline_;
};

#endif
