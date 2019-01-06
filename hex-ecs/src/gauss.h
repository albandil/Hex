//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2016, Jakub Benda, Charles University in Prague                    //
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
         * @param points Gauss-Legendre points half-count. If too low/high, the return value
         *               will contain the (used) nearest implemented value.
         * @param vx     On return, the Gauss-Legendre nodes (nonnegative half of them).
         * @param vw     On return, the corresponding Gauss-Legendre weights.
         */
        static void gauss_nodes_and_weights (int points, const Real*& vx, const Real*& vw);

        /**
         * @brief Precalculate nodes and weights so that the retrieval is fast and thread-safe.
         * 
         * Uses the GSL routine \c gsl_integration_glfixed_table_alloc.
         */
        static void precompute_nodes_and_weights (int points);

    private:

        // precomputed nodes and weights, common to all instances
        static std::vector<std::pair<Real*,Real*>> data_;
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

        /**
         * @brief Get Gauss-Legendre quadrature points and weights in interval.
         * 
         * Map Gauss-Legendre points provided by @ref gauss_nodes_and_weights
         * to a complex interval @f$ (x_1,x_2) @f$.
         * 
         */
        void scaled_nodes_and_weights (int points, Complex x1, Complex x2, Complex* xs, Complex* ws) const;

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
            Functor f, Bspline const & bspline,
            int points, int iknot, Complex x1, Complex x2,
            Data... data
        ) const
        {
            // check boundaries
            if (x1.real() < bspline.t(iknot).real() or bspline.t(iknot+1).real() < x1.real() or
                x2.real() < bspline.t(iknot).real() or bspline.t(iknot+1).real() < x2.real())
            {
                HexException("[quad] Error: boundaries not for this iknot!");
            }

            // get evaluation points and weights
            std::vector<Complex> xs(points), ws(points), fs(points);
            scaled_nodes_and_weights(points, x1, x2, xs.data(), ws.data());

            // evaluate the function
            f(points, xs.data(), fs.data(), data...);

            // sum the results
            Complex result = 0.;
            for (int ipt = 0; ipt < points; ipt++)
                result += ws[ipt] * fs[ipt];

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
            ClassPtr ptr, Functor f, Bspline const & bspline,
            int points, int iknot, Complex x1, Complex x2,
            Data... data
        ) const
        {
            // check boundaries
            if (x1.real() < bspline.t(iknot).real() or bspline.t(iknot+1).real() < x1.real() or
                x2.real() < bspline.t(iknot).real() or bspline.t(iknot+1).real() < x2.real())
            {
                HexException("[quad] Error: boundaries not for this iknot!");
            }

            // get evaluation points and weights
            std::vector<Complex> xs(points), ws(points), fs(points);
            scaled_nodes_and_weights(points, x1, x2, xs.data(), ws.data());

            // evaluate the function
            (ptr->*f)(points, xs.data(), fs.data(), data...);

            // sum the results
            Complex result = 0.;
            for (int ipt = 0; ipt < points; ipt++)
                result += ws[ipt] * fs[ipt];

            return result;
        }
};

#endif
