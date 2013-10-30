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

#ifndef HEX_GAUSS
#define HEX_GAUSS

#include <map>
#include <vector>

#include "bspline.h"

class GaussLegendre
{
    public:
        
        // constructor
        GaussLegendre(Bspline const & bspline)
            : bspline_(bspline) {}
        
        /**
         * Get pointer to precomputed values of Gauss-Legendre points.
         * @param points Gauss-Legendre points half-count. If too low/high, the return value
         *               will contain the (used) nearest implemented value.
         * @param vx     On return, the Gauss-Legendre nodes (nonnegative half of them).
         * @param vw     On return, the corresponding Gauss-Legendre weights.
         */
        int gauss_nodes_and_weights (int points, const double* &vx, const double* &vw) const;

        /**
         * Get Gauss-Legendre points in complex interval.
         * @param points Number of Gauss-Legendre points.
         * @param x1 Left boundary of the interval.
         * @param x2 Right boundary of the interval.
         */
        cArray p_points (int points, Complex x1, Complex x2) const;

        /**
         * Get Gauss-Legendre weights in complex interval.
         * @param points Number of Gauss-Legendtre points.
         * @param x1 Left boundary of the interval.
         * @param x2 Right boundary of the interval.
         */
        cArray p_weights (int points, Complex x1, Complex x2) const;

        /** Gauss-Legendre integrator
         * 
         * Integrate function.
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
        template <class Functor, class... Data> Complex quad (
            Functor f,
            int points, int iknot, Complex x1, Complex x2,
            Data... data
        ) const {
            // check boundaries
            if (x1.real() < bspline_.t(iknot).real() or bspline_.t(iknot+1).real() < x1.real() or
                x2.real() < bspline_.t(iknot).real() or bspline_.t(iknot+1).real() < x2.real())
            {
                throw exception ("[quad] Error: boundaries not for this iknot!");
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
        template <class ClassPtr, class Functor, class... Data> Complex quadMFP (
            ClassPtr ptr, Functor f, int points, int iknot, Complex x1, Complex x2, Data... data
        ) const {
            // check boundaries
            if (x1.real() < bspline_.t(iknot).real() or bspline_.t(iknot+1).real() < x1.real() or
                x2.real() < bspline_.t(iknot).real() or bspline_.t(iknot+1).real() < x2.real())
            {
                throw exception ("[quad] Error: boundaries not for this iknot!");
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
        
        // precomputed nodes and weights, common to all instances
        static std::vector<std::pair<double*,double*>> data_;
};

#endif
