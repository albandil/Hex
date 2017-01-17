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

#ifndef HEX_BSPLINE
#define HEX_BSPLINE

#include "hex-arrays.h"
#include "hex-misc.h"

/**
 * @brief B-spline environment.
 * 
 * This class is used to manage B-spline computations within the context
 * of the exterior complex scaling. That is, it allows (upon construction)
 * specification of the real and complex knot sequences
 * @f[
 *     t_1^{R}, t_2^{R}, ..., t_{\mathrm{Nreknot}}^{R} \ ,
 * @f]
 * @f[
 *     t_{\mathrm{Nreknot}+1}^{C}, ..., t_{\mathrm{Nknot}}^{C} \ ,
 * @f]
 * which then serve as the knot data of the B-spline set. See @ref bspline
 * and @ref dspline on the definition formulae.
 */
class Bspline
{
    public:
    
        /**
         * @brief Constructor.
         * 
         * Setup the knot sequence, which will consist of two parts.
         * @param order  B-spline order.
         * @param rknots Real knot array (usually including R₀).
         * @param th     ECS angle in radians.
         * @param cknots To-be-complex knot array (including R₀ and Rmax).
         */
        Bspline (int order, rArrayView const & rknots, Real th, rArrayView const & cknots);
        
        // move constructor
        Bspline (Bspline && other);
        
        // destructor
        ~Bspline ();
        
        /**
         * @brief Evaluate B-spline.
         * 
         * Function evaluates a B-spline in one point. The recursive Cox de Boor formula is
         * used:
         * @f[
         *     B_i^{(0)}(x) = \cases {
         *              1 & t_i \le x \le t_{i+1} \cr
         *              0 & \mbox{otherwise}
         *     }\ ,
         * @f]
         * @f[
         *     B_i^{(k)}(x) = \frac{x - t_i}{t_{i+k-1} - t_i} B_i^{(k-1)}(x)
         *                  + \frac{t_{i+k} - x}{t_{i+k} - t_{i+1}} B_{i+1}^{(k-1)}(x) \ .
         * @f]
         * @param i       Index of the B-spline.
         * @param iknot   Index of the left knot.
         * @param k       B-spline order.
         * @param r       Coordinate (independent variable).
         */
        Complex bspline (int i, int iknot, int k, Complex r) const;
        
        /**
         * @brief Evaluate derivative of a B-spline.
         * 
         * This function evaluated the derivative of a B-spline, using the formula
         * 
         * 
         * 
         * @param i       Index of the B-spline.
         * @param iknot   Index of the left knot.
         * @param k       B-spline order.
         * @param r       Coordinate (independent variable).
         */
        Complex dspline (int i, int iknot, int k, Complex r) const;
        
        /** 
         * @brief B-spline.
         * 
         * Evaluates single B-spline at several points at one. This results in
         * increased speed, because some data need to be computed only once.
         * Also, some operation have been reorganized to a form that allows
         * SIMD auto-vectorization.
         * 
         * Note that all evaluation points must be in the interval bounded by
         * 'iknot'-th and 'iknot+1'-th knot.
         * 
         * For parameters description see @ref bspline.
         */
        void B (int i, int iknot, int n, const Complex* x, Complex* y) const;
        
        /**
         * @brief Derivative of a B-spline.
         * 
         * For parameters description see @ref dspline.
         * 
         * This routine simply calls @ref dspline for every evaluation point.
         */
        void dB (int i, int iknot, int n, const Complex* x, Complex* y) const;
        
        /**
         * @brief Apply the ECS transformation.
         * 
         * @param r Real coordinate.
         */
        inline Complex rotate (Real r) const
        {
            return (r <= R0_) ? r : R0_ + (r - R0_) * rotation_;
        };
        
        /**
         * @brief Apply the inverse ECS transformation.
         * 
         * @param z Complex coordinate.
         */
        inline Real unrotate (Complex z) const
        {
            return (z.imag() == 0.) ? z.real() : R0_ + ((z - R0_) / rotation_).real();
        };
        
        /**
         * @brief Zip 1D expansion.
         * 
         * Evaluate 1D function given as a B-spline expansion over a grid.
         * @param coeff Expansion coefficients.
         * @param grid Evaluation points (unrotated). They are assumed to be sorted.
         * The length of \c coeff must be at least equal to the spline count and it is
         * these first \c Nspline coefficients that are used in evaluation.
         */
        cArray zip (const cArrayView coeff, const rArrayView grid) const;
        
        /**
         * @brief Zip 2D expansion.
         * 
         * Evaluate 2D function given as a B-spline expansion over a carthesian product
         * of two 1D grids.
         * @param coeff Expansion coefficients.
         * @param xgrid Evaluation points at x-axis (unrotated). They are assumed to be sorted.
         * @param ygrid Evaluation points at y-axis (unrotated). They are assumed to be sorted.
         * The length of \c coeff must be at least equal to the spline count squared and it is
         * these first \c Nspline**2 coefficients that are used in evaluation.
         */
        cArray zip (const cArrayView coeff, const rArrayView xgrid, const rArrayView ygrid) const;
        static cArray zip (Bspline const & bx, Bspline const & by, const cArrayView coeff, const rArrayView xgrid, const rArrayView ygrid);
        
        /**
         * @brief Get knot index for coordinate.
         * 
         * Finds knot (interval @f$ \left< t_i,t_{i+1} \right> @f$) for 'x'.
         * @param x Complex coordinate.
         */
        int knot (Complex x) const;
        
        /**
         * @brief Evaluate 1D B-spline expansion.
         * 
         * Evaluates a B-spline expansion in a single point x.
         * This function is faster than calling
         * @code
         *     Complex z = zip (coeff, rArray({x}));
         * @endcode
         * though the results are the same.
         */
        Complex eval (const cArrayView coeff, Real x) const;
        
        /**
         * @brief Evaluate 2D B-spline expansion.
         * 
         * Evaluates a double B-spline expansion in a single point (x,y).
         * This function is faster than calling
         * @code
         *     Complex z = zip (coeff, rArray({x}), rArray({y}));
         * @endcode
         * though the results are the same.
         */
        Complex eval (const cArrayView coeff, Real x, Real y) const;
        
        // getters
        
        /// B-spline knot sequence.
        inline cArrayView t () const { return cArrayView(Nknot_, t_); }
        
        /// B-spline knot sequence.
        inline Complex const & t (int i) const { return *(t_ + i); }
        
        /// Number of B-splines.
        inline int Nspline () const { return Nspline_; }
        
        /// Number of knots.
        inline int Nknot () const { return Nknot_; }
        
        /// Number of real knots.
        inline int Nreknot () const { return Nreknot_; }
        
        /// B-spline order.
        inline int order () const { return order_; }
        
        /// End of real grid.
        inline Real R0 () const { return R0_; };
        
        /// End of complex grid (real, unrotated).
        inline Real Rmax () const { return Rmax_; };
        
        /// ECS rotation angle.
        inline Real ECStheta () const { return theta_; }
        
        /// real knots
        inline rArray const & rknots () const { return rknots_; }
        
        /// complex knots
        inline rArray const & cknots () const { return cknots_; }
        
        /**
         * @brief Return (almost) unique identification for this B-spline object.
         */
        std::size_t hash () const;
        
    private:
    
        /// real knots
        rArray rknots_;
        
        /// complex knots
        rArray cknots_;
        
        /// knot sequence
        Complex * restrict t_;
        
        /// rotation angle (rad)
        Real theta_;
        
        /// ECS rotation factor
        Complex rotation_;
        
        /// end of real grid
        Real R0_;
        
        /// end of complex grid
        Real Rmax_;
        
        /// knot count
        int Nknot_;
        
        /// real knot count
        int Nreknot_;
        
        /// B-spline count
        int Nspline_;
        
        /// inter-knot interval count
        int Nintval_;
        
        /// B-spline order
        int order_;
        
        /// Work arrays (one per thread)
        mutable rArrays work_;
        
        /// Work max size
        static const std::size_t work_size_;
        
        /// Disallow implicit bitwise copy.
        Bspline (Bspline const &) = delete;
        Bspline const & operator= (Bspline const &) = delete;
};

#endif
