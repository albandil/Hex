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

#ifndef HEX_ECS_BSPLINE
#define HEX_ECS_BSPLINE

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"
#include "hex-misc.h"

// --------------------------------------------------------------------------------- //

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
         * @param order   B-spline order.
         * @param th      ECS angle in radians.
         * @param cknots1 Leading complex-rotated knots.
         * @param rknots  Real knot array.
         * @param cknots2 Trailing complex-rotated knots.
         */
        Bspline
        (
            int order,
            Real th,
            rArrayView cknots1,
            rArrayView rknots,
            rArrayView cknots2
        );
        
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
        Complex rotate (Real r) const
        {
            // leading complex part
            if (r < R1_)
                return R1_ + (r - R1_) * rotation_;
            
            // trailing complex part
            if (r > R2_)
                return R2_ + (r - R2_) * rotation_;
            
            // real part
            return r;
        };
        
        /**
         * @brief Apply the inverse ECS transformation.
         * 
         * @param z Complex coordinate.
         */
        Real unrotate (Complex z) const
        {
            // leading complex part
            if (z.imag() < 0)
                return R1_ + (z.real() - R1_) / rotation_.real();
            
            // trailing complex part
            if (z.imag() > 0)
                return R2_ + (z.real() - R2_) / rotation_.real();
            
            // real part
            return z.real();
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
        
        /**
         * @nrief Zip 2D expansion.
         * 
         * Evaluate 2D function given as a B-spline expansion over a carthesian product
         * of two 1D grids. The user may supply the B-spline evaluation functions for the
         * individual axes, thus allowing evaluation of partial derivatives besides the
         * plain B-spline expansion.
         */
        static cArray zip
        (
            Bspline const & bx,
            Bspline const & by,
            const cArrayView coeff,
            const rArrayView xgrid,
            const rArrayView ygrid,
            Complex (Bspline::* evalXBSpline) (int,int,int,Complex) const = &Bspline::bspline,
            Complex (Bspline::* evalYBSpline) (int,int,int,Complex) const = &Bspline::bspline
        );
        
        /**
         * @brief Get knot index for coordinate.
         * 
         * Finds knot (interval @f$ \langle t_i,t_{i+1} \rangle @f$) for 'x'.
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
        cArrayView t () const { return t_; }
        
        /// B-spline knot sequence.
        Complex const & t (int i) const { return t_[i]; }
        
        /// Number of B-splines.
        int Nspline () const { return Nspline_; }
        
        /// Number of knots.
        int Nknot () const { return Nknot_; }
        
        /// Number of real knots.
        int Nreknot () const { return Nreknot_; }
        
        /// B-spline order.
        int order () const { return order_; }
        
        /**
         * @brief Index of knot between leading complex grid and real grid.
         * 
         * This is also the index of the first possibly purely real B-spline.
         */
        int iR1 () const { return cknots1_.empty() ? 0 : cknots1_.size() - 1; }
        
        /**
         * @brief Index of knot between real grid and trailing complex grid.
         * 
         * This is the index of the trailing knot of the last possibly purely real B-spline.
         * The index of that B-spline is iR2 - order - 1.
         */
        int iR2 () const { return cknots2_.empty() ? t_.size() - 1 : t_.size() - cknots2_.size(); }
        
        /// Beginning of the real grid.
        Real R1 () const { return R1_; }
        
        /// End of the real grid.
        Real R2 () const { return R2_; }
        
        /// Beginning of the grid (unrotated).
        Real Rmin () const { return Rmin_; }
        
        /// End of the grid (unrotated).
        Real Rmax () const { return Rmax_; };
        
        /// ECS rotation angle.
        Real ECStheta () const { return theta_; }
        
        /// real knots
        rArray const & rknots () const { return rknots_; }
        
        /// complex knots
        rArray const & cknots1 () const { return cknots1_; }
        rArray const & cknots2 () const { return cknots2_; }
        
        /**
         * @brief Return (almost) unique identification for this B-spline object.
         * 
         * The identifier is calculated as a simple hash of the data contained
         * in the object.
         */
        std::size_t hash () const;
        
    private:
    
        /// real knots
        rArray rknots_;
        
        /// complex knots
        rArray cknots1_, cknots2_;
        
        /// knot sequence (rotated)
        cArray t_;
        
        /// rotation angle (rad)
        Real theta_;
        
        /// ECS rotation factor
        Complex rotation_;
        
        /// beginning of the real grid
        Real R1_;
        
        /// end of the real grid
        Real R2_;
        
        /// beginning of the full grid
        Real Rmin_;
        
        /// end of the full grid
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
};

// --------------------------------------------------------------------------------- //

#endif // HEX_ECS_BSPLINE
