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

#ifndef HEX_BSPLINE
#define HEX_BSPLINE

#include "arrays.h"
#include "complex.h"

/**
 * @brief B-spline environment.
 * 
 * This class is used to manage B-spline computations within the context
 * of the exterior complex scaling. That is, it allows (upon construction)
 * specification of the real and complex knot sequences, which then serve
 * as the knot data of the B-spline set.
 */
class Bspline
{
    public:
    
        /**
        * Setup the knot sequence, which will consist of two parts.
        * @param order  B-spline order.
        * @param rknots Real knot array (usually including R₀).
        * @param th     ECS angle in radians.
        * @param cknots To-be-complex knot array (including R₀ and Rmax).
        */
        Bspline (int order, rArrayView const & rknots, double th, rArrayView const & cknots);
        
        /**
         * Function evaluates a B-spline in one point
         * @param i       Index of the B-spline.
         * @param iknot   Index of the left knot.
         * @param k       B-spline order.
         * @param r       Coordinate (independent variable).
         */
        Complex bspline(int i, int iknot, int k, Complex r) const;
        
        /** Derivative of a B-spline
         * @param i       Index of the B-spline.
         * @param iknot   Index of the left knot.
         * @param k       B-spline order.
         * @param r       Coordinate (independent variable).
         */
        Complex dspline(int i, int iknot, int k, Complex r) const;
        
        /** 
         * Vectorized interface to the _bspline function.
         * For parameters description see _bspline.
         */
        void B(int i, int iknot, int n, const Complex* x, Complex* y) const;
        
        /** 
         * Vectorized interface to the _dspline function.
         * For parameters description see _bspline.
         */
        void dB(int i, int iknot, int n, const Complex* x, Complex* y) const;
        
        /**
         * Apply the ECS transformation.
         * @param r Real coordinate.
         */
        inline Complex rotate (double r) const
        {
            return (r <= R0_) ? r : R0_ + (r - R0_) * rotation_;
        };
        
        /**
         * Apply the inverse ECS transformation.
         * @param z Complex coordinate.
         */
        inline double unrotate (Complex z) const
        {
            return (z.imag() == 0.) ? z.real() : R0_ + ((z - R0_) / rotation_).real();
        };
        
        /**
         * Evaluate 1D function given as a B-spline expansion over a grid.
         * @param coeff Expansion coefficients.
         * @param grid Evaluation points (unrotated). They are assumed to be sorted.
         * The length of \c coeff must be at least equal to the spline count and it is
         * these first \c Nspline coefficients that are used in evaluation.
         */
        cArray zip (cArrayView const & coeff, rArrayView const & grid) const;
        
        /**
         * Evaluate 2D function given as a B-spline expansion over a carthesian product
         * of two 1D grids
         * @param coeff Expansion coefficients.
         * @param xgrid Evaluation points at x-axis (unrotated). They are assumed to be sorted.
         * @param ygrid Evaluation points at y-axis (unrotated). They are assumed to be sorted.
         * The length of \c coeff must be at least equal to the spline count squared and it is
         * these first \c Nspline**2 coefficients that are used in evaluation.
         */
        cArray zip (cArrayView const & coeff, rArrayView const & xgrid, rArrayView const & ygrid) const;
        
        /**
         * Finds knot for 'x'.
         * @param x Complex coordinate.
         */
        int knot (Complex x) const;
        
        /**
         * Evaluates a B-spline expansion in a single point x.
         * This function is faster than calling
         * @code
         *     Complex z = zip (coeff, rArray({x}));
         * @endcode
         * though the results are the same.
         */
        Complex eval(cArrayView const & coeff, double x) const;
        
        /**
         * Evaluates a double B-spline expansion in a single point (x,y).
         * This function is faster than calling
         * @code
         *     Complex z = zip (coeff, rArray({x}), rArray({y}));
         * @endcode
         * though the results are the same.
         */
        Complex eval(cArrayView const & coeff, double x, double y) const;
        
        // getters
        
        inline Complex const & t (int i) const { return *(t_ + i); }
        inline int Nspline() const { return Nspline_; }
        inline int Nknot() const { return Nknot_; }
        inline int Nreknot() const { return Nreknot_; }
        inline int order() const { return order_; }
        inline double R0() const { return R0_; };
        inline double Rmax() const { return Rmax_; };
    
    private:
    
        /// knot sequence
        Complex * restrict t_;
        
        /// ECS rotation factor
        Complex rotation_;
        
        /// end of real grid
        double R0_;
        
        /// end of complex grid
        double Rmax_;
        
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
};

#endif
