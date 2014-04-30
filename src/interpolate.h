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

#ifndef HEX_INTERPOLATE
#define HEX_INTERPOLATE

#include <algorithm>

#include <gsl/gsl_interp.h>

#include "arrays.h"

/**
 * @brief Return linearly interpolated values.
 * 
 * Returns an array of interpolates of the array y0 for every value of x.
 * @param x0 X-values for the discrete samples.
 * @param y0 Discrete samples
 * @param x  Evaluation (interpolation) points.
 */
template <typename T> NumberArray<T> interpolate (rArray const & x0, NumberArray<T> const & y0, rArray const & x)
{
    assert(x0.size() == y0.size());
    
    // return 0 if the array to interpolate is empty
    if (x0.size() == 0)
        return NumberArray<T>(x.size(), T(0));
    
    // return the source element if the array to interpolate contains just the one element
    if (x0.size() == 1)
        return NumberArray<T>(x.size(), y0[0]);
    
    // output array
    NumberArray<T> y(x.size());
    
    for (size_t i = 0; i < x.size(); i++)
    {
        // right guardian
        rArray::const_iterator right = std::upper_bound(x0.begin(), x0.end(), x[i]);
        
        if (right == x0.end() or right == x0.begin())
        {
            y[i] = T(0);
            continue;
        }
        
        // neighbours
        double x0_left = *(right - 1);
        double x0_right = *right;
        T y0_left = y0[right - 1 - x0.begin()];
        T y0_right = y0[right - x0.begin()];
        
        // compute linear average
        y[i] = (y0_left * (x0_right - x[i]) + y0_right * (x[i] - x0_left)) / (x0_right - x0_left);
    }
    
    // FIXME implement better schemes
    
    return y;
}

/**
 * @brief Return values interpolated by O₂scl.
 * 
 * Returns an array of interpolates of the array y0 for every value of x.
 * @param x0 X-values for the discrete samples.
 * @param y0 Discrete samples
 * @param x  Evaluation (interpolation) points.
 * @param interpolation Interpolation type.
 *  - gsl_interp_linear : Linear interpolation. This interpolation method
 *    does not require any additional memory. 
 *  - gsl_interp_polynomial : Polynomial interpolation. This method should only
 *    be used for interpolating small numbers of points because polynomial
 *    interpolation introduces large oscillations, even for well-behaved datasets.
 *    The number of terms in the interpolating polynomial is equal to the number of points. 
 *  - gsl_interp_cspline : Cubic spline with natural boundary conditions.
 *    The resulting curve is piecewise cubic on each interval, with matching
 *    first and second derivatives at the supplied data-points. The second
 *    derivative is chosen to be zero at the first point and last point. 
 *  - gsl_interp_cspline_periodic : Cubic spline with periodic boundary 
 *    conditions. The resulting curve is piecewise cubic on each interval, 
 *    with matching first and second derivatives at the supplied data-points. 
 *    The derivatives at the first and last points are also matched. Note 
 *    that the last point in the data must have the same y-value as the 
 *    first point, otherwise the resulting periodic interpolation will have 
 *    a discontinuity at the boundary.
 *  - gsl_interp_akima : Non-rounded Akima spline with natural boundary 
 *    conditions. This method uses the non-rounded corner algorithm of Wodicka. 
 *  - gsl_interp_akima_periodic : Non-rounded Akima spline with periodic 
 *    boundary conditions. This method uses the non-rounded corner algorithm of Wodicka.
 */
inline rArray interpolate_real (rArray const & x0, rArray const & y0, rArray const & x, const gsl_interp_type * interpolation)
{
    assert(x0.size() == y0.size());
    
    // return 0 if the array to interpolate is empty
    if (x0.size() == 0)
        return rArray(x.size(), 0.);
    
    // return the source element if the array to interpolate contains just the one element
    if (x0.size() == 1)
        return rArray(x.size(), y0[0]);
    
    // setup the interpolator
    gsl_interp *itp = gsl_interp_alloc (interpolation, x0.size());
    gsl_interp_init (itp, x0.data(), y0.data(), x0.size());
    gsl_interp_accel *accel = gsl_interp_accel_alloc ();
    
    // interpolate
    rArray y(x.size());
    for (size_t i = 0; i < x.size(); i++)
    {
        // check that we are not extrapolating
        if (x0.front() <= x[i] and x[i] <= x0.back())
            y[i] = gsl_interp_eval (itp, x0.data(), y0.data(), x[i], accel);
        else
            y[i] = 0.;
    }
    
    // release memory
    gsl_interp_accel_free (accel);
    gsl_interp_free (itp);
    
    return y;
}

#endif
