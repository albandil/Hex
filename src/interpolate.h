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

#ifndef HEX_INTERPOLATE
#define HEX_INTERPOLATE

#include <algorithm>

#include <o2scl/interp.h>

#include "arrays.h"

/**
 * \brief Return linearly interpolated values.
 * 
 * Returns an array of interpolates of the array y0 for every value of x.
 * \param x0 X-values for the discrete samples.
 * \param y0 Discrete samples
 * \param x  Evaluation (interpolation) points.
 */
template <typename T> NumberArray<T> interpolate (rArray const & x0, NumberArray<T> const & y0, rArray const & x)
{
//     if (x0.size() == 0)
//         throw exception ("Nothing to interpolate.\n");
    
    if (x0.size() < 2)
        return y0;
    
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
 * \brief Return values interpolated by O₂scl.
 * 
 * Returns an array of interpolates of the array y0 for every value of x.
 * \param x0 X-values for the discrete samples.
 * \param y0 Discrete samples
 * \param x  Evaluation (interpolation) points.
 * \param interpolation Interpolation type.
 * \code
   enum {
     // Linear
     itp_linear=0,
     // Cubic spline for natural boundary conditions
     itp_cspline=1,
     // Cubic spline for periodic boundary conditions
     itp_cspline_peri=2,
     // Akima spline for natural boundary conditions
     itp_akima=3,
     // Akima spline for periodic boundary conditions
     itp_akima_peri=4
   };
 * \endcode
 */
inline rArray interpolate_real (rArray const & x0, rArray const & y0, rArray const & x, int interpolation)
{
//     if (x0.size() == 0)
//         throw exception ("Nothing to interpolate.\n");
    
    if (x0.size() < 2)
        return y0;
    
    // setup the interpolator
    o2scl::interp_o2scl_vec<const double*> itp (
        x0.size(),
        x0.data(), y0.data(),
        interpolation
    );
    
    // interpolate
    rArray y(x.size());
    for (size_t i = 0; i < x.size(); i++)
    {
        // check that we are not extrapolating
        if (x0.front() <= x[i] and x[i] <= x0.back())
            y[i] = itp(x[i]);
        else
            y[i] = 0.;
    }
    
    return y;
}

#endif
