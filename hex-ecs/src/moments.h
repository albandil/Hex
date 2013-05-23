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

#ifndef _MOMTS_H_
#define _MOMTS_H_

#include <vector>
#include <complex>

#include "bspline.h"

/**
 * Compute derivative overlap of B-splines \f$ B_i \f$ and \f$ B_j \f$
 * over the knot "iknot", using Gauss-Legendre integration.
 * \param i      B-spline index.
 * \param j      B-spline index.
 * \param iknot  Interval index.
 */
Complex computeD_iknot(int i, int j, int iknot);

/**
 * Compute derivative overlap for B-splines \f$ B_i \f$ and \f$ B_j \f$.
 * \param i B-spline index.
 * \param j B-spline index.
 * \param maxknot Right-most knot of any integration.
 */
Complex computeD(int i, int j, int maxknot = Bspline::ECS().Nknot() - 1);

/**
 * Compute integral moment of coordinate power between the B-splines
 * \f$ B_i \f$ and \f$ B_j \f$
 * over the knot "iknot", using Gauss-Legendre integration.
 * \param a      Exponent.
 * \param i      B-spline index.
 * \param j      B-spline index.
 * \param iknot  Interval index.
 * \param R      Potential truncation point.
 */
Complex computeM_iknot(int a, int i, int j, int iknot, Complex R);

/**
 * Compute integral moment of coordinate power between the B-splines
 * \f$ B_i \f$ and \f$ B_j \f$
 * \param a Exponent.
 * \param i B-spline index.
 * \param j B-spline index.
 * \param maxknot Right-most knot of any integration.
 */
Complex computeM(int a, int i, int j, int maxknot = 0);

/**
 * Compute integral moment of degree "a" for every B-spline pair and every
 * interknot sub-interval. Store in 1-D array of shape
 * \code
 * [ Nspline × Nspline × Nintval ]
 * \endcode
 * Zero entries are stores as well, to allow for better caching.
 * 
 * \param a Moment degree.
 */
cArray computeMi(int a, int iknotmax = 0);

#endif
