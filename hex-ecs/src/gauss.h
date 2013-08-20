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

#ifndef HEX_QUAD
#define HEX_QUAD

/**
 * Get pointer to precomputed values of Gauss-Legendre points.
 * \param points Gauss-Legendre points half-count. If too low/high, the return value
 *               will contain the (used) nearest implemented value.
 * \param vx     On return, the Gauss-Legendre nodes (nonnegative half of them).
 * \param vw     On return, the corresponding Gauss-Legendre weights.
 */
int gauss_nodes_and_weights(int points, const double* &vx, const double* &vw);

/**
 * Get Gauss-Legendre points in complex interval.
 * \param points Number of Gauss-Legendre points. If too low/high, it will be
 *               changed to nearest implemented value.
 * \param x1 Left boundary of the interval.
 * \param x2 Right boundary of the interval.
 */
cArray p_points(int& points, Complex x1, Complex x2);

/**
 * Get Gauss-Legendre weights in complex interval.
 * \param points Number of Gauss-Legendtre points. If too low/high, it will be
 *               changed to nearest implemented value.
 * \param x1 Left boundary of the interval.
 * \param x2 Right boundary of the interval.
 */
cArray p_weights(int& points, Complex x1, Complex x2);

/** Gauss-Legendre integrator
 * 
 * Integrate function.
 * 
 * \param f Function of type (void (*) (unsigned, complex*, complex*, void*)),
 *          where the first parameter is number of points to evaluate, the second parameter
 *          is pointer to an array of double complex values at which to evaluate, the third
 *          is pointer to an array into which the data will be written (responsibility for
 *          memory management is on callers side) and finally the fourth parameter is 
 *          any other data to be suplied to the function.
 * \param data Data to pass to the function.
 * \param points Gauss-Legendre points count.
 * \param iknot Knot index.
 * \param x1 Left integration boundary.
 * \param x2 Right integration boundary.
 * 
 * It must be
 * \code
 *    t[iknot] <= x1 <= t[iknot+1]
 *    t[iknot] <= x2 <= t[iknot+1]
 * \endcode
 */
Complex quad(
	void (*f)(int, Complex*, Complex*, void*), void *data,
	int points, int iknot, Complex x1, Complex x2
);

#endif
