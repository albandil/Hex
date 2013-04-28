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

#include <complex>
#include <vector>

#include "arrays.h"
#include "complex.h"

extern Complex* t;			// knot sequence
extern Complex rotation;	// ECS rotation factor
extern double R0;			// ECS edge
extern double Rmax;			// grid end
extern int Nknot;			// knot count
extern int Nreknot;			// real knot count
extern int Nspline;			// B-spline count
extern int Nintval;			// interval count
extern int order;			// B-spline order

/** Function evaluates a B-spline in one point
 * \param i       Index of the B-spline.
 * \param iknot   Index of the left knot.
 * \param k       B-spline order.
 * \param r       Coordinate (independent variable).
 */
Complex _bspline(int i, int iknot, int k, Complex r);

/** Derivative of a B-spline
 * \param i       Index of the B-spline.
 * \param iknot   Index of the left knot.
 * \param k       B-spline order.
 * \param r       Coordinate (independent variable).
 */
Complex _dspline(int i, int iknot, int k, Complex r);

/** Vectorized interface to the _bspline function.
 * For parameters description see _bspline.
 */
void B(int i, int iknot, int n, const Complex* x, Complex* y);

/** Vectorized interface to the _dspline function.
 * For parameters description see _bspline.
 */
void dB(int i, int iknot, int n, const Complex* x, Complex* y);

/**
 * Apply the ECS transformation.
 * \param r Real coordinate.
 */
inline Complex ECS_rotate(double r)
{
	return (r <= R0) ? r : R0 + (r - R0) * rotation;
};

/**
 * Apply the inverse ECS transformation.
 * \param z Complex coordinate.
 */
inline double ECS_unrotate(Complex z)
{
	return (z.imag() == 0.) ? z.real() : R0 + ((z - R0) / rotation).real();
};

/**
 * Evaluate 1D function given as a B-spline expansion over a grid.
 * \param coeff Expansion coefficients.
 * \param grid Evaluation points (unrotated). They are assumed to be sorted.
 * The length of \c coeff must be at least equal to the spline count and it is
 * these first \c Nspline coefficients that are used in evaluation.
 */
cArray zip(
	cArray const & coeff,
	rArray const & grid
);

/**
 * Evaluate 2D function given as a B-spline expansion over a carthesian product
 * of two 1D grids
 * \param coeff Expansion coefficients.
 * \param xgrid Evaluation points at x-axis (unrotated). They are assumed to be sorted.
 * \param ygrid Evaluation points at y-axis (unrotated). They are assumed to be sorted.
 * The length of \c coeff must be at least equal to the spline count squared and it is
 * these first \c Nspline**2 coefficients that are used in evaluation.
 */
cArray zip(
	cArray const & coeff,
	rArray const & xgrid,
	rArray const & ygrid
);

/**
 * Setup the knot sequence, which will consist of two parts.
 * \param order  B-spline order.
 * \param rknots Real knot array (usually including R₀).
 * \param R0     ECS turning point.
 * \param th     ECS angle in radians.
 * \param cknots To-be-complex knot array (including R₀ and Rmax).
 * \param Rmax   Last knot.
 */
void setup_knot_sequence(int _order_, rArray rknots, double _R0_, double th, rArray cknots, double _Rmax_);

#endif
