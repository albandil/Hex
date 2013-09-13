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

#ifndef HEX_SLATER_H
#define HEX_SLATER_H

#include <complex>

#include "arrays.h"

/**
 * Compute the two-electron (Slater-type) four-B-spline integral.
 * 
 * @param lambda Multipole degree.
 * @param i First (x-dependent) B-spline index.
 * @param j Second (y-dependent) B-spline index.
 * @param k Third (x-dependent) B-spline index.
 * @param l Fourth (y-dependent) B-spline index.
 * @param Mtr_L As above, for R₀-truncated moments.
 * @param Mtr_mLm1 As above, for R₀-truncated moments.
 * 
 * Given the R-type integral symmetry, following calls will produce identical results:
 * 
 * @code
 * computeR(lambda, i, j, k, l);
 * computeR(lambda, j, i, l, k);
 * computeR(lambda, k, j, i, l);
 * computeR(lambda, i, l, k, j);
 * computeR(lambda, k, l, i, j);
 * @endcode
 * 
 * @return Value of Rtr.
 */
Complex computeR (
    int lambda,
    int i, int j, int k, int l,
    cArray const & Mtr_L,
    cArray const & Mtr_mLm1
);

void allSymmetries (
    int i, int j, int k, int l,
    Complex Rijkl_tr,
    NumberArray<long> & R_tr_i,
    NumberArray<long> & R_tr_j,
    NumberArray<Complex> & R_tr_v
);

#endif
