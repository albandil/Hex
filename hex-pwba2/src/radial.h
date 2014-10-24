//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2014, Jakub Benda, Charles University in Prague                    //
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

#ifndef HEX_PWBA2_RADIAL_H
#define HEX_PWBA2_RADIAL_H

#include <iostream>

#include "arrays.h"

/**
 * @brief Calculate bound-bound intermediate multipole potential.
 * 
 * Evaluate the bound-bound multipole intermediate potential at grid 'x'.
 * The potential is defined in the following way:
 * @f[
 *     V_{ab}^\lambda(r) = \int_0^\infty P_a(r')
 *     \left(
 *         \frac{r_<^\lambda}{r_>^{\lambda+1}}
 *          - \frac{\delta_{\lambda 0}}{r}
 *     \right) P_b(r') \mathrm{d}r' \,,
 * @f]
 * where @f$ P(r) = rR(r) @f$ is the bound hydrogen radial orbital.
 * 
 * The potential is evaluated exactly as the result is sum of gamma-functions.
 * 
 * @param x Evaluation grid.
 * @param lambda Multipole (= angular momentum transfer).
 * @param Na Left state principal quantum number.
 * @param La Left state orbital quantum number.
 * @param Nb Right state principal quantum number.
 * @param Lb Right state orbital quantum number.
 */
rArray interpolate_bound_bound_potential
(
    rArray const & x,
    int lambda,
    int Na, int La, int Nb, int Lb
);

/**
 * @brief Calculate bound-bound intermediate multipole potential.
 * 
 * Evaluate the bound-bound multipole intermediate potential at grid 'x'.
 * The potential is defined in the following way:
 * @f[
 *     V_{ab}^\lambda(r) = \int_0^\infty P_a(r')
 *     \left(
 *         \frac{r_<^\lambda}{r_>^{\lambda+1}}
 *          - \frac{\delta_{\lambda 0}}{r}
 *     \right) F_b(r') \mathrm{d}r' \,,
 * @f]
 * where @f$ P(r) = rR(r) @f$ is the bound hydrogen radial orbital
 * and @f$ F(r) @f$ is the radial part of the continuum hydrogen function
 * defined in @ref Hydrogen::F.
 * 
 * The straightforward evaluation of this potential is a time consuming process
 * due to difficulties with evaluation of the Coulomb functions (see @ref Hydrogen::F).
 * To speed up the calculation, the product @f$ P_a(r') F_b(r') @f$ is evaluated
 * on the grid 'x' and replaced by a cubic spline approximation. The subsequent
 * adaptive integration uses that spline instead of exact evaluation of the product.
 * 
 * @param x Evaluation grid.
 * @param lambda Multipole (= angular momentum transfer).
 * @param Na Left state principal quantum number.
 * @param La Left state orbital quantum number.
 * @param Nb Right state principal quantum number.
 * @param Lb Right state orbital quantum number.
 */
rArray interpolate_bound_free_potential
(
    rArray const & x,
    int lambda,
    int Na, int La, double Kb, int Lb
);

rArray interpolate_free_bound_potential
(
    rArray const & x,
    int lambda,
    double Ka, int La, int Nb, int Lb
);

rArray interpolate_riccati_bessel_j
(
    rArray const & x,
    int l, double k
);

rArray interpolate_riccati_bessel_y
(
    rArray const & x,
    int l, double k
);

rArray interpolate_riccati_bessel_iscaled
(
    rArray const & x,
    int l, double k
);

rArray interpolate_riccati_bessel_kscaled
(
    rArray const & x,
    int l, double k
);

Complex Idir_allowed
(
    rArray const & grid,
    rArray const & jf, rArray const & Vfn,
    rArray const & jn, rArray const & yn_ji, rArray const & yn_jf,
    rArray const & ji, rArray const & Vni
);

double Idir_forbidden
(
    rArray const & grid,
    rArray const & jf, rArray const & Vfn,
    rArray const & iscaled_n, rArray const & kscaled_n,
    rArray const & ji, rArray const & Vni
);

Complex Idir_nBound_allowed
(
    rArray const & grid, int L,
    int Nf, int Lf, double kf, int lf,
    int Nn, int Ln, double kn, int ln,
    int Ni, int Li, double ki, int li,
    std::ostream & log = std::cout
);

double Idir_nBound_forbidden
(
    rArray const & grid, int L,
    int Nf, int Lf, double kf, int lf,
    int Nn, int Ln, double kappan, int ln,
    int Ni, int Li, double ki, int li,
    std::ostream & log = std::cout
);

Complex Idir_nFree_allowed
(
    rArray const & grid, int L,
    int Nf, int Lf, double kf, int lf,
    double Kn, int Ln, double kn, int ln,
    int Ni, int Li, double ki, int li,
    std::ostream & log = std::cout
);

double Idir_nFree_forbidden
(
    rArray const & grid, int L,
    int Nf, int Lf, double kf, int lf,
    double Kn, int Ln, double kappan, int ln,
    int Ni, int Li, double ki, int li,
    std::ostream & log = std::cout
);

#endif
