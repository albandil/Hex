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

#ifndef HEX_PWBA2_RADIAL_H
#define HEX_PWBA2_RADIAL_H

#include "arrays.h"

rArray interpolate_bound_bound_potential
(
    rArray const & x,
    int lambda,
    int Na, int La, int Nb, int Lb
);

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
    rArray const & jn, rArray const & yn,
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
    int Ni, int Li, double ki, int li
);

Complex Idir_nBound_forbidden
(
    rArray const & grid, int L,
    int Nf, int Lf, double kf, int lf,
    int Nn, int Ln, double kappan, int ln,
    int Ni, int Li, double ki, int li
);

Complex Idir_nFree_allowed
(
    rArray const & grid, int L,
    int Nf, int Lf, double kf, int lf,
    double Kn, int Ln, double kn, int ln,
    int Ni, int Li, double ki, int li
);

Complex Idir_nFree_forbidden
(
    rArray const & grid, int L,
    int Nf, int Lf, double kf, int lf,
    double Kn, int Ln, double kappan, int ln,
    int Ni, int Li, double ki, int li
);

#endif
