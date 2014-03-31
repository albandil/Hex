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

#ifndef HEX_PWBA2_RADIAL
#define HEX_PWBA2_RADIAL

Complex Idir_nBound
(
    int lambdaf, int lambdai,
    int Nf, int Lf, int lf, double kf,
    int Nn, int Ln, int ln, double kn,
    int Ni, int Li, int li, double ki
);

Complex Idir_nFree
(
    int lambdaf, int lambdai,
    int Nf, int Lf, int lf, double kf,
    double Kn, int Ln, int ln, double kn,
    int Ni, int Li, int li, double ki
);

#endif
