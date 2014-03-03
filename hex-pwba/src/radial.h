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

#ifndef HEX_PWBA_RADIAL
#define HEX_PWBA_RADIAL

double compute_Idir
(
    int li, int lf, int lambda,
    int Ni, int Li, double ki,
    int Nf, int Lf, double kf
);

double compute_Iexc
(
    int li, int lf, int lambda,
    int Ni, int Li, double ki,
    int Nf, int Lf, double kf
);

#endif /* HEX_PWBA_RADIAL */
