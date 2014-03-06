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

#ifndef HEX_PWBA
#define HEX_PWBA

#include "arrays.h"

/**
 * @brief Plane-wave Born approximation.
 * 
 * This function computes all contributions to the T-matrix for the specified
 * total angular momentum L in the Born approximation of the first-order.
 * The contributions depend on all such "li" and "lf" that
 * the following conditions must be satisfied:
 *    |li - Li| <= L <= li + Li
 *    |lf - Lf| <= L <= lf + Lf
 * This results in the following bounds on "li" and "lf":
 *    |Li - L| <= li <= Li + L
 *    |Lf - L| <= lf <= Lf + L
 * 
 * @param Ni Initial atomic state.
 * @param Li Initial atomic state.
 * @param Nf Final atomic state.
 * @param Lf Final atomic state.
 * @param ki Initial projectile momentum.
 * @param kf Final projectile momentum.
 * @param L Total angular momentum.
 * @param Tdir List of direct scattering T-matrices for all allowed "lf".
 * @param Texc List of exchange scattering T-matrices for all allowed "lf".
 * @param direct Whether to compute direct scattering contributions.
 * @param exchange Whether to compute exhange scattering contributions.
 */
void pwba
(
    int Ni, int Li, double ki,
    int Nf, int Lf, double kf,
    int L,
    cArrays & Tdir, cArrays & Texc,
    bool direct = true, bool exchange = true
);

#endif /* HEX_PWBA */
