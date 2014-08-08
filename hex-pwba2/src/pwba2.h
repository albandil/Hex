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

#ifndef HEX_PWBA2_H
#define HEX_PWBA2_H

#include "arrays.h"

namespace PWBA2
{
    
    cArrays PartialWave_direct
    (
        rArray grid,
        int L, int Pi,
        int Ni, int Li, double ki,
        int Nf, int Lf, double kf,
        int nL, int maxNn, double Enmax,
        bool integrate_allowed, bool integrate_forbidden
    );
    
    cArrays FullTMatrix_direct
    (
        rArray grid,
        int Ni, int Li, double ki,
        int Nf, int Lf, double kf,
        int maxNn, int maxLn, double maxEn,
        bool integrate_allowed, bool integrate_forbidden
    );
    
}; /* namespace PWBA2 */

#endif
