//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2015, Jakub Benda, Charles University in Prague                    //
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

namespace PWBA1
{

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

}; /* namespace PWBA1 */

#endif /* HEX_PWBA */
