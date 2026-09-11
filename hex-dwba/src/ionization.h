//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2026, Jakub Benda, Charles University in Prague                    //
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

#ifndef HEX_DWBA_IONIZATION
#define HEX_DWBA_IONIZATION

#include "hex-arrays.h"

/**
 * @brief Partial ionization cross sections in the plane-wave Born approximation.
 *
 * Computes the contribution of every total angular momentum @f$ L @f$ to the cross
 * section of
 * @f[
 *     e(k_i) + \mathrm{H}(1s) \rightarrow e(k_f) + p + e(\kappa) \,,
 * @f]
 * in the first Born approximation with plane waves and without exchange. The result
 * is integrated over the energy of the ejected electron and summed over its angular
 * momentum, so that the sum over @f$ L @f$ reproduces the total cross section written
 * by @c hex-fullborn @c --ionization.
 *
 * The projectile is undistorted, so the initial partial wave is pinned to
 * @f$ \ell_i = L @f$ by the sphericity of the initial state, and for the same reason a
 * given angular momentum @f$ \ell_a @f$ of the ejected electron is fed by the single
 * multipole @f$ \lambda = \ell_a @f$. What is left is the sum over the partial wave
 * @f$ \ell_f @f$ of the scattered projectile, which runs over the triangle of
 * @f$ (\ell_a, \ell_f, L) @f$ with @f$ L + \ell_a + \ell_f @f$ even.
 *
 * @param Ei Impact energy (Ry).
 * @param Lmin First total angular momentum.
 * @param Lmax Last total angular momentum.
 * @param sigma Output cross sections (a_0^2), one entry per L from Lmin to Lmax.
 * @param lamax Output: the highest angular momentum of the ejected electron that was
 *              needed to reach the requested accuracy.
 * @param eps Relative accuracy at which the sum over the ejected angular momentum is
 *            truncated.
 */
void pwba_ionization
(
    double Ei,
    int Lmin, int Lmax,
    rArray & sigma,
    int & lamax,
    double eps = 1e-6
);

#endif
