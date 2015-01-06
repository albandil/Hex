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

#include <cmath>
#include <cstdlib>
#include <complex>
#include <iostream>
#include <vector>
#include <tuple>

#ifndef NO_OPENCL
#include "opencl.h"
#endif

#include "arrays.h"
#include "bspline.h"
#include "gauss.h"
#include "radial.h"

/**
 * Potential suppressing factor. 
 * @param y Radial coordinate of some electron.
 * @param x Radial coordinate of the other electron.
 * @param R Truncation radius.
 */
inline double damp (Complex y, Complex x, Complex R)
{
    // compute hyperradius
    double r = std::hypot(x.real(), y.real());
    
    // if sufficiently far, return clean zero
    if (r > R.real())
        return 0.;
    
    // else damp using tanh(x) distribution
    return std::tanh(0.125 * (R.real() - r));
}

void RadialIntegrals::R_inner_integrand (int n, Complex* in, Complex* out, int i, int j, int L, int iknot, int iknotmax, Complex x) const
{
    Complex R = bspline_.t(iknotmax);
    
    // evaluate B-splines
    Complex values_i[n], values_j[n];
    bspline_.B(i, iknot, n, in, values_i);
    bspline_.B(j, iknot, n, in, values_j);
    
    // fill output array
    for (int k = 0; k < n; k++)
        out[k] = values_i[k] * values_j[k] * std::pow(in[k]/x,L) * damp(in[k], 0, R);
}

void RadialIntegrals::R_outer_integrand (int n, Complex* in, Complex* out, int i, int j, int k, int l, int L, int iknot, int iknotmax) const
{
    // extract data
    Complex R = bspline_.t(iknotmax);
    
    // evaluate B-splines
    Complex values_i[n], values_j[n];
    bspline_.B(i, iknot, n, in, values_i);
    bspline_.B(j, iknot, n, in, values_j);
    
    // use at least 2nd order
    int points2 = std::max(2, bspline_.order() + L + 1);
    
    // evaluate inner integral, fill output array
    for (int u = 0; u < n; u++)
    {
        out[u] = values_i[u] * values_j[u] / in[u] * damp(0., in[u], R) * g_.quadMFP
        (
            this, &RadialIntegrals::R_inner_integrand,                    // integrand pointers
            points2, iknot, bspline_.t(iknot), in[u],                     // integrator parameters
            k, l, L, iknot, iknotmax, in[u]                               // integrand data
        );
    }
}

Complex RadialIntegrals::computeRtri (int L, int k, int l, int m, int n, int iknot, int iknotmax) const
{
    // integration points (the integrand is a poly of order equal to the four times the order of B-splines)
    int points = 2 * bspline_.order() + 1; 
    
    // integrate
    return g_.quadMFP
    (
        this, &RadialIntegrals::R_outer_integrand,                // integrand pointers
        points, iknot, bspline_.t(iknot), bspline_.t(iknot+1),    // integrator parameters
        k, l, m, n, L, iknot, iknotmax                            // integrand data
    );
}

Complex RadialIntegrals::computeRdiag (int L, int a, int b, int c, int d, int iknot, int iknotmax) const
{
    int order = bspline_.order();
    
    // throw away if any B-spline identically zero here
    if (iknot < a or a + order < iknot or
        iknot < b or b + order < iknot or
        iknot < c or c + order < iknot or
        iknot < d or d + order < iknot)
        return 0;
    
    // throw away zero length intervals as well
    if (bspline_.t(iknot) == bspline_.t(iknot + 1))
        return 0.;
    
    return computeRtri(L,b,d,a,c,iknot,iknotmax) + computeRtri(L,a,c,b,d,iknot,iknotmax);
}

Complex RadialIntegrals::computeR
(
    int lambda,
    int a, int b, int c, int d,
    cArray const & Mtr_L, cArray const & Mtr_mLm1
) const
{
    int order = bspline_.order();
    int Nreknot = bspline_.Nreknot();
    
    // check overlaps
    if (std::abs(a - c) > order or std::abs(b - d) > order)
        return 0.;
    
    // diagonal part
    Complex Rtr_Labcd_diag = 0;
    
    // off-diagonal part
    Complex Rtr_Labcd_offdiag = 0;
    
    // sum the diagonal (iknot_x = iknot_y = iknot) contributions
    for (int iknot = mmax(a,b,c,d); iknot <= mmin(a,b,c,d) + order and iknot < Nreknot - 1; iknot++)
        Rtr_Labcd_diag += computeRdiag(lambda,a,b,c,d,iknot,Nreknot-1);
    
    // Further parts are a bit cryptical, because we are using precomputed
    // (partial, per knot) integral moments, which are quite compactly stored
    // in arrays M_L and M_mLm1 of shape [Nspline × (2*order+1) × (order+1)],
    // but the aim is straightforward: Just to sum the offdiagonal elements,
    // i.e. the products of two two-spline integrals, when ix ≠ iy.
    
    // shorthands
    Complex const * const restrict Mtr_L_ac    = Mtr_L.data()    + (a * (2*order+1) + c - (a-order)) * (order+1);
    Complex const * const restrict Mtr_mLm1_ac = Mtr_mLm1.data() + (a * (2*order+1) + c - (a-order)) * (order+1);
    Complex const * const restrict Mtr_L_bd    = Mtr_L.data()    + (b * (2*order+1) + d - (b-order)) * (order+1);
    Complex const * const restrict Mtr_mLm1_bd = Mtr_mLm1.data() + (b * (2*order+1) + d - (b-order)) * (order+1);
    
    // sum the off-diagonal (iknot_x ≠ iknot_y) contributions for R_tr
    
    // ix < iy
    for (int ix = a; ix <= a + order and ix < Nreknot - 1; ix++)
    for (int iy = std::max(b,ix+1); iy <= b + order and iy < Nreknot - 1; iy++)
    {
        Complex lg = Mtr_L_ac[ix - a] + Mtr_mLm1_bd[iy - b];
        if (std::isfinite(lg.imag()))
            Rtr_Labcd_offdiag += std::exp(lg);
    }
    
    // ix > iy (by renaming the ix,iy indices)
    for (int ix = b; ix <= b + order and ix < Nreknot - 1; ix++)
    for (int iy = std::max(a,ix+1); iy <= a + order and iy < Nreknot - 1; iy++)
    {
        Complex lg = Mtr_L_bd[ix - b] + Mtr_mLm1_ac[iy - a];
        if (std::isfinite(lg.imag()))
            Rtr_Labcd_offdiag += std::exp(lg);
    }
    
    // sum the diagonal and offdiagonal contributions
    return Rtr_Labcd_diag + Rtr_Labcd_offdiag;
}

void RadialIntegrals::allSymmetries
(
    int i, int j, int k, int l,
    Complex Rijkl_tr,
    SymDiaMatrix & R_tr_dia
) const
{
    // shorthand
    int Nspline = bspline_.Nspline();
    
    {
        // store the integral
        R_tr_dia(i * Nspline + j, k * Nspline + l) = Rijkl_tr;
    }
    
    if (i != k) // i.e. i < k
    {
        // swap i <-> k (symmetry 1)
        R_tr_dia(k * Nspline + j, i * Nspline + l) = Rijkl_tr;
    }
    
    if (j != l) // i.e. j < l
    {
        // swap j <-> l (symmetry 2)
        R_tr_dia(i * Nspline + l, k * Nspline + j) = Rijkl_tr;
    }
    
    if (i != j or k != l) // i.e. i < j or k < l
    {
        // swap i <-> j and k <-> l (symmetry 3)
        R_tr_dia(j * Nspline + i, l * Nspline + k) = Rijkl_tr;
    }
    
    if (i != k and (i != j or k != l)) // i.e. i < k and (i < j or k < l)
    {
        // swap i <-> k (symmetry 1) and i <-> j and k <-> l (symmetry 3)
        R_tr_dia(l * Nspline + i, j * Nspline + k) = Rijkl_tr;
    }
    
    if (j != l and (i != j or k != l)) // i.e. j < l and (i < j or k < l)
    {
        // swap j <-> l (symmetry 2) and i <-> j and k <-> l (symmetry 3)
        R_tr_dia(j * Nspline + k, l * Nspline + i) = Rijkl_tr;
    }
    
    // NOTE : There are two more symmetries, (1)+(2) and (1)+(2)+(3),
    // but we don't need those for the construction of a symmetrical
    // matrix.
}
