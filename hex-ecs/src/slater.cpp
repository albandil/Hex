//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2017, Jakub Benda, Charles University in Prague                    //
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

// --------------------------------------------------------------------------------- //

#include <gsl/gsl_sf.h>

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"
#include "hex-special.h"

// --------------------------------------------------------------------------------- //

#include "bspline.h"
#include "gauss.h"
#include "radial.h"

// --------------------------------------------------------------------------------- //

void RadialIntegrals::R_inner_integrand
(
    int n, Complex* in, Complex* out,
    int i, int j,
    int L, int iknot, int iknotmax, Complex x
) const
{
    //Complex R = bspline_full_y_.t(iknotmax);
    
    // evaluate B-splines
    cArray values_i(n), values_j(n);
    bspline_y_.B(i, iknot, n, in, values_i.data());
    bspline_y_.B(j, iknot, n, in, values_j.data());
    
    // fill output array
    for (int k = 0; k < n; k++)
        out[k] = values_i[k] * values_j[k] * special::pow_int<Complex>(in[k]/x,L) /* * damp(in[k], R) */;
}

void RadialIntegrals::R_outer_integrand
(
    int n, Complex* in, Complex* out,
    int i, int j,
    int k, int l,
    int L, int iknot, int iknotmax
) const
{
    // extract data
    //Complex R = bspline_full_x_.t(iknotmax);
    
    // evaluate B-splines
    cArray values_i(n), values_j(n);
    bspline_x_.B(i, iknot, n, in, values_i.data());
    bspline_x_.B(j, iknot, n, in, values_j.data());
    
    // use at least 2nd order
    int points2 = std::max(2, bspline_x_.order() + L + 1);
    
    // evaluate inner integral, fill output array
    for (int u = 0; u < n; u++)
    {
        out[u] = values_i[u] * values_j[u] / in[u] /** damp(in[u], R)*/ * g_y_.quadMFP
        (
            this, &RadialIntegrals::R_inner_integrand,      // integrand pointers
            bspline_y_, points2, iknot, bspline_y_.t(iknot), in[u],     // integrator parameters
            k, l, L, iknot, iknotmax, in[u]     // integrand data
        );
    }
}

Complex RadialIntegrals::computeRtri
(
    int L,
    int k, int l,
    int m, int n,
    int iknot, int iknotmax
) const
{
    // compute integral of Bk(1) Bl(1) V(1,2) Bm(2) Bn(2)
    
    // integration points (the integrand is a poly of order equal to the four times the order of B-splines)
    int points = 2 * bspline_x_.order() + 1; 
    
    // integrate
    return g_x_.quadMFP
    (
        this, &RadialIntegrals::R_outer_integrand,                      // integrand pointers
        bspline_x_, points, iknot, bspline_x_.t(iknot), bspline_x_.t(iknot+1),      // integrator parameters
        k, l, m, n, L, iknot, iknotmax    // integrand data
    );
}

Complex RadialIntegrals::computeRdiag (int L, int a, int b, int c, int d, int iknot, int iknotmax) const
{
    // shorthands
    int order = bspline_x_.order();
    
    // throw away if any B-spline identically zero here
    if (iknot < a or a + order < iknot or
        iknot < b or b + order < iknot or
        iknot < c or c + order < iknot or
        iknot < d or d + order < iknot)
        return 0.;
    
    // throw away zero length intervals as well
    if (bspline_x_.t(iknot).real() == bspline_x_.t(iknot + 1).real())
        return 0.;
    
    // sum the two triangle integrals
    return computeRtri(L,b,d,a,c,iknot,iknotmax)
         + computeRtri(L,a,c,b,d,iknot,iknotmax);
}

cArray RadialIntegrals::diagonalR (int lambda) const
{
    // assume bspline_atom == bspline_proj
    int order = bspline_x_.order();
    int Nspline = bspline_x_.Nspline();
    int Nreknot = bspline_x_.Nreknot();
    
    // allocate space
    std::size_t O = order + 1;
    cArray R (Nspline * O * O * O, 0.);
    
    // calculate elements
    if (bspline_x_.hash() == bspline_y_.hash())
    {
        # pragma omp parallel for
        for (int a = 0; a < Nspline; a++)
        for (int b = a; b <= a + order; b++)
        for (int c = a; c <= a + order; c++)
        for (int d = a; d <= a + order; d++)
        {
            for (int iknot = mmax(a,b,c,d); iknot <= mmin(a,b,c,d) + order and iknot < Nreknot - 1; iknot++)
                R[((a * O + (b-a)) * O + (c-a)) * O + (d-a)] += computeRdiag(lambda, a, b, c, d, iknot, Nreknot - 1);
        }
    }
    
    return R;
}

Complex RadialIntegrals::computeR (int lambda, int a, int b, int c, int d) const
{
    //
    // Preparations.
    //
    
        // shorthands
        int order = bspline_x_.order();
        int iR1x = bspline_x_.iR1(), iR2x = bspline_x_.iR2();
        int iR1y = bspline_y_.iR1(), iR2y = bspline_y_.iR2();
        
        // leading and trailing knots of the B-splines
        Real ta1 = bspline_x_.t(a).real(), ta2 = bspline_x_.t(a + order + 1).real();
        Real tb1 = bspline_y_.t(b).real(), tb2 = bspline_y_.t(b + order + 1).real();
        Real tc1 = bspline_x_.t(c).real(), tc2 = bspline_x_.t(c + order + 1).real();
        Real td1 = bspline_y_.t(d).real(), td2 = bspline_y_.t(d + order + 1).real();
    
    //
    // Return clean zero if there are no pair overlaps.
    //
    
        if (std::max(ta1,tc1) > std::min(ta2,tc2) or std::max(tb1,td1) > std::min(tb2,td2))
            return 0.;
    
    // If the overlaps of the same-coordinate B-splines do not overlap
    // with each other, the four-B-spline integral decouples into a product
    // of two two-B-spline integrals. These have been already calculated.
    
        // the overlap of Bb,Bd precedes the overlap of Ba,Bc on the radial axis
        if (std::min(tb2,td2) <= std::max(ta1,tc1))
        {
            Real t_ac = std::min(ta2,tc2);
            Real t_bd = std::min(tb2,td2);
            Real scale = gsl_sf_pow_int(t_bd / t_ac, lambda) / t_ac;
            return scale * Mtr_mLm1_x_[lambda](a,c) * Mtr_L_y_[lambda](b,d);
        }
        
        // the overlap of Ba,Bc precedes the overlap of Bb,Bd on the radial axis
        if (std::min(ta2,tc2) <= std::max(tb1,td1))
        {
            Real t_ac = std::min(ta2,tc2);
            Real t_bd = std::min(tb2,td2);
            Real scale = gsl_sf_pow_int(t_ac / t_bd, lambda) / t_bd;
            return scale * Mtr_L_x_[lambda](a,c) * Mtr_mLm1_y_[lambda](b,d);
        }
    
    // The four-B-spline integrals that do not decouple into a product of two one-dimensional integrals need to be carefully calculated 
    // on a Cartesian product of the knot sequences, which transforms the full integral into a sum of integrals from several two-dimensional
    // cells. The integral does decouple, again, in such two-dimensional cells, where either r1 > r2 or r1 < r2 everywhere. This is called
    // the "off-diagonal" contribution to the R-integral. On the contrary, the integration cells where the notion of larger and smaller
    // coordinate changes from place to place contribute to the "diagonal" part. The diagonal part has been precomputed and is stored
    // in the array R_tr_dia_diag_ for all real B-splines. It is assumed that the complex B-splines lie so far from the origin
    // that the integral can be approximated by the arithmetic average of its two possible decoupled forms.
    
        // handle overlap of complex B-splines in much the same way as a decoupled integral above
        if
        (
            std::min(a,c)             < iR1x or std::min(b,d)             < iR1y or
            std::max(a,c) + order + 1 > iR2x or std::max(b,d) + order + 1 > iR2y
        )
        {
            Real t_ac = std::min(ta2,tc2);
            Real t_bd = std::min(tb2,td2);
            
            Real scale1 = gsl_sf_pow_int(t_bd / t_ac, lambda) / t_ac;
            Real scale2 = gsl_sf_pow_int(t_ac / t_bd, lambda) / t_bd;
            
            // return average of the two possible decoupled forms
            return 0.5 *
            (
                scale1 * Mtr_mLm1_x_[lambda](a,c) * Mtr_L_y_[lambda](b,d) +
                scale2 * Mtr_L_x_[lambda](a,c) * Mtr_mLm1_y_[lambda](b,d)
            );
        }
        
        // diagonal part (only between real B-splines on conforming knot sequences)
        Complex Rtr_Labcd_diag = 0;
        
        // locate the diagonal contribution in R_tr_dia_diag_
        std::size_t idx;
        if (a <= b and a <= c and a <= d)
            idx = ((a * (order + 1) + (b - a)) * (order + 1) + (c - a)) * (order + 1) + (d - a);
        else if (b <= a and b <= c and b <= d)
            idx = ((b * (order + 1) + (a - b)) * (order + 1) + (d - b)) * (order + 1) + (c - b);
        else if (c <= a and c <= b and c <= d)
            idx = ((c * (order + 1) + (b - c)) * (order + 1) + (a - c)) * (order + 1) + (d - c);
        else // (d <= a and d <= b and d <= c)
            idx = ((d * (order + 1) + (a - d)) * (order + 1) + (b - d)) * (order + 1) + (c - d);
        
        // retrieve the diagonal contribution (if available)
        if (idx < R_tr_dia_diag_[lambda].size())
            Rtr_Labcd_diag = R_tr_dia_diag_[lambda][idx];
    
    // Further parts are a bit cryptical, because we are using precomputed
    // (partial, per knot) integral moments, which are quite compactly stored
    // in arrays M_L and M_mLm1 of shape [Nspline × (2*order+1) × (order+1)],
    // but the aim is straightforward: Just to sum the offdiagonal elements,
    // i.e. the products of two two-spline integrals, when ix ≠ iy.
    
        // off-diagonal part
        Complex Rtr_Labcd_offdiag = 0;
        
        // shorthands
        Complex const * const restrict Mitr_L_ac    = Mitr_L_x(lambda).data()    + (a * (2*order+1) + c - (a-order)) * (order+1);
        Complex const * const restrict Mitr_mLm1_ac = Mitr_mLm1_x(lambda).data() + (a * (2*order+1) + c - (a-order)) * (order+1);
        Complex const * const restrict Mitr_L_bd    = Mitr_L_y(lambda).data()    + (b * (2*order+1) + d - (b-order)) * (order+1);
        Complex const * const restrict Mitr_mLm1_bd = Mitr_mLm1_y(lambda).data() + (b * (2*order+1) + d - (b-order)) * (order+1);
        
        // sum the off-diagonal (iknot_x ≠ iknot_y) contributions for R_tr
        
        // ta[ix] < tp[iy]
        for (int ix = a;                   ix < std::min(a + order + 1, iR2x); ix++) if (bspline_x_.t(ix+1).real() != 0.)
        for (int iy = std::max(b, ix + 1); iy < std::min(b + order + 1, iR2y); iy++)
        {
            // calculate scale factor
            Real tx = bspline_x_.t(ix+1).real();
            Real ty = bspline_y_.t(iy+1).real();
            Real scale = gsl_sf_pow_int(tx / ty, lambda) / ty;
            
            // calculate contribution to the integral
            Rtr_Labcd_offdiag += Mitr_L_ac[ix-a] * Mitr_mLm1_bd[iy-b] * scale;
        }
        
        // ta[ix] > tp[iy] (by swapping (a,c) and (b,d) multi-indices)
        for (int ix = b;                   ix < std::min(b + order + 1, iR2x); ix++) if (bspline_x_.t(ix+1).real() != 0.)
        for (int iy = std::max(a, ix + 1); iy < std::min(a + order + 1, iR2y); iy++)
        {
            // calculate scale factor
            Real tx = bspline_x_.t(ix+1).real();
            Real ty = bspline_y_.t(iy+1).real();
            Real scale = gsl_sf_pow_int(tx / ty, lambda) / ty;
            
            // calculate contribution to the integral
            Rtr_Labcd_offdiag += Mitr_L_bd[ix-b] * Mitr_mLm1_ac[iy-a] * scale;
        }
    
    // sum the diagonal and offdiagonal contributions
    return Rtr_Labcd_diag + Rtr_Labcd_offdiag;
}
