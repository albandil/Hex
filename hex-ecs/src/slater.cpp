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

#include <gsl/gsl_sf.h>

#include "arrays.h"
#include "bspline.h"
#include "gauss.h"
#include "radial.h"
#include "special.h"

void RadialIntegrals::R_inner_integrand
(
    int n, Complex* in, Complex* out,
    Bspline const & bspline_ij, int i, int j,
    int L, int iknot, int iknotmax, Complex x
) const
{
    Complex R = bspline_ij.t(iknotmax);
    
    // evaluate B-splines
    cArray values_i(n), values_j(n);
    bspline_ij.B(i, iknot, n, in, values_i.data());
    bspline_ij.B(j, iknot, n, in, values_j.data());
    
    // fill output array
    for (int k = 0; k < n; k++)
        out[k] = values_i[k] * values_j[k] * special::pow_int<Complex>(in[k]/x,L) * damp(in[k], 0, R);
}

void RadialIntegrals::R_outer_integrand
(
    int n, Complex* in, Complex* out,
    Bspline const & bspline_ij, int i, int j,
    Bspline const & bspline_kl, GaussLegendre const & g_kl, int k, int l,
    int L, int iknot, int iknotmax
) const
{
    // extract data
    Complex R = bspline_ij.t(iknotmax);
    
    // evaluate B-splines
    cArray values_i(n), values_j(n);
    bspline_ij.B(i, iknot, n, in, values_i.data());
    bspline_ij.B(j, iknot, n, in, values_j.data());
    
    // use at least 2nd order
    int points2 = std::max(2, bspline_kl.order() + L + 1);
    
    // evaluate inner integral, fill output array
    for (int u = 0; u < n; u++)
    {
        out[u] = values_i[u] * values_j[u] / in[u] * damp(0., in[u], R) * g_kl.quadMFP
        (
            this, &RadialIntegrals::R_inner_integrand,      // integrand pointers
            points2, iknot, bspline_kl.t(iknot), in[u],     // integrator parameters
            bspline_kl, k, l, L, iknot, iknotmax, in[u]     // integrand data
        );
    }
}

Complex RadialIntegrals::computeRtri
(
    int L,
    Bspline const & bspline_kl, GaussLegendre const & g_kl, int k, int l,
    Bspline const & bspline_mn, GaussLegendre const & g_mn, int m, int n,
    int iknot, int iknotmax
) const
{
    // compute integral of Bk(1) Bl(1) V(1,2) Bm(2) Bn(2)
    
    // integration points (the integrand is a poly of order equal to the four times the order of B-splines)
    int points = 2 * bspline_kl.order() + 1; 
    
    // integrate
    return g_kl.quadMFP
    (
        this, &RadialIntegrals::R_outer_integrand,                      // integrand pointers
        points, iknot, bspline_kl.t(iknot), bspline_kl.t(iknot+1),      // integrator parameters
        bspline_kl, k, l, bspline_mn, g_mn, m, n, L, iknot, iknotmax    // integrand data
    );
}

Complex RadialIntegrals::computeRdiag (int L, int a, int b, int c, int d, int iknot, int iknotmax) const
{
    int order = bspline_atom_.order();
    
    // throw away if any B-spline identically zero here
    if (iknot < a or a + order < iknot or
        iknot < b or b + order < iknot or
        iknot < c or c + order < iknot or
        iknot < d or d + order < iknot)
        return 0;
    
    // throw away zero length intervals as well
    if (bspline_atom_.t(iknot) == bspline_atom_.t(iknot + 1))
        return 0.;
    
    // sum the two triangle integrals
    return computeRtri(L,bspline_atom_,g_atom_,b,d,bspline_atom_,g_atom_,a,c,iknot,iknotmax)
         + computeRtri(L,bspline_atom_,g_atom_,a,c,bspline_atom_,g_atom_,b,d,iknot,iknotmax);
}

cArray RadialIntegrals::diagonalR (int lambda) const
{
    // assume bspline_atom == bspline_proj
    int order = bspline_atom_.order();
    int Nspline = bspline_atom_.Nspline();
    int Nreknot = bspline_atom_.Nreknot();
    
    // allocate space
    std::size_t O = order + 1;
    cArray R (Nspline * O * O * O, 0.);
    
    // calculate elements
    # pragma omp parallel for
    for (int a = 0; a < Nspline; a++)
    for (int b = a; b <= a + order; b++)
    for (int c = a; c <= a + order; c++)
    for (int d = a; d <= a + order; d++)
    {
        for (int iknot = mmax(a,b,c,d); iknot <= mmin(a,b,c,d) + order and iknot < Nreknot - 1; iknot++)
            R[((a * O + (b-a)) * O + (c-a)) * O + (d-a)] += computeRdiag(lambda, a, b, c, d, iknot, Nreknot - 1);
    }
    
    return R;
}

Complex RadialIntegrals::computeR
(
    int lambda,
    int a, int b, int c, int d,
    bool simple
) const
{
    int order = bspline_atom_.order();
    int Nspline_atom = bspline_atom_.Nspline();
    int Nreknot_atom = bspline_atom_.Nreknot();
    int Nreknot_proj = bspline_proj_.Nreknot();
    
    double ta1 = bspline_atom_.t(a).real(), ta2 = bspline_atom_.t(a + order + 1).real();
    double tb1 = bspline_proj_.t(b).real(), tb2 = bspline_proj_.t(b + order + 1).real();
    double tc1 = bspline_atom_.t(c).real(), tc2 = bspline_atom_.t(c + order + 1).real();
    double td1 = bspline_proj_.t(d).real(), td2 = bspline_proj_.t(d + order + 1).real();
    
    // check overlaps
    if (ta2 <= tc1 or tc2 <= ta1 or tb2 <= td1 or td2 <= tb1)
        return 0.;
    
    // Are the integral moments completely decoupled, i.e. there is there no overlap between Ba, Bb, Bc and Bd?
    // In such cases we can compute the off-diagonal contribution just as a product of the two
    // integral moments of order "lambda" and "-lambda-1", respectively. Moreover, in such
    // case there is no diagonal contribution, because there is no overlap between the _four_
    // participating B-splines.
    
    // (b,d) << (a,c)
    if (std::min(tb2,td2) <= std::max(ta1,tc1))
    {
        double t_ac = std::min(ta2,tc2);
        double t_bd = std::min(tb2,td2);
        double scale = gsl_sf_pow_int(t_bd / t_ac, lambda) / t_ac;
        return scale * Mtr_mLm1_atom_[lambda](a,c) * Mtr_L_proj_[lambda](b,d);
    }
    
    // (a,c) << (b,d)   [this is used by non-first panels]
    if (std::min(ta2,tc2) <= std::max(tb1,td1))
    {
        double t_ac = std::min(ta2,tc2);
        double t_bd = std::min(tb2,td2);
        double scale = gsl_sf_pow_int(t_ac / t_bd, lambda) / t_bd;
        return scale * Mtr_L_atom_[lambda](a,c) * Mtr_mLm1_proj_[lambda](b,d);
    }
    
    // The rest allows overlap of the four B-splines and is used only by the first
    // (origin) panel.
    
    // diagonal part
    Complex Rtr_Labcd_diag = 0;
    
    // off-diagonal part
    Complex Rtr_Labcd_offdiag = 0;
    
    // sum the diagonal (iknot_x = iknot_y = iknot) contributions
    if (not simple)
    {
        // work in atomic basis index space -> translate indices
        int Ta = a;
        int Tb = b + pshift_;
        int Tc = c;
        int Td = d + pshift_;
        
        // check that all indices are valid and that there is a four-B-spline overlap
        if (Tb < Nspline_atom and Td < Nspline_atom and mmax(Ta,Tb,Tc,Td) + order < mmin(Ta,Tb,Tc,Td))
        {
            // auxiliary parameter
            std::size_t O = order + 1;
            
            // retrieve diagonal contribution
            if (Ta <= Tb and Ta <= Tc and Ta <= Td)
                Rtr_Labcd_diag = R_tr_dia_diag_[lambda][((Ta * O + (Tb-Ta)) * O + (Tc-Ta)) * O + (Td-Ta)];
            else if (Tb <= Ta and Tb <= Tc and Tb <= Td)
                Rtr_Labcd_diag = R_tr_dia_diag_[lambda][((Tb * O + (Ta-Tb)) * O + (Td-Tb)) * O + (Tc-Tb)];
            else if (Tc <= Ta and Tc <= Tb and Tc <= Td)
                Rtr_Labcd_diag = R_tr_dia_diag_[lambda][((Tc * O + (Tb-Tc)) * O + (Ta-Tc)) * O + (Td-Tc)];
            else // (Td <= Ta and Td <= Tb and Td <= Tc)
                Rtr_Labcd_diag = R_tr_dia_diag_[lambda][((Td * O + (Ta-Td)) * O + (Tb-Td)) * O + (Tc-Td)];
        }
    }
    
    // Further parts are a bit cryptical, because we are using precomputed
    // (partial, per knot) integral moments, which are quite compactly stored
    // in arrays M_L and M_mLm1 of shape [Nspline × (2*order+1) × (order+1)],
    // but the aim is straightforward: Just to sum the offdiagonal elements,
    // i.e. the products of two two-spline integrals, when ix ≠ iy.
    
    // shorthands
    Complex const * const restrict Mitr_L_ac    = Mitr_L_atom(lambda).data()    + (a * (2*order+1) + c - (a-order)) * (order+1);
    Complex const * const restrict Mitr_mLm1_ac = Mitr_mLm1_atom(lambda).data() + (a * (2*order+1) + c - (a-order)) * (order+1);
    Complex const * const restrict Mitr_L_bd    = Mitr_L_proj(lambda).data()    + (b * (2*order+1) + d - (b-order)) * (order+1);
    Complex const * const restrict Mitr_mLm1_bd = Mitr_mLm1_proj(lambda).data() + (b * (2*order+1) + d - (b-order)) * (order+1);
    
    // sum the off-diagonal (iknot_x ≠ iknot_y) contributions for R_tr
    
    // ix < iy
    for (int ix = a;                   ix < std::min(a + order + 1, Nreknot_atom - 1); ix++) if (bspline_atom_.t(ix+1).real() > 0)
    for (int iy = std::max(b, ix + 1); iy < std::min(b + order + 1, Nreknot_proj - 1); iy++)
    {
        // calculate scale factor
        double tx = bspline_atom_.t(ix+1).real();
        double ty = bspline_proj_.t(iy+1).real();
        double scale = gsl_sf_pow_int(tx / ty, lambda) / ty;
        
        // calculate contribution to the integral
        Rtr_Labcd_offdiag += Mitr_L_ac[ix-a] * Mitr_mLm1_bd[iy-b] * scale;
    }
    
    // ix > iy (by swapping (a,c) and (b,d) multi-indices)
    for (int ix = b;                   ix < std::min(b + order + 1, Nreknot_proj - 1); ix++) if (bspline_proj_.t(ix+1).real() > 0)
    for (int iy = std::max(a, ix + 1); iy < std::min(a + order + 1, Nreknot_atom - 1); iy++)
    {
        // calculate scale factor
        double tx = bspline_proj_.t(ix+1).real();
        double ty = bspline_atom_.t(iy+1).real();
        double scale = gsl_sf_pow_int(tx / ty, lambda) / ty;
        
        // calculate contribution to the integral
        Rtr_Labcd_offdiag += Mitr_L_bd[ix-b] * Mitr_mLm1_ac[iy-a] * scale;
    }
    
    // sum the diagonal and offdiagonal contributions
    return Rtr_Labcd_diag + Rtr_Labcd_offdiag;
}
