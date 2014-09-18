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

#include <cmath>
#include <cstdlib>
#include <complex>
#include <iostream>
#include <vector>
#include <tuple>

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
        out[k] = values_i[k] * values_j[k] * pow(in[k]/x,L) * damp(in[k], 0, R);
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
    // raise point count - this is not a poly!
    // TODO Estimate the order.
    int points = bspline_.order() + L + 10; 
    
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
    for (int iknot = 0; iknot < (int)Nreknot - 1; iknot++)
        Rtr_Labcd_diag += computeRdiag(lambda,a,b,c,d,iknot,Nreknot-1);
    
    // Further parts are a bit cryptical, because we are using precomputed
    // (partial, per knot) integral moments, which are quite compactly stored
    // in arrays M_L and M_mLm1 of shape [Nspline × (2*order+1) × (order+1)],
    // but the aim is straightforward: Just to sum the offdiagonal elements,
    // i.e. the products of two two-spline integrals, when ix ≠ iy.
    
    // shorthands
    Complex const * const restrict Mtr_L_ac    = Mtr_L.data()    + (a * (2*order+1) + c - (a - order)) * (order+1);
    Complex const * const restrict Mtr_mLm1_ac = Mtr_mLm1.data() + (a * (2*order+1) + c - (a - order)) * (order+1);
    Complex const * const restrict Mtr_L_bd    = Mtr_L.data()    + (b * (2*order+1) + d - (b - order)) * (order+1);
    Complex const * const restrict Mtr_mLm1_bd = Mtr_mLm1.data() + (b * (2*order+1) + d - (b - order)) * (order+1);
    
    // sum the off-diagonal (iknot_x ≠ iknot_y) contributions for R_tr
    for (int ix = 0; ix < (int)Nreknot - 1; ix++)
    {
        for (int iy = ix + 1; iy < (int)Nreknot - 1; iy++)
        {
            // ix < iy
            if (a <= ix and ix <= a + (int)order and
                b <= iy and iy <= b + (int)order)
            {
                Complex lg = Mtr_L_ac[ix - a] + Mtr_mLm1_bd[iy - b];
                if (std::isfinite(lg.imag()))
                    Rtr_Labcd_offdiag += std::exp(lg);
            }
            
            // ix > iy (by renaming the ix,iy indices)
            if (b <= ix and ix <= b + (int)order and
                a <= iy and iy <= a + (int)order)
            {
                Complex lg = Mtr_L_bd[ix - b] + Mtr_mLm1_ac[iy - a];
                if (std::isfinite(lg.imag()))
                    Rtr_Labcd_offdiag += std::exp(lg);
            }
        }
    }
    
    // sum the diagonal and offdiagonal contributions
    return Rtr_Labcd_diag + Rtr_Labcd_offdiag;
}

void RadialIntegrals::allSymmetries
(
    int i, int j, int k, int l,
    Complex Rijkl_tr,
    NumberArray<long> & R_tr_i,
    NumberArray<long> & R_tr_j,
    NumberArray<Complex> & R_tr_v
) const
{
    // shorthand
    int Nspline = bspline_.Nspline();
    
    {
        // store the integral
        R_tr_i.push_back(i * Nspline + j);
        R_tr_j.push_back(k * Nspline + l);
        R_tr_v.push_back(Rijkl_tr);
    }
    
    if (i != k) // i.e. i < k
    {
        // swap i <-> k (symmetry 1)
        R_tr_i.push_back(k * Nspline + j);
        R_tr_j.push_back(i * Nspline + l);
        R_tr_v.push_back(Rijkl_tr);
    }
    
    if (j != l) // i.e. j < l
    {
        // swap j <-> l (symmetry 2)
        R_tr_i.push_back(i * Nspline + l);
        R_tr_j.push_back(k * Nspline + j);
        R_tr_v.push_back(Rijkl_tr);
    }
    
    if (i != j or k != l) // i.e. i < j or k < l
    {
        // swap i <-> j and k <-> l (symmetry 3)
        R_tr_i.push_back(j * Nspline + i);
        R_tr_j.push_back(l * Nspline + k);
        R_tr_v.push_back(Rijkl_tr);
    }
    
    if (i != k and (i != j or k != l)) // i.e. i < k and (i < j or k < l)
    {
        // swap i <-> k (symmetry 1) and i <-> j and k <-> l (symmetry 3)
        R_tr_i.push_back(l * Nspline + i);
        R_tr_j.push_back(j * Nspline + k);
        R_tr_v.push_back(Rijkl_tr);
    }
    
    if (j != l and (i != j or k != l)) // i.e. j < l and (i < j or k < l)
    {
        // swap j <-> l (symmetry 2) and i <-> j and k <-> l (symmetry 3)
        R_tr_i.push_back(j * Nspline + k);
        R_tr_j.push_back(l * Nspline + i);
        R_tr_v.push_back(Rijkl_tr);
    }
    
    // NOTE there are two more symmetries, (1)+(2) and (1)+(2)+(3),
    // but we don't need those for the construction of a symmetrical
    // matrix, so we will not store them at all
}
