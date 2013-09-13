/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2013                              *
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
#include "moments.h"
#include "gauss.h"

/**
 * Potential suppressing factor. 
 * @param y Radial coordinate of some electron.
 * @param x Radial coordinate of the other electron.
 * @param R Truncation radius.
 */
inline double damp(Complex y, Complex x, Complex R)
{
    // compute hyperradius
    double r = hypot(x.real(), y.real());
    
    // if sufficiently far, return clean zero
    if (r > R.real())
        return 0.;
    
    // else damp using tanh(x) distribution
    return tanh(0.125 * (R.real() - r));
}

void R_inner_integrand(int n, Complex* in, Complex* out, void* data)
{
    // extract data
    auto params = *(std::tuple<int,int,int,int,int,Complex>*)data;
    int i = std::get<0>(params);
    int j = std::get<1>(params);
    int L = std::get<2>(params);
    int iknot = std::get<3>(params);
    int iknotmax = std::get<4>(params);
    Complex x = std::get<5>(params);
    Complex R = Bspline::ECS().t(iknotmax);
    
    // evaluate B-splines
    Complex values_i[n], values_j[n];
    Bspline::ECS().B(i, iknot, n, in, values_i);
    Bspline::ECS().B(j, iknot, n, in, values_j);
    
    // fill output array
    for (int k = 0; k < n; k++)
        out[k] = values_i[k] * values_j[k] * pow(in[k]/x,L) * damp(in[k], 0, R);
}


void R_outer_integrand(int n, Complex* in, Complex* out, void* data)
{
    // extract data
    int i = ((int*)data)[0];
    int j = ((int*)data)[1];
    int k = ((int*)data)[2];
    int l = ((int*)data)[3];
    int L = ((int*)data)[4];
    int iknot = ((int*)data)[5];
    int iknotmax = ((int*)data)[6];
    Complex R = Bspline::ECS().t(iknotmax);
    
    // evaluate B-splines
    Complex values_i[n], values_j[n];
    Bspline::ECS().B(i, iknot, n, in, values_i);
    Bspline::ECS().B(j, iknot, n, in, values_j);
    
    int points2 = Bspline::ECS().order() + L + 1;
    
    // evaluate inner integral, fill output array
    for (int u = 0; u < n; u++)
    {
        auto data2 = std::make_tuple(k,l,L,iknot,iknotmax,in[u]);
        out[u] = values_i[u] * values_j[u] / in[u] * damp(0, in[u], R)
             * quad(&R_inner_integrand, &data2, points2, iknot, Bspline::ECS().t(iknot), in[u]);
    }
}


/**
 * Triangle integral
 */
Complex computeRtri(int L, int k, int l, int m, int n, int iknot, int iknotmax)
{
    // data for the outer integral
    int data[7] = {k, l, m, n, L, iknot, iknotmax};
    
    // raise point count - this is not a poly!
    // TODO Estimate the order.
    int points = Bspline::ECS().order() + L + 10; 
    
    // integrate
    return quad(&R_outer_integrand, data, points, iknot, Bspline::ECS().t(iknot), Bspline::ECS().t(iknot+1));
}


/**
 * Diagonal part of Slater integral
 */
Complex computeRdiag(int L, int a, int b, int c, int d, int iknot, int iknotmax)
{
    int order = Bspline::ECS().order();
    
    // throw away if any B-spline identically zero here
    if (iknot < a or a + order < iknot or
        iknot < b or b + order < iknot or
        iknot < c or c + order < iknot or
        iknot < d or d + order < iknot)
        return 0;
    
    // throw away zero length intervals as well
    if (Bspline::ECS().t(iknot) == Bspline::ECS().t(iknot + 1))
        return 0.;
    
    return computeRtri(L,b,d,a,c,iknot,iknotmax) + computeRtri(L,a,c,b,d,iknot,iknotmax);
}


Complex computeR(int lambda,
                 int a, int b, int c, int d,
                 const cArray& Mtr_L, const cArray& Mtr_mLm1
                 
){
    int order = Bspline::ECS().order();
    int Nreknot = Bspline::ECS().Nreknot();
    
    // check overlaps
    if (abs(a-c) > order or abs(b-d) > order)
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
    const Complex* const Mtr_L_ac    = Mtr_L.data()    + (a * (2*order+1) + c - (a - order)) * (order+1);
    const Complex* const Mtr_mLm1_ac = Mtr_mLm1.data() + (a * (2*order+1) + c - (a - order)) * (order+1);
    const Complex* const Mtr_L_bd    = Mtr_L.data()    + (b * (2*order+1) + d - (b - order)) * (order+1);
    const Complex* const Mtr_mLm1_bd = Mtr_mLm1.data() + (b * (2*order+1) + d - (b - order)) * (order+1);
    
    // sum the off-diagonal (iknot_x ≠ iknot_y) contributions for R_tr
    double B_re = 0., B_im = 0.;
    for (int ix = 0; ix < (int)Nreknot - 1; ix++)
    {
        for (int iy = ix + 1; iy < (int)Nreknot - 1; iy++)
        {
            Complex B (0.,0.);
            
            // ix < iy
            if (a <= ix and ix <= a + (int)order and
                b <= iy and iy <= b + (int)order)
                B += Mtr_L_ac[ix - a] * Mtr_mLm1_bd[iy - b];
            
            // ix > iy (by renaming the ix,iy indices)
            if (b <= ix and ix <= b + (int)order and
                a <= iy and iy <= a + (int)order)
                B += Mtr_L_bd[ix - b] * Mtr_mLm1_ac[iy - a];
            
            // store full integral contribution
            B_re += B.real();
            B_im += B.imag();
        }
    }
    
    // store offdiagonal R-type integral contributions
    Rtr_Labcd_offdiag = Complex(B_re, B_im);

    // sum the diagonal and offdiagonal contributions
    return Rtr_Labcd_diag + Rtr_Labcd_offdiag;
}

void allSymmetries (
    int i, int j, int k, int l,
    Complex Rijkl_tr,
    NumberArray<long> & R_tr_i,
    NumberArray<long> & R_tr_j,
    NumberArray<Complex> & R_tr_v
){
    // shorthand
    int Nspline = Bspline::ECS().Nspline();
    
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
