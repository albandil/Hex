/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2012                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

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
 * \param r Radius.
 * \param R Truncation radius.
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
	Complex R = t[iknotmax];
	
	// evaluate B-splines
	Complex values_i[n], values_j[n];
	B(i, iknot, n, in, values_i);
	B(j, iknot, n, in, values_j);
	
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
	Complex R = t[iknotmax];	
	
	// evaluate B-splines
	Complex values_i[n], values_j[n];
	B(i, iknot, n, in, values_i);
	B(j, iknot, n, in, values_j);
	
	int points2 = order + L + 1;
	
	// evaluate inner integral, fill output array
	for (int u = 0; u < n; u++)
	{
		auto data2 = std::make_tuple(k,l,L,iknot,iknotmax,in[u]);
		out[u] = values_i[u] * values_j[u] / in[u] * damp(0, in[u], R)
 					* quad(&R_inner_integrand, &data2, points2, iknot, t[iknot], in[u]);
	}
}


/**
 * Triangle integral
 */
Complex computeRtri(int L, int k, int l, int m, int n, int iknot, int iknotmax)
{
	// data for the outer integral
	int data[7] = {k, l, m, n, L, iknot, iknotmax};
	int points = 20; //15; // raise point count - this is not a poly! TODO ???
	
	return quad(&R_outer_integrand, data, points, iknot, t[iknot], t[iknot+1]);
}


/**
 * Diagonal part of Slater integral
 */
Complex computeRdiag(int L, int a, int b, int c, int d, int iknot, int iknotmax)
{
	// throw away if any B-spline identically zero here
	if (iknot < a or a + order < iknot or
		iknot < b or b + order < iknot or
		iknot < c or c + order < iknot or
		iknot < d or d + order < iknot)
		return 0;
	
	// throw away zero length intervals as well
	if (t[iknot] == t[iknot + 1])
		return 0.;
	
	return computeRtri(L,b,d,a,c,iknot,iknotmax) + computeRtri(L,a,c,b,d,iknot,iknotmax);
}


Complex computeR(int lambda,
				 int a, int b, int c, int d,
				 const cArray& Mtr_L, const cArray& Mtr_mLm1
				 
){
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
