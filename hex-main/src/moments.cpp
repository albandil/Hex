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

#include <algorithm>
#include <complex>

#include "arrays.h"
#include "bspline.h"
#include "gauss.h"

/**
 * Potential suppressing factor. 
 * \param r Radius.
 * \param R Truncation radius.
 */
inline double damp(Complex r, Complex R)
{
	// if sufficiently far, return clean zero
	if (r.real() > R.real())
		return 0.;
	
	// else damp using tanh(x) distribution
	return tanh(0.125 * (R.real() - r.real()));
}

// ----------------------------------------------------------------------- //
//  Computation of B-spline moments                                        //
// ----------------------------------------------------------------------- //

void M_integrand(int n, Complex *in, Complex *out, void *data)
{
	// extract data
	int i = ((int*)data)[0];
	int j = ((int*)data)[1];
	int a = ((int*)data)[2];
	int iknot = ((int*)data)[3];
	int iknotmax = ((int*)data)[4];
	Complex R = t[iknotmax];
	
	// evaluate B-splines
	Complex values_i[n], values_j[n];
	B(i, iknot, n, in, values_i);
	B(j, iknot, n, in, values_j);
	
	// fill output array
	int k;
	if (R != 0.)
	{
		for (k = 0; k < n; k++)
			out[k] = values_i[k] * values_j[k] * pow(in[k],a) * damp(in[k],R);
	}
	else
	{
		for (k = 0; k < n; k++)
			out[k] = values_i[k] * values_j[k] * pow(in[k],a);
	}
}

/** 
 * Compute two-spline integral moment of degree "a"
 * \param a Moment degree.
 * \param iknotmax Truncation knot.
 */
cArray computeMi(int a, int iknotmax)
{
	int i, j, iknot;
	size_t size = Nspline * (2 * order + 1) * (order + 1);
	cArray m(size);
	
	// for all B-splines
	for (i = 0; i < (int)Nspline; i++)
	{
		// for all B-splines (use symmetry)
		for (j = i; j <= i + (int)order and j < (int)Nspline; j++)
		{
			// determine relevant knots
			int ileft = (i < j) ? j : i;
			int iright = ((i < j) ? i : j) + order + 1;
			
			// "right" has to be smaller than "Rmax"
			if (iright > iknotmax)
				iright = iknotmax;
			
			// for all relevant knots
			for (iknot = ileft; iknot < iright; iknot++)
			{
				// get integration boundaries
				Complex xa = t[iknot];
				Complex xb = t[iknot+1];
				
				// throw away zero length intervals
				if (xa == xb)
					continue;
				
				// compute the moment
				int points = 20;//order + abs(a) + 1;
				int data[5] = {i, j, a, iknot, iknotmax};
				Complex integral = quad(&M_integrand, data, points, iknot, xa, xb);
				
				// get the coordinates in m-matrix
				int x_1 = i;				// reference spline is i-th
				int y_1 = j - (i - order);
				int z_1 = iknot - i;
				
				// get the coordinates in m-matrix of the symmetric case
				int x_2 = j;				// reference spline is j-th
				int y_2 = i - (j - order);
				int z_2 = iknot - j;
				
				// save to m-matrix
				m[(x_1 * (2 * order + 1) + y_1) * (order + 1) + z_1] = integral;
				m[(x_2 * (2 * order + 1) + y_2) * (order + 1) + z_2] = integral;
			}
		}
	}
	
	return m;
}


Complex computeD_iknot(int i, int j, int iknot)
{
	// get interval boundaries
	Complex x1 = t[iknot];
	Complex x2 = t[iknot + 1];
	
	// throw away zero-length intervals
	if (x1 == x2)
		return 0;
	
	// get Gauss-Legendre nodes and weights for the interval [-1, 1]
	const double *vx = 0, *vw = 0;
	int points = p_roots(20/*order*/, vx, vw);
	
	// actual stored data count (only non-negative values)
	int Ndata = (points + 1) / 2;
	
	// mirror evaluation nodes and weights
	Complex xs[points], ws[points];
	for (int dat = 0; dat < Ndata; dat++)
	{
		// stretch nodes
		xs[dat]          = x1 + 0.5 * (1.0 - vx[Ndata-1-dat]) * (x2 - x1);
		xs[points-1-dat] = x1 + 0.5 * (1.0 + vx[Ndata-1-dat]) * (x2 - x1);
		
		// scale weights
		ws[dat]          = 0.5 * vw[Ndata-1-dat] * (x2 - x1);
		ws[points-1-dat] = 0.5 * vw[Ndata-1-dat] * (x2 - x1);
	}
	
	// evaluate B-splines at Gauss-Legendre nodes
	Complex values_i[points], values_j[points];
	dB(i, iknot, points, xs, values_i);
	dB(j, iknot, points, xs, values_j);
	
	// result
	Complex res = 0;
	
	// accumulate the result
	for (int k = 0; k < points; k++)
		res += values_i[k] * values_j[k] * ws[k];
	
	return res;
}

Complex computeD(int i, int j, int maxknot)
{
	// get boundary iknots
	int left = std::max(i, j);
	int right = std::min(i, j) + order;
	
	// cut at maxknot
	if (right > maxknot)
		right = maxknot;
	
	// the result
	Complex res = 0;
	
	// undergo integration on sub-intervals
	for (int iknot = left; iknot <= right; iknot++)
		res += computeD_iknot(i, j, iknot);
		
	return res;
}

Complex computeM_iknot(int a, int i, int j, int iknot, Complex R)
{
	// get interval boundaries
	Complex x1 = t[iknot];
	Complex x2 = t[iknot + 1];
	
	// throw away zero-length intervals
	if (x1 == x2)
		return 0;
	
	// get Gauss-Legendre nodes and weights for the interval [-1, 1]
	const double *vx = 0, *vw = 0;
	int points = p_roots(20/*order + 10*/, vx, vw);
	
	// actual stored data count (only non-negative values)
	int Ndata = (points + 1) / 2;
	
	// mirror evaluation nodes and weights
	Complex xs[points], ws[points];
	for (int dat = 0; dat < Ndata; dat++)
	{
		// stretch nodes
		xs[dat]          = x1 + 0.5 * (1.0 - vx[Ndata-1-dat]) * (x2 - x1);
		xs[points-1-dat] = x1 + 0.5 * (1.0 + vx[Ndata-1-dat]) * (x2 - x1);
		
		// scale weights
		ws[dat]          = 0.5 * vw[Ndata-1-dat] * (x2 - x1);
		ws[points-1-dat] = 0.5 * vw[Ndata-1-dat] * (x2 - x1);
	}
	
	// evaluate B-splines at Gauss-Legendre nodes
	Complex values_i[points], values_j[points];
	B(i, iknot, points, xs, values_i);
	B(j, iknot, points, xs, values_j);
	
	// result
	Complex res = 0;
	
	// accumulate the (damped) result
	if (R != 0.)
	{
		for (int k = 0; k < points; k++)
			res += values_i[k] * values_j[k] * pow(xs[k],a) * ws[k] * damp(xs[k],R);
	}
	else
	{
		for (int k = 0; k < points; k++)
			res += values_i[k] * values_j[k] * pow(xs[k],a) * ws[k];
	}
	
	return res;
}

Complex computeM(int a, int i, int j, int maxknot)
{
	// get boundary iknots
	int left = std::max(i, j);
	int right = std::min(i, j) + order;
	
	// cut at maxknot
	if (maxknot != 0 and right > maxknot)
		right = maxknot;
	
	// the result
	Complex res = 0;
	
	// undergo integration on sub-intervals
	for (int iknot = left; iknot <= right; iknot++)
		res += computeM_iknot(a, i, j, iknot, t[maxknot]);
	
	return res;
}
