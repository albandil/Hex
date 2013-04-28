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
#include <vector>

#include "arrays.h"
#include "bspline.h"

// ----------------------------------------------------------------------- //
//  Global data                                                            //
// ----------------------------------------------------------------------- //


Complex* t = 0;	// knot sequence
Complex rotation;	// ECS rotation factor
double R0;			// ECS edge
double Rmax;		// grid end
int Nknot;		// knot count
int Nreknot;	// real knot count
int Nspline;	// B-spline count
int Nintval;	// interval count (just shorthand for Nknot - 1)
int order;		// B-spline order


// ----------------------------------------------------------------------- //
//  Recursive single-point B-spline evaluation                             //
// ----------------------------------------------------------------------- //


Complex _bspline(int i, int iknot, int k, Complex r)
{
	if (k == 0)
		return  (i == iknot) ? 1. : 0.;
	
	Complex A = _bspline(i,   iknot, k-1, r);
	Complex B = _bspline(i+1, iknot, k-1, r);
	
	Complex S1 = (t[i+k]   != t[i])   ? (r - t[i]) / (t[i+k] - t[i]) : 0.;
	Complex S2 = (t[i+k+1] != t[i+1]) ? (t[i+k+1] - r) / (t[i+k+1] - t[i+1]) : 0.;
	
	return A * S1 + B * S2;
}

Complex _dspline(int i, int iknot, int k, Complex r)
{
	if (k == 0)
		return 0.;
	
	Complex A = _bspline(i,   iknot, k-1, r);
	Complex B = _bspline(i+1, iknot, k-1, r);
	
	Complex S1 = (t[i+k]   != t[i])   ? double(k) / (t[i+k] - t[i]) : 0.;
	Complex S2 = (t[i+k+1] != t[i+1]) ? double(k) / (t[i+k+1] - t[i+1]) : 0.;
	
	return A * S1 - B * S2;
}


// ----------------------------------------------------------------------- //
//  Evaluation of B-splines on a grid                                      //
// ----------------------------------------------------------------------- //


void B(int i, int iknot, int n, const Complex* x, Complex* y)
{
	for (int j = 0; j < n; j++)
		y[j] = _bspline(i, iknot, order, x[j]);
}

void dB(int i, int iknot, int n, const Complex* x, Complex* y)
{
	for (int j = 0; j < n; j++)
		y[j] = _dspline(i, iknot, order, x[j]);
}


// ----------------------------------------------------------------------- //
//  Evaluation of B-spline expansion on a grid                             //
// ----------------------------------------------------------------------- //


cArray zip(
	cArray const & coeff,
	rArray const & grid
) {
	// evaluated function
	cArray f(grid.size());
	
	// rotate the grid
	cArray x(grid.size());
	for (size_t i = 0; i < grid.size(); i++)
		x[i] = ECS_rotate(grid[i]);
	
	// iterators pointing to the left and right boundary of evaluation
	// points subset belonging to some interknot interval
	cArray::iterator left = x.begin();
	cArray::iterator right = x.begin();
	
	// comparator of points' distance from origin along the contour
	auto comp = [ = ](const Complex& a,
					  const Complex& b) -> bool {
						  return a.real() < b.real();
					};
	
	// for all intervals
	for (int iknot = 0; iknot < Nknot - 1; iknot++)
	{
		// skip zero-length intervals
		if (t[iknot] == t[iknot+1])
			continue;
		
		// get subset of "x", that belongs to the interval
		// t[iknot] .. t[iknot + 1]
		left = std::lower_bound(right, x.end(), t[iknot], comp);
		right = std::upper_bound(left, x.end(), t[iknot + 1], comp);
		
		// for all points in range
		for (auto ix = left; ix != right; ix++)
		{
			// evaluated B-splines
			Complex yn = 0;
			
			// for all splines
			for (int ispline = std::max((int)iknot-(int)order,0); ispline <= iknot and ispline < Nspline; ispline++)
			{
				Complex y1;
				B(ispline, iknot, 1, &*ix, &y1);
				yn += coeff[ispline] * y1;
			}
			
			// add evaluated point
			f[ix-x.begin()] = yn;
		}
	}
	
	return f;
}

cArray zip(
	cArray const & coeff,
	rArray const & xgrid,
	rArray const & ygrid
){
	// evaluated function
	cArray f(xgrid.size() * ygrid.size());
	
	// rotate the grid
	cArray x(xgrid.size()), y(ygrid.size());
	for (size_t i = 0; i < xgrid.size(); i++)
		x[i] = ECS_rotate(xgrid[i]);
	for (size_t i = 0; i < ygrid.size(); i++)
		y[i] = ECS_rotate(ygrid[i]);
	
	// iterators pointing to the left and right boundary of evaluation
	// points subset belonging to some interknot interval
	cArray::iterator xleft, xright, yleft, yright;
	
	// comparator of points' distance from origin along the contour
	auto comp = [ = ](const Complex& a, const Complex& b) -> bool {
		return a.real() < b.real();
	};
	
	xleft = x.begin();
	xright = x.begin();
	
	// for all x-intervals
	for (int ixknot = 0; ixknot < Nknot - 1; ixknot++)
	{
		// skip zero-length intervals
		if (t[ixknot] == t[ixknot+1])
			continue;
		
		// get subset of "x", that belongs to the interval
		// t[iknot] .. t[iknot + 1]
		xleft = std::lower_bound(xright, x.end(), t[ixknot], comp);
		xright = std::upper_bound(xleft, x.end(), t[ixknot + 1], comp);
		
		yleft = y.begin();
		yright = y.begin();
		
		// for all y-intervals
		for (int iyknot = 0; iyknot < Nknot - 1; iyknot++)
		{
			// skip zero-length intervals
			if (t[iyknot] == t[iyknot+1])
				continue;
			
			// get subset of "y", that belongs to the interval
			// t[iknot] .. t[iknot + 1]
			yleft = std::lower_bound(yright, y.end(), t[iyknot], comp);
			yright = std::upper_bound(yleft, y.end(), t[iyknot + 1], comp);
			
			// for all points in x-range
			for (auto ix = xleft; ix != xright; ix++)
			{
				// for all points in y-range
				for (auto iy = yleft; iy != yright; iy++)
				{
					// evaluated B-splines
					Complex zn = 0;
					
					// for all x-splines contributing at ixknot
					for (int ixspline = std::max((int)ixknot-(int)order,0);
						 ixspline <= ixknot and ixspline < Nspline; ixspline++)
					{
						Complex z1;
						B(ixspline, ixknot, 1, &*ix, &z1);
						
						// for all y-splines contributing at iyknot
						for (int iyspline = std::max((int)iyknot-(int)order,0);
						     iyspline <= iyknot and iyspline < Nspline; iyspline++)
						{
							Complex z2;
							B(iyspline, iyknot, 1, &*iy, &z2);
							
							zn += coeff[ixspline*Nspline+iyspline] * z1 * z2;
						}
					}
					
					// add evaluated point
					f[(ix-x.begin())*ygrid.size()+(iy-y.begin())] = zn;
				}
			}
		}
	}
	
	return f;
}


// ----------------------------------------------------------------------- //
//  B-spline environment setup                                             //
// ----------------------------------------------------------------------- //


void setup_knot_sequence(int _order_, rArray rknots, double _R0_, double th, rArray cknots, double _Rmax_)
{
	// globalize
	order = _order_;
	R0 = _R0_;
	Rmax = _Rmax_;
	
	// real and complex knot counts; both include the knot R₀
	int rknots_len = rknots.size();
	int cknots_len = cknots.size();
	
	// the complex rotation factor
	rotation = Complex(cos(th),sin(th));
	
	// join the sequences; rotate the complex one
	t = new Complex [rknots_len + cknots_len - 1];
	for (int i = 0; i < rknots_len; i++)
		t[i] = rknots[i];
	for (int i = 1; i < cknots_len; i++)
		t[i + rknots_len - 1] = ECS_rotate(cknots[i]);
	
	// set knot count and other global variables
	Nreknot = rknots_len;
	Nknot = rknots_len + cknots_len - 1;
	Nintval = Nknot - 1;
	Nspline = Nknot - order - 1;
}
