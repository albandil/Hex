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
//  Recursive single-point B-spline evaluation                             //
// ----------------------------------------------------------------------- //


Complex Bspline::bspline(int i, int iknot, int k, Complex r) const
{
	if (k == 0)
		return  (i == iknot) ? 1. : 0.;
	
	Complex A = bspline(i,   iknot, k-1, r);
	Complex B = bspline(i+1, iknot, k-1, r);
	
	Complex S1 = (t_[i+k]   != t_[i])   ? (r - t_[i]) / (t_[i+k] - t_[i]) : 0.;
	Complex S2 = (t_[i+k+1] != t_[i+1]) ? (t_[i+k+1] - r) / (t_[i+k+1] - t_[i+1]) : 0.;
	
	return A * S1 + B * S2;
}

Complex Bspline::dspline(int i, int iknot, int k, Complex r) const
{
	if (k == 0)
		return 0.;
	
	Complex A = bspline(i,   iknot, k-1, r);
	Complex B = bspline(i+1, iknot, k-1, r);
	
	Complex S1 = (t_[i+k]   != t_[i])   ? double(k) / (t_[i+k] - t_[i]) : 0.;
	Complex S2 = (t_[i+k+1] != t_[i+1]) ? double(k) / (t_[i+k+1] - t_[i+1]) : 0.;
	
	return A * S1 - B * S2;
}


// ----------------------------------------------------------------------- //
//  Evaluation of B-splines on a grid                                      //
// ----------------------------------------------------------------------- //


void Bspline::B(int i, int iknot, int n, const Complex* x, Complex* y) const
{
	for (int j = 0; j < n; j++)
		y[j] = bspline(i, iknot, order_, x[j]);
}

void Bspline::dB(int i, int iknot, int n, const Complex* x, Complex* y) const
{
	for (int j = 0; j < n; j++)
		y[j] = dspline(i, iknot, order_, x[j]);
}


// ----------------------------------------------------------------------- //
//  Evaluation of B-spline expansion on a grid                             //
// ----------------------------------------------------------------------- //


cArray Bspline::zip (cArray const & coeff, rArray const & grid) const
{
	// evaluated function
	cArray f(grid.size());
	
	// rotate the grid
	cArray x(grid.size());
	for (size_t i = 0; i < grid.size(); i++)
		x[i] = rotate(grid[i]);
	
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
	for (int iknot = 0; iknot < Nknot_ - 1; iknot++)
	{
		// skip zero-length intervals
		if (t_[iknot] == t_[iknot+1])
			continue;
		
		// get subset of "x", that belongs to the interval
		// t[iknot] .. t[iknot + 1]
		left = std::lower_bound(right, x.end(), t_[iknot], comp);
		right = std::upper_bound(left, x.end(), t_[iknot + 1], comp);
		
		// for all points in range
		for (auto ix = left; ix != right; ix++)
		{
			// evaluated B-splines
			Complex yn = 0;
			
			// for all splines
			for (int ispline = std::max(iknot-order_,0); ispline <= iknot and ispline < Nspline_; ispline++)
			{
				Complex y1;
				B (ispline, iknot, 1, &*ix, &y1);
				yn += coeff[ispline] * y1;
			}
			
			// add evaluated point
			f[ix-x.begin()] = yn;
		}
	}
	
	return f;
}

cArray Bspline::zip (
	cArray const & coeff,
	rArray const & xgrid,
	rArray const & ygrid
) const {
	// evaluated function
	cArray f(xgrid.size() * ygrid.size());
	
	// rotate the grid
	cArray x(xgrid.size()), y(ygrid.size());
	for (size_t i = 0; i < xgrid.size(); i++)
		x[i] = rotate(xgrid[i]);
	for (size_t i = 0; i < ygrid.size(); i++)
		y[i] = rotate(ygrid[i]);
	
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
	for (int ixknot = 0; ixknot < Nknot_ - 1; ixknot++)
	{
		// skip zero-length intervals
		if (t_[ixknot] == t_[ixknot+1])
			continue;
		
		// get subset of "x", that belongs to the interval
		// t[iknot] .. t[iknot + 1]
		xleft = std::lower_bound(xright, x.end(), t_[ixknot], comp);
		xright = std::upper_bound(xleft, x.end(), t_[ixknot + 1], comp);
		
		yleft = y.begin();
		yright = y.begin();
		
		// for all y-intervals
		for (int iyknot = 0; iyknot < Nknot_ - 1; iyknot++)
		{
			// skip zero-length intervals
			if (t_[iyknot] == t_[iyknot+1])
				continue;
			
			// get subset of "y", that belongs to the interval
			// t[iknot] .. t[iknot + 1]
			yleft = std::lower_bound(yright, y.end(), t_[iyknot], comp);
			yright = std::upper_bound(yleft, y.end(), t_[iyknot + 1], comp);
			
			// for all points in x-range
			for (auto ix = xleft; ix != xright; ix++)
			{
				// for all points in y-range
				for (auto iy = yleft; iy != yright; iy++)
				{
					// evaluated B-splines
					Complex zn = 0;
					
					// for all x-splines contributing at ixknot
					for (int ixspline = std::max(ixknot-order_,0);
						 ixspline <= ixknot and ixspline < Nspline_; ixspline++)
					{
						Complex z1;
						B(ixspline, ixknot, 1, &*ix, &z1);
						
						// for all y-splines contributing at iyknot
						for (int iyspline = std::max(iyknot-order_,0);
						     iyspline <= iyknot and iyspline < Nspline_; iyspline++)
						{
							Complex z2;
							B(iyspline, iyknot, 1, &*iy, &z2);
							
							zn += coeff[ixspline*Nspline_+iyspline] * z1 * z2;
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


void Bspline::init (
	int order, rArray rknots, double R0, double th, rArray cknots, double Rmax
){
	// globalize
	order_ = order;
	R0_ = R0;
	Rmax_ = Rmax;
	
	// real and complex knot counts; both include the knot R₀
	int rknots_len = rknots.size();
	int cknots_len = cknots.size();
	
	// the complex rotation factor
	rotation_ = Complex(cos(th),sin(th));
	
	// join the sequences; rotate the complex one
	t_ = new Complex [rknots_len + cknots_len - 1];
	for (int i = 0; i < rknots_len; i++)
		t_[i] = rknots[i];
	for (int i = 1; i < cknots_len; i++)
		t_[i + rknots_len - 1] = rotate(cknots[i]);
	
	// set knot count and other global variables
	Nreknot_ = rknots_len;
	Nknot_ = rknots_len + cknots_len - 1;
	Nintval_ = Nknot_ - 1;
	Nspline_ = Nknot_ - order - 1;
}
