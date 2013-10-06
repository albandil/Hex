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

#include <algorithm>
#include <vector>

#include "arrays.h"
#include "bspline.h"
#include "complex.h"

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


void Bspline::B (int i, int iknot, int n, Complex const * const restrict x, Complex * const restrict y) const
{
    for (int j = 0; j < n; j++)
        y[j] = bspline(i, iknot, order_, x[j]);
}

void Bspline::dB(int i, int iknot, int n, Complex const * const restrict x, Complex * const restrict y) const
{
    for (int j = 0; j < n; j++)
        y[j] = dspline(i, iknot, order_, x[j]);
}


// ----------------------------------------------------------------------- //
//  Evaluation of B-spline expansion on a grid                             //
// ----------------------------------------------------------------------- //


cArray Bspline::zip (cArrayView const & coeff, rArrayView const & grid) const
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
    
    // for all intervals
    for (int iknot = 0; iknot < Nknot_ - 1; iknot++)
    {
        // skip zero-length intervals
        if (t_[iknot] == t_[iknot+1])
            continue;
        
        // get subset of "x", that belongs to the interval
        // t[iknot] .. t[iknot + 1]
        left  = std::lower_bound(right, x.end(), t_[iknot],     Complex_realpart_less);
        right = std::upper_bound(left,  x.end(), t_[iknot + 1], Complex_realpart_less);
        
        // for all points in range
        for (auto ix = left; ix != right; ix++)
        {
            // evaluated B-splines
            Complex yn = 0;
            
            // for all splines
            for (int ispline = std::max(iknot-order_,0); ispline <= iknot and ispline < Nspline_; ispline++)
            {
                if (coeff[ispline] == 0.)
                    continue;
                
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
    cArrayView const & coeff,
    rArrayView const & xgrid,
    rArrayView const & ygrid
) const {
    
    // evaluated function
    cArray f(xgrid.size() * ygrid.size());
    
    // rotate the grids
    cArray x(xgrid.size()), y(ygrid.size());
    for (size_t i = 0; i < xgrid.size(); i++)
        x[i] = rotate(xgrid[i]);
    for (size_t i = 0; i < ygrid.size(); i++)
        y[i] = rotate(ygrid[i]);
    
    // evaluate B-splines on x-grid
    cArrays evBx(Nspline_);
    cArray::const_iterator xleft = x.begin(), xright = x.begin();
    for (int ispline = 0; ispline < Nspline_; ispline++)
    {
        // get relevant subset of x[]
        xleft = std::upper_bound (
            x.cbegin(), x.cend(),
            t_[ispline],
            Complex_realpart_less
        );
        xright = std::lower_bound (
            xleft, x.cend(),
            t_[ispline + order_ + 1],
            Complex_realpart_less
        );
        
        // setup evaluation vector
        evBx[ispline].resize(xright-xleft);
        
        // evaluation interval
        int iknot = ispline;
        
        // evaluate at x[]
        for (cArray::const_iterator ix = xleft; ix != xright; ix++)
        {
            // increment knot
            while ( ix->real() > t_[iknot+1].real() )
                if (++iknot >= Nknot_)
                    throw exception("Some evaluation points are outside of grid.");
            
            // evaluate this spline at *ix
            evBx[ispline][ix-xleft] = bspline(ispline, iknot, order_, *ix);
        }
    }
    
    // evaluate B-splines on y-grid
    cArrays evBy(Nspline_);
    cArray::const_iterator yleft = y.begin(), yright = y.begin();
    for (int ispline = 0; ispline < Nspline_; ispline++)
    {
        // get relevant subset of y[]
        yleft = std::upper_bound (
            y.cbegin(), y.cend(),
            t_[ispline],
            Complex_realpart_less
        );
        yright = std::lower_bound (
            yleft, y.cend(),
            t_[ispline + order_ + 1],
            Complex_realpart_less
        );
        
        // setup evaluation vector
        evBy[ispline].resize(yright-yleft);
        
        // evaluation interval
        int iknot = ispline;
        
        // evaluate at y[]
        for (cArray::const_iterator iy = yleft; iy != yright; iy++)
        {
            // increment knot
            while ( iy->real() > t_[iknot+1].real() )
                if (++iknot >= Nknot_)
                    throw exception("Some evaluation points are outside of grid.");
            
            // evaluate this spline at *ix
            evBy[ispline][iy-yleft] = bspline(ispline, iknot, order_, *iy);
        }
    }
    
    // zip double expansion
    xleft = xright = x.begin();
    for (int ixspline = 0; ixspline < Nspline_; ixspline++)
    {
        // get relevant subset of x[]
        xleft = std::upper_bound (
            x.cbegin(), x.cend(),
            t_[ixspline],
            Complex_realpart_less
        );
        xright = std::lower_bound (
            xleft, x.cend(),
            t_[ixspline + order_ + 1],
            Complex_realpart_less
        );
        
        // get relevant evaluations
        cArray const & Bx_row = evBx[ixspline];
        
        yleft = yright = y.begin();
        for (int iyspline = 0; iyspline < Nspline_; iyspline++)
        {
            // get relevant subset of y[]
            yleft = std::upper_bound (
                y.cbegin(), y.cend(),
                t_[iyspline],
                Complex_realpart_less
            );
            yright = std::lower_bound (
                yleft, y.cend(),
                t_[iyspline + order_ + 1],
                Complex_realpart_less
            );
            
            // get relevant evaluations
            cArray const & By_row = evBy[iyspline];
            
            // get coefficient of the expansion
            Complex C = coeff[ixspline * Nspline_ + iyspline];
            
            // loop over relevant points
            for (cArray::const_iterator ix = xleft; ix != xright; ix++)
            {
                Complex CBx = C * Bx_row[ix-xleft];
                Complex *fx = &f[0] + (ix-x.begin()) * y.size();
                
                for (cArray::const_iterator iy = yleft; iy != yright; iy++)
                    fx[iy-y.begin()] += CBx * By_row[iy-yleft];
            }
        }
    }
    
    return f;
}


// ----------------------------------------------------------------------- //
//  B-spline environment setup                                             //
// ----------------------------------------------------------------------- //


Bspline::Bspline (int order, rArrayView const & rknots, double th, rArrayView const & cknots)
{
    // globalize
    order_ = order;
    R0_ = rknots.back();
    Rmax_ = cknots.back();
    
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

// ----------------------------------------------------------------------- //
//  Others                                                                 //
// ----------------------------------------------------------------------- //

int Bspline::knot(Complex x) const
{
    // get "lower" bound by bisection (will return the first equal or greater element)
    Complex* iknot_notless = std::lower_bound (
        t_,             // search from here ...
        t_ + Nknot_,    // ... to here (exclusively)
        x,              // and compare with respect to this item
        // comparator using the real parts
        [](Complex const & a, Complex const & x) -> bool {
            return a.real() < x.real();
        }
    );
    
    // check if this is a valid knot
    if (iknot_notless != t_ /* x > 0 */ and iknot_notless != t_ + Nknot_ /* x < Rmax */)
        return iknot_notless - t_ - 1;
    else
        return -1;
}

Complex Bspline::eval(cArrayView const & coeff, double x) const
{
    Complex z = rotate(x);
    
    // get knot index
    int iknot = knot(z);
    
    // get bounding B-splines
    int leftspline = iknot-order_;
    int rightspline = iknot;
    
    // evaluate B-splines
    cArray evB(Nspline_);
    for (int ispline = leftspline; ispline <= rightspline; ispline++)
        evB[ispline] = bspline(ispline,iknot,order_,z);
    
    // sum expansion
    Complex result = 0.;
    for (int ispline = leftspline; ispline <= rightspline; ispline++)
        result += coeff[ispline] * evB[ispline];
    return result;
}

Complex Bspline::eval(cArrayView const & coeff, double x, double y) const
{
    Complex w = rotate(x);
    Complex z = rotate(y);
    
    // get knot indices
    int xknot = knot(w);
    int yknot = knot(z);
    
    // get bounding B-splines
    int leftxspline = xknot-order_;
    int rightxspline = xknot;
    int leftyspline = yknot-order_;
    int rightyspline = yknot;
    
    // evaluate B-splines
    cArray evBx(Nspline_), evBy(Nspline_);
    for (int ixspline = leftxspline; ixspline <= rightxspline; ixspline++)
        evBx[ixspline] = bspline(ixspline,xknot,order_,w);
    for (int iyspline = leftyspline; iyspline <= rightyspline; iyspline++)
        evBx[iyspline] = bspline(iyspline,yknot,order_,z);
    
    // sum the expansion
    Complex result = 0.;
    for (int ixspline = leftxspline; ixspline <= rightxspline; ixspline++)
    for (int iyspline = leftyspline; iyspline <= rightyspline; iyspline++)
        result += coeff[ixspline * Nspline_ + iyspline] * evBx[ixspline] * evBy[iyspline];
    return result;
}
