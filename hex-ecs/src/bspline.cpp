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

#include <algorithm>
#include <vector>

#include "arrays.h"
#include "bspline.h"
#include "memory.h"


// ----------------------------------------------------------------------- //
//  Recursive single-point B-spline evaluation                             //
// ----------------------------------------------------------------------- //


Complex Bspline::bspline (int i, int iknot, int k, Complex r) const
{
    // NOTE: The following bounds check is the caller's responsibility:
    // if (r <= t_[i] or t_[i + k] <= r)
    //    return 0.;
    
    // value of the parent B-splines of the requested B-spline
    Complex b[k + 1];
    
    // initialize zero-order B-splines
    for (int n = 0; n <= k; n++)
        b[n] = (i + n == iknot ? 1. : 0.);
    
    // calculate higher orders
    for (int ord = 1; ord <= k; ord++)
    {
        // update splines
        for (int n = 0; n <= k-ord; n++)
        {
            b[n] = (t_[i+ord+n]   == t_[i+n]   ? 0. : b[n]   * (r    -    t_[i+n]) / (t_[i+ord+n]   -   t_[i+n]))
                 + (t_[i+ord+n+1] == t_[i+n+1] ? 0. : b[n+1] * (t_[i+ord+n+1] - r) / (t_[i+ord+n+1] - t_[i+n+1]));
        }
    }
    
    // return the collected value of the requested B-spline
    return b[0];
}

Complex Bspline::dspline (int i, int iknot, int k, Complex r) const
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

void Bspline::B (int i, int iknot, int M, Complex const * const restrict x, Complex * const restrict y) const
{
    // NOTE: The caller's responsibility is to check that all 'x' lie in the interval (t[iknot],t[iknot+1])
    //       and that the i-th B-spline is defined there.
    
    //
    // a) Use real arithmetic when using only real knots. Allow SIMD auto-vectorization.
    //
    
    if (i + order_ + 1 < Nreknot_)
    {
        // number of needed SIMD vectors for M doubles
        int nvec = (M + simd_double_vec_size - 1) / simd_double_vec_size;
        
        // copy real parts of the evaluation points, pad by zeros
        simd_double_vec_t rx[nvec];
        for (int m = 0; m < nvec * (int)simd_double_vec_size; m++)
            rx[m / simd_double_vec_size][m % simd_double_vec_size] = (m < M ? x[m].real() : 0);
        
        // evaluations of the parent B-splines of the wanted B-spline
        simd_double_vec_t b[order_ + 1][nvec];
        
        // initialize all ancestral zero-order B-splines
        for (int n = 0; n <= order_; n++)
        {
            // determine value of the zero-order B-spline B_n on interval (t[iknot],t[iknot+1])
            double val = (n + i == iknot ? 1. : 0.);
            
            // store the value at all points
            for (int m = 0; m < nvec; m++)
            for (int v = 0; v < (int)simd_double_vec_size; v++) // <-- likely to auto-vectorize
                b[n][m][v] = val;
        }
        
        // precomputed denominators (used later)
        double invden[order_ + 1];
        
        // real knots restricted pointer (for fast access)
        double const * const restrict rknots = rknots_.data();
        
        // calculate B-splines of higher orders
        for (int ord = 1; ord <= order_; ord++)
        {
            // precompute denominators
            for (int n = 0; n <= order_ - ord + 1; n++)
                invden[n] = (rknots[i+ord+n] == rknots[i+n] ? 0. : 1. / (rknots[i+ord+n] - rknots_[i+n]));
            
            // evaluate B-splines from lower orders
            for (int n = 0; n <= order_ - ord; n++)
            for (int m = 0; m < nvec; m++)
            for (int v = 0; v < (int)simd_double_vec_size; v++) // <-- likely to auto-vectorize
                b[n][m][v] = b[n][m][v] * (rx[m][v] - rknots[i+n]) * invden[n] + b[n+1][m][v] * (rknots[i+ord+n+1] - rx[m][v]) * invden[n+1];
        }
        
        // return the collected value of the requested B-spline
        for (int m = 0; m < M; m++)
            y[m] = b[0][m / simd_double_vec_size][m % simd_double_vec_size];
    }
    
    //
    // b) Use complex arithmetic otherwise.
    //
    
    else
    {
        // value of the parent B-splines of the requested B-spline
        Complex b[M][order_ + 1];
        
        // initialize zero-order B-splines
        for (int m = 0; m < M; m++)
        for (int n = 0; n <= order_; n++)
            b[m][n] = (i + n == iknot ? 1. : 0.);
        
        // calculate higher orders
        for (int ord = 1; ord <= order_; ord++)
        {
            // update splines
            for (int n = 0; n <= order_ - ord; n++)
            {
                Complex invden1 = (t_[i+ord+n]   == t_[i+n]   ? 0. : 1. / (t_[i+ord+n]   - t_[i+n]));
                Complex invden2 = (t_[i+ord+n+1] == t_[i+n+1] ? 0. : 1. / (t_[i+ord+n+1] - t_[i+n+1]));
                
                // for all evaluation points
                for (int m = 0; m < M; m++)
                    b[m][n] = b[m][n] * (x[m] - t_[i+n]) * invden1 + b[m][n+1] * (t_[i+ord+n+1] - x[m]) * invden2;
            }
        }
        
        // return the collected value of the requested B-spline
        for (int m = 0; m < M; m++)
            y[m] = b[m][0];
    }
}

void Bspline::dB (int i, int iknot, int n, Complex const * const restrict x, Complex * const restrict y) const
{
    for (int j = 0; j < n; j++)
        y[j] = dspline(i, iknot, order_, x[j]);
}


// ----------------------------------------------------------------------- //
//  Evaluation of B-spline expansion on a grid                             //
// ----------------------------------------------------------------------- //


cArray Bspline::zip (const cArrayView coeff, const rArrayView grid) const
{
    // evaluated function
    cArray f(grid.size());
    
    // rotate the grid
    cArray x(grid.size());
    for (std::size_t i = 0; i < grid.size(); i++)
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

cArray Bspline::zip (Bspline const & bx, Bspline const & by, const cArrayView coeff, const rArrayView xgrid, const rArrayView ygrid)
{
    // evaluated function
    cArray f (xgrid.size() * ygrid.size());
    
    // rotate the grids
    cArray x (xgrid.size()), y (ygrid.size());
    for (std::size_t i = 0; i < xgrid.size(); i++)
        x[i] = bx.rotate(xgrid[i]);
    for (std::size_t i = 0; i < ygrid.size(); i++)
        y[i] = by.rotate(ygrid[i]);
    
    // evaluate B-splines (from "bx") on x-grid
    cArrays evBx (bx.Nspline_);
    cArray::const_iterator xleft = x.begin(), xright;
    for (int ispline = 0; ispline < bx.Nspline_; ispline++)
    {
        // get first x-coordinate within this spline's support
        while (xleft < x.end() and xleft->real() < bx.t_[ispline].real())
            xleft++;
        
        // get (one beyond the) last x-coordinate within this spline's support
        xright = xleft;
        while (xright < x.end() and xright->real() < bx.t_[ispline + bx.order_ + 1].real())
            xright++;
        
        // setup evaluation vector
        evBx[ispline].resize(xright - xleft);
        
        // evaluation interval
        int iknot = ispline;
        
        // evaluate at x[]
        for (cArray::const_iterator ix = xleft; ix != xright; ix++)
        {
            // increment knot to match current evaluation point
            while ( ix->real() > bx.t_[iknot + 1].real() )
                if (++iknot >= bx.Nknot_)
                    HexException("Some evaluation points are outside of grid.");
            
            // evaluate this spline at the current evaluation point
            evBx[ispline][ix - xleft] = bx.bspline(ispline, iknot, bx.order_, *ix);
        }
    }
    
    // evaluate B-splines (from "by") on y-grid
    cArrays evBy (by.Nspline_);
    cArray::const_iterator yleft = y.begin(), yright;
    for (int ispline = 0; ispline < by.Nspline_; ispline++)
    {
        // get first y-coordinate within this spline's support
        while (yleft < y.end() and yleft->real() < by.t_[ispline].real())
            yleft++;
        
        // get (one beyond the) last y-coordinate within this spline's support
        yright = yleft;
        while (yright < y.end() and yright->real() < by.t_[ispline + by.order_ + 1].real())
            yright++;
        
        // setup evaluation vector
        evBy[ispline].resize(yright - yleft);
        
        // evaluation interval
        int iknot = ispline;
        
        // evaluate at y[]
        for (cArray::const_iterator iy = yleft; iy != yright; iy++)
        {
            // increment knot to match current evaluation point
            while ( iy->real() > by.t_[iknot + 1].real() )
                if (++iknot >= by.Nknot_)
                    HexException("Some evaluation points are outside of grid.");
            
            // evaluate this spline at the current evaluation point
            evBy[ispline][iy - yleft] = by.bspline(ispline, iknot, by.order_, *iy);
        }
    }
    
    // zip double expansion
    xleft = x.begin();
    for (int ixspline = 0; ixspline < bx.Nspline_; ixspline++)
    {
        // get first x-coordinate within this spline's support
        while (xleft < x.end() and xleft->real() < bx.t_[ixspline].real())
            xleft++;
        
        // get (one beyond the) last x-coordinate within this spline's support
        xright = xleft;
        while (xright < x.end() and xright->real() < bx.t_[ixspline + bx.order_ + 1].real())
            xright++;
        
        // get relevant evaluations
        cArray const & Bx_row = evBx[ixspline];
        
        // for all y-splines
        yleft = y.begin();
        for (int iyspline = 0; iyspline < by.Nspline_; iyspline++)
        {
            // get first y-coordinate within this spline's support
            while (yleft < y.end() and yleft->real() < by.t_[iyspline].real())
                yleft++;
            
            // get (one beyond the) last y-coordinate within this spline's support
            yright = yleft;
            while (yright < y.end() and yright->real() < by.t_[iyspline + by.order_ + 1].real())
                yright++;
            
            // get relevant evaluations
            cArray const & By_row = evBy[iyspline];
            
            // get coefficient of the expansion
            Complex C = coeff[ixspline * by.Nspline_ + iyspline];
            
            // loop over relevant points
            for (cArray::const_iterator ix = xleft; ix != xright; ix++)
            for (cArray::const_iterator iy = yleft; iy != yright; iy++)
            {
                // get position indices within various arrays
                std::size_t idx = ix - x.begin(), evx = ix - xleft;
                std::size_t idy = iy - y.begin(), evy = iy - yleft;
                
                // update the evaluation
                f[idy * x.size() + idx] += C * Bx_row[evx] * By_row[evy];
            }
        }
    }
    
    return f;
}


cArray Bspline::zip (const cArrayView coeff, const rArrayView xgrid, const rArrayView ygrid) const
{
    return zip (*this, *this, coeff, xgrid, ygrid);
}

// ----------------------------------------------------------------------- //
//  B-spline environment setup                                             //
// ----------------------------------------------------------------------- //


Bspline::Bspline (int order, rArrayView const & rknots, double th, rArrayView const & cknots)
    : rknots_(rknots), cknots_(cknots), theta_(th), rotation_(Complex(cos(th),sin(th))),
      R0_(rknots_.back()), Rmax_(cknots_.back()), Nknot_(rknots.size() + cknots.size() - 1),
      Nreknot_(rknots.size()), Nspline_(Nknot_ - order - 1), Nintval_(Nknot_ - 1),
      order_(order)
{
    // real and complex knot counts; both include the knot R₀
    int rknots_len = rknots.size();
    int cknots_len = cknots.size();
    
    // join the sequences; rotate the complex one
    t_ = new Complex [rknots_len + cknots_len - 1];
    for (int i = 0; i < rknots_len; i++)
        t_[i] = rknots[i];
    for (int i = 1; i < cknots_len; i++)
        t_[i + rknots_len - 1] = rotate(cknots[i]);
}


// ----------------------------------------------------------------------- //
//  Others                                                                 //
// ----------------------------------------------------------------------- //


int Bspline::knot (Complex x) const
{
    // get "lower" bound by bisection (will return the first equal or greater element)
    Complex* iknot_notless = std::lower_bound
    (
        t_,                     // search from here ...
        t_ + Nknot_,            // ... to here (exclusively)
        x,                      // and compare with respect to this item
        Complex_realpart_less   // comparator using the real parts
    );
    
    // check if this is a valid knot
    if (iknot_notless != t_ /* x > 0 */ and iknot_notless != t_ + Nknot_ /* x < Rmax */)
        return iknot_notless - t_ - 1;
    else
        return -1;
}

Complex Bspline::eval(const cArrayView coeff, double x) const
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

Complex Bspline::eval(const cArrayView coeff, double x, double y) const
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
