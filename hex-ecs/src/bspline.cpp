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

#include <algorithm>
#include <vector>

// --------------------------------------------------------------------------------- //

#ifdef _OPENMP
#include <omp.h>
#endif

// --------------------------------------------------------------------------------- //


#include "hex-arrays.h"
#include "hex-memory.h"
#include "hex-special.h"

// --------------------------------------------------------------------------------- //

#include "bspline.h"

// --------------------------------------------------------------------------------- //


const std::size_t Bspline::work_size_ = 1024;

// ----------------------------------------------------------------------- //
//  Recursive single-point B-spline evaluation                             //
// ----------------------------------------------------------------------- //


Complex Bspline::bspline (int i, int iknot, int k, Complex r) const
{
    // NOTE: The following bounds check is the caller's responsibility:
    // if (r <= t_[i] or t_[i + k] <= r)
    //    return 0.;

#ifdef _OPENMP
    unsigned ithread = omp_get_thread_num();
#else
    unsigned ithread = 0;
#endif

    // value of the parent B-splines of the requested B-spline
    Complex * const restrict b = reinterpret_cast<Complex*>(work_[ithread].data());

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

    Complex S1 = (t_[i+k]   != t_[i])   ? Real(k) / (t_[i+k] - t_[i]) : 0.;
    Complex S2 = (t_[i+k+1] != t_[i+1]) ? Real(k) / (t_[i+k+1] - t_[i+1]) : 0.;

    return A * S1 - B * S2;
}


// ----------------------------------------------------------------------- //
//  Evaluation of B-splines on a grid                                      //
// ----------------------------------------------------------------------- //

void Bspline::B (int i, int iknot, int M, Complex const * const restrict x, Complex * const restrict y) const
{
    // NOTE: The caller's responsibility is to check that all 'x' lie in the interval (t[iknot],t[iknot+1])
    //       and that the i-th B-spline is defined there.

#ifdef _OPENMP
    unsigned ithread = omp_get_thread_num();
#else
    unsigned ithread = 0;
#endif

    // calculate higher orders using real or complex path
    if (t_[i].imag() == 0 and t_[i + order_ + 1].imag() == 0)
    {
        // value of the parent B-splines of the requested B-spline
        Real * const restrict b = reinterpret_cast<Real*>(work_[ithread].data());

        // initialize zero-order B-splines
        for (int n = 0; n <= order_; n++)
        for (int m = 0; m < M; m++)
            b[n * M + m] = (i + n == iknot ? 1. : 0.);

        // evaluating real B-spline
        for (int ord = 1; ord <= order_; ord++)
        {
            // update splines
            for (int n = 0; n <= order_ - ord; n++)
            {
                Real invden1 = (t_[i+ord+n].real()   == t_[i+n].real()   ? 0.0_r : 1.0_r / (t_[i+ord+n].real()   - t_[i+n].real()));
                Real invden2 = (t_[i+ord+n+1].real() == t_[i+n+1].real() ? 0.0_r : 1.0_r / (t_[i+ord+n+1].real() - t_[i+n+1].real()));

                Real * const restrict pb  = b +  n      * M;
                Real * const restrict pbn = b + (n + 1) * M;

                // for all evaluation points
                # pragma omp simd
                for (int m = 0; m < M; m++)
                    pb[m] = pb[m] * (x[m].real() - t_[i+n].real()) * invden1 + pbn[m] * (t_[i+ord+n+1].real() - x[m].real()) * invden2;
            }
        }

        // return the collected value of the requested B-spline
        for (int m = 0; m < M; m++)
            y[m] = b[m];
    }
    else
    {
        // value of the parent B-splines of the requested B-spline
        Complex * const restrict b = reinterpret_cast<Complex*>(work_[ithread].data());

        // initialize zero-order B-splines
        for (int n = 0; n <= order_; n++)
        for (int m = 0; m < M; m++)
            b[n * M + m] = (i + n == iknot ? 1. : 0.);

        // general complex variant
        for (int ord = 1; ord <= order_; ord++)
        {
            // update splines
            for (int n = 0; n <= order_ - ord; n++)
            {
                Complex invden1 = (t_[i+ord+n]   == t_[i+n]   ? 0.0_r : 1.0_r / (t_[i+ord+n]   - t_[i+n]));
                Complex invden2 = (t_[i+ord+n+1] == t_[i+n+1] ? 0.0_r : 1.0_r / (t_[i+ord+n+1] - t_[i+n+1]));

                Complex * const restrict pb  = b +  n      * M;
                Complex * const restrict pbn = b + (n + 1) * M;

                // for all evaluation points
                for (int m = 0; m < M; m++)
                    pb[m] = pb[m] * (x[m] - t_[i+n]) * invden1 + pbn[m] * (t_[i+ord+n+1] - x[m]) * invden2;
            }
        }

        // return the collected value of the requested B-spline
        for (int m = 0; m < M; m++)
            y[m] = b[m];
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
                if (coeff[ispline] == 0.0_r)
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

cArray Bspline::zip
(
    Bspline const & bx,
    Bspline const & by,
    const cArrayView coeff,
    const rArrayView xgrid,
    const rArrayView ygrid,
    Complex (Bspline::* evalXBSpline) (int,int,int,Complex) const,
    Complex (Bspline::* evalYBSpline) (int,int,int,Complex) const
)
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
            evBx[ispline][ix - xleft] = (bx.*evalXBSpline)(ispline, iknot, bx.order_, *ix);
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
            evBy[ispline][iy - yleft] = (by.*evalYBSpline)(ispline, iknot, by.order_, *iy);
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

cArray Bspline::zip
(
    Bspline const & bx,
    Bspline const & by,
    Bspline const & bz,
    const cArrayView coeff,
    const rArrayView xgrid,
    const rArrayView ygrid,
    const rArrayView zgrid,
    Complex (Bspline::* evalXBSpline) (int,int,int,Complex) const,
    Complex (Bspline::* evalYBSpline) (int,int,int,Complex) const,
    Complex (Bspline::* evalZBSpline) (int,int,int,Complex) const
)
{
    // evaluated function
    cArray f (xgrid.size() * ygrid.size() * zgrid.size());

    // rotate the grids
    cArray x (xgrid.size()), y (ygrid.size()), z (zgrid.size());
    for (std::size_t i = 0; i < xgrid.size(); i++)
        x[i] = bx.rotate(xgrid[i]);
    for (std::size_t i = 0; i < ygrid.size(); i++)
        y[i] = by.rotate(ygrid[i]);
    for (std::size_t i = 0; i < zgrid.size(); i++)
        z[i] = bz.rotate(zgrid[i]);

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
            evBx[ispline][ix - xleft] = (bx.*evalXBSpline)(ispline, iknot, bx.order_, *ix);
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
            evBy[ispline][iy - yleft] = (by.*evalYBSpline)(ispline, iknot, by.order_, *iy);
        }
    }

    // evaluate B-splines (from "bz") on z-grid
    cArrays evBz (bz.Nspline_);
    cArray::const_iterator zleft = z.begin(), zright;
    for (int ispline = 0; ispline < bz.Nspline_; ispline++)
    {
        // get first z-coordinate within this spline's support
        while (zleft < z.end() and zleft->real() < bz.t_[ispline].real())
            zleft++;

        // get (one beyond the) last y-coordinate within this spline's support
        zright = zleft;
        while (zright < z.end() and zright->real() < bz.t_[ispline + bz.order_ + 1].real())
            zright++;

        // setup evaluation vector
        evBz[ispline].resize(zright - zleft);

        // evaluation interval
        int iknot = ispline;

        // evaluate at y[]
        for (cArray::const_iterator iz = zleft; iz != zright; iz++)
        {
            // increment knot to match current evaluation point
            while ( iz->real() > bz.t_[iknot + 1].real() )
                if (++iknot >= bz.Nknot_)
                    HexException("Some evaluation points are outside of grid.");

            // evaluate this spline at the current evaluation point
            evBz[ispline][iz - zleft] = (bz.*evalZBSpline)(ispline, iknot, bz.order_, *iz);
        }
    }

    // zip triplet expansion
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

            // for all z-splines
            zleft = z.begin();
            for (int izspline = 0; izspline < bz.Nspline_; izspline++)
            {
                // get first z-coordinate within this spline's support
                while (zleft < z.end() and zleft->real() < bz.t_[izspline].real())
                    zleft++;

                // get (one beyond the) last z-coordinate within this spline's support
                zright = zleft;
                while (zright < z.end() and zright->real() < bz.t_[izspline + bz.order_ + 1].real())
                    zright++;

                // get relevant evaluations
                cArray const & Bz_row = evBz[izspline];

                // get coefficient of the expansion
                Complex C = coeff[(ixspline * by.Nspline_ + iyspline) * by.Nspline_ + izspline];

                // loop over relevant points
                for (cArray::const_iterator ix = xleft; ix != xright; ix++)
                for (cArray::const_iterator iy = yleft; iy != yright; iy++)
                for (cArray::const_iterator iz = zleft; iz != zright; iz++)
                {
                    // get position indices within various arrays
                    std::size_t idx = ix - x.begin(), evx = ix - xleft;
                    std::size_t idy = iy - y.begin(), evy = iy - yleft;
                    std::size_t idz = iz - z.begin(), evz = iz - zleft;

                    // update the evaluation
                    f[(idz * y.size() + idy) * x.size() + idx] += C * Bx_row[evx] * By_row[evy] * Bz_row[evz];
                }
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


Bspline::Bspline
(
    int order,
    Real th,
    rArrayView cknots1,
    rArrayView rknots,
    rArrayView cknots2
)   : rknots_(rknots),
      cknots1_(cknots1),
      cknots2_(cknots2),
      t_(),
      theta_(th),
      rotation_(Complex(cos(th),sin(th))),
      R1_(rknots_.empty() ? 0 : rknots_.front()),
      R2_(rknots_.empty() ? 0 : rknots_.back()),
      Rmin_(cknots1_.empty() ? R1_ : cknots1_.front()),
      Rmax_(cknots2_.empty() ? R2_ : cknots2_.back()),
      Nknot_(rknots.size() + (cknots1.empty() ? 0 : cknots1.size() - 1) + (cknots2.empty() ? 0 : cknots2.size() - 1)),
      Nreknot_(rknots.size()),
      Nspline_(Nknot_ > order + 1 ? Nknot_ - order - 1 : 0),
      Nintval_(Nknot_ > 1 ? Nknot_ - 1 : 0),
      order_(order)
{
    // check connection of knots
    if (not cknots1.empty() and not rknots.empty() and cknots1.back() != rknots.front())
        HexException("Last knot of leading complex knots must be equal to the first real knot.");
    if (not cknots2.empty() and not rknots.empty() and cknots2.front() != rknots.back())
        HexException("First knot of trailing complex knots must be equal to the last real knot.");
    if (rknots.empty() and not cknots1.empty() and not cknots2.empty() and cknots1.back() != cknots2.front())
        HexException("First knot of trailing complex knots must be equal to the last knot of leading complex knots.");

    // merge knots
    rArray knots;
    knots.append(cknots1);
    if (not cknots1.empty()) knots.pop_back();
    knots.append(rknots);
    if (not rknots.empty()) knots.pop_back();
    knots.append(cknots2);

    // rotate knots
    for (Real t : knots)
        t_.push_back(rotate(t));

    // allocate workspace for all threads
    unsigned nthreads = 1;
#ifdef _OPENMP
    # pragma omp parallel
    # pragma omp master
    nthreads = omp_get_num_threads();
#endif
    while (work_.size() <= nthreads)
        work_.push_back(rArray(work_size_));
}


// ----------------------------------------------------------------------- //
//  Others                                                                 //
// ----------------------------------------------------------------------- //


int Bspline::knot (Complex x) const
{
    // empty grid contains no knots
    if (t_.empty())
        return -1;

    // boundary search first so that we don't need to consider it below
    if (x.real() == t_.front().real()) return 0;
    if (x.real() == t_.back().real())  return Nknot_ - 1;

    // is the point inside the grid at all?
    if (x.real() < t_.front().real() or t_.back().real() < x.real())
        return -1;

    // the point 'x' is now strictly inside the knot sequence 't_'
    // - let's find the first greater knot by bisection
    // - the wanted knot index is then one position to the left
    return std::upper_bound
    (
        t_.begin(),             // search from here ...
        t_.end(),               // ... to here (exclusively)
        x,                      // and compare with respect to this item
        Complex_realpart_less   // comparator using the real parts
    ) - t_.begin() - 1;
}

Complex Bspline::eval (const cArrayView coeff, Real x) const
{
    // map real coordinate to complex contour
    Complex z = rotate(x);

    // get knot index
    int iknot = knot(z);

    // evaluate the B-spline expansion
    Complex result = 0.;
    for (int ispline = std::max(0, iknot - order_); ispline <= std::min(iknot, Nspline_ - 1); ispline++)
        result += coeff[ispline] * bspline(ispline, iknot, order_, z);

    return result;
}

Complex Bspline::eval (const cArrayView coeff, Real x, Real y) const
{
    Complex w = rotate(x);
    Complex z = rotate(y);

    // get knot indices
    int xknot = knot(w);
    int yknot = knot(z);

    // get bounding B-splines
    int leftxspline = xknot - order_;
    int rightxspline = xknot;
    int leftyspline = yknot - order_;
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

std::size_t Bspline::hash () const
{
    std::size_t seed = 0;

    for (auto & i : cknots1_)
        seed ^= *reinterpret_cast<typeinfo<Real>::inttype const*>(&i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);

    for (auto & i : rknots_)
        seed ^= *reinterpret_cast<typeinfo<Real>::inttype const*>(&i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);

    for (auto & i : cknots2_)
        seed ^= *reinterpret_cast<typeinfo<Real>::inttype const*>(&i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);

    for (unsigned i = 0; i < sizeof(theta_); i++)
        seed ^= *(i + (char*)&theta_) + 0x9e3779b9 + (seed << 6) + (seed >> 2);

    return seed % 65536;
}

Complex Bspline::clamp (Complex z, Real a, Real b) const
{
    return rotate
    (
        special::clamp
        (
            unrotate(z),
            a, b
        )
    );
}
