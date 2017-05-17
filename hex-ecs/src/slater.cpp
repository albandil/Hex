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

#include <cmath>
#include <complex>

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"
#include "hex-special.h"

// --------------------------------------------------------------------------------- //

#include "bspline.h"
#include "gauss.h"
#include "radial.h"

// --------------------------------------------------------------------------------- //

Complex RadialIntegrals::computeRtri
(
    int lambda,
    Bspline const & bsplinex, GaussLegendre const & gx,
    int k, int l, int ixknot, Complex xmin, Complex xmax,
    Bspline const & bspliney, GaussLegendre const & gy,
    int m, int n, int iyknot, Complex ymin, Complex ymax
) const
{
    // Compute integral of Bk(1) Bl(1) V(1,2) Bm(2) Bn(2) over two-dimensional interval
    //    xmin <= x <= xmax
    //    ymin <= y <= x
    
    // number of quadrature points
    int points = std::max(2, bspline().order() + lambda + 1);
    
    // allocate workspace at once
    cArray workspace (8 * points);
    
    // get workspace pointers for individual data
    Complex * xs  = &workspace[0] + points * 0;
    Complex * wxs = &workspace[0] + points * 1;
    Complex * Bxk = &workspace[0] + points * 2;
    Complex * Bxl = &workspace[0] + points * 3;
    Complex * ys  = &workspace[0] + points * 4;
    Complex * wys = &workspace[0] + points * 5;
    Complex * Bym = &workspace[0] + points * 6;
    Complex * Byn = &workspace[0] + points * 7;
    
    // restrict the x-integration interval to a subset, where y < x is possible
    xmin = bsplinex.rotate(std::max(bsplinex.unrotate(xmin), bspliney.unrotate(ymin)));
    
    // x evaluation points and evaluated B-splines
    gx.scaled_nodes_and_weights(points, xmin, xmax, xs, wxs);
    bsplinex.B(k, ixknot, points, xs, Bxk);
    bsplinex.B(l, ixknot, points, xs, Bxl);
    
    // resulting integral
    Complex result = 0;
    
    // for all x points
    for (int ixpoint = 0; ixpoint < points; ixpoint++)
    {
        // The integration over 'y' goes from 'ymin' to 'x'. But, generally,
        // 'y' and 'x' lie on distinct ECS contours. So we need, first,
        // to transform the current 'x' to the 'y' contour.
        Complex ytop = bspliney.rotate(std::min(bsplinex.unrotate(xs[ixpoint]), bspliney.unrotate(ymax)));
        
        // get y evaluation points satisfyng ymin < y < top
        gy.scaled_nodes_and_weights(points, ymin, ytop, ys, wys);
        bspliney.B(m, iyknot, points, ys, Bym);
        bspliney.B(n, iyknot, points, ys, Byn);
        
        // restrict effective x-radius
        Complex x = bsplinex.clamp(xs[ixpoint], rxmin_, rxmax_);
        
        // for all y points
        for (int iypoint = 0; iypoint < points; iypoint++)
        {
            // restrict the effective y-radius
            Complex y = bspliney.clamp(ys[iypoint], rymin_, rymax_);
            
            // evalute the potential
            Complex V = special::pow_int(y / x, lambda) / x;
            
            // evalutate the integral
            result += Bxk[ixpoint] * Bxl[ixpoint] * Bym[iypoint] * Byn[iypoint] * wxs[ixpoint] * wys[iypoint] * V;
        }
    }
    
    // done
    return result;
}

cArray RadialIntegrals::diagonalR (int lambda) const
{
    // assume bspline_atom == bspline_proj
    int order   = bspline_x_.order();
    int Nspline = bspline_x_.Nspline();
    
    // allocate space
    std::size_t O = order + 1;
    cArray R (Nspline * O * O * O, 0.);
    
    // calculate elements
    # pragma omp parallel for
    for (int a = 0; a < Nspline; a++)
    for (int b = a; b < Nspline and b <= a + order; b++)
    for (int c = a; c < Nspline and c <= a + order; c++)
    for (int d = a; d < Nspline and d <= a + order; d++)
    {
        R[((a * O + (b-a)) * O + (c-a)) * O + (d-a)] = computeR(lambda, a, b, c, d);
    }
    
    return R;
}

Complex RadialIntegrals::computeR (int lambda, int a, int b, int c, int d) const
{
    //
    // Preparations.
    //
    
        // shorthands
        int order = bspline().order();
        std::size_t O = order + 1;
        int e = mmin(a,b,c,d);
    
    //
    // This integral might be already precomputed.
    //
    
        if
        (
            lambda < (int)R_tr_dia_diag_.size() and
            not R_tr_dia_diag_[lambda].empty() and
            a - e <= order and b - e <= order and
            c - e <= order and d - e <= order
        )
        {
            if (a <= b and a <= c and a <= d)
                return R_tr_dia_diag_[lambda][((a * O + (b-a)) * O + (c-a)) * O + (d-a)];
            else if (b <= a and b <= c and b <= d)
                return R_tr_dia_diag_[lambda][((b * O + (a-b)) * O + (d-b)) * O + (c-b)];
            else if (c <= a and c <= b and c <= d)
                return R_tr_dia_diag_[lambda][((c * O + (b-c)) * O + (a-c)) * O + (d-c)];
            else // (d <= a and d <= b and d <= c)
                return R_tr_dia_diag_[lambda][((d * O + (a-d)) * O + (b-d)) * O + (c-d)];
        }
    
    //
    // Get bounding knots.
    //
    
        // leading knots of the B-splines
        Real ta1 = bspline_x_.unrotate(bspline_x_.t(a));
        Real tb1 = bspline_y_.unrotate(bspline_y_.t(b));
        Real tc1 = bspline_x_.unrotate(bspline_x_.t(c));
        Real td1 = bspline_y_.unrotate(bspline_y_.t(d));
        
        // trailing knots of the B-splines
        Real ta2 = bspline_x_.unrotate(bspline_x_.t(a + order + 1));
        Real tb2 = bspline_y_.unrotate(bspline_y_.t(b + order + 1));
        Real tc2 = bspline_x_.unrotate(bspline_x_.t(c + order + 1));
        Real td2 = bspline_y_.unrotate(bspline_y_.t(d + order + 1));
    
    //
    // Return clean zero if there are no pair overlaps or they have null measure.
    //
    
        if (std::max(ta1,tc1) >= std::min(ta2,tc2) or std::max(tb1,td1) >= std::min(tb2,td2))
            return 0.;
    
    // If the overlaps of the same-coordinate B-splines do not overlap
    // with each other, the four-B-spline integral decouples into a product
    // of two two-B-spline integrals. These have been already calculated.
    
        // the overlap of Bb,Bd precedes the overlap of Ba,Bc on the radial axis
        if (std::min(tb2,td2) <= std::max(ta1,tc1))
        {
            Real t_ac = special::clamp(std::min(ta2,tc2), rxmin_, rxmax_);
            Real t_bd = special::clamp(std::min(tb2,td2), rymin_, rymax_);
            Real scale = special::pow_int(t_bd / t_ac, lambda) / t_ac;
            return scale * Mtr_mLm1_x_[lambda](a,c) * Mtr_L_y_[lambda](b,d);
        }
        
        // the overlap of Ba,Bc precedes the overlap of Bb,Bd on the radial axis
        if (std::min(ta2,tc2) <= std::max(tb1,td1))
        {
            Real t_ac = special::clamp(std::min(ta2,tc2), rxmin_, rxmax_);
            Real t_bd = special::clamp(std::min(tb2,td2), rymin_, rymax_);
            Real scale = special::pow_int(t_ac / t_bd, lambda) / t_bd;
            return scale * Mtr_L_x_[lambda](a,c) * Mtr_mLm1_y_[lambda](b,d);
        }
    
    // The four-B-spline integrals that do not decouple into a product of two one-dimensional integrals need to be carefully calculated 
    // on a Cartesian product of the knot sequences, which transforms the full integral into a sum of integrals from several two-dimensional
    // cells. The integral does decouple, again, in such two-dimensional cells, where either r1 > r2 or r1 < r2 everywhere. This is called
    // the "off-diagonal" contribution to the R-integral. On the contrary, the integration cells where the notion of larger and smaller
    // coordinate changes from place to place contribute to the "diagonal" part. Computation of the diagonal part is computationally
    // more expensive than of the off-diagonal part.
    
        // result
        Complex R = 0;
        
        // shorthands for per-knot reduced integral moments
        Complex const * const restrict Mi_L_ac    = Mitr_L_x(lambda).data()    + (a * (2*order+1) + c - (a-order)) * (order+1);
        Complex const * const restrict Mi_mLm1_ac = Mitr_mLm1_x(lambda).data() + (a * (2*order+1) + c - (a-order)) * (order+1);
        Complex const * const restrict Mi_L_bd    = Mitr_L_y(lambda).data()    + (b * (2*order+1) + d - (b-order)) * (order+1);
        Complex const * const restrict Mi_mLm1_bd = Mitr_mLm1_y(lambda).data() + (b * (2*order+1) + d - (b-order)) * (order+1);
        
        // for all non-empty sub-intervals of the B-spline overlaps
        for (int ix = std::max(a,c); ix <= std::min(a,c) + order; ix++) if (bspline_x_.t(ix) != bspline_x_.t(ix + 1))
        for (int iy = std::max(b,d); iy <= std::min(b,d) + order; iy++) if (bspline_y_.t(iy) != bspline_y_.t(iy + 1))
        {
            // get integration rectangle
            Complex tx1 = bspline_x_.t(ix    );  Real rx1 = bspline_x_.unrotate(tx1);
            Complex ty1 = bspline_y_.t(iy    );  Real ry1 = bspline_y_.unrotate(ty1);
            Complex tx2 = bspline_x_.t(ix + 1);  Real rx2 = bspline_x_.unrotate(tx2);
            Complex ty2 = bspline_y_.t(iy + 1);  Real ry2 = bspline_y_.unrotate(ty2);
            
            // avoid negative radii (happens only in multi-domain calculation for some complex absorption grids)
            if (rx1 < 0 or rx2 < 0 or ry1 < 0 or ry2 < 0)
                continue;
            
            // restrict effective radius to potential truncation bounds
            Real rx = special::clamp(rx2, rxmin_, rxmax_);
            Real ry = special::clamp(ry2, rymin_, rymax_);
            
            // interval fully decoupled (x < y) ...?
            if (rx2 <= ry1)
            {
                Real scale = special::pow_int(rx / ry, lambda) / ry;
                R += Mi_L_ac[ix - a] * Mi_mLm1_bd[iy - b] * scale;
            }
            
            // inteval fully decoupled (y < x) ...?
            else if (ry2 <= rx1)
            {
                Real scale = special::pow_int(ry / rx, lambda) / rx;
                R += Mi_L_bd[iy - b] * Mi_mLm1_ac[ix - a] * scale;
            }
            
            // interval partially coupled (x ~ y)
            else
            {
                R += computeRtri(lambda, bspline_x_, g_x_, a, c, ix, tx1, tx2, bspline_y_, g_y_, b, d, iy, ty1, ty2);
                R += computeRtri(lambda, bspline_y_, g_y_, b, d, iy, ty1, ty2, bspline_x_, g_x_, a, c, ix, tx1, tx2);
            }
        }
    
    if (sqrabs(R) > 1e+6)
    {
        std::cout << "Warning: Radial integral too large!" << std::endl;
        std::cout << "         R[" << lambda << "](" << a << "," << b << "," << c << "," << d << ") = " << R << std::endl;
        std::cout << "         If this happens on panel in the DOM preconditioner, most likely the domains are too small." << std::endl;
    }
    
    return R;
}
