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
    int k, int l, int ixknot, Real xmin, Real xmax,
    Bspline const & bspliney, GaussLegendre const & gy,
    int m, int n, int iyknot, Real ymin, Real ymax
) const
{
    // Compute integral of Bk(x) Bl(x) V(x,y) Bm(y) Bn(y) over two-dimensional interval
    //    xmin <= x <= xmax
    //    ymin <= y <= x
    // Assume that all coordinates are real.
    
    // number of quadrature points
    int points = std::max(2, bspline().order() + lambda + 1);
    
    // allocate workspace at once
    cArray workspace (8 * points);
    
    // get workspace pointers for individual data chunks
    Complex * xs  = &workspace[0] + points * 0;
    Complex * wxs = &workspace[0] + points * 1;
    Complex * Bxk = &workspace[0] + points * 2;
    Complex * Bxl = &workspace[0] + points * 3;
    Complex * ys  = &workspace[0] + points * 4;
    Complex * wys = &workspace[0] + points * 5;
    Complex * Bym = &workspace[0] + points * 6;
    Complex * Byn = &workspace[0] + points * 7;
    
    // restrict the x-integration interval to a subset, where y < x is possible
    xmin = std::max(xmin, ymin);
    
    // x evaluation points and evaluated B-splines
    gx.scaled_nodes_and_weights(points, xmin, xmax, xs, wxs);
    bsplinex.B(k, ixknot, points, xs, Bxk);
    bsplinex.B(l, ixknot, points, xs, Bxl);
    
    // resulting integral
    Complex result = 0;
    
    // for all x points
    for (int ixpoint = 0; ixpoint < points; ixpoint++)
    {
        Real x = xs[ixpoint].real();
        Real ytop = std::min(x, ymax);
        
        // get y evaluation points satisfying ymin < y < top
        gy.scaled_nodes_and_weights(points, ymin, ytop, ys, wys);
        bspliney.B(m, iyknot, points, ys, Bym);
        bspliney.B(n, iyknot, points, ys, Byn);
        
        // for all y points
        for (int iypoint = 0; iypoint < points; iypoint++)
        {
            // evalute the potential
            Real y = ys[iypoint].real();
            Complex V = special::pow_int(y / x, lambda) / x;
            
            // evalutate the integral
            result += Bxk[ixpoint] * Bxl[ixpoint] * Bym[iypoint] * Byn[iypoint] * wxs[ixpoint] * wys[iypoint] * V;
        }
    }
    
    // done
    return result;
}

void RadialIntegrals::coupledR (int lambda)
{
    std::size_t order    = bspline_x_.order();
    std::size_t Nxspline = bspline_x_.Nspline();
    std::size_t Nyspline = bspline_y_.Nspline();
    
    CooMatrix<LU_int_t,Complex> ucR (Nxspline * (order + 1), Nyspline * (order + 1));
    
    # pragma omp parallel for collapse (2)
    for (std::size_t a = 0; a < Nxspline; a++)
    for (std::size_t b = 0; b < Nyspline; b++)
    {
        for (std::size_t c = a; c < Nxspline and c <= a + order; c++)
        for (std::size_t d = b; d < Nyspline and d <= b + order; d++)
        {
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
            
            // get bounding knots of pair overlap
            Real t_xmin = special::clamp(std::max(ta1,tc1), rxmin_, rxmax_);
            Real t_xmax = special::clamp(std::min(ta2,tc2), rxmin_, rxmax_);
            Real t_ymin = special::clamp(std::max(tb1,td1), rymin_, rymax_);
            Real t_ymax = special::clamp(std::min(tb2,td2), rymin_, rymax_);
            
            // get bounding knots of four-overlap
            Real t_min = std::max(t_xmin, t_ymin);
            Real t_max = std::min(t_xmax, t_ymax);
            
            // is the integral coupled?
            if (t_min < t_max)
            {
                // resulting integral
                Complex R = 0;
                
                // integrate per cell
                for (unsigned ix = std::max(a,c); ix <= std::min(a,c) + order; ix++)
                for (unsigned iy = std::max(b,d); iy <= std::min(b,d) + order; iy++)
                {
                    // get integration cell
                    Complex tx1 = bspline_x_.t(ix    );  Real rx1 = tx1.real();
                    Complex ty1 = bspline_y_.t(iy    );  Real ry1 = ty1.real();
                    Complex tx2 = bspline_x_.t(ix + 1);  Real rx2 = tx2.real();
                    Complex ty2 = bspline_y_.t(iy + 1);  Real ry2 = ty2.real();
                    
                    // skip zero-length intervals
                    if (rx1 == rx2 or ry1 == ry2)
                        continue;
                    
                    // restrict effective radius to potential truncation bounds
                    Real ex1 = special::clamp(rx1, bspline_x_.R1(), bspline_x_.R2());
                    Real ey1 = special::clamp(ry1, bspline_y_.R1(), bspline_y_.R2());
                    Real ex2 = special::clamp(rx2, bspline_x_.R1(), bspline_x_.R2());
                    Real ey2 = special::clamp(ry2, bspline_y_.R1(), bspline_y_.R2());
                    
                    // shorthands for per-knot reduced integral moments
                    Complex const * const restrict Mi_L_ac    = Mitr_L_x(lambda).data()    + (a * (2*order+1) + c - (a-order)) * (order+1);
                    Complex const * const restrict Mi_mLm1_ac = Mitr_mLm1_x(lambda).data() + (a * (2*order+1) + c - (a-order)) * (order+1);
                    Complex const * const restrict Mi_L_bd    = Mitr_L_y(lambda).data()    + (b * (2*order+1) + d - (b-order)) * (order+1);
                    Complex const * const restrict Mi_mLm1_bd = Mitr_mLm1_y(lambda).data() + (b * (2*order+1) + d - (b-order)) * (order+1);
                    
                    // interval fully decoupled (x < y) ...?
                    if (ex2 <= ey1)
                    {
                        Real scale = special::pow_int(ex2 / ey2, lambda) / ey2;
                        R += Mi_L_ac[ix - a] * Mi_mLm1_bd[iy - b] * scale;
                    }
                    
                    // inteval fully decoupled (y < x) ...?
                    else if (ey2 <= ex1)
                    {
                        Real scale = special::pow_int(ey2 / ex2, lambda) / ex2;
                        R += Mi_L_bd[iy - b] * Mi_mLm1_ac[ix - a] * scale;
                    }
                    
                    // interval partially coupled (x ~ y) -> split to triangles
                    else
                    {
                        Complex tri1 = computeRtri(lambda, bspline_x_, g_x_, a, c, ix, rx1, rx2, bspline_y_, g_y_, b, d, iy, ry1, ry2);
                        Complex tri2 = computeRtri(lambda, bspline_y_, g_y_, b, d, iy, ry1, ry2, bspline_x_, g_x_, a, c, ix, rx1, rx2);
                        R += tri1 + tri2;
                    }
                }
                
                // add the calculated coupled integral
                # pragma omp critical
                ucR.add(a * (order + 1) + (c - a), b * (order + 1) + (d - b), R);
            }
        }
    }
    
    R_coupled_[lambda] = ucR.tocsr();
}

Complex RadialIntegrals::computeR (int lambda, int a, int b, int c, int d) const
{
    // shorthands
    int order = bspline_x_.order();
    
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
    
    // return clean zero if there are no pair overlaps
    if (std::max(ta1,tc1) >= std::min(ta2,tc2) or std::max(tb1,td1) >= std::min(tb2,td2))
        return 0.;
    
    // effective boundaries
    Real t_xmin = special::clamp(std::max(ta1,tc1), rxmin_, rxmax_);
    Real t_xmax = special::clamp(std::min(ta2,tc2), rxmin_, rxmax_);
    Real t_ymin = special::clamp(std::max(tb1,td1), rymin_, rymax_);
    Real t_ymax = special::clamp(std::min(tb2,td2), rymin_, rymax_);
    
    // decoupled y < x
    if (t_ymax <= t_xmin)
    {
        Real scale = special::pow_int(t_ymax / t_xmax, lambda) / t_xmax;
        return scale * Mtr_mLm1_x_[lambda](a,c) * Mtr_L_y_[lambda](b,d);
    }
    
    // decoupled x < y
    if (t_xmax <= t_ymin)
    {
        Real scale = special::pow_int(t_xmax / t_ymax, lambda) / t_ymax;
        return scale * Mtr_L_x_[lambda](a,c) * Mtr_mLm1_y_[lambda](b,d);
    }
    
    // coupled
    return R_coupled_[lambda]
    (
        std::min(a, c) * (order + 1) + std::abs(a - c),
        std::min(b, d) * (order + 1) + std::abs(b - d)
    );
}
