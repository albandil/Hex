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
    
    // x evaluation points and evaluated B-splines
    gx.scaled_nodes_and_weights(points, xmin, xmax, xs, wxs);
    bsplinex.B(k, ixknot, points, xs, Bxk);
    bsplinex.B(l, ixknot, points, xs, Bxl);
    
    // resulting integral
    Complex result = 0;
    
//     std::cout << "        xmin = " << xmin << ", xmax = " << xmax << ", x = " << cArrayView(points, xs) << std::endl;
//     std::cout << "        ymin = " << ymin << ", ymax = " << ymax << std::endl;
    
    // for all x points
    for (int ixpoint = 0; ixpoint < points; ixpoint++)
    {
        // The integration over 'y' goes from 'ymin' to 'x'. But, generally,
        // 'y' and 'x' lie on distinct ECS contours. So we need, first,
        // to transform the current 'x' to the 'y' contour.
        ymax = bspliney.rotate(bsplinex.unrotate(xs[ixpoint]));
        
        // get y evaluation points satisfyng ymin < y < ytop
        gy.scaled_nodes_and_weights(points, ymin, ymax, ys, wys);
        bspliney.B(m, iyknot, points, ys, Bym);
        bspliney.B(n, iyknot, points, ys, Byn);
        
//         std::cout << "        ymin = " << ymin << ", ymax = " << ymax << ", y = " << cArrayView(points, ys) << std::endl;
        
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

Complex RadialIntegrals::computeR (int lambda, int a, int b, int c, int d) const
{
    //
    // Preparations.
    //
    
        // shorthands
        int order = bspline().order();
        
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
        
        // shorthands
        Complex const * const restrict Mi_L_ac    = Mitr_L_x(lambda).data()    + (a * (2*order+1) + c - (a-order)) * (order+1);
        Complex const * const restrict Mi_mLm1_ac = Mitr_mLm1_x(lambda).data() + (a * (2*order+1) + c - (a-order)) * (order+1);
        Complex const * const restrict Mi_L_bd    = Mitr_L_y(lambda).data()    + (b * (2*order+1) + d - (b-order)) * (order+1);
        Complex const * const restrict Mi_mLm1_bd = Mitr_mLm1_y(lambda).data() + (b * (2*order+1) + d - (b-order)) * (order+1);
        
//         std::cout << "R[" << lambda << "](" << a << "," << b << "," << c << "," << d << ")" << std::endl;
        
        // for all knot pairs of the Cartesian product of same-coordinate overlaps
        for (int ix = std::max(a,c); ix <= std::min(a,c) + order; ix++) if (bspline_x_.t(ix) != bspline_x_.t(ix + 1))
        for (int iy = std::max(b,d); iy <= std::min(b,d) + order; iy++) if (bspline_y_.t(iy) != bspline_y_.t(iy + 1))
        {
            // get position of the rectangle to integrate
            Complex tx1 = bspline_x_.t(ix    );
            Complex ty1 = bspline_y_.t(iy    );
            Complex tx2 = bspline_x_.t(ix + 1);
            Complex ty2 = bspline_y_.t(iy + 1);
            
            // get effective scaling distances
            Real rx = special::clamp(bspline_x_.unrotate(tx2), rxmin_, rxmax_);
            Real ry = special::clamp(bspline_y_.unrotate(ty2), rymin_, rymax_);
            
            // decoupled [a,c] << [b,d]
            if (bspline_x_.unrotate(tx2) <= bspline_y_.unrotate(ty1))
            {
                Real scale = special::pow_int(rx / ry, lambda) / ry;
                R += Mi_L_ac[ix - a] * Mi_mLm1_bd[iy - b] * scale;
//                 std::cout << "    1 : " << ix << " " << iy << " " << R << std::endl;
            }
            
            // decoupled [b,d] << [a,c]
            else if (bspline_y_.unrotate(ty2) <= bspline_x_.unrotate(tx1))
            {
                Real scale = special::pow_int(ry / rx, lambda) / rx;
                R += Mi_L_bd[iy - b] * Mi_mLm1_ac[ix - a] * scale;
//                 std::cout << "    2 : " << ix << " " << iy << " " << R << std::endl;
            }
            
            // coupled [a,c] ~ [b,d]
            else
            {
                R += computeRtri(lambda, bspline_x_, g_x_, a, c, ix, tx1, tx2, bspline_y_, g_y_, b, d, iy, ty1, ty2);
//                 std::cout << "    3a : " << ix << " " << iy << " " << R << std::endl;
                R += computeRtri(lambda, bspline_y_, g_y_, b, d, iy, ty1, ty2, bspline_x_, g_x_, a, c, ix, tx1, tx2);
//                 std::cout << "    3b : " << ix << " " << iy << " " << R << std::endl;
            }
        }
    
    if (std::abs(R) > 1000)
    {
        std::cout << "Abort!" << std::endl;
        std::exit(0);
    }
    return R;
}
