//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2016, Jakub Benda, Charles University in Prague                    //
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
#include <complex>
#include <cmath>
#include <cstdio>
#include <cstring>

#include "hex-arrays.h"
#include "hex-special.h"
#include "hex-symbandmatrix.h"

#include "bspline.h"
#include "gauss.h"
#include "parallel.h"
#include "radial.h"

RadialIntegrals::RadialIntegrals
(
    Bspline const & bspline_inner,
    Bspline const & bspline_full,
    int Nlambdas
) : RadialIntegrals
(
    bspline_inner, bspline_full,
    bspline_inner, bspline_full,
    Nlambdas
)
{
    // nothing
}

RadialIntegrals::RadialIntegrals
(
    Bspline const & bspline_inner_x,
    Bspline const & bspline_full_x,
    Bspline const & bspline_inner_y,
    Bspline const & bspline_full_y,
    int Nlambdas
)
  : bspline_inner_x_(bspline_inner_x), bspline_full_x_ (bspline_full_x),
    bspline_inner_y_(bspline_inner_y), bspline_full_y_ (bspline_full_y),
    D_inner_x_     (bspline_inner_x.Nspline(), bspline_inner_x.order() + 1),
    S_inner_x_     (bspline_inner_x.Nspline(), bspline_inner_x.order() + 1),
    Mm1_inner_x_   (bspline_inner_x.Nspline(), bspline_inner_x.order() + 1),
    Mm1_tr_inner_x_(bspline_inner_x.Nspline(), bspline_inner_x.order() + 1),
    Mm2_inner_x_   (bspline_inner_x.Nspline(), bspline_inner_x.order() + 1),
    D_inner_y_     (bspline_inner_y.Nspline(), bspline_inner_y.order() + 1),
    S_inner_y_     (bspline_inner_y.Nspline(), bspline_inner_y.order() + 1),
    Mm1_inner_y_   (bspline_inner_y.Nspline(), bspline_inner_y.order() + 1),
    Mm1_tr_inner_y_(bspline_inner_y.Nspline(), bspline_inner_y.order() + 1),
    Mm2_inner_y_   (bspline_inner_y.Nspline(), bspline_inner_y.order() + 1),
    D_full_x_      (bspline_full_x .Nspline(), bspline_full_x .order() + 1),
    S_full_x_      (bspline_full_x .Nspline(), bspline_full_x .order() + 1),
    Mm1_full_x_    (bspline_full_x .Nspline(), bspline_full_x .order() + 1),
    Mm1_tr_full_x_ (bspline_full_x .Nspline(), bspline_full_x .order() + 1),
    Mm2_full_x_    (bspline_full_x .Nspline(), bspline_full_x .order() + 1),
    D_full_y_      (bspline_full_y .Nspline(), bspline_full_y .order() + 1),
    S_full_y_      (bspline_full_y .Nspline(), bspline_full_y .order() + 1),
    Mm1_full_y_    (bspline_full_y .Nspline(), bspline_full_y .order() + 1),
    Mm1_tr_full_y_ (bspline_full_y .Nspline(), bspline_full_y .order() + 1),
    Mm2_full_y_    (bspline_full_y .Nspline(), bspline_full_y .order() + 1),
    verbose_(true),
    Nlambdas_(Nlambdas)
{
    // maximal number of evaluation points (quadrature rule)
    int npts = std::max(EXPANSION_QUADRATURE_POINTS, bspline_inner_x_.order() + Nlambdas + 1);
    
    // precompute Gaussian weights
    g_inner_x_.precompute_nodes_and_weights(npts);
    g_inner_y_.precompute_nodes_and_weights(npts);
    g_full_x_ .precompute_nodes_and_weights(npts);
    g_full_y_ .precompute_nodes_and_weights(npts);
}

cArray RadialIntegrals::computeMi (Bspline const & bspline, GaussLegendre const & g, int a) const
{
    int Nspline = bspline.Nspline();
    int order = bspline.order();
    int Nknot = bspline.Nknot();
    
    // partial integral moments
    cArray m (Nspline * (2 * order + 1) * (order + 1));
    
    // quadrature point count (use at least 2nd order)
    int points = std::max(2, order + std::abs(a) + 1);
    
    // get quadrature points and weights
    cArray xs ((Nknot - 1) * points), ws ((Nknot - 1) * points);
    # pragma omp parallel for
    for (int iknot = 0; iknot < Nknot - 1; iknot++)
    if (bspline.t(iknot).real() != bspline.t(iknot + 1).real())
    {
        g.scaled_nodes_and_weights
        (
            points,
            bspline.t(iknot),
            bspline.t(iknot + 1),
            xs.data() + iknot * points,
            ws.data() + iknot * points
        );
    }
    
    // calculate (scaled) powers
    cArray x_a ((Nknot - 1) * points);
    # pragma omp parallel for
    for (int iknot = 0; iknot < Nknot - 1; iknot++)
    if (bspline.t(iknot).real() != bspline.t(iknot + 1).real())
    for (int ipoint = 0; ipoint < points; ipoint++)
    {
        x_a[iknot * points + ipoint] = special::pow_int(xs[iknot * points + ipoint] / bspline.t(iknot + 1).real(), a);
    }
    
    // damp weights
    # pragma omp parallel for
    for (std::size_t ix = 0; ix < xs.size(); ix++)
    {
        if (xs[ix].real() < bspline.R1() or bspline.R2() < xs[ix].real())
            ws[ix] = 0;
        else
            ws[ix] *= std::tanh(0.125 * (bspline.R2() - xs[ix].real()));
    }
    
    // evaluate all B-splines in all points
    cArray B_x (Nspline * (order + 1) * points);
    # pragma omp parallel for
    for (int ispline = 0; ispline < Nspline; ispline++)
    for (int iknot = ispline; iknot <= ispline + order; iknot++)
    if (bspline.t(iknot).real() != bspline.t(iknot + 1).real())
    {
        bspline.B
        (
            ispline, iknot, points,
            xs.data() + iknot * points,
            B_x.data() + (ispline * (order + 1) + (iknot - ispline)) * points
        );
    }
    
    # pragma omp parallel for
    for (int i = 0; i < Nspline; i++)
    for (int j = i; j <= i + order and j < Nspline; j++)
    for (int iknot = j; iknot <= i + order; iknot++)
    if (bspline.t(iknot).real() != bspline.t(iknot + 1).real())
    {
        // results of the quadrature
        Complex integral = 0;
        
        // sum the quadrature rule
        for (int ipoint = 0; ipoint < points; ipoint++)
        {
            Complex Bi = B_x[(i * (order + 1) + (iknot - i)) * points + ipoint];
            Complex Bj = B_x[(j * (order + 1) + (iknot - j)) * points + ipoint];
            Complex xa = x_a[iknot * points + ipoint];
            Complex wx = ws[iknot * points + ipoint];
            
            integral += Bi * xa * Bj * wx;
        }
        
        // get the coordinates in m-matrix
        int x_1 = i;                // reference spline is i-th
        int y_1 = j - (i - order);
        int z_1 = iknot - i;
        
        // get the coordinates in m-matrix of the symmetric case
        int x_2 = j;                // reference spline is j-th
        int y_2 = i - (j - order);
        int z_2 = iknot - j;
        
        // save to m-matrix
        m[(x_1 * (2 * order + 1) + y_1) * (order + 1) + z_1] = integral;
        m[(x_2 * (2 * order + 1) + y_2) * (order + 1) + z_2] = integral;
    }
    
    return m;
}

Complex RadialIntegrals::computeD_iknot
(
    Bspline const & bspline,
    GaussLegendre const & g,
    int i, int j, int iknot
) const
{
    if (iknot < 0)
        iknot = bspline.Nknot() - 1;
    
    // get interval boundaries
    Complex x1 = bspline.t(iknot);
    Complex x2 = bspline.t(iknot + 1);
    
    // throw away zero-length intervals
    if (x1 == x2)
        return 0;
    
    // get Gauss-Legendre nodes and weights for the interval [-1, 1]
    // - use at least 2nd order
    int points = std::max(2, bspline.order());
    cArray xs(points), ws(points);
    g.scaled_nodes_and_weights(points, x1, x2, xs.data(), ws.data());
    
    // evaluate B-splines at Gauss-Legendre nodes
    cArray values_i(points), values_j(points);
    bspline.dB(i, iknot, points, xs.data(), values_i.data());
    bspline.dB(j, iknot, points, xs.data(), values_j.data());
    
    // result
    Complex res = 0;
    
    // accumulate the result
    for (int k = 0; k < points; k++)
        res += values_i[k] * values_j[k] * ws[k];
    
    return res;
}

Complex RadialIntegrals::computeD
(
    Bspline const & bspline,
    GaussLegendre const & g,
    int i, int j
) const
{
    // get boundary iknots
    int left = std::max(i, j);
    int right = std::min(i, j) + bspline.order();
    
    // the result
    Complex res = 0;
    
    // undergo integration on sub-intervals
    for (int iknot = left; iknot <= right; iknot++)
        res += computeD_iknot(bspline, g, i, j, iknot);
        
    return res;
}

Complex RadialIntegrals::computeOverlapMatrixElement_iknot
(
    Bspline const & bspline, GaussLegendre const & g,
    int i, int j,
    std::function<Complex(Complex)> weight,
    int iknot
) const
{
    // get interval boundaries
    Complex x1 = bspline.t(iknot);
    Complex x2 = bspline.t(iknot + 1);
    
    // throw away zero-length intervals
    if (x1 == x2)
        return 0;
    
    // get Gauss-Legendre nodes and weights for the interval [-1, 1]
    // - use at least 2nd order
    int points = EXPANSION_QUADRATURE_POINTS;
    cArray xs(points), ws(points);
    g.scaled_nodes_and_weights(points, x1, x2, xs.data(), ws.data());
    
    // evaluate B-splines at Gauss-Legendre nodes
    cArray values_i(points), values_j(points);
    bspline.B(i, iknot, points, xs.data(), values_i.data());
    bspline.B(j, iknot, points, xs.data(), values_j.data());
    
    // result
    Complex res = 0;
    
    // accumulate the result
    for (int k = 0; k < points; k++)
        res += values_i[k] * values_j[k] * weight(xs[k]) * ws[k];
    
    return res;
}

Complex RadialIntegrals::computeOverlapMatrixElement
(
    Bspline const & bspline, GaussLegendre const & g,
    int i, int j,
    std::function<Complex(Complex)> weight
) const
{
    // get boundary iknots
    int left = std::max(i, j);
    int right = std::min(i, j) + bspline.order() + 1;
    
    // the result
    Complex res = 0;
    
    // undergo integration on sub-intervals
    for (int iknot = left; iknot < right; iknot++)
        res += computeOverlapMatrixElement_iknot(bspline, g, i, j, weight, iknot);
    
    return res;
}

SymBandMatrix<Complex> RadialIntegrals::computeOverlapMatrix
(
    Bspline const & bspline,
    std::function<Complex(Complex)> weight
) const
{
    GaussLegendre g;
    g.precompute_nodes_and_weights(EXPANSION_QUADRATURE_POINTS);
    
    return SymBandMatrix<Complex>(bspline.Nspline(), bspline.order() + 1).populate
    (
        [&](int i, int j)
        {
            return computeOverlapMatrixElement(bspline, g, i, j, weight);
        }
    );
}

Complex RadialIntegrals::computeM_iknot
(
    Bspline const & bspline, GaussLegendre const & g,
    int a, int i, int j,
    int iknot, CCFunction weight, Real scale
) const
{
    // get interval boundaries
    Complex x1 = bspline.t(iknot);
    Complex x2 = bspline.t(iknot + 1);
    
    // throw away zero-length intervals
    if (x1 == x2)
        return 0;
    
    // get Gauss-Legendre nodes and weights for the interval [-1, 1]
    // - use at least 2nd order
    int points = std::max(2, bspline.order() + std::abs(a) + 1);
    cArray xs(points), ws(points);
    g.scaled_nodes_and_weights(points, x1, x2, xs.data(), ws.data());
    
    // evaluate B-splines at Gauss-Legendre nodes
    cArray values_i(points), values_j(points);
    bspline.B(i, iknot, points, xs.data(), values_i.data());
    bspline.B(j, iknot, points, xs.data(), values_j.data());
    
    // result
    Complex res = 0;
    
    // accumulate the (weighted) result
    for (int k = 0; k < points; k++)
        res += weight(xs[k]) * values_i[k] * values_j[k] * special::pow_int(scale * xs[k], a) * ws[k];
    
    return res;
}

Complex RadialIntegrals::computeM
(
    Bspline const & bspline, GaussLegendre const & g,
    int a, int i, int j,
    bool truncate, bool scale
) const
{
    // get boundary iknots
    int left = std::max(i, j);
    int right = std::min(i, j) + bspline.order() + 1;
    
    // truncating function
    auto weight = [&](Complex x) -> Complex
    {
        if (not truncate)
            return 1.0_z;
        else if (x.real() < bspline.R1() or bspline.R2() < x.real())
            return 0.0_z;
        else
            return std::tanh(0.125 * (bspline.R2() - x.real()));
        
    };
    
    // calculate scaling factor
    Real scalefactor = (scale ? 1.0_r / bspline.t(right).real() : 1.0_r);
    
    // the result
    Complex res = 0;
    
    // undergo integration on sub-intervals
    for (int iknot = left; iknot < right; iknot++)
        res += computeM_iknot(bspline, g, a, i, j, iknot, weight, scalefactor);
    
    return res;
}

void RadialIntegrals::setupOneElectronIntegrals (Parallel const & par, CommandLine const & cmd)
{
    setupOneElectronIntegrals(cmd.shared_scratch, par.IamMaster());
}

void RadialIntegrals::setupOneElectronIntegrals (bool shared_scratch, bool IamMaster)
{
    // inner basis one-electron integrals
    #define SetupOneElectronIntegrals(REG,AXIS) \
        std::size_t hash_##REG##_##AXIS = bspline_##REG##_##AXIS##_.hash(); \
        \
        D_##REG##_##AXIS##_     .hdflink(format("rad-D-%.4lx.hdf",      hash_##REG##_##AXIS)); \
        S_##REG##_##AXIS##_     .hdflink(format("rad-S-%.4lx.hdf",      hash_##REG##_##AXIS)); \
        Mm1_##REG##_##AXIS##_   .hdflink(format("rad-Mm1-%.4lx.hdf",    hash_##REG##_##AXIS)); \
        Mm1_tr_##REG##_##AXIS##_.hdflink(format("rad-Mm1_tr-%.4lx.hdf", hash_##REG##_##AXIS)); \
        Mm2_##REG##_##AXIS##_   .hdflink(format("rad-Mm2-%.4lx.hdf",    hash_##REG##_##AXIS)); \
        \
        D_##REG##_##AXIS##_     .populate([=](int m, int n) -> Complex { return computeD(bspline_##REG##_##AXIS##_, g_##REG##_##AXIS##_,     m, n              ); }); \
        S_##REG##_##AXIS##_     .populate([=](int m, int n) -> Complex { return computeM(bspline_##REG##_##AXIS##_, g_##REG##_##AXIS##_,  0, m, n, false, false); }); \
        Mm1_##REG##_##AXIS##_   .populate([=](int m, int n) -> Complex { return computeM(bspline_##REG##_##AXIS##_, g_##REG##_##AXIS##_, -1, m, n, false, false); }); \
        Mm1_tr_##REG##_##AXIS##_.populate([=](int m, int n) -> Complex { return computeM(bspline_##REG##_##AXIS##_, g_##REG##_##AXIS##_, -1, m, n, true,  false); }); \
        Mm2_##REG##_##AXIS##_   .populate([=](int m, int n) -> Complex { return computeM(bspline_##REG##_##AXIS##_, g_##REG##_##AXIS##_, -2, m, n, false, false); }); \
        \
        if (not shared_scratch or IamMaster) \
        { \
            D_##REG##_##AXIS##_     .hdfsave(); \
            S_##REG##_##AXIS##_     .hdfsave(); \
            Mm1_##REG##_##AXIS##_   .hdfsave(); \
            Mm1_tr_##REG##_##AXIS##_.hdfsave(); \
            Mm2_##REG##_##AXIS##_   .hdfsave(); \
        }
    
    if (verbose_)
        std::cout << "Precomputing one-electron integrals ... " << std::flush;
    
    Timer t;
    
    SetupOneElectronIntegrals(inner,x)
    SetupOneElectronIntegrals(inner,y)
    SetupOneElectronIntegrals(full,x)
    SetupOneElectronIntegrals(full,y)
    
    if (verbose_)
        std::cout << "done in " << t.nice_time() << std::endl << std::endl;
}

void RadialIntegrals::setupTwoElectronIntegrals (Parallel const & par, CommandLine const & cmd)
{
    // shorthands
    int order = bspline_inner_x_.order();
    int Nspline_inner_x = bspline_inner_x_.Nspline();
    int Nspline_inner_y = bspline_inner_y_.Nspline();
    int Nspline_full_x  = bspline_full_x_ .Nspline();
    int Nspline_full_y  = bspline_full_y_ .Nspline();
    
    // set number of two-electron integrals
    R_tr_dia_.resize(Nlambdas_);
    
    if (verbose_)
        std::cout << "Precomputing partial integral moments (lambda = 0 .. " << Nlambdas_ - 1 << ")" << std::endl;
    
    // partial moments sizes
    std::size_t mi_size_inner_x = Nspline_inner_x * (2 * order + 1) * (order + 1);
    std::size_t mi_size_inner_y = Nspline_inner_y * (2 * order + 1) * (order + 1);
    std::size_t mi_size_full_x  = Nspline_full_x  * (2 * order + 1) * (order + 1);
    std::size_t mi_size_full_y  = Nspline_full_y  * (2 * order + 1) * (order + 1);
    
    // resize partial moments arrays
    Mitr_L_inner_x_   .resize(Nlambdas_ * mi_size_inner_x);  Mitr_L_full_x_   .resize(Nlambdas_ * mi_size_full_x); 
    Mitr_L_inner_y_   .resize(Nlambdas_ * mi_size_inner_y);  Mitr_L_full_y_   .resize(Nlambdas_ * mi_size_full_y); 
    Mitr_mLm1_inner_x_.resize(Nlambdas_ * mi_size_inner_x);  Mitr_mLm1_full_x_.resize(Nlambdas_ * mi_size_full_x);
    Mitr_mLm1_inner_y_.resize(Nlambdas_ * mi_size_inner_y);  Mitr_mLm1_full_y_.resize(Nlambdas_ * mi_size_full_y);
    
    // resize full moments arrays, initialize with empty matrices
    Mtr_L_inner_x_ .resize(Nlambdas_, SymBandMatrix<Complex>(Nspline_inner_x, order + 1));
    Mtr_L_inner_y_ .resize(Nlambdas_, SymBandMatrix<Complex>(Nspline_inner_y, order + 1));
    Mtr_L_full_x_  .resize(Nlambdas_, SymBandMatrix<Complex>(Nspline_full_x, order + 1));
    Mtr_L_full_y_  .resize(Nlambdas_, SymBandMatrix<Complex>(Nspline_full_y, order + 1));
    Mtr_mLm1_inner_x_ .resize(Nlambdas_, SymBandMatrix<Complex>(Nspline_inner_x, order + 1));
    Mtr_mLm1_inner_y_ .resize(Nlambdas_, SymBandMatrix<Complex>(Nspline_inner_y, order + 1));
    Mtr_mLm1_full_x_  .resize(Nlambdas_, SymBandMatrix<Complex>(Nspline_full_x, order + 1));
    Mtr_mLm1_full_y_  .resize(Nlambdas_, SymBandMatrix<Complex>(Nspline_full_y, order + 1));
    
    // resize two-electron integrals array
    R_tr_dia_diag_.resize(Nlambdas_);
    
    // for all multipole moments
    for (int lambda = 0; lambda < Nlambdas_; lambda++)
    {
        // calculate inner basis partial integral moments
        cArrayView(Mitr_L_inner_x_,    lambda * mi_size_inner_x, mi_size_inner_x) = computeMi(bspline_inner_x_, g_inner_x_,   lambda  );
        cArrayView(Mitr_L_inner_y_,    lambda * mi_size_inner_y, mi_size_inner_y) = computeMi(bspline_inner_y_, g_inner_y_,   lambda  );
        cArrayView(Mitr_mLm1_inner_x_, lambda * mi_size_inner_x, mi_size_inner_x) = computeMi(bspline_inner_x_, g_inner_x_,  -lambda-1);
        cArrayView(Mitr_mLm1_inner_y_, lambda * mi_size_inner_y, mi_size_inner_y) = computeMi(bspline_inner_y_, g_inner_y_,  -lambda-1);
        
        // calculate inner basis full integral moments
        Mtr_L_inner_x_[lambda]   .populate([&](int i, int j) -> Complex { return computeM(bspline_inner_x_, g_inner_x_,  lambda,   i, j, true, true); });
        Mtr_L_inner_y_[lambda]   .populate([&](int i, int j) -> Complex { return computeM(bspline_inner_y_, g_inner_y_,  lambda,   i, j, true, true); });
        Mtr_mLm1_inner_x_[lambda].populate([&](int i, int j) -> Complex { return computeM(bspline_inner_x_, g_inner_x_, -lambda-1, i, j, true, true); });
        Mtr_mLm1_inner_y_[lambda].populate([&](int i, int j) -> Complex { return computeM(bspline_inner_y_, g_inner_y_, -lambda-1, i, j, true, true); });
        
        // calculate full basis partial integral moments
        cArrayView(Mitr_L_full_x_,    lambda * mi_size_full_x, mi_size_full_x) = computeMi(bspline_full_x_, g_full_x_,  lambda  );
        cArrayView(Mitr_L_full_y_,    lambda * mi_size_full_y, mi_size_full_y) = computeMi(bspline_full_y_, g_full_y_,  lambda  );
        cArrayView(Mitr_mLm1_full_x_, lambda * mi_size_full_x, mi_size_full_x) = computeMi(bspline_full_x_, g_full_x_, -lambda-1);
        cArrayView(Mitr_mLm1_full_y_, lambda * mi_size_full_y, mi_size_full_y) = computeMi(bspline_full_y_, g_full_y_, -lambda-1);
        
        // calculate full basis full integral moments
        Mtr_L_full_x_[lambda]   .populate([&](int i, int j) -> Complex { return computeM(bspline_full_x_, g_full_x_,  lambda,   i, j, true, true); });
        Mtr_L_full_y_[lambda]   .populate([&](int i, int j) -> Complex { return computeM(bspline_full_y_, g_full_y_,  lambda,   i, j, true, true); });
        Mtr_mLm1_full_x_[lambda].populate([&](int i, int j) -> Complex { return computeM(bspline_full_x_, g_full_x_, -lambda-1, i, j, true, true); });
        Mtr_mLm1_full_y_[lambda].populate([&](int i, int j) -> Complex { return computeM(bspline_full_y_, g_full_y_, -lambda-1, i, j, true, true); });
        
        // diagonal contributions to two-electron integrals
        std::string filename = format("rad-R_tr_dia_diag_%d-%.4lx-%.4lx.hdf", lambda, bspline_inner_x_.hash(), bspline_inner_y_.hash());
        if (not cmd.shared_scratch or par.isMyWork(lambda))
        {
            if (R_tr_dia_diag_[lambda].hdfload(filename))
            {
                if (verbose_)
                    std::cout << "\t- integrals for lambda = " << lambda << " loaded from \"" << filename << "\"" << std::endl;
            }
            else
            {
                Timer t;
                R_tr_dia_diag_[lambda] = diagonalR(lambda);
                
                if (verbose_)
                    std::cout << "\t- integrals for lambda = " << lambda << " computed after " << t.nice_time() << std::endl;
                
                R_tr_dia_diag_[lambda].hdfsave(filename);
            }
        }
    }
    
    if (cmd.shared_scratch)
    {
        // wait for all processes, so that the work sharing works
        par.wait();
        
        // this process will load all necessary data that were calculated by other processes
        for (int lambda = 0; lambda < Nlambdas_; lambda++) if (not par.isMyWork(lambda))
        {
            std::string filename = format("rad-R_tr_dia_diag_%d-%.4lx-%lx..hdf", lambda, bspline_inner_x_.hash(), bspline_inner_y_.hash());
            R_tr_dia_diag_[lambda].hdfload(filename);
            std::cout << "\t- integrals for lambda = " << lambda << " loaded from \"" << filename << "\"" << std::endl;
        }
    }
    
    if (verbose_)
        std::cout << std::endl;
    
    // abandon their computation, if not necessary
    if (cmd.lightweight_radial_cache)
        return;
    
    // allocate storage and associate names
    for (int lambda = 0; lambda < Nlambdas_; lambda++)
    {
        // keep data in memory?
        bool keep_in_memory = ((par.isMyWork(lambda) and cmd.cache_own_radint) or cmd.cache_all_radint);
        
        // compose the file name
        std::string filename = format("rad-R_tr_dia_%d-%.4lx-%.4lx.hdf", lambda, bspline_inner_x_.hash(), bspline_inner_y_.hash());
        
        // create the block matrix for radial integrals of this multipole
        R_tr_dia_[lambda] = BlockSymBandMatrix<Complex>
        (
            Nspline_full_x,     // block count
            order + 1,          // block structure half-bandwidth
            Nspline_full_y,     // block size
            order + 1,          // block half-bandwidth
            keep_in_memory,     // whether to keep in memory
            filename            // HDF scratch disk file name
        );
    }
    
    // print information
    if (verbose_)
        std::cout << "Precomputing multipole integrals (lambda = 0 .. " << Nlambdas_ - 1 << ")." << std::endl;
    
    // for all multipoles : compute / load
    for (int lambda = 0; lambda < Nlambdas_; lambda++)
    {
        // if the radial integrals are shared, this process will only compute the owned subset of radial integrals
        if (cmd.shared_scratch and not par.isMyWork(lambda))
            continue;
        
        // look for precomputed data on disk
        if (R_tr_dia_[lambda].hdfcheck())
        {
            if (not par.isMyWork(lambda) or not cmd.cache_own_radint)
            {
                if (verbose_)
                    std::cout << "\t- integrals for lambda = " << lambda << " present in \"" << R_tr_dia_[lambda].hdfname() << "\"" << std::endl;
                continue;
            }
            
            if (R_tr_dia_[lambda].hdfload())
            {
                if (verbose_)
                    std::cout << "\t- integrals for lambda = " << lambda << " loaded from \"" << R_tr_dia_[lambda].hdfname() << "\"" << std::endl;
                continue;
            }
        }
        else
        {
            R_tr_dia_[lambda].hdfinit();
        }
        
        Timer t;
        
        # pragma omp parallel firstprivate (lambda)
        {
            // for all blocks of the radial matrix
            # pragma omp for schedule (dynamic,1)
            for (int i = 0; i < Nspline_full_x; i++)
            for (int d = 0; d <= order; d++)
            if (i + d < Nspline_full_x)
            {
                // calculate the block
                SymBandMatrix<Complex> block = calc_R_tr_dia_block(lambda, i, i + d);
                
                // write the finished block to disk
                # pragma omp critical
                R_tr_dia_[lambda].setBlock(i * (order + 1) + d, block.data());
            }
        }
        
        if (verbose_)
            std::cout << "\t- integrals for lambda = " << lambda << " computed after " << t.nice_time() << std::endl;
        
        // save to disk even if the integrals are to be cached
        if (R_tr_dia_[lambda].inmemory())
            R_tr_dia_[lambda].hdfsave();
    }
    
    
    // wait for completition of all processes
    par.wait();
    
    // if this process skipped some lambda-s due to scratch sharing and still it needs them in memory, load them
    if (cmd.shared_scratch and cmd.cache_all_radint)
    {
        for (int lambda = 0; lambda < Nlambdas_; lambda++)
        {
            // skip own data (already loaded since calculation)
            if (par.isMyWork(lambda))
                continue;
            
            // load radial integrals
            if (R_tr_dia_[lambda].hdfload())
            {
                if (verbose_)
                    std::cout << "\t- integrals for lambda = " << lambda << " loaded from shared file \"" << R_tr_dia_[lambda].hdfname() << "\"" << std::endl;
            }
            else
                HexException("Can't read shared radial integral file \"%s\".", R_tr_dia_[lambda].hdfname().c_str());
        }
    }
    
    if (verbose_)
        std::cout << std::endl;
}

SymBandMatrix<Complex> RadialIntegrals::calc_R_tr_dia_block (unsigned int lambda, int i, int k, bool simple) const
{
    // shorthands
    int Nspline_full_y  = bspline_full_y_.Nspline();
    int Nspline_inner_y = bspline_inner_y_.Nspline(); // WARNING : Computing only inner-region integrals.
    int order = bspline_full_y_.order();
    
    // (i,k)-block data
    SymBandMatrix<Complex> block_ik (Nspline_full_y, order + 1);
    
    // for all elements in the symmetrical block : evaluate 2-D integral of Bi(1)Bj(2)V(1,2)Bk(1)Bl(2)
    for (int j = 0; j < Nspline_inner_y; j++)
    for (int l = j; l < Nspline_full_y and l <= j + order; l++)
        block_ik(j,l) = computeR(lambda, i, j, k, l, simple);
    
    return block_ik;
}

void RadialIntegrals::apply_R_matrix
(
    unsigned lambda,
    Complex a, const cArrayView src,
    Complex b,       cArrayView dst,
    bool simple
) const
{
    // shorthands
    std::size_t Nxspline = bspline_inner_x_.Nspline(); // WARNING : Computing only inner-region integrals.
    std::size_t Nyspline = bspline_inner_y_.Nspline(); // WARNING : Computing only inner-region integrals.
    std::size_t order = bspline_full_x_.order();
    
    // update destination vector
    # pragma omp simd
    for (std::size_t j = 0; j < dst.size(); j++)
        dst[j] *= b;
    
    // workspace
    cArray prod (Nyspline);
    
    // for all blocks of the radial matrix
    # pragma omp parallel for firstprivate (prod) schedule (dynamic, 1)
    for (std::size_t i = 0; i < Nxspline; i++)
    for (std::size_t k = i; k < Nxspline and k <= i + order; k++)
    {
        // (i,k)-block data (= concatenated non-zero upper diagonals)
        SymBandMatrix<Complex> block_ik = std::move ( calc_R_tr_dia_block(lambda, i, k, simple) );
        
        // multiply source vector by this block
        block_ik.dot(1., cArrayView(src, k * Nyspline, Nyspline), 0., prod);
        
        // update destination vector
        # pragma omp critical
        {
            # pragma omp simd
            for (std::size_t j  = 0; j < Nyspline; j++)
                dst[i * Nyspline + j] += a * prod[j];
        }
        
        // also handle symmetric case
        if (i != k)
        {
            // multiply source vector by this block
            block_ik.dot(1., cArrayView(src, i * Nyspline, Nyspline), 0., prod);
            
            // update destination vector
            # pragma omp critical
            {
                # pragma omp simd
                for (std::size_t j  = 0; j < Nyspline; j++)
                    dst[k * Nyspline + j] += a * prod[j];
            }
        }
    }
}

cArray RadialIntegrals::overlap (Bspline const & bspline, GaussLegendre const & g, std::function<Complex(Complex)> funct) const
{
    // result
    cArray res (bspline.Nspline());
    
    // per interval
    int points = EXPANSION_QUADRATURE_POINTS;
    
    // evaluated B-spline and function
    cArray evalB (points), evalF (points);
    
    // quadrature nodes and weights
    cArray xs (points), ws (points);
    
    // for all knots
    for (int iknot = 0; iknot < bspline.Nknot() - 1; iknot++)
    {
        // skip zero length intervals
        if (bspline.t(iknot) == bspline.t(iknot+1))
            continue;
        
        // which points are to be used here?
        g.scaled_nodes_and_weights(points, bspline.t(iknot), bspline.t(iknot+1), &xs[0], &ws[0]);
        
        // evaluate the function
        std::transform(xs.begin(), xs.end(), evalF.begin(), funct);
        
        // for all relevant B-splines
        for (int ispline = std::max(iknot - bspline.order(),0); ispline < bspline.Nspline() and ispline <= iknot; ispline++)
        {
            // evaluate the B-spline
            bspline.B(ispline, iknot, points, xs.data(), evalB.data());
            
            // sum with weights
            Complex sum = 0.;
            for (int ipoint = 0; ipoint < points; ipoint++)
                sum += ws[ipoint] * evalF[ipoint] * evalB[ipoint];
            
            // store the overlap
            res[ispline] += sum;
        }
    }
    
    return res;
}

cArray RadialIntegrals::overlapP (Bspline const & bspline, GaussLegendre const & g, int n, int l) const
{
    // result
    cArray res (bspline.Nspline());
    
    // per interval
    int points = EXPANSION_QUADRATURE_POINTS;
    
    // evaluated B-spline and hydrogenic functions
    cArray evalB (points), evalP (points);
    
    // quadrature nodes and weights
    cArray xs (points), ws (points);
    
    // for all knots
    # pragma omp parallel for firstprivate (points,xs,ws,evalB,evalP)
    for (int iknot = 0; iknot < bspline.Nknot() - 1; iknot++)
    {
        // skip zero length intervals
        if (bspline.t(iknot) == bspline.t(iknot+1))
            continue;
        
        // which points are to be used here?
        g.scaled_nodes_and_weights(points, bspline.t(iknot), bspline.t(iknot+1), &xs[0], &ws[0]);
        
        // evaluate the hydrogenic function
        std::transform
        (
            xs.begin(), xs.end(), evalP.begin(),
            [ = ](Complex x) -> Complex
            {
                gsl_sf_result R;
                if (gsl_sf_hydrogenicR_e(n, l, 1., x.real(), &R) == GSL_EUNDRFLW)
                    return 0.;
                else
                    return /* weightf(x) * */ x * Real(R.val);
            }
        );
        
        // for all relevant B-splines
        for (int ispline = std::max(iknot-bspline.order(),0); ispline < bspline.Nspline() and ispline <= iknot; ispline++)
        {
            // evaluate the B-spline
            bspline.B(ispline, iknot, points, xs.data(), evalB.data());
            
            // sum with weights
            Complex sum = 0.;
            for (int ipoint = 0; ipoint < points; ipoint++)
                sum += ws[ipoint] * evalP[ipoint] * evalB[ipoint];
            
            // store the overlap
            # pragma omp critical
            res[ispline] += sum;
        }
    }
    
    return res;
}

cArray RadialIntegrals::overlapj (Bspline const & bspline, GaussLegendre const & g, int maxell, const rArrayView vk, bool fast_bessel) const
{
    // shorthands
    int Nenergy = vk.size();
    int Nspline = bspline.Nspline();
    int Nknot = bspline.Nknot();
    int order = bspline.order();
    
    // reserve space for the output array
    cArray res (Nspline * Nenergy * (maxell + 1));
    
    // per interval
    int points = EXPANSION_QUADRATURE_POINTS;
    
    // quadrature weights and nodes
    cArray xs (points), ws (points);
    
    // other auxiliary variables
    cArray evalB ((order + 1) * points), evalj (points * (maxell + 1));
    
    // for all knots
    # pragma omp parallel for firstprivate (points,xs,ws,evalB,evalj)
    for (int iknot = 0; iknot < Nknot - 1; iknot++)
    {
        // skip zero length intervals
        if (bspline.t(iknot) == bspline.t(iknot+1))
            continue;
        
        // which points are to be used here?
        g.scaled_nodes_and_weights(points, bspline.t(iknot), bspline.t(iknot+1), &xs[0], &ws[0]);
        
        // evaluate relevant B-splines on this knot
        for (int ispline = std::max(iknot-order,0); ispline < Nspline and ispline <= iknot; ispline++)
            bspline.B(ispline, iknot, points, xs.data(), &evalB[(iknot - ispline) * points]);
        
        // for all linear momenta (= energies)
        for (int ie = 0; ie < Nenergy; ie++)
        {
            // evaluate the Riccati-Bessel function for this knot and energy and for all angular momenta
            for (int ipoint = 0; ipoint < points; ipoint++)
            {
                // compute the damping factor
                Real damp = 1;
                if (bspline.t(iknot).imag() != 0 or bspline.t(iknot + 1).imag() != 0)
                {
                    // total damping in complex region
                    damp = 0;
                }
                else
                {
                    // smooth damping in vicinity of the complex rotation point
                    damp = 0.5 * (1.0 + std::tanh(0.5 * (bspline.t(bspline.iR2() - order).real() - bspline.t(iknot + 1).real())));
                }
                
                // if the factor is numerical zero, do not evaluate the function at all, just fill zeros
                if (damp == 0.0_r)
                {
                    evalj.fill(0.0_z);
                    continue;
                }
                
                // which Bessel function evaluator to use?
                std::function<int(int,double,double*)> jv = (fast_bessel ? gsl_sf_bessel_jl_array : gsl_sf_bessel_jl_steed_array);
                
                // evaluate all Riccati-Bessel functions in point
                cArrayView(evalj, ipoint * (maxell + 1), maxell + 1) = damp * special::ric_jv(maxell, vk[ie] * xs[ipoint], jv);
                
                // clear all possible NaN entries (these may occur for far radii, where should be zero)
                for (int l = 0; l <= maxell; l++) if (not Complex_finite(evalj[ipoint * (maxell + 1) + l]))
                    evalj[ipoint * (maxell + 1) + l] = 0.0_z;
            }
            
            // for all angular momenta
            for (int l = 0; l <= maxell; l++)
            {
                // for all relevant B-splines
                for (int ispline = std::max(iknot-order,0); ispline < Nspline and ispline <= iknot; ispline++)
                {
                    // sum with weights
                    Complex sum = 0.;
                    for (int ipoint = 0; ipoint < points; ipoint++)
                        sum += ws[ipoint] * evalj[ipoint * (maxell + 1) + l] * evalB[(iknot - ispline) * points + ipoint];
                    
                    // store the overlap; keep the shape Nmomenta × Nspline × (maxl+1)
                    # pragma omp critical
                    res[(ie * (maxell + 1) + l) * Nspline + ispline] += sum;
                }
            }
        }
    }
    
    return res;
}
