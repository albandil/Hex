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
#include <complex>
#include <cmath>
#include <cstdio>
#include <cstring>

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"
#include "hex-csrmatrix.h"
#include "hex-hdffile.h"
#include "hex-special.h"
#include "hex-symbandmatrix.h"

// --------------------------------------------------------------------------------- //

#include "bspline.h"
#include "gauss.h"
#include "parallel.h"
#include "radial.h"

// --------------------------------------------------------------------------------- //

RadialIntegrals::RadialIntegrals
(
    Bspline const & bspline_x,
    Bspline const & bspline_y,
    int Nlambdas
)
  : bspline_x_(bspline_x),
    bspline_y_(bspline_y),
    rxmin_(bspline_x.Rmin() >= 0 ? bspline_x.Rmin() : 1.0_r),
    rymin_(bspline_y.Rmin() >= 0 ? bspline_y.Rmin() : 1.0_r),
    rxmax_(bspline_x.R2()),
    rymax_(bspline_y.R2()),
    D_x_  (bspline_x.Nspline(), bspline_x.order() + 1),
    S_x_  (bspline_x.Nspline(), bspline_x.order() + 1),
    Mm1_x_(bspline_x.Nspline(), bspline_x.order() + 1),
    Mm2_x_(bspline_x.Nspline(), bspline_x.order() + 1),
    D_y_  (bspline_y.Nspline(), bspline_y.order() + 1),
    S_y_  (bspline_y.Nspline(), bspline_y.order() + 1),
    Mm1_y_(bspline_y.Nspline(), bspline_y.order() + 1),
    Mm2_y_(bspline_y.Nspline(), bspline_y.order() + 1),
    verbose_(true),
    Nlambdas_(Nlambdas)
{
    // maximal number of evaluation points (quadrature rule)
    int npts = std::max(EXPANSION_QUADRATURE_POINTS, bspline_x_.order() + Nlambdas + 2);

    // precompute Gaussian weights
    g_x_.precompute_nodes_and_weights(npts);
    g_y_.precompute_nodes_and_weights(npts);
}

cArray RadialIntegrals::computeMi (Bspline const & bspline, GaussLegendre const & g, int a, Real rmin, Real rmax) const
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
    {
        // get restricted scale factor
        Real scale = special::clamp(bspline.unrotate(bspline.t(iknot + 1)), rmin, rmax);

        for (int ipoint = 0; ipoint < points; ipoint++)
        {
            // get restricted radius
            Complex r = bspline.clamp(xs[iknot * points + ipoint], rmin, rmax);

            // calculate effective power
            x_a[iknot * points + ipoint] = special::pow_int(r / scale, a);
        }
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

/*SymBandMatrix<Complex> RadialIntegrals::computeOverlapMatrix (Bspline const & bspline, CCFunction func) const
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
}*/

CsrMatrix<LU_int_t, Complex> RadialIntegrals::computeOverlapMatrix
(
    Bspline const & bspline1,
    Bspline const & bspline2,
    Real R_left,
    Real R_right
)
{
    CooMatrix<LU_int_t, Complex> coo (bspline1.Nspline(), bspline2.Nspline());

    int points = (bspline1.order() + bspline2.order()) / 2 + 1;

    GaussLegendre g;
    g.precompute_nodes_and_weights(points);

    cArray xs(points), ws(points), B1(points), B2(points);

    for (int i = 0; i < bspline1.Nspline(); i++)
    for (int j = 0; j < bspline2.Nspline(); j++)
    {
        // Check that the two B-splines overlap. This happens when
        // both B-splines start before the end of the other B-spline.
        if
        (
            bspline1.t(i).real() < bspline2.t(j + bspline2.order() + 1).real() and
            bspline2.t(j).real() < bspline1.t(i + bspline1.order() + 1).real()
        )
        {
            // overlap integral
            Complex v = 0;

            // starting knot of the other B-spline
            int iknot = std::max(i, bspline1.knot(bspline2.t(j)));
            int jknot = std::max(j, bspline2.knot(bspline1.t(i)));

            // integration bounds
            Complex left = (bspline1.t(iknot).real() > bspline2.t(jknot).real() ? bspline1.t(iknot) : bspline2.t(jknot));
            Complex right;

            // integrate on all knots
            while (iknot <= i + bspline1.order() and jknot <= j + bspline2.order())
            {
                bool ishift = false, jshift = false;

                // get nearest right knot
                if (bspline1.t(iknot + 1).real() == bspline2.t(jknot + 1).real())
                {
                    right = bspline1.t(iknot + 1);
                    ishift = true;
                    jshift = true;
                }
                else if (bspline1.t(iknot + 1).real() < bspline2.t(jknot + 1).real())
                {
                    right = bspline1.t(iknot + 1);
                    ishift = true;
                }
                else /* (bspline1.t(iknot + 1).real() > bspline2.t(jknot + 1).real()) */
                {
                    right = bspline2.t(jknot + 1);
                    jshift = true;
                }

                // integrate the product of the two B-splines on the interval (left, right)
                if (R_left <= left.real() and right.real() <= R_right)
                {
                    g.scaled_nodes_and_weights(points, left, right, &xs[0], &ws[0]);
                    bspline1.B(i, iknot, points, &xs[0], &B1[0]);
                    bspline2.B(j, jknot, points, &xs[0], &B2[0]);
                    for (int ipt = 0; ipt < points; ipt++)
                        v += B1[ipt] * B2[ipt] * ws[ipt];
                }

                // move on to the next interval
                if (ishift) iknot++;
                if (jshift) jknot++;
                left = right;
            }

            // store the element
            coo.add(i, j, v);
        }
    }

    return coo.tocsr();
}

Complex RadialIntegrals::computeM_iknot
(
    Bspline const & bspline, GaussLegendre const & g,
    int a, int i, int j,
    int iknot, Real rmin, Real rmax, Real scale
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

    // restrict radial dependence to (rxmin_,rxmax_) or (rymin_,rymax_)
    for (int k = 0; k < points; k++)
        xs[k] = bspline.clamp(xs[k], rmin, rmax);

    // accumulate the (weighted) result
    for (int k = 0; k < points; k++)
        res += values_i[k] * values_j[k] * special::pow_int(scale * xs[k], a) * ws[k];

    return res;
}

Complex RadialIntegrals::computeM
(
    Bspline const & bspline, GaussLegendre const & g,
    int a, int i, int j,
    Real rmin, Real rmax, bool scale
) const
{
    // get indices of knots that enclose the overlap of the two B-splines
    int left  = std::max(i, j);
    int right = std::min(i, j) + bspline.order() + 1;

    return computeM(bspline, g, a, i, j, left, right, rmin, rmax, scale);
}

Complex RadialIntegrals::computeM
(
    Bspline const & bspline, GaussLegendre const & g,
    int a, int i, int j,
    int left, int right,
    Real rmin, Real rmax,
    bool scale
) const
{
    // the result
    Complex res = 0;

    // calculate scale factor
    Real scalefactor = 1;
    if (scale)
    {
        scalefactor /= special::clamp
        (
            bspline.unrotate(bspline.t(std::min(i,j) + bspline.order() + 1)),
            rmin, rmax
        );
    }

    // integration on sub-intervals
    for (int iknot = left; iknot < right; iknot++)
        res += computeM_iknot(bspline, g, a, i, j, iknot, rmin, rmax, scalefactor);

    return res;
}

void RadialIntegrals::setupOneElectronIntegrals (Parallel const & par, CommandLine const & cmd)
{
    setupOneElectronIntegrals(cmd.shared_scratch, par.IamMaster());

    // Precompute also diagonal parts of the two-electron integrals. These are more computationally intensive
    // than the rest and the possibility to reuse them leads to a great speedup.

    R_coupled_.resize(Nlambdas_);

    if (verbose_) std::cout << "Precomputing coupled two-electron integrals ..." << std::endl;

    for (int lambda = 0; lambda < Nlambdas_; lambda++)
    {
        std::string filename = format("rad-R-diag-%d-%.4lx-%.4lx.hdf", lambda, bspline_x_.hash(), bspline_y_.hash());

        if (not cmd.shared_scratch or par.isMyWork(lambda))
        {
            // try to load the data from disk
            if (R_coupled_[lambda].hdfload(filename))
            {
                if (verbose_) std::cout << "\t- integrals for lambda = " << lambda << " loaded from \"" << filename << "\"" << std::endl;
            }

            // calculate new integrals
            else
            {
                Timer t;
                coupledR(lambda);
                if (verbose_) std::cout << "\t- integrals for lambda = " << lambda << " computed after " << t.nice_time() << std::endl;

                R_coupled_[lambda].hdfsave(filename);
            }
        }
    }

    par.wait();

    // Load the integrals from disk in case of distributed shared-scratch mode.

    for (int lambda = 0; lambda < Nlambdas_; lambda++)
    {
        std::string filename = format("rad-R-diag-%d-%.4lx-%.4lx.hdf", lambda, bspline_x_.hash(), bspline_y_.hash());

        if (cmd.shared_scratch and not par.isMyWork(lambda))
        {
            if (R_coupled_[lambda].hdfload(filename))
            {
                if (verbose_) std::cout << "\t- integrals for lambda = " << lambda << " loaded from \"" << filename << "\"" << std::endl;
            }
            else
            {
                HexException("Failed to load radial integrals calculated by other process.");
            }
        }
    }

    if (verbose_) std::cout << std::endl;
}

void RadialIntegrals::setupOneElectronIntegrals (bool shared_scratch, bool IamMaster)
{
    // inner basis one-electron integrals
    #define SetupOneElectronIntegrals(AXIS) \
        std::size_t hash_##AXIS = bspline_##AXIS##_.hash(); \
        \
        D_##AXIS##_     .hdflink(format("rad-D-%.4lx.hdf",      hash_##AXIS)); \
        S_##AXIS##_     .hdflink(format("rad-S-%.4lx.hdf",      hash_##AXIS)); \
        Mm1_##AXIS##_   .hdflink(format("rad-Mm1-%.4lx.hdf",    hash_##AXIS)); \
        Mm2_##AXIS##_   .hdflink(format("rad-Mm2-%.4lx.hdf",    hash_##AXIS)); \
        \
        D_##AXIS##_     .populate([=](int m, int n) -> Complex { return computeD(bspline_##AXIS##_, g_##AXIS##_,     m, n); }); \
        S_##AXIS##_     .populate([=](int m, int n) -> Complex { return computeM(bspline_##AXIS##_, g_##AXIS##_,  0, m, n, bspline_##AXIS##_.Rmin(), bspline_##AXIS##_.Rmax(), false); }); \
        Mm1_##AXIS##_   .populate([=](int m, int n) -> Complex { return computeM(bspline_##AXIS##_, g_##AXIS##_, -1, m, n, r##AXIS##min_,            r##AXIS##max_,            false); }); \
        Mm2_##AXIS##_   .populate([=](int m, int n) -> Complex { return computeM(bspline_##AXIS##_, g_##AXIS##_, -2, m, n, bspline_##AXIS##_.Rmin(), bspline_##AXIS##_.Rmax(), false); }); \
        \
        if (not shared_scratch or IamMaster) \
        { \
            D_##AXIS##_     .hdfsave(); \
            S_##AXIS##_     .hdfsave(); \
            Mm1_##AXIS##_   .hdfsave(); \
            Mm2_##AXIS##_   .hdfsave(); \
        }

    if (verbose_)
        std::cout << "Precomputing one-electron integrals ... " << std::flush;

    Timer t;

    SetupOneElectronIntegrals(x)
    SetupOneElectronIntegrals(y)

    if (verbose_)
        std::cout << "done in " << t.nice_time() << std::endl << std::endl;

    int order = bspline().order();
    t.reset();

    if (verbose_)
        std::cout << "Precomputing partial integral moments (lambda = 0 .. " << Nlambdas_ - 1 << ") ... " << std::flush;

    // partial moments sizes
    std::size_t mi_size_x = bspline_x_ .Nspline() * (2 * order + 1) * (order + 1);
    std::size_t mi_size_y = bspline_y_ .Nspline() * (2 * order + 1) * (order + 1);

    // resize partial moments arrays
    Mitr_L_x_   .resize(Nlambdas_ * mi_size_x);
    Mitr_L_y_   .resize(Nlambdas_ * mi_size_y);
    Mitr_mLm1_x_.resize(Nlambdas_ * mi_size_x);
    Mitr_mLm1_y_.resize(Nlambdas_ * mi_size_y);

    // resize full moments arrays, initialize with empty matrices
    Mtr_L_x_   .resize(Nlambdas_, SymBandMatrix<Complex>(bspline_x_ .Nspline(), order + 1));
    Mtr_L_y_   .resize(Nlambdas_, SymBandMatrix<Complex>(bspline_y_ .Nspline(), order + 1));
    Mtr_mLm1_x_.resize(Nlambdas_, SymBandMatrix<Complex>(bspline_x_ .Nspline(), order + 1));
    Mtr_mLm1_y_.resize(Nlambdas_, SymBandMatrix<Complex>(bspline_y_ .Nspline(), order + 1));

    // for all multipole moments
    for (int lambda = 0; lambda < Nlambdas_; lambda++)
    {
        // calculate partial integral moments
        cArrayView(Mitr_L_x_,    lambda * mi_size_x, mi_size_x) = computeMi(bspline_x_, g_x_,   lambda  , rxmin_, rxmax_);
        cArrayView(Mitr_L_y_,    lambda * mi_size_y, mi_size_y) = computeMi(bspline_y_, g_y_,   lambda  , rymin_, rymax_);
        cArrayView(Mitr_mLm1_x_, lambda * mi_size_x, mi_size_x) = computeMi(bspline_x_, g_x_,  -lambda-1, rxmin_, rxmax_);
        cArrayView(Mitr_mLm1_y_, lambda * mi_size_y, mi_size_y) = computeMi(bspline_y_, g_y_,  -lambda-1, rymin_, rymax_);

        // calculate full (scaled) integral moments
        Mtr_L_x_[lambda]   .populate([&](int i, int j) -> Complex { return computeM(bspline_x_, g_x_,  lambda,   i, j, rxmin_, rxmax_, true); });
        Mtr_L_y_[lambda]   .populate([&](int i, int j) -> Complex { return computeM(bspline_y_, g_y_,  lambda,   i, j, rymin_, rymax_, true); });
        Mtr_mLm1_x_[lambda].populate([&](int i, int j) -> Complex { return computeM(bspline_x_, g_x_, -lambda-1, i, j, rxmin_, rxmax_, true); });
        Mtr_mLm1_y_[lambda].populate([&](int i, int j) -> Complex { return computeM(bspline_y_, g_y_, -lambda-1, i, j, rymin_, rymax_, true); });
    }

    if (verbose_)
        std::cout << "done in " << t.nice_time() << std::endl << std::endl;
}

void RadialIntegrals::setupTwoElectronIntegrals (Parallel const & par, CommandLine const & cmd)
{
    // shorthands
    int order = bspline_x_.order();
    int Nspline_x  = bspline_x_ .Nspline();
    int Nspline_y  = bspline_y_ .Nspline();

    // set number of two-electron integrals
    R_tr_dia_.resize(Nlambdas_);

    // abandon their computation, if not necessary
    if (cmd.lightweight_radial_cache)
        return;

    // Precompute the full two-electron integrals. These matrices span a lot of memory and it may be
    // advantageous to use --lightweight-radial-cache to avoid their computation altogether. Or one
    // can use --no-radial-cache/--own-radial-cache to keep them on disk, though that may slow down
    // the calculation unnecessarily.

    // allocate storage and associate names
    for (int lambda = 0; lambda < Nlambdas_; lambda++)
    {
        // keep data in memory?
        bool keep_in_memory = ((par.isMyWork(lambda) and cmd.cache_own_radint) or cmd.cache_all_radint);

        // compose the file name
        std::string filename = format("rad-R-%d-%.4lx-%.4lx.hdf", lambda, bspline_x_.hash(), bspline_y_.hash());

        // create the block matrix for radial integrals of this multipole
        R_tr_dia_[lambda] = BlockSymBandMatrix<Complex>
        (
            Nspline_x,          // block count
            order + 1,          // block structure half-bandwidth
            Nspline_y,          // block size
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

        // for all blocks of the radial matrix
        # pragma omp parallel for
        for (int i = 0; i < Nspline_x; i++)
        for (int k = i; k < Nspline_x and k <= i + order; k++)
        {
            // calculate the block
            SymBandMatrix<Complex> block = calc_R_tr_dia_block(lambda, i, k);

            // update the finished block
            R_tr_dia_[lambda].setBlock(i * (order + 1) + k - i, block.data());
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

SymBandMatrix<Complex> RadialIntegrals::calc_R_tr_dia_block (unsigned int lambda, int i, int k) const
{
    // shorthands
    int Nspline_y = bspline_y_.Nspline();
    int order = bspline_y_.order();

    // (i,k)-block data
    SymBandMatrix<Complex> block_ik (Nspline_y, order + 1);

    // for all elements in the symmetrical block
    for (int j = 0; j < Nspline_y; j++)
    for (int l = j; l < Nspline_y and l <= j + order; l++)
        block_ik(j,l) = computeR(lambda, i, j, k, l);

    return block_ik;
}

void RadialIntegrals::apply_R_matrix
(
    unsigned lambda,
    Complex a, const cArrayView src,
    Complex b,       cArrayView dst,
    int x_row_limit,
    int y_row_limit
) const
{
    // shorthands
    std::size_t Nxspline = bspline_x_.Nspline();
    std::size_t Nyspline = bspline_y_.Nspline();
    std::size_t order = bspline_x_.order();

    // get expected number of blocks/rows
    std::size_t x_rows = x_row_limit >= 0 ? std::min<std::size_t>(x_row_limit, Nxspline) : Nxspline;
    std::size_t y_rows = y_row_limit >= 0 ? std::min<std::size_t>(y_row_limit, Nyspline) : Nyspline;

    // check compatibility
    if (x_rows * y_rows != src.size())
        HexException("Incompatible sizes %ldx%ld != %ld (src) passed to apply_R_matrix.", x_rows, y_rows, src.size());
    if (x_rows * y_rows != dst.size())
        HexException("Incompatible sizes %ldx%ld != %ld (dst) passed to apply_R_matrix.", x_rows, y_rows, dst.size());

    // update destination vector
    # pragma omp simd
    for (std::size_t j = 0; j < dst.size(); j++)
        dst[j] *= b;

    // workspace
    cArray prod (y_rows);

    // for all blocks of the radial matrix
    # pragma omp parallel for firstprivate (prod) schedule (dynamic, 1)
    for (std::size_t i = 0; i < x_rows; i++)
    for (std::size_t k = i; k < x_rows and k <= i + order; k++)
    {
        // (i,k)-block data (= concatenated non-zero upper diagonals)
        SymBandMatrix<Complex> block_ik
        (
            y_rows,
            order + 1,
            calc_R_tr_dia_block(lambda, i, k).data().slice(0, y_rows * (order + 1))
        );

        // multiply source vector by this block
        block_ik.dot(1., cArrayView(src, k * y_rows, y_rows), 0., prod);

        // update destination vector
        # pragma omp critical
        {
            # pragma omp simd
            for (std::size_t j  = 0; j < y_rows; j++)
                dst[i * y_rows + j] += a * prod[j];
        }

        // also handle symmetric case
        if (i != k)
        {
            // multiply source vector by this block
            block_ik.dot(1., cArrayView(src, i * y_rows, y_rows), 0., prod);

            // update destination vector
            # pragma omp critical
            {
                # pragma omp simd
                for (std::size_t j  = 0; j < y_rows; j++)
                    dst[k * y_rows + j] += a * prod[j];
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

cArray RadialIntegrals::overlap
(
    Bspline const & bspline1, GaussLegendre const & g1,
    Bspline const & bspline2, GaussLegendre const & g2,
    C2CFunction funct, int points
) const
{
    // resulting overlap
    cArray res (bspline1.Nspline() * bspline2.Nspline());

    // evaluate quadrature nodes and weights
    cArray xs ((bspline1.Nknot() - 1) * points), wxs ((bspline1.Nknot() - 1) * points);
    for (int ixknot = 0; ixknot < bspline1.Nknot() - 1; ixknot++)
        g1.scaled_nodes_and_weights(points, bspline1.t(ixknot), bspline1.t(ixknot+1), xs.data() + ixknot * points, wxs.data() + ixknot * points);
    cArray ys ((bspline2.Nknot() - 1) * points), wys ((bspline2.Nknot() - 1) * points);
    for (int iyknot = 0; iyknot < bspline2.Nknot() - 1; iyknot++)
        g2.scaled_nodes_and_weights(points, bspline2.t(iyknot), bspline2.t(iyknot+1), ys.data() + iyknot * points, wys.data() + iyknot * points);

    // evaluate B-splines
    cArray Bx (bspline1.Nspline() * (bspline1.order() + 1) * points);
    for (int ixspline = 0; ixspline < bspline1.Nspline(); ixspline++)
    for (int ixknot = ixspline; ixknot <= ixspline + bspline1.order() and ixknot < bspline1.Nknot() - 1; ixknot++)
        bspline1.B(ixspline, ixknot, points, xs.data() + ixknot * points, Bx.data() + (ixspline * (bspline1.order() + 1) + ixknot - ixspline) * points);
    cArray By (bspline2.Nspline() * (bspline2.order() + 1) * points);
    for (int iyspline = 0; iyspline < bspline2.Nspline(); iyspline++)
    for (int iyknot = iyspline; iyknot <= iyspline + bspline2.order() and iyknot < bspline2.Nknot() - 1; iyknot++)
        bspline2.B(iyspline, iyknot, points, ys.data() + iyknot * points, By.data() + (iyspline * (bspline2.order() + 1) + iyknot - iyspline) * points);

    // evaluated function
    cArray evalF (points * points);

    // for all knots
    # pragma omp parallel for firstprivate(evalF)
    for (int ixknot = 0; ixknot < bspline1.Nknot() - 1; ixknot++)
    for (int iyknot = 0; iyknot < bspline2.Nknot() - 1; iyknot++)
    {
        // skip zero length intervals
        if (bspline1.t(ixknot) == bspline1.t(ixknot+1) or bspline2.t(iyknot) == bspline2.t(iyknot+1))
            continue;

        // evaluate the function (with weights) in quadrature points
        for (int ixpoint = 0; ixpoint < points; ixpoint++)
        for (int iypoint = 0; iypoint < points; iypoint++)
        {
            evalF[ixpoint * points + iypoint] =
                funct(xs[ixknot * points + ixpoint], ys[iyknot * points + iypoint])
                * wxs[ixknot * points + ixpoint]
                * wys[ixknot * points + ixpoint];
        }

        // for all relevant B-splines
        for (int ixspline = std::max(ixknot - bspline1.order(), 0); ixspline < bspline1.Nspline() and ixspline <= ixknot; ixspline++)
        for (int iyspline = std::max(iyknot - bspline2.order(), 0); iyspline < bspline2.Nspline() and iyspline <= iyknot; iyspline++)
        {
            // get pointer to the B-spline evaluations
            Complex const * const Bxp = Bx.data() + (ixspline * (bspline1.order() + 1) + ixknot - ixspline) * points;
            Complex const * const Byp = By.data() + (iyspline * (bspline2.order() + 1) + iyknot - iyspline) * points;

            // sum with weights
            Complex sum = 0.;
            for (int ixpoint = 0; ixpoint < points; ixpoint++)
            for (int iypoint = 0; iypoint < points; iypoint++)
                sum += evalF[ixpoint * points + iypoint] * Bxp[ixpoint] * Byp[iypoint];

            // store the overlap
            # pragma omp critical
            res[ixspline * bspline2.Nspline() + iyspline] += sum;
        }
    }

    return res;
}

cArray RadialIntegrals::overlapP (Bspline const & bspline, GaussLegendre const & g, Real Z, int n, int l)
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
                if (gsl_sf_hydrogenicR_e(n, l, Z, x.real(), &R) == GSL_EUNDRFLW)
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

cArray RadialIntegrals::overlapj (Bspline const & bspline, GaussLegendre const & g, int maxell, const rArrayView vk, bool fast_bessel)
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
                // which Bessel function evaluator to use?
                std::function<int(int,double,double*)> jv = (fast_bessel ? gsl_sf_bessel_jl_array : gsl_sf_bessel_jl_steed_array);

                // evaluate all Riccati-Bessel functions in point
                cArrayView(evalj, ipoint * (maxell + 1), maxell + 1) = special::ric_jv(maxell, vk[ie] * xs[ipoint].real(), jv);

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
