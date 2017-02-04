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
#include "hex-special.h"

// --------------------------------------------------------------------------------- //

#include "bspline.h"
#include "gauss.h"
#include "parallel.h"
#include "radial.h"

// --------------------------------------------------------------------------------- //

RadialIntegrals::RadialIntegrals
(
    Bspline const & bspline_inner,
    Bspline const & bspline_full,
    int Nlambdas
)
  : bspline_inner_(bspline_inner),
    bspline_full_ (bspline_full),
    D_inner_     (bspline_inner.Nspline(), bspline_inner.order() + 1),
    S_inner_     (bspline_inner.Nspline(), bspline_inner.order() + 1),
    Mm1_inner_   (bspline_inner.Nspline(), bspline_inner.order() + 1),
    Mm1_tr_inner_(bspline_inner.Nspline(), bspline_inner.order() + 1),
    Mm2_inner_   (bspline_inner.Nspline(), bspline_inner.order() + 1),
    D_full_      (bspline_full .Nspline(), bspline_full .order() + 1),
    S_full_      (bspline_full .Nspline(), bspline_full .order() + 1),
    Mm1_full_    (bspline_full .Nspline(), bspline_full .order() + 1),
    Mm1_tr_full_ (bspline_full .Nspline(), bspline_full .order() + 1),
    Mm2_full_    (bspline_full .Nspline(), bspline_full .order() + 1),
    verbose_(true),
    Nlambdas_(Nlambdas)
{
    // determine the needed order of Gauss-Legendre quadrature
    unsigned pts = std::max(EXPANSION_QUADRATURE_POINTS, bspline_inner.order() + Nlambdas + 2);
    
    // precompute the nodes and weights
    g_inner_.precompute_nodes_and_weights(pts);
    g_full_ .precompute_nodes_and_weights(pts);
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
        ws[ix] *= damp(0, xs[ix], bspline.R0());
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
    Bspline const & bspline, GaussLegendre const & g,
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
    Bspline const & bspline, GaussLegendre const & g,
    int i, int j, int maxknot
) const
{
    // get boundary iknots
    int left = std::max(i, j);
    int right = std::min(i, j) + bspline.order();
    
    // cut at maxknot
    if (right > maxknot)
        right = maxknot;
    
    // the result
    Complex res = 0;
    
    // undergo integration on sub-intervals
    for (int iknot = left; iknot <= right; iknot++)
        res += computeD_iknot(bspline, g, i, j, iknot);
        
    return res;
}

Complex RadialIntegrals::computeM_iknot
(
    Bspline const & bspline, GaussLegendre const & g,
    int a, int i, int j,
    int iknot, Complex R, Real scale
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
    
    // accumulate the (damped) result
    if (R != 0.0_r)
    {
        for (int k = 0; k < points; k++)
            res += values_i[k] * values_j[k] * special::pow_int(scale * xs[k], a) * ws[k] * damp(xs[k], 0., R);
    }
    else
    {
        for (int k = 0; k < points; k++)
            res += values_i[k] * values_j[k] * special::pow_int(scale * xs[k], a) * ws[k];
    }
    
    return res;
}

Complex RadialIntegrals::computeM
(
    Bspline const & bspline, GaussLegendre const & g,
    int a, int i, int j,
    int maxknot, bool doscale
) const
{
    // get boundary iknots
    int left = std::max(i, j);
    int right = std::min(i, j) + bspline.order() + 1;
    
    // calculate scaling factor
    double scale = (doscale ? 1 / bspline.t(right).real() : 1);
    
    // cut at maxknot
    if (maxknot != -1 and right > maxknot)
        right = maxknot;
    
    // the result
    Complex res = 0;
    
    // undergo integration on sub-intervals
    for (int iknot = left; iknot < right; iknot++)
        res += computeM_iknot(bspline, g, a, i, j, iknot, maxknot >= 0 ? bspline.t(maxknot) : 0., scale);
    
    return res;
}

Complex RadialIntegrals::computeS12
(
    GaussLegendre const & g,
    Bspline const & bspline1,
    Bspline const & bspline2,
    int i, int j
) const
{
    // skip non-overlapping B-splines
    if (bspline1.t(i).real() > bspline2.t(j + bspline2.order() + 1).real() or
        bspline2.t(j).real() > bspline1.t(i + bspline1.order() + 1).real())
        return 0.;
    
    // get all knots that define the B-splines
    std::vector<Complex> knots;
    for (int k = 0; k <= bspline1.order(); k++) knots.push_back(bspline1.t(i + k));
    for (int k = 0; k <= bspline2.order(); k++) knots.push_back(bspline2.t(j + k));
    
    // sort the knots, remove duplicates
    std::sort(knots.begin(), knots.end(), Complex_realpart_less);
    knots.erase(std::unique(knots.begin(), knots.end()), knots.end());
    
    // quadrature degree
    int points = std::max(bspline1.order() + 1, bspline2.order() + 1);
    
    // auxiliary arrays
    cArray eval1 (points), eval2 (points), nodes (points), weights (points);
    
    // the result
    Complex res = 0;
    
    // undergo integration on sub-intervals
    for (unsigned n = 1; n < knots.size(); n++)
    {
        // integrate on knots[n-1] ... knots[n]; check that both B-splines are defined there
        if (bspline1.t(i).real() <= knots[n-1].real() and knots[n].real() <= bspline1.t(i + bspline1.order() + 1).real() and
            bspline2.t(j).real() <= knots[n-1].real() and knots[n].real() <= bspline2.t(j + bspline2.order() + 1).real())
        {
            // get middle of the integration interval
            Complex mid = 0.5_r * (knots[n-1] + knots[n]);
            
            // find to which interval this point belongs in both bases
            int k1 = bspline1.knot(mid);
            int k2 = bspline2.knot(mid);
            
            // get quadrature data
            g.scaled_nodes_and_weights(points, knots[n-1], knots[n], nodes.data(), weights.data());
            
            // evaluate both B-splines on this interval
            bspline1.B(i, k1, points, nodes.data(), eval1.data());
            bspline2.B(j, k2, points, nodes.data(), eval2.data());
            
            // sum the evaluations
            res += sum(eval1 * eval2 * weights);
        }
    }
    
    return res;
}

void RadialIntegrals::setupOneElectronIntegrals (Parallel const & par, CommandLine const & cmd)
{
    // inner basis one-electron integrals
    #define SetupOneElectronIntegrals(REG) \
        int Nknot_##REG   = bspline_##REG##_.Nknot(); \
        int Nreknot_##REG = bspline_##REG##_.Nreknot(); \
        \
        std::size_t hash_##REG = bspline_##REG##_.hash(); \
        \
        D_##REG##_     .hdflink(format("rad-D-%.4lx.hdf",      hash_##REG)); \
        S_##REG##_     .hdflink(format("rad-S-%.4lx.hdf",      hash_##REG)); \
        Mm1_##REG##_   .hdflink(format("rad-Mm1-%.4lx.hdf",    hash_##REG)); \
        Mm1_tr_##REG##_.hdflink(format("rad-Mm1_tr-%.4lx.hdf", hash_##REG)); \
        Mm2_##REG##_   .hdflink(format("rad-Mm2-%.4lx.hdf",    hash_##REG)); \
        \
        D_##REG##_     .populate([=](int m, int n) -> Complex { return computeD(bspline_##REG##_, g_##REG##_,     m, n, Nknot_##REG - 1  ); }); \
        S_##REG##_     .populate([=](int m, int n) -> Complex { return computeM(bspline_##REG##_, g_##REG##_,  0, m, n                   ); }); \
        Mm1_##REG##_   .populate([=](int m, int n) -> Complex { return computeM(bspline_##REG##_, g_##REG##_, -1, m, n                   ); }); \
        Mm1_tr_##REG##_.populate([=](int m, int n) -> Complex { return computeM(bspline_##REG##_, g_##REG##_, -1, m, n, Nreknot_##REG - 1); }); \
        Mm2_##REG##_   .populate([=](int m, int n) -> Complex { return computeM(bspline_##REG##_, g_##REG##_, -2, m, n                   ); }); \
        \
        if (not cmd.shared_scratch or par.IamMaster()) \
        { \
            D_##REG##_     .hdfsave(); \
            S_##REG##_     .hdfsave(); \
            Mm1_##REG##_   .hdfsave(); \
            Mm1_tr_##REG##_.hdfsave(); \
            Mm2_##REG##_   .hdfsave(); \
        }
    
    Timer t;
    if (verbose_)
        std::cout << "Precomputing one-electron integrals ... " << std::flush;
    
    SetupOneElectronIntegrals(inner)
    SetupOneElectronIntegrals(full)
    
    if (verbose_)
        std::cout << "done in " << t.nice_time() << std::endl << std::endl;
    
//     // compute inter-basis overlaps
//     if (not cmd.map_solution.empty())
//     {
//         S12_ = RowMatrix<Complex>(Nspline_atom, Nspline_proj);
//         S21_ = RowMatrix<Complex>(Nspline_proj, Nspline_atom);
//         for (int i = 0; i < Nspline_atom; i++)
//         for (int j = 0; j < Nspline_proj; j++)
//             S12_(i,j) = S21_(j,i) = computeS12(g_atom_, bspline_atom_, bspline_proj_, i, j);
//     }
}

void RadialIntegrals::setupTwoElectronIntegrals (Parallel const & par, CommandLine const & cmd)
{
    // shorthands
    int order = bspline_inner_.order();
    int Nspline_inner = bspline_inner_.Nspline(), Nreknot_inner = bspline_inner_.Nreknot();
    int Nspline_full  = bspline_full_ .Nspline(), Nreknot_full  = bspline_full_ .Nreknot();
    
    // set number of two-electron integrals
    R_tr_dia_.resize(Nlambdas_);
    
    if (verbose_)
        std::cout << "Precomputing partial integral moments (lambda = 0 .. " << Nlambdas_ - 1 << ")" << std::endl;
    
    // partial moments
    std::size_t mi_size_inner = Nspline_inner * (2 * order + 1) * (order + 1);
    std::size_t mi_size_full  = Nspline_full  * (2 * order + 1) * (order + 1);
    Mitr_L_inner_   .resize(Nlambdas_ * mi_size_inner);  Mitr_L_full_   .resize(Nlambdas_ * mi_size_full); 
    Mitr_mLm1_inner_.resize(Nlambdas_ * mi_size_inner);  Mitr_mLm1_full_.resize(Nlambdas_ * mi_size_full);
    Mtr_L_inner_    .resize(Nlambdas_);                  Mtr_L_full_    .resize(Nlambdas_);
    Mtr_mLm1_inner_ .resize(Nlambdas_);                  Mtr_mLm1_full_ .resize(Nlambdas_);
    
    // resize vector of two-electron integrals
    R_tr_dia_diag_ .resize(Nlambdas_);
    
    // for all multipole moments
    for (int lambda = 0; lambda < Nlambdas_; lambda++)
    {
        // inner basis
        cArrayView(Mitr_L_inner_,    lambda * mi_size_inner, mi_size_inner) = computeMi(bspline_inner_, g_inner_,   lambda  );
        cArrayView(Mitr_mLm1_inner_, lambda * mi_size_inner, mi_size_inner) = computeMi(bspline_inner_, g_inner_,  -lambda-1);
        Mtr_L_inner_[lambda]    = SymBandMatrix<Complex>(Nspline_inner, order + 1).populate([&](int i, int j) -> Complex { return computeM(bspline_inner_, g_inner_,  lambda,   i, j, Nreknot_inner - 1, true); });
        Mtr_mLm1_inner_[lambda] = SymBandMatrix<Complex>(Nspline_inner, order + 1).populate([&](int i, int j) -> Complex { return computeM(bspline_inner_, g_inner_, -lambda-1, i, j, Nreknot_inner - 1, true); });
        
        // full basis
        cArrayView(Mitr_L_full_,    lambda * mi_size_full, mi_size_full) = computeMi(bspline_full_, g_full_,  lambda  );
        cArrayView(Mitr_mLm1_full_, lambda * mi_size_full, mi_size_full) = computeMi(bspline_full_, g_full_, -lambda-1);
        Mtr_L_full_[lambda]    = SymBandMatrix<Complex>(Nspline_full, order + 1).populate([&](int i, int j) -> Complex { return computeM(bspline_full_, g_full_,  lambda,   i, j, Nreknot_full - 1, true); });
        Mtr_mLm1_full_[lambda] = SymBandMatrix<Complex>(Nspline_full, order + 1).populate([&](int i, int j) -> Complex { return computeM(bspline_full_, g_full_, -lambda-1, i, j, Nreknot_full - 1, true); });
        
        // no need to do anything else for non-identical B-spline bases
//         if (bspline_atom_.hash() != bspline_proj_.hash())
//             continue;
        
        // diagonal contributions to two-electron integrals
        std::string filename = format("rad-R_tr_dia_diag_%d-%.4lx.hdf", lambda, bspline_inner_.hash());
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
                if (verbose_) std::cout << "\t- integrals for lambda = " << lambda << " computed after " << t.nice_time() << std::endl;
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
            std::string filename = format("rad-R_tr_dia_diag_%d-%.4lx.hdf", lambda, bspline_inner_.hash());
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
        std::string filename = format("rad-R_tr_dia_%d-%.4lx.hdf", lambda, bspline_inner_.hash());
        
        // create the block matrix for radial integrals of this multipole
        R_tr_dia_[lambda] = BlockSymBandMatrix<Complex>
        (
            Nspline_full,      // block count
            order + 1,          // block structure half-bandwidth
            Nspline_full,      // block size
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
            for (int i = 0; i < Nspline_full; i++)
            for (int d = 0; d <= order; d++)
            if (i + d < Nspline_full)
            {
                // calculate the block
                SymBandMatrix<Complex> block = calc_R_tr_dia_block(lambda, i, i + d, false);
                
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

SymBandMatrix<Complex> RadialIntegrals::calc_R_tr_dia_block (unsigned int lambda, int i, int k, bool inner_only, bool simple) const
{
    // shorthands
    int Nspline = inner_only ? bspline_inner_.Nspline() : bspline_full_.Nspline();
    int order = bspline_full_.order();
    
    // (i,k)-block data
    SymBandMatrix<Complex> block_ik (Nspline, order + 1);
    
    // for all elements in the symmetrical block : evaluate 2-D integral of Bi(1)Bj(2)V(1,2)Bk(1)Bl(2)
    for (int j = 0; j < Nspline; j++)
    for (int l = j; l < Nspline and l <= j + order; l++)
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
    std::size_t Nspline = bspline_inner_.Nspline(); // WARNING : Computing only inner-region integrals.
    std::size_t order = bspline_full_.order();
    
    // update destination vector
    # pragma omp simd
    for (std::size_t j = 0; j < dst.size(); j++)
        dst[j] *= b;
    
    // workspace
    cArray prod (Nspline);
    
    // for all blocks of the radial matrix
    # pragma omp parallel for firstprivate (prod) schedule (dynamic, 1)
    for (unsigned i = 0; i < Nspline; i++)
    for (std::size_t k = i; k < Nspline and k <= i + order; k++)
    {
        // (i,k)-block data (= concatenated non-zero upper diagonals)
        SymBandMatrix<Complex> block_ik = std::move ( calc_R_tr_dia_block(lambda, i, k, true, simple) );
        
        // multiply source vector by this block
        block_ik.dot(1., cArrayView(src, k * Nspline, Nspline), 0., prod);
        
        // update destination vector
        # pragma omp critical
        {
            # pragma omp simd
            for (std::size_t j  = 0; j < Nspline; j++)
                dst[i * Nspline + j] += a * prod[j];
        }
        
        // also handle symmetric case
        if (i != k)
        {
            // multiply source vector by this block
            block_ik.dot(1., cArrayView(src, i * Nspline, Nspline), 0., prod);
            
            // update destination vector
            # pragma omp critical
            {
                # pragma omp simd
                for (std::size_t j  = 0; j < Nspline; j++)
                    dst[k * Nspline + j] += a * prod[j];
            }
        }
    }
}

cArray RadialIntegrals::overlap (Bspline const & bspline, GaussLegendre const & g, std::function<Complex(Complex)> funct, std::function<Real(Complex)> weightf) const
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

cArray RadialIntegrals::overlapP (Bspline const & bspline, GaussLegendre const & g, int n, int l, std::function<Real(Complex)> weightf) const
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
                    return weightf(x) * x * Real(R.val);
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

cArray RadialIntegrals::overlapj (Bspline const & bspline, GaussLegendre const & g, int maxell, const rArrayView vk, std::function<Real(Complex)> weightf, bool fast_bessel) const
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
                Real damp = weightf(xs[ipoint]);
                
                // if the factor is numerical zero, do not evaluate the function, just fill zeros
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
