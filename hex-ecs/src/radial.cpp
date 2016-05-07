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

#include "bspline.h"
#include "gauss.h"
#include "parallel.h"
#include "radial.h"

RadialIntegrals::RadialIntegrals
(
    Bspline const & bspline_atom,
    Bspline const & bspline_proj,
    Bspline const & bspline_proj_full,
    int Nlambdas
)
  : bspline_atom_(bspline_atom), bspline_proj_(bspline_proj), bspline_proj_full_(bspline_proj_full),
    g_atom_(bspline_atom), g_proj_(bspline_proj),
    D_atom_(bspline_atom.Nspline(),bspline_atom.order()+1),
    S_atom_(bspline_atom.Nspline(),bspline_atom.order()+1),
    Mm1_atom_(bspline_atom.Nspline(),bspline_atom.order()+1),
    Mm1_tr_atom_(bspline_atom.Nspline(),bspline_atom.order()+1),
    Mm2_atom_(bspline_atom.Nspline(),bspline_atom.order()+1),
    D_proj_(bspline_proj.Nspline(),bspline_proj.order()+1),
    S_proj_(bspline_proj.Nspline(),bspline_proj.order()+1),
    Mm1_proj_(bspline_proj.Nspline(),bspline_proj.order()+1),
    Mm1_tr_proj_(bspline_proj.Nspline(),bspline_proj.order()+1),
    Mm2_proj_(bspline_proj.Nspline(),bspline_proj.order()+1),
    verbose_(true), Nlambdas_(Nlambdas)
{
    // maximal number of evaluation points (quadrature rule)
    int npts = std::max(EXPANSION_QUADRATURE_POINTS, bspline_atom_.order() + Nlambdas + 1);
    
    // get projectile basis shift
    proj_basis_shift_ = std::find(bspline_atom_.t().begin(), bspline_atom_.t().end(), bspline_proj_.t().front()) - bspline_atom_.t().begin();
    
    // precompute Gaussian weights
    g_atom_.precompute_nodes_and_weights(npts);
    g_proj_.precompute_nodes_and_weights(npts);
}

void RadialIntegrals::Mi_integrand
(
    int n,
    Complex * const restrict in,
    Complex * const restrict out,
    Bspline const & bspline_ij, int i, int j, int a,
    int iknot, int iknotmax
) const
{
    // extract data
    Complex R = bspline_ij.t(iknotmax);
    
    // evaluate B-splines
    cArray values_i(n), values_j(n);
    bspline_ij.B(i, iknot, n, in, values_i.data());
    bspline_ij.B(j, iknot, n, in, values_j.data());
    
    // get upper bound
    Real t = bspline_ij.t(iknot + 1).real();
    
    // scale factor for the multipole
    Real scalef = 1 / t;
    
    // fill output array
    if (R != 0.0_r)
    {
        // use damping
        for (int k = 0; k < n; k++)
            out[k] = values_i[k] * values_j[k] * special::pow_int<Complex>(scalef*in[k],a) * damp(in[k],0.,R);
    }
    else
    {
        // do not use damping
        for (int k = 0; k < n; k++)
            out[k] = values_i[k] * values_j[k] * special::pow_int<Complex>(scalef*in[k],a);
    }
}

cArray RadialIntegrals::computeMi (Bspline const & bspline, GaussLegendre const & g, int a, int iknotmax) const
{
    int Nspline = bspline.Nspline();
    int order = bspline.order();
    
    // partial integral moments
    cArray m (Nspline * (2 * order + 1) * (order + 1));
    
    // for all B-splines
    # pragma omp parallel
    for (int i = 0; i < Nspline; i++)
    {
        // for all B-splines (use symmetry)
        for (int j = i; j <= i + order and j < Nspline; j++)
        {
            // determine relevant knots
            int ileft = j;
            int iright = i + order + 1;
            
            // "right" has to be smaller than "Rmax"
            if (iright > iknotmax)
                iright = iknotmax;
            
            // for all relevant knots
            for (int iknot = ileft; iknot < iright; iknot++)
            {
                // get integration boundaries
                Complex xa = bspline.t(iknot);
                Complex xb = bspline.t(iknot+1);
                
                // throw away zero length intervals
                if (xa == xb)
                    continue;
                
                // results of the quadrature
                Complex integral;
                
                // use at least 2nd order
                int points = std::max(2, order + std::abs(a) + 1);
                
                // integrate
                integral = g.quadMFP
                (
                    this, &RadialIntegrals::Mi_integrand,      // integrand pointer
                    points, iknot, xa, xb,                     // integration parameters
                    bspline, i, j, a, iknot, iknotmax          // data to pass to the integrator
                );
                
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
        }
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
    // shorthands
    int Nspline_atom = bspline_atom_.Nspline(); int Nspline_proj = bspline_proj_.Nspline();
    int Nknot_atom   = bspline_atom_.Nknot();   int Nknot_proj   = bspline_proj_.Nknot();
    int Nreknot_atom = bspline_atom_.Nreknot(); int Nreknot_proj = bspline_proj_.Nreknot();
    
    // create file names for these radial integrals
    D_atom_     .hdflink(format("rad-D-%.4lx.hdf",      bspline_atom_.hash()));  D_proj_     .hdflink(format("rad-D-%.4lx.hdf",      bspline_proj_.hash()));
    S_atom_     .hdflink(format("rad-S-%.4lx.hdf",      bspline_atom_.hash()));  S_proj_     .hdflink(format("rad-S-%.4lx.hdf",      bspline_proj_.hash()));
    Mm1_atom_   .hdflink(format("rad-Mm1-%.4lx.hdf",    bspline_atom_.hash()));  Mm1_proj_   .hdflink(format("rad-Mm1-%.4lx.hdf",    bspline_proj_.hash()));
    Mm1_tr_atom_.hdflink(format("rad-Mm1_tr-%.4lx.hdf", bspline_atom_.hash()));  Mm1_tr_proj_.hdflink(format("rad-Mm1_tr-%.4lx.hdf", bspline_proj_.hash()));
    Mm2_atom_   .hdflink(format("rad-Mm2-%.4lx.hdf",    bspline_atom_.hash()));  Mm2_proj_   .hdflink(format("rad-Mm2-%.4lx.hdf",    bspline_proj_.hash()));
    
    if (verbose_)
        std::cout << "Precomputing one-electron integrals ... " << std::flush;
    
    // compute one-electron matrices (atom basis)
    D_atom_     .populate([=](int m, int n) -> Complex { return computeD(bspline_atom_, g_atom_,     m, n, Nknot_atom - 1  ); });
    S_atom_     .populate([=](int m, int n) -> Complex { return computeM(bspline_atom_, g_atom_,  0, m, n                  ); });
    Mm1_atom_   .populate([=](int m, int n) -> Complex { return computeM(bspline_atom_, g_atom_, -1, m, n                  ); });
    Mm1_tr_atom_.populate([=](int m, int n) -> Complex { return computeM(bspline_atom_, g_atom_, -1, m, n, Nreknot_atom - 1); });
    Mm2_atom_   .populate([=](int m, int n) -> Complex { return computeM(bspline_atom_, g_atom_, -2, m, n                  ); });
    
    // compute one-electron matrices (projectile basis)
    D_proj_     .populate([=](int m, int n) -> Complex { return computeD(bspline_proj_, g_proj_,     m, n, Nknot_proj - 1  ); });
    S_proj_     .populate([=](int m, int n) -> Complex { return computeM(bspline_proj_, g_proj_,  0, m, n                  ); });
    Mm1_proj_   .populate([=](int m, int n) -> Complex { return computeM(bspline_proj_, g_proj_, -1, m, n                  ); });
    Mm1_tr_proj_.populate([=](int m, int n) -> Complex { return computeM(bspline_proj_, g_proj_, -1, m, n, Nreknot_proj - 1); });
    Mm2_proj_   .populate([=](int m, int n) -> Complex { return computeM(bspline_proj_, g_proj_, -2, m, n                  ); });
    
    // compute inter-basis overlaps
    if (not cmd.map_solution.empty())
    {
        S12_ = RowMatrix<Complex>(Nspline_atom, Nspline_proj);
        S21_ = RowMatrix<Complex>(Nspline_proj, Nspline_atom);
        for (int i = 0; i < Nspline_atom; i++)
        for (int j = 0; j < Nspline_proj; j++)
            S12_(i,j) = S21_(j,i) = computeS12(g_atom_, bspline_atom_, bspline_proj_, i, j);
    }
    
    // save the matrices to disk
    if (not cmd.shared_scratch or par.IamMaster())
    {
        D_atom_     .hdfsave();  D_proj_     .hdfsave();
        S_atom_     .hdfsave();  S_proj_     .hdfsave();
        Mm1_atom_   .hdfsave();  Mm1_proj_   .hdfsave();
        Mm1_tr_atom_.hdfsave();  Mm1_tr_proj_.hdfsave();
        Mm2_atom_   .hdfsave();  Mm2_proj_   .hdfsave();
    }
    
    if (verbose_)
        std::cout << "ok" << std::endl << std::endl;
}

void RadialIntegrals::setupTwoElectronIntegrals (Parallel const & par, CommandLine const & cmd)
{
    // shorthands
    int order = bspline_atom_.order();
    int Nspline_atom = bspline_atom_.Nspline();
    int Nreknot_atom = bspline_atom_.Nreknot();
    int Nspline_proj = bspline_proj_.Nspline();
    int Nreknot_proj = bspline_proj_.Nreknot();
    
    // get knot that terminates (x^lambda)-scaled region
    for (int i = 0; i < bspline_atom_.Nknot(); i++)
    {
        if (bspline_atom_.t(i).real() < 1)
            lastscaled_ = i;
        else
            break;
    }
    
    // set number of two-electron integrals
    R_tr_dia_.resize(Nlambdas_);
    
    if (verbose_)
        std::cout << "Precomputing partial integral moments (lambda = 0 .. " << Nlambdas_ - 1 << ")" << std::endl;
    
    // compute partial moments
    std::size_t mi_size_atom = Nspline_atom * (2 * order + 1) * (order + 1);
    std::size_t mi_size_proj = Nspline_proj * (2 * order + 1) * (order + 1);
    Mitr_L_atom_   .resize(Nlambdas_ * mi_size_atom);    Mitr_L_proj_   .resize(Nlambdas_ * mi_size_proj);
    Mitr_mLm1_atom_.resize(Nlambdas_ * mi_size_atom);    Mitr_mLm1_proj_.resize(Nlambdas_ * mi_size_proj);
    Mtr_L_atom_    .resize(Nlambdas_);                   Mtr_L_proj_    .resize(Nlambdas_);
    Mtr_mLm1_atom_ .resize(Nlambdas_);                   Mtr_mLm1_proj_ .resize(Nlambdas_);
    R_tr_dia_diag_ .resize(Nlambdas_);
    for (int lambda = 0; lambda < (int)Nlambdas_; lambda++)
    {
        // atomic basis
        cArrayView(Mitr_L_atom_,    lambda * mi_size_atom, mi_size_atom) = computeMi(bspline_atom_, g_atom_,   lambda,   Nreknot_atom - 1);
        cArrayView(Mitr_mLm1_atom_, lambda * mi_size_atom, mi_size_atom) = computeMi(bspline_atom_, g_atom_,  -lambda-1, Nreknot_atom - 1);
        Mtr_L_atom_[lambda]    = SymBandMatrix<Complex>(Nspline_atom, order + 1).populate([&](int i, int j) -> Complex { return computeM(bspline_atom_, g_atom_,  lambda,   i, j, Nreknot_atom - 1, true); });
        Mtr_mLm1_atom_[lambda] = SymBandMatrix<Complex>(Nspline_atom, order + 1).populate([&](int i, int j) -> Complex { return computeM(bspline_atom_, g_atom_, -lambda-1, i, j, Nreknot_atom - 1, true); });
        
        // projectile basis
        cArrayView(Mitr_L_proj_,    lambda * mi_size_proj, mi_size_proj) = computeMi(bspline_proj_, g_proj_,  lambda,   Nreknot_proj - 1);
        cArrayView(Mitr_mLm1_proj_, lambda * mi_size_proj, mi_size_proj) = computeMi(bspline_proj_, g_proj_, -lambda-1, Nreknot_proj - 1);
        Mtr_L_proj_[lambda]    = SymBandMatrix<Complex>(Nspline_proj, order + 1).populate([&](int i, int j) -> Complex { return computeM(bspline_proj_, g_proj_,  lambda,   i, j, Nreknot_proj - 1, true); });
        Mtr_mLm1_proj_[lambda] = SymBandMatrix<Complex>(Nspline_proj, order + 1).populate([&](int i, int j) -> Complex { return computeM(bspline_proj_, g_proj_, -lambda-1, i, j, Nreknot_proj - 1, true); });
        
        // no need to do anything else for non-identical B-spline bases
//         if (bspline_atom_.hash() != bspline_proj_.hash())
//             continue;
        
        // diagonal contributions to two-electron integrals
        std::string filename = format("rad-R_tr_dia_diag_%d-%.4lx.hdf", lambda, bspline_atom_.hash());
        if (not cmd.shared_scratch or par.isMyWork(lambda))
        {
            if (R_tr_dia_diag_[lambda].hdfload(filename))
            {
                if (verbose_)
                    std::cout << "\t- integrals for lambda = " << lambda << " loaded from \"" << filename << "\"" << std::endl;
            }
            else
            {
                R_tr_dia_diag_[lambda] = diagonalR(lambda);
                R_tr_dia_diag_[lambda].hdfsave(filename);
                if (verbose_)
                    std::cout << "\t- integrals for lambda = " << lambda << " computed" << std::endl;
            }
        }
    }
    
    if (cmd.shared_scratch)
    {
        // wait for all processes, so that the work sharing works
        par.wait();
        
        // this process will load all necessary data that were calculated by other processes
        for (int lambda = 0; lambda < (int)Nlambdas_; lambda++) if (not par.isMyWork(lambda))
        {
            std::string filename = format("rad-R_tr_dia_diag_%d-%.4lx.hdf", lambda, bspline_atom_.hash());
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
        std::string filename = format("rad-R_tr_dia_%d-%.4lx-%.4lx.hdf", lambda, bspline_atom_.hash(), bspline_proj_.hash());
        
        // create the block matrix for radial integrals of this multipole
        R_tr_dia_[lambda] = BlockSymBandMatrix<Complex>
        (
            Nspline_atom,       // block count
            order + 1,          // block structure half-bandwidth
            Nspline_proj,       // block size
            order + 1,          // block half-bandwidth
            keep_in_memory,     // whether to keep in memory
            filename            // HDF scratch disk file name
        );
    }
    
    // print information
    if (verbose_)
        std::cout << "Precomputing multipole integrals (lambda = 0 .. " << Nlambdas_ - 1 << ")." << std::endl;
    
    // for all multipoles : compute / load
    for (int lambda = 0; lambda < (int)Nlambdas_; lambda++)
    {
        // if the radial integrals are shared, this process will only compute the owned subset of radial integrals
        if (cmd.shared_scratch and not par.isMyWork(lambda))
            continue;
        
        // look for precomputed data on disk
        if (R_tr_dia_[lambda].hdfcheck())
        {
            if (not par.isMyWork(lambda) and not cmd.cache_own_radint)
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
        
        # pragma omp parallel firstprivate (lambda)
        {
            // for all blocks of the radial matrix
            # pragma omp for schedule (dynamic,1)
            for (int i = 0; i < Nspline_atom; i++)
            for (int d = 0; d <= order; d++)
            if (i + d < Nspline_atom)
            {
                // calculate the block
                SymBandMatrix<Complex> block = calc_R_tr_dia_block(lambda, i, i + d);
                
                // write the finished block to disk
                # pragma omp critical
                R_tr_dia_[lambda].setBlock(i * (order + 1) + d, block.data());
            }
        }
        
        if (verbose_)
            std::cout << "\t- integrals for lambda = " << lambda << " computed" << std::endl;
        
        // save to disk even if the integrals are to be cached
        if (R_tr_dia_[lambda].inmemory())
            R_tr_dia_[lambda].hdfsave();
    }
    
    
    // wait for completition of all processes
    par.wait();
    
    // if this process skipped some lambda-s due to scratch sharing and still it needs them in memory, load them
    if (cmd.shared_scratch and cmd.cache_all_radint)
    {
        for (int lambda = 0; lambda < (int)Nlambdas_; lambda++)
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
    int Nspline_proj = bspline_proj_.Nspline();
    int order = bspline_proj_.order();
    
    // (i,k)-block data
    SymBandMatrix<Complex> block_ik (Nspline_proj, order + 1);
    
    // for all elements in the symmetrical block : evaluate 2-D integral of Bi(1)Bj(2)V(1,2)Bk(1)Bl(2)
    for (int j = 0; j < Nspline_proj; j++)
    for (int l = j; l < Nspline_proj and l <= j + order; l++)
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
    std::size_t Nspline_atom = bspline_atom_.Nspline();
    std::size_t Nspline_proj = bspline_proj_.Nspline();
    std::size_t order = bspline_atom_.order();
    
    // update destination vector
    # pragma omp simd
    for (std::size_t j = 0; j < dst.size(); j++)
        dst[j] *= b;
    
    // workspace
    cArray prod (Nspline_proj);
    
    // for all blocks of the radial matrix
    # pragma omp parallel for firstprivate (prod) schedule (dynamic, 1)
    for (unsigned i = 0; i < Nspline_atom; i++)
    for (std::size_t k = i; k < Nspline_atom and k <= i + order; k++)
    {
        // (i,k)-block data (= concatenated non-zero upper diagonals)
        SymBandMatrix<Complex> block_ik = std::move ( calc_R_tr_dia_block(lambda, i, k, simple) );
        
        // multiply source vector by this block
        block_ik.dot(1., cArrayView(src, k * Nspline_proj, Nspline_proj), 0., prod);
        
        // update destination vector
        # pragma omp critical
        {
            # pragma omp simd
            for (std::size_t j  = 0; j < Nspline_proj; j++)
                dst[i * Nspline_proj + j] += a * prod[j];
        }
        
        // also handle symmetric case
        if (i != k)
        {
            // multiply source vector by this block
            block_ik.dot(1., cArrayView(src, i * Nspline_proj, Nspline_proj), 0., prod);
            
            // update destination vector
            # pragma omp critical
            {
                # pragma omp simd
                for (std::size_t j  = 0; j < Nspline_proj; j++)
                    dst[k * Nspline_proj + j] += a * prod[j];
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
            res[ispline] += sum;
        }
    }
    
    return res;
}

cArray RadialIntegrals::overlapj (Bspline const & bspline, GaussLegendre const & g, int maxell, const rArrayView vk, std::function<Real(Complex)> weightf) const
{
    // shorthands
    int Nenergy = vk.size();
    int Nspline = bspline.Nspline();
    int Nknot = bspline.Nknot();
    int order = bspline.order();
    
    // reserve space for the output array
    std::size_t size = Nspline * Nenergy * (maxell + 1);
    cArray res (size);
    
    // per interval
    int points = EXPANSION_QUADRATURE_POINTS;
    
    // quadrature weights and nodes
    cArray xs (points), ws (points);
    
    // for all knots
    # pragma omp parallel for firstprivate (xs,ws)
    for (int iknot = 0; iknot < Nknot - 1; iknot++)
    {
        // skip zero length intervals
        if (bspline.t(iknot) == bspline.t(iknot+1))
            continue;
        
        // which points are to be used here?
        g.scaled_nodes_and_weights(points, bspline.t(iknot), bspline.t(iknot+1), &xs[0], &ws[0]);
        
        // evaluate relevant B-splines on this knot
        cArrays evalB(Nspline);
        for (int ispline = std::max(iknot-order,0); ispline < Nspline and ispline <= iknot; ispline++)
        {
            evalB[ispline] = cArray(points);
            bspline.B(ispline, iknot, points, xs.data(), evalB[ispline].data());
        }
        
        // for all linear momenta (= energies)
        for (int ie = 0; ie < Nenergy; ie++)
        {
            // evaluate the Riccati-Bessel function for this knot and energy and for all angular momenta
            cArrays evalj(points);
            for (int ipoint = 0; ipoint < points; ipoint++)
            {
                // compute the damping factor
                Real damp = weightf(xs[ipoint]);
                
                // if the factor is numerical zero, do not evaluate the function, just allocate zeros
                if (damp == 0)
                {
                    evalj[ipoint] = cArray(maxell + 1);
                    continue;
                }
                
                // evaluate all Riccati-Bessel functions in point
                evalj[ipoint] = damp * special::ric_jv(maxell, vk[ie] * xs[ipoint]);
                
                // clear all possible NaN entries (these may occur for far radii, where should be zero)
                for (int l = 0; l <= maxell; l++) if (not Complex_finite(evalj[ipoint][l]))
                    evalj[ipoint][l] = 0.;
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
                        sum += ws[ipoint] * evalj[ipoint][l] * evalB[ispline][ipoint];
                    
                    // store the overlap; keep the shape Nmomenta × Nspline × (maxl+1)
                    res[(ie * (maxell + 1) + l) * Nspline + ispline] += sum;
                }
            }
        }
    }
    
    return res;
}
