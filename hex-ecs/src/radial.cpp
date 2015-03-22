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
#include <complex>
#include <cmath>
#include <cstdio>
#include <cstring>

#include "arrays.h"
#include "bspline.h"
#include "gauss.h"
#include "parallel.h"
#include "radial.h"

void RadialIntegrals::Mi_integrand
(
    int n,
    Complex * const restrict in,
    Complex * const restrict out,
    int i, int j, int a,
    int iknot, int iknotmax
) const
{
    // extract data
    Complex R = bspline_.t(iknotmax);
    
    // evaluate B-splines
    Complex values_i[n], values_j[n];
    bspline_.B(i, iknot, n, in, values_i);
    bspline_.B(j, iknot, n, in, values_j);
    
    // get upper bound
    double t = bspline_.t()[iknot + 1].real();
    
    // scale factor for the multipole
    double scalef = (t < 1 ? 1/t : 1);
    
    // fill output array
    if (R != 0.)
    {
        // use damping
        for (int k = 0; k < n; k++)
            out[k] = values_i[k] * values_j[k] * pow(scalef*in[k],a) * damp(in[k],0.,R);
    }
    else
    {
        // do not use damping
        for (int k = 0; k < n; k++)
            out[k] = values_i[k] * values_j[k] * pow(scalef*in[k],a);
    }
}

cArray RadialIntegrals::computeMi (int a, int iknotmax) const
{
    int Nspline = bspline_.Nspline();
    int order = bspline_.order();
    
    // partial integral moments
    cArray m (Nspline * (2 * order + 1) * (order + 1));
    
    // for all B-splines
    for (int i = 0; i < (int)Nspline; i++)
    {
        // for all B-splines (use symmetry)
        for (int j = i; j <= i + (int)order and j < (int)Nspline; j++)
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
                Complex xa = bspline_.t(iknot);
                Complex xb = bspline_.t(iknot+1);
                
                // throw away zero length intervals
                if (xa == xb)
                    continue;
                
                // results of the quadrature
                Complex integral;
                
                // use at least 2nd order
                int points = std::max(2, order + std::abs(a) + 1);
                
                // integrate
                integral = g_.quadMFP
                (
                    this, &RadialIntegrals::Mi_integrand,      // integrand pointer
                    points, iknot, xa, xb,                     // integration parameters
                    i, j, a, iknot, iknotmax                   // data to pass to the integrator
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


Complex RadialIntegrals::computeD_iknot (int i, int j, int iknot) const
{
    if (iknot < 0)
        iknot = bspline_.Nknot() - 1;
    
    // get interval boundaries
    Complex x1 = bspline_.t(iknot);
    Complex x2 = bspline_.t(iknot + 1);
    
    // throw away zero-length intervals
    if (x1 == x2)
        return 0;
    
    // get Gauss-Legendre nodes and weights for the interval [-1, 1]
    // - use at least 2nd order
    int points = std::max (2, bspline_.order());
    Complex xs[points], ws[points];
    g_.scaled_nodes_and_weights(points, x1, x2, xs, ws);
    
    // evaluate B-splines at Gauss-Legendre nodes
    Complex values_i[points], values_j[points];
    bspline_.dB(i, iknot, points, xs, values_i);
    bspline_.dB(j, iknot, points, xs, values_j);
    
    // result
    Complex res = 0;
    
    // accumulate the result
    for (int k = 0; k < points; k++)
        res += values_i[k] * values_j[k] * ws[k];
    
    return res;
}

Complex RadialIntegrals::computeD (int i, int j, int maxknot) const
{
    // get boundary iknots
    int left = std::max(i, j);
    int right = std::min(i, j) + bspline_.order();
    
    // cut at maxknot
    if (right > maxknot)
        right = maxknot;
    
    // the result
    Complex res = 0;
    
    // undergo integration on sub-intervals
    for (int iknot = left; iknot <= right; iknot++)
        res += computeD_iknot(i, j, iknot);
        
    return res;
}

Complex RadialIntegrals::computeM_iknot (int a, int i, int j, int iknot, Complex R) const
{
    // get interval boundaries
    Complex x1 = bspline_.t(iknot);
    Complex x2 = bspline_.t(iknot + 1);
    
    // throw away zero-length intervals
    if (x1 == x2)
        return 0;
    
    // get Gauss-Legendre nodes and weights for the interval [-1, 1]
    // - use at least 2nd order
    int points = std::max (2, bspline_.order() + std::abs(a) + 1);
    Complex xs[points], ws[points];
    g_.scaled_nodes_and_weights(points, x1, x2, xs, ws);
    
    // evaluate B-splines at Gauss-Legendre nodes
    Complex values_i[points], values_j[points];
    bspline_.B(i, iknot, points, xs, values_i);
    bspline_.B(j, iknot, points, xs, values_j);
    
    // result
    Complex res = 0;
    
    // accumulate the (damped) result
    if (R != 0.)
    {
        for (int k = 0; k < points; k++)
            res += values_i[k] * values_j[k] * pow(xs[k],a) * ws[k] * damp(xs[k],0.,R);
    }
    else
    {
        for (int k = 0; k < points; k++)
            res += values_i[k] * values_j[k] * pow(xs[k],a) * ws[k];
    }
    
    return res;
}

Complex RadialIntegrals::computeM (int a, int i, int j, int maxknot) const
{
    // get boundary iknots
    int left = std::max(i, j);
    int right = std::min(i, j) + bspline_.order();
    
    // cut at maxknot
    if (maxknot != 0 and right > maxknot)
        right = maxknot;
    
    // the result
    Complex res = 0;
    
    // undergo integration on sub-intervals
    for (int iknot = left; iknot <= right; iknot++)
        res += computeM_iknot(a, i, j, iknot, bspline_.t(maxknot));
    
    return res;
}

void RadialIntegrals::setupOneElectronIntegrals (Parallel const & par, CommandLine const & cmd)
{
    // shorthands
    int order   = bspline_.order();
    int Nknot   = bspline_.Nknot();
    int Nreknot = bspline_.Nreknot();
    
    // create file names for this radial integrals
    D_.hdflink(format("%d-D.hdf", order));
    S_.hdflink(format("%d-S.hdf", order));
    Mm1_.hdflink(format("%d-Mm1.hdf", order));
    Mm1_tr_.hdflink(format("%d-Mm1_tr.hdf", order));
    Mm2_.hdflink(format("%d-Mm2.hdf", order));
    
    std::cout << "Precomputing one-electron matrices... " << std::flush;
    
    // compute one-electron matrices
    D_.populate([=](int m, int n) -> Complex { return computeD(m, n, Nknot - 1); });
    S_.populate([=](int m, int n) -> Complex { return computeM(0, m, n); });
    Mm1_.populate([=](int m, int n) -> Complex { return computeM(-1, m, n); });
    Mm1_tr_.populate([=](int m, int n) -> Complex { return computeM(-1, m, n, Nreknot - 1);});
    Mm2_.populate([=](int m, int n) -> Complex { return computeM(-2, m, n); });
    
    // save the matrices to disk
    if (not cmd.shared_scratch or par.IamMaster())
    {
        D_.hdfsave();
        S_.hdfsave();
        Mm1_.hdfsave();
        Mm1_tr_.hdfsave();
        Mm2_.hdfsave();
    }
    
    std::cout << "ok" << std::endl << std::endl;
}

void RadialIntegrals::setupTwoElectronIntegrals (Parallel const & par, CommandLine const & cmd, Array<bool> const & lambdas)
{
    // shorthands
    int Nspline = bspline_.Nspline();
    int order = bspline_.order();
    
    // set number of two-electron integrals
    R_tr_dia_.resize(lambdas.size());
    
    // print information
    std::cout << "Precomputing partial integral moments." << std::endl << std::endl;
    
    // compute partial moments
    std::size_t mi_size = Nspline * (2 * order + 1) * (order + 1);
    Mitr_L_.resize((lambdas.size() + 1) * mi_size);
    Mitr_mLm1_.resize((lambdas.size() + 1) * mi_size);
    for (int lambda = 0; lambda < (int)lambdas.size(); lambda++)
    {
        cArrayView(Mitr_L_, lambda * mi_size, mi_size) = computeMi(lambda, bspline_.Nreknot() - 1);
        cArrayView(Mitr_mLm1_, lambda * mi_size, mi_size) = computeMi(-lambda-1, bspline_.Nreknot() - 1);
    }
    
    // abandon their computation, if not necessary
    if (cmd.lightweight_radial_cache)
        return;
    
    // allocate storage and associate names
    for (unsigned lambda = 0; lambda < lambdas.size(); lambda++)
    {
        bool keep_in_memory = ((par.isMyWork(lambda) and cmd.cache_all_radint) or cmd.cache_all_radint);
        
        R_tr_dia_[lambda] = BlockSymBandMatrix
        (
            Nspline,            // block count (and size)
            order + 1,          // half-bandwidth
            keep_in_memory,     // whether to keep in memory
            format("%d-R_tr_dia_%d.hdf", bspline_.order(), lambda) // HDF scratch disk file name
        );
    }
    
    // print information
    std::cout << "Precomputing multipole integrals (lambda = 0 .. " << lambdas.size() - 1 << ")." << std::endl;
    
    // for all multipoles : compute / load
    for (int lambda = 0; lambda < (int)lambdas.size(); lambda++)
    {
        // if the radial integrals are shared, this process will only compute the owned subset of radial integrals
        if (cmd.shared_scratch and not par.isMyWork(lambda))
            continue;
        
        // look for precomputed data on disk
        if (R_tr_dia_[lambda].hdfcheck())
        {
            if (/*not par.isMyWork(lambda) or*/ not cmd.cache_own_radint)
            {
                std::cout << "\t- integrals for lambda = " << lambda << " present in \"" << R_tr_dia_[lambda].hdfname() << "\"\n";
                continue;
            }
            
            if (R_tr_dia_[lambda].hdfload())
            {
                std::cout << "\t- integrals for lambda = " << lambda << " loaded from \"" << R_tr_dia_[lambda].hdfname() << "\"\n";
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
            for (int i = 0; i < Nspline; i++)
            for (int d = 0; d <= order; d++)
            if (i + d < Nspline)
            {
                // calculate the block
                SymBandMatrix block = calc_R_tr_dia_block(lambda, i, i + d);
                
                // write the finished block to disk
                # pragma omp critical
                R_tr_dia_[lambda].setBlock(i * (order + 1) + d, block.data());
            }
        }
        
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
        for (int lambda = 0; lambda < (int)lambdas.size(); lambda++)
        {
            // skip own data (already loaded since calculation)
            if (par.isMyWork(lambda))
                continue;
            
            // load radial integrals
            if (R_tr_dia_[lambda].hdfload())
                std::cout << "\t- integrals for lambda = " << lambda << " loaded from shared file \"" << R_tr_dia_[lambda].hdfname() << "\"\n";
            else
                HexException("Can't read shared radial integral file \"%s\".", R_tr_dia_[lambda].hdfname().c_str());
        }
    }
    
    std::cout << std::endl;
}

SymBandMatrix RadialIntegrals::calc_R_tr_dia_block (unsigned int lambda, int i, int k, bool simple) const
{
    // shorthands
    int Nspline = bspline_.Nspline();
    int order = bspline_.order();
    
    // (i,k)-block data
    SymBandMatrix block_ik (Nspline, order + 1);
    
    // for all elements in the symmetrical block : evaluate 2-D integral of Bi(1)Bj(2)V(1,2)Bk(1)Bl(2)
    for (int j = 0; j < Nspline; j++)
    for (int l = j; l < Nspline and l <= j + order; l++)
        block_ik(j,l) = computeR(lambda, i, j, k, l, simple);
    
    return block_ik;
}

cArray RadialIntegrals::apply_R_matrix (unsigned lambda, cArray const & src, bool simple) const
{
    // number of source vectors
    int N = src.size();
    
    // shorthands
    std::size_t Nspline = bspline_.Nspline();
    std::size_t order = bspline_.order();
    
    // output array (the product)
    cArray dst(N);
    
    // for all blocks of the radial matrix
    # pragma omp parallel for firstprivate (lambda) schedule (dynamic, 1)
    for (unsigned i = 0; i < Nspline; i++)
    for (std::size_t k = i; k < Nspline and k <= i + order; k++)
    {
        // (i,k)-block data (= concatenated non-zero upper diagonals)
        SymBandMatrix block_ik = calc_R_tr_dia_block(lambda, i, k, simple);
        
        // multiply source vector by this block
        # pragma omp critical
        {
            // source and destination segment
            cArrayView(dst, i * Nspline, Nspline) += block_ik.dot(cArrayView(src, k * Nspline, Nspline));
            
            // take care also of symmetric position of the off-diagonal block
            if (i != k)
            cArrayView(dst, k * Nspline, Nspline) += block_ik.dot(cArrayView(src, i * Nspline, Nspline));
        }
    }
    
    // return result
    return dst;
}