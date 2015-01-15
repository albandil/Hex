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

inline double damp (Complex r, Complex R)
{
    // if sufficiently far, return clean zero
    if (r.real() > R.real())
        return 0.;
    
    // else damp using tanh(x) distribution
    return std::tanh(0.125 * (R.real() - r.real()));
}

rArray RadialIntegrals::computeScale (int lambda, int iknotmax) const
{
    // get last knot (end of last interval)
    if (iknotmax == 0)
        iknotmax = bspline_.Nknot() - 1;
    
    // quadrature order
    // NOTE : must match that in RadialIntegrals::computeMi !
    int Npts = std::max(2, bspline_.order() + std::abs(lambda) + 1);
    
    // output arrays
    rArray data (iknotmax);
    
    // for all knots
    for (int iknot = 0; iknot < iknotmax; iknot++)
    {
        // skip zero-length intervals
        if (bspline_.t(iknot) == bspline_.t(iknot + 1))
            continue;
        
        // get (real part of) the left-most Gauss-Legendre point
        double rho = g_.p_points(Npts, bspline_.t(iknot), bspline_.t(iknot+1))[0].real();
        
        // compute logarithms of the scale factors
        data[iknot] = log(rho);
    }
    
    return data;
}

void RadialIntegrals::M_integrand
(
    int n,
    Complex * const restrict in,
    Complex * const restrict out,
    int i, int j, int a,
    int iknot, int iknotmax,
    double & logscale
) const
{
    // extract data
    Complex R = bspline_.t(iknotmax);
    
    // evaluate B-splines
    Complex values_i[n], values_j[n];
    bspline_.B(i, iknot, n, in, values_i);
    bspline_.B(j, iknot, n, in, values_j);
    
    // all evaluations produced finite results
    bool all_finite = true;
    
    // fill output array
    if (R != 0.)
    {
        for (int k = 0; k < n; k++)
        {
            out[k] = values_i[k] * values_j[k] * pow(in[k],a) * damp(in[k],R);
            
            if (not (all_finite = Complex_finite(out[k])))
                break;
        }
    }
    else
    {
        for (int k = 0; k < n; k++)
        {
            out[k] = values_i[k] * values_j[k] * pow(in[k],a);
            
            if (not (all_finite = Complex_finite(out[k])))
                break;
        }
    }
    
    //
    // check that all elements are finite
    //
    
    if (not all_finite)
    {
        // compute logarithms of the integrand
        if (R != 0.)
        {
            for (int k = 0; k < n; k++)
            {
                out[k] = std::log(values_i[k]) + std::log(values_j[k]) + double(a) * std::log(in[k]) + std::log(damp(in[k],R));
                
                // use newly computed value as scale, if larger than the current value
                if (out[k].real() > logscale)
                    logscale = out[k].real();
            }
        }
        else
        {
            for (int k = 0; k < n; k++)
            {
                out[k] = std::log(values_i[k]) + std::log(values_j[k]) + double(a) * std::log(in[k]);
                
                // use newly computed value as scale, if larger than the current value
                if (out[k].real() > logscale)
                    logscale = out[k].real();
            }
        }
        
        // scale values by subtracting logarithms and exponentialize, so that they can be integrated by weighed summing
        for (int k = 0; k < n; k++)
            out[k] = std::exp(Complex(out[k].real() - logscale, out[k].imag()));
    }
}

cArray RadialIntegrals::computeMi (int a, int iknotmax) const
{
    int Nspline = bspline_.Nspline();
    int order = bspline_.order();
    
    // (logarithms of) partial integral moments
    cArray m
    (
        Nspline * (2 * order + 1) * (order + 1),
        Complex(0.,special::constant::Inf)
    );
    
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
                double logscale = 0.; // logarithm of the scale
                
                // use at least 2nd order
                int points = std::max(2, order + std::abs(a) + 1);
                
                // integrate
                integral = g_.quadMFP
                (
                    this, &RadialIntegrals::M_integrand,       // integrand pointer
                    points, iknot, xa, xb,                     // integration parameters
                    i, j, a, iknot, iknotmax, logscale         // data to pass to the integrator
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
                if (integral != 0.)
                {
                    Complex lg = std::log(integral) + logscale;
                    m[(x_1 * (2 * order + 1) + y_1) * (order + 1) + z_1] = lg;
                    m[(x_2 * (2 * order + 1) + y_2) * (order + 1) + z_2] = lg;
                }
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
    cArray xs = g_.p_points(points, x1, x2);
    cArray ws = g_.p_weights(points, x1, x2);
    
    // evaluate B-splines at Gauss-Legendre nodes
    Complex values_i[points], values_j[points];
    bspline_.dB(i, iknot, points, xs.data(), values_i);
    bspline_.dB(j, iknot, points, xs.data(), values_j);
    
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
    cArray xs = g_.p_points(points, x1, x2);
    cArray ws = g_.p_weights(points, x1, x2);
    
    // evaluate B-splines at Gauss-Legendre nodes
    Complex values_i[points], values_j[points];
    bspline_.B(i, iknot, points, xs.data(), values_i);
    bspline_.B(j, iknot, points, xs.data(), values_j);
    
    // result
    Complex res = 0;
    
    // accumulate the (damped) result
    if (R != 0.)
    {
        for (int k = 0; k < points; k++)
            res += values_i[k] * values_j[k] * pow(xs[k],a) * ws[k] * damp(xs[k],R);
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

void RadialIntegrals::setupOneElectronIntegrals ()
{
    // create file names for this radial integrals
    char D_name[20], S_name[20], Mm1_name[20], Mm1_tr_name[20], Mm2_name[20];
    std::snprintf(D_name,      sizeof(D_name),      "%d-D.hdf",      bspline_.order());
    std::snprintf(S_name,      sizeof(S_name),      "%d-S.hdf",      bspline_.order());
    std::snprintf(Mm1_name,    sizeof(Mm1_name),    "%d-Mm1.hdf",    bspline_.order());
    std::snprintf(Mm1_tr_name, sizeof(Mm1_tr_name), "%d-Mm1_tr.hdf", bspline_.order());
    std::snprintf(Mm2_name,    sizeof(Mm2_name),    "%d-Mm2.hdf",    bspline_.order());
    
    // load/compute derivative overlaps
    std::cout << "Loading/precomputing derivative overlaps... " << std::flush;
    D_.hdfload(D_name) or D_.populate (
        bspline_.order(), [=](int i, int j) -> Complex { return computeD(i, j, bspline_.Nknot() - 1); }
    ).hdfsave(D_name); D_d_ = D_.torow();
    
    // load/compute integral moments
    std::cout << "ok\n\nLoading/precomputing integral moments... " << std::flush;
    S_.hdfload(S_name) or S_.populate (
        bspline_.order(), [=](int m, int n) -> Complex { return computeM(0, m, n); }
    ).hdfsave(S_name); S_d_ = S_.torow();
    Mm1_.hdfload(Mm1_name) or Mm1_.populate (
        bspline_.order(), [=](int m, int n) -> Complex { return computeM(-1, m, n); }
    ).hdfsave(Mm1_name); Mm1_d_ = Mm1_.torow();
    Mm1_tr_.hdfload(Mm1_tr_name) or Mm1_tr_.populate (
        bspline_.order(),    [=](int m, int n) -> Complex { return computeM(-1, m, n, bspline_.Nreknot() - 1);}
    ).hdfsave(Mm1_tr_name); Mm1_tr_d_ = Mm1_tr_.torow();
    Mm2_.hdfload(Mm2_name) or Mm2_.populate (
        bspline_.order(), [=](int m, int n) -> Complex { return computeM(-2, m, n); }
    ).hdfsave(Mm2_name); Mm2_d_ = Mm2_.torow();
    std::cout << "ok\n\n";
}

void RadialIntegrals::setupTwoElectronIntegrals (Parallel const & par, CommandLine const & cmd, Array<bool> const & lambdas)
{
    // set number of two-electron integrals
    R_tr_dia_.resize(lambdas.size());
    
    // abandon their computation, if not necessary
    if (cmd.lightweight)
        return;
    
    // shorthands
    int Nspline = bspline_.Nspline();
    
    // allocate storage and associate names
    for (unsigned lambda = 0; lambda < lambdas.size(); lambda++)
    {
        bool keep_in_memory = ((par.isMyWork(lambda) and cmd.cache_all_radint) or cmd.cache_all_radint);
        
        R_tr_dia_[lambda] = BlockSymDiaMatrix
        (
            Nspline,            // block count (and size)
            S_.nzpattern(),     // block structure
            S_.diag(),          // non-zero diagonals
            keep_in_memory,     // whether to keep in memory
            format("%d-R_tr_dia_%d.hdf", bspline_.order(), lambda) // HDF scratch disk file name
        );
    }
    
    if (!cmd.gpu_slater)
    {
        // print information
        std::cout << "Precomputing multipole integrals (lambda = 0 .. " << lambdas.size() - 1 << ")." << std::endl;
    }
    
    // for all multipoles : compute / load
    for (int lambda = 0; lambda < (int)lambdas.size(); lambda++)
    {
        // this process will only compute a subset of radial integrals
//         if (not par.isMyWork(lambda))
//             continue;
        
        // look for precomputed data on disk
        if (R_tr_dia_[lambda].hdfcheck())
        {
            if (not par.isMyWork(lambda) or not cmd.cache_own_radint)
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
        
        // logarithms of partial integral moments
        cArray Mtr_L, Mtr_mLm1;
        
        // compute partial moments
        Mtr_L    = std::move(computeMi( lambda,   bspline_.Nreknot() - 1));
        Mtr_mLm1 = std::move(computeMi(-lambda-1, bspline_.Nreknot() - 1));
        
        # pragma omp parallel firstprivate (lambda, Mtr_L, Mtr_mLm1)
        {
            // get recursive structure
            std::vector<std::pair<int,int>> const & structure = R_tr_dia_[lambda].structure();
            
            // for all blocks of the radial matrix
            # pragma omp for schedule (dynamic,1)
            for (unsigned iblock = 0; iblock < structure.size(); iblock++)
            {
                // block indices
                int i = structure[iblock].first;
                int k = structure[iblock].second;
                
                // create a new block of the radial integral matrix
                cArray block (structure.size());
                
                // for all elements in the symmetrical block
                for (unsigned n = 0; n < structure.size(); n++)
                {
                    // element indices
                    int j = structure[n].first;
                    int l = structure[n].second;
                    
                    // evaluate 2-D integral of Bi(1)Bj(2)V(1,2)Bk(1)Bl(2)
                    block[n] = computeR(lambda, i, j, k, l, Mtr_L, Mtr_mLm1);
                }
                
                // write the finished block to disk
                # pragma omp critical
                R_tr_dia_[lambda].setBlock(iblock, block);
            }
        }
        
        std::cout << "\t- integrals for lambda = " << lambda << " computed" << std::endl;
    }

    std::cout << std::endl;
}

cArrays RadialIntegrals::apply_R_matrix (unsigned lambda, cArrays const & src) const
{
    // logarithms of partial integral moments
    cArray Mtr_L    = std::move(computeMi( lambda,   bspline_.Nreknot() - 1));
    cArray Mtr_mLm1 = std::move(computeMi(-lambda-1, bspline_.Nreknot() - 1));
    
    // number of source vectors
    int N = src.size();
    
    // size of a block
    std::size_t chunk = bspline_.Nspline();
    
    // output array (the product)
    cArrays dst(N);
    for (int n = 0; n < N; n++)
        dst[n] = cArray(src[n].size());
    
    // get recursive structure
    std::vector<std::pair<int,int>> structure = S_.nzpattern();
    
    // for all blocks of the radial matrix
    # pragma omp parallel for firstprivate (structure, lambda, Mtr_L, Mtr_mLm1) schedule (dynamic, 1)
    for (unsigned iblock = 0; iblock < structure.size(); iblock++)
    {
        // block indices
        int i = structure[iblock].first;
        int k = structure[iblock].second;
        
        // (i,k)-block data (= concatenated non-zero upper diagonals)
        cArray block_ik (structure.size());
        
        // for all elements in the symmetrical block
        for (unsigned n = 0; n < structure.size(); n++)
        {
            // element indices
            int j = structure[n].first;
            int l = structure[n].second;
            
            // evaluate 2-D integral of Bi(1)Bj(2)V(1,2)Bk(1)Bl(2)
            block_ik[n] = computeR(lambda, i, j, k, l, Mtr_L, Mtr_mLm1);
        }
        
        // multiply all source vectors by this block
        # pragma omp critical
        for (int isrc = 0; isrc < N; isrc++)
        {
            // source and destination segment
            dst[isrc].slice(i * chunk, (i + 1) * chunk) += SymDiaMatrix::sym_dia_dot
            (
                chunk,
                S_.diag(),
                block_ik.data(),
                src[isrc].data() + k * chunk
            );
            
            // take care of symmetric position of the off-diagonal block
            if (i != k)
            {
                dst[isrc].slice(k * chunk, (k + 1) * chunk) += SymDiaMatrix::sym_dia_dot
                (
                    chunk,
                    S_.diag(),
                    block_ik.data(),
                    src[isrc].data() + i * chunk
                );
            }
        }
    }
    
    // return result
    return dst;
}
