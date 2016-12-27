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

#include <iostream>

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"
#include "hex-csrmatrix.h"
#include "hex-misc.h"
#include "hex-openmp.h"
#include "hex-symbandmatrix.h"

// --------------------------------------------------------------------------------- //

#include "gauss.h"
#include "parallel.h"
#include "radial.h"

// --------------------------------------------------------------------------------- //

#include "NoPreconditioner.h"
#include "CGPreconditioner.h"

// --------------------------------------------------------------------------------- //

NoPreconditioner::NoPreconditioner ()
  : PreconditionerBase(),
    E_(0), cmd_(nullptr), par_(nullptr), inp_(nullptr), ang_(nullptr), rad_(nullptr)
{
    // nothing to do
}

NoPreconditioner::NoPreconditioner
(
    CommandLine  const & cmd,
    InputFile    const & inp,
    Parallel     const & par,
    AngularBasis const & ang,
    Bspline const & bspline_x_inner,
    Bspline const & bspline_x_full,
    Bspline const & bspline_y_inner,
    Bspline const & bspline_y_full
) : PreconditionerBase(),
    E_(0), cmd_(&cmd), par_(&par), inp_(&inp), ang_(&ang),
    A_blocks_ (ang.states().size() * ang.states().size()),
    B1_blocks_(ang.states().size() * ang.states().size()),
    B2_blocks_(ang.states().size() * ang.states().size()),
    Cu_blocks_(ang.states().size() * ang.states().size()),
    Cl_blocks_(ang.states().size() * ang.states().size()),
    rad_
    (
        new RadialIntegrals
        (
            bspline_x_inner, bspline_x_full,
            bspline_y_inner, bspline_y_full,
            ang.maxlambda() + 1
        )
    )
{
    // nothing to do
}

NoPreconditioner::~NoPreconditioner ()
{
    if (rad_)
        delete rad_;
}

std::string NoPreconditioner::description () const
{
    return "\"Preconditioning\" by the identity matrix.";
}

void NoPreconditioner::setup ()
{
    rad_->setupOneElectronIntegrals(*par_, *cmd_);
    rad_->setupTwoElectronIntegrals(*par_, *cmd_);
}

BlockSymBandMatrix<Complex> NoPreconditioner::calc_A_block (int ill, int illp, bool twoel) const
{
    // B-spline order
    int order = inp_->order;
    
    // B-spline count
    int Nspline_inner_x = rad_->bspline_inner_x().Nspline();
    int Nspline_inner_y = rad_->bspline_inner_y().Nspline();
    
    // other shorthands
    Real Zp = inp_->Zp;
    
    // angular momenta
    int l1 = ang_->states()[ill].first;
    int l2 = ang_->states()[ill].second;
    
    // inner region matrix
    BlockSymBandMatrix<Complex> A
    (
        Nspline_inner_x,            // block count
        order + 1,                  // block structure half-bandwidth
        Nspline_inner_y,            // block size
        order + 1,                  // block half-bandwidth
        !cmd_->outofcore,           // keep in memory?
        format("blk-A-%d-%d.ooc", ill, illp)  // scratch disk file name
    );
    
    // sub-block of the inner region matrix
    SymBandMatrix<Complex> subblock (Nspline_inner_y, order + 1);
    
    // for all sub-blocks
    # pragma omp parallel for firstprivate (subblock) if (!cmd_->outofcore)
    for (int i = 0; i < Nspline_inner_x; i++)
    for (int d = 0; d <= order; d++)
    if (i + d < Nspline_inner_x)
    {
        int k = i + d;
        
        subblock.data().fill(0.0_z);
        
        Complex Sik = rad_->S_full_x(i,k);
        Complex Dik = rad_->D_full_x(i,k);
        Complex Mm1ik = rad_->Mm1_tr_full_x(i,k);
        Complex Mm2ik = rad_->Mm2_full_x(i,k);
        
        // one-electron part
        if (ill == illp)
        {
            subblock.populate
            (
                [&](int j, int l)
                {
                    Complex Sjl = rad_->S_full_y(j,l);
                    Complex Djl = rad_->D_full_y(j,l);
                    Complex Mm1jl = rad_->Mm1_tr_full_y(j,l);
                    Complex Mm2jl = rad_->Mm2_full_y(j,l);
                    
                    return E_ * Sik * Sjl
                            - (0.5_z * Dik - Mm1ik) * Sjl
                            - 0.5_r * l1 * (l1 + 1) * Mm2ik * Sjl
                            - Sik * (0.5_z * Djl + Zp * Mm1jl)
                            - 0.5_r * l2 * (l2 + 1) * Sik * Mm2jl;
                }
            );
        }
        
        // two-electron part
        if(twoel)
        for (int lambda = 0; lambda <= rad_->maxlambda(); lambda++) if (ang_->f(ill,illp,lambda) != 0)
        {
            // calculate two-electron term
            if (not cmd_->lightweight_radial_cache)
            {
                // use precomputed block from scratch file or from memory
                subblock.data() += Zp * ang_->f(ill,illp,lambda) * rad_->R_tr_dia(lambda).getBlock(i * (order + 1) + d).slice(0, Nspline_inner_y * (order + 1));
            }
            else
            {
                // compute the data anew
                subblock.data() += Zp * ang_->f(ill,illp,lambda) * rad_->calc_R_tr_dia_block(lambda, i, k).data().slice(0, Nspline_inner_y * (order + 1));
            }
        }
        
        // save block
        A.setBlock(i * (order + 1) + d, subblock.data());
    }
    
    if (cmd_->outofcore)
        A.drop();
    
    return A;
}

void NoPreconditioner::update (Real E)
{
    // shorthands
    Real Zp = inp_->Zp;
    int order = inp_->order;
    int Nang = ang_->states().size();
    int Nspline_inner_x = rad_->bspline_inner_x().Nspline();
    int Nspline_inner_y = rad_->bspline_inner_y().Nspline();
    int Nspline_full_x  = rad_->bspline_full_x().Nspline();
    int Nspline_full_y  = rad_->bspline_full_y().Nspline();
    int Nspline_outer_x = Nspline_full_x - Nspline_inner_x;
    int Nspline_outer_y = Nspline_full_y - Nspline_inner_y;
    std::size_t A_size = std::size_t(Nspline_inner_x) * std::size_t(Nspline_inner_y);
    
    // update energy
    E_ = E;
    
    // get maximal asymptotic principal quantum number
    max_n_ = (E_ >= 0 ? 0 : 1.0 / std::sqrt(-2 * E_));
    
    // update number of asymptotic channels
    Nchan_.clear();
    for (int ill = 0; ill < Nang; ill++)
    {
        int l1 = ang_->states()[ill].first;
        int l2 = ang_->states()[ill].second;
        
        // number of channels when r1 -> inf (i.e. second electron is bound)
        int Nchan1 = (inp_->Zp > 0 ? 0 : std::max(0, max_n_ - l2));
        
        // number of channels when r2 -> inf (i.e. first electron is bound)
        int Nchan2 = std::max(0, max_n_ - l1);
        
        Nchan_.push_back(std::make_pair(Nchan1, Nchan2));
    }
    
    std::cout << "\tPrecompute matrix blocks ... " << std::flush;
    Timer t;
    
    // LU-factorize the overlap matrices
    CsrMatrix<LU_int_t,Complex> csr_Sx = rad_->S_full_x().tocoo<LU_int_t>().tocsr();
    CsrMatrix<LU_int_t,Complex> csr_Sy = rad_->S_full_y().tocoo<LU_int_t>().tocsr();
    std::shared_ptr<LUft> lu_Sx, lu_Sy;
    lu_Sx.reset(LUft::Choose("lapack"));
    lu_Sy.reset(LUft::Choose("lapack"));
    lu_Sx->factorize(csr_Sx);
    lu_Sy->factorize(csr_Sy);
    
    // outer one-electron overlap matrix
    SymBandMatrix<Complex> Sx_outer (Nspline_outer_x, order + 1);
    SymBandMatrix<Complex> Sy_outer (Nspline_outer_y, order + 1);
    Sx_outer.populate([&](int m, int n) { return rad_->S_full_x()(Nspline_inner_x + m, Nspline_inner_x + n); });
    Sy_outer.populate([&](int m, int n) { return rad_->S_full_y()(Nspline_inner_y + m, Nspline_inner_y + n); });
    
    // outer one-electron derivative matrix
    SymBandMatrix<Complex> Dx_outer (Nspline_outer_x, order + 1);
    SymBandMatrix<Complex> Dy_outer (Nspline_outer_y, order + 1);
    Dx_outer.populate([&](int m, int n) { return rad_->D_full_x()(Nspline_inner_x + m, Nspline_inner_x + n); });
    Dy_outer.populate([&](int m, int n) { return rad_->D_full_y()(Nspline_inner_y + m, Nspline_inner_y + n); });
    
    // outer one-electron centrifugal moment matrix
    SymBandMatrix<Complex> Mm2x_outer (Nspline_outer_x, order + 1);
    SymBandMatrix<Complex> Mm2y_outer (Nspline_outer_y, order + 1);
    Mm2x_outer.populate([&](int m, int n) { return rad_->Mm2_full_x()(Nspline_inner_x + m, Nspline_inner_x + n); });
    Mm2y_outer.populate([&](int m, int n) { return rad_->Mm2_full_y()(Nspline_inner_y + m, Nspline_inner_y + n); });
    
    // outer one-electron multipole moment matrices
    std::vector<SymBandMatrix<Complex>> Mtr_mLm1_outer_x, Mtr_mLm1_outer_y;
    for (int lambda = 0; lambda <= rad_->maxlambda(); lambda++)
    {
        Mtr_mLm1_outer_x.push_back(SymBandMatrix<Complex>(Nspline_outer_x, order + 1));
        Mtr_mLm1_outer_y.push_back(SymBandMatrix<Complex>(Nspline_outer_y, order + 1));
        
        Mtr_mLm1_outer_x.back().populate
        (
            [ & ] (int m, int n)
            {
                return rad_->Mtr_mLm1_full_x(lambda)(Nspline_inner_x + m, Nspline_inner_x + n)
                     * special::pow_int(rad_->bspline_full_x().t(Nspline_inner_x + std::min(m,n) + order + 1).real(), -lambda-1);
            }
        );
        Mtr_mLm1_outer_y.back().populate
        (
            [ & ] (int m, int n)
            {
                return rad_->Mtr_mLm1_full_y(lambda)(Nspline_inner_y + m, Nspline_inner_y + n)
                     * special::pow_int(rad_->bspline_full_y().t(Nspline_inner_y + std::min(m,n) + order + 1).real(), -lambda-1);
            }
        );
    }
    
    // setup blocks
    for (int ill = 0; ill < Nang; ill++) if (par_->isMyGroupWork(ill))
    for (int illp = 0; illp < Nang; illp++)
    {
        // angular momenta
        int l1 = ang_->states()[ill].first;
        int l2 = ang_->states()[ill].second;
        int l1p = ang_->states()[illp].first;
        int l2p = ang_->states()[illp].second;
        
        // get number of asymptotic bound channels
        int Nchan1 = Nchan_[ill].first;     // # r1 -> inf, l2 bound
        int Nchan2 = Nchan_[ill].second;    // # r2 -> inf, l1 bound
        int Nchan1p = Nchan_[illp].first;   // # r1 -> inf, l2p bound
        int Nchan2p = Nchan_[illp].second;  // # r2 -> inf, l1p bound
        
        // calculate all missing hydrogen overlaps
        while (Spx.size() <= (unsigned)std::max(l1,l1p)) Spx.push_back(cArrays());
        while (Spy.size() <= (unsigned)std::max(l2,l2p)) Spy.push_back(cArrays());
        while (Xpx.size() <= (unsigned)std::max(l1,l1p)) Xpx.push_back(cArrays());
        while (Xpy.size() <= (unsigned)std::max(l2,l2p)) Xpy.push_back(cArrays());
        while (Spx[l1] .size() < (unsigned)Nchan2 ) Spx[l1] .push_back(cArray());
        while (Spy[l2] .size() < (unsigned)Nchan1 ) Spy[l2] .push_back(cArray());
        while (Spx[l1p].size() < (unsigned)Nchan2p) Spx[l1p].push_back(cArray());
        while (Spy[l2p].size() < (unsigned)Nchan1p) Spy[l2p].push_back(cArray());
        while (Xpx[l1] .size() < (unsigned)Nchan2 ) Xpx[l1] .push_back(cArray());
        while (Xpy[l2] .size() < (unsigned)Nchan1 ) Xpy[l2] .push_back(cArray());
        while (Xpx[l1p].size() < (unsigned)Nchan2p) Xpx[l1p].push_back(cArray());
        while (Xpy[l2p].size() < (unsigned)Nchan1p) Xpy[l2p].push_back(cArray());
        for (int n = 0; n < Nchan2; n++) if (Spx[l1][n].empty())
        {
            Spx[l1][n] = rad_->overlapP(rad_->bspline_full_x(), rad_->gaussleg_full_x(), l1 + n + 1, l1);
            Xpx[l1][n] = lu_Sx->solve(Spx[l1][n]);
        }
        for (int n = 0; n < Nchan1; n++) if (Spy[l2][n].empty())
        {
            Spy[l2][n] = rad_->overlapP(rad_->bspline_full_y(), rad_->gaussleg_full_y(), l2 + n + 1, l2);
            Xpy[l2][n] = lu_Sy->solve(Spy[l2][n]);
        }
        for (int n = 0; n < Nchan2p; n++) if (Spx[l1p][n].empty())
        {
            Spx[l1p][n] = rad_->overlapP(rad_->bspline_full_x(), rad_->gaussleg_full_x(), l1p + n + 1, l1p);
            Xpx[l1p][n] = lu_Sx->solve(Spx[l1p][n]);
        }
        for (int n = 0; n < Nchan1p; n++) if (Spy[l2p][n].empty())
        {
            Spy[l2p][n] = rad_->overlapP(rad_->bspline_full_y(), rad_->gaussleg_full_y(), l2p + n + 1, l2p);
            Xpy[l2p][n] = lu_Sy->solve(Spy[l2p][n]);
        }
        
        // initialize diagonal block of the inner problem
        // - do not precompute off-diagonal blocks in lightweight mode
        if (not cmd_->lightweight_full or ill == illp)
            A_blocks_[ill * Nang + illp] = calc_A_block(ill, illp);
        
        // create inner-outer coupling blocks
        Cu_blocks_[ill * Nang + illp] = CooMatrix<LU_int_t,Complex>
        (
            A_size + Nchan1  * Nspline_outer_x + Nchan2  * Nspline_outer_y,
            A_size + Nchan1p * Nspline_outer_x + Nchan2p * Nspline_outer_y
        );
        Cl_blocks_[ill * Nang + illp] = CooMatrix<LU_int_t,Complex>
        (
            A_size + Nchan1  * Nspline_outer_x + Nchan2  * Nspline_outer_y,
            A_size + Nchan1p * Nspline_outer_x + Nchan2p * Nspline_outer_y
        );
        
        // setup extended inner-outer problem
        if (not inp_->inner_only)
        {
            // outer problem matrix : r2 -> inf, l1 bound
            B2_blocks_[ill * Nang + illp].resize(Nchan2 * Nchan2p);
            for (int m = 0; m < Nchan2; m++)
            for (int n = 0; n < Nchan2p; n++)
            {
                SymBandMatrix<Complex> subblock (Nspline_outer_y, order + 1);
                
                // channel-diagonal contribution
                if (ill == illp and m == n)
                {
                    subblock += (E_ + 1.0_z / (2.0_z * (l1 + m + 1.0_r) * (l1 + m + 1.0_r))) * Sy_outer
                             - 0.5_z * Dy_outer
                             - 0.5_z * (l2 * (l2 + 1.0_r)) * Mm2y_outer;
                }
                
                // channel-offdiagonal contribution
                for (int lambda = 1; lambda <= rad_->maxlambda(); lambda++) if (ang_->f(ill,illp,lambda) != 0.0_r)
                    subblock += Complex(Zp * ang_->f(ill,illp,lambda) * special::hydro_rho(l1 + m + 1, l1, l1p + n + 1, l1p, lambda)) * Mtr_mLm1_outer_y[lambda];
                
                // use the block
                B2_blocks_[ill * Nang + illp][m * Nchan2p + n].hdflink(format("blk-B2-%d-%d-%d-%d.ooc", ill, illp, m, n));
                B2_blocks_[ill * Nang + illp][m * Nchan2p + n] = std::move(subblock);
                if (cmd_->outofcore)
                {
                    B2_blocks_[ill * Nang + illp][m * Nchan2p + n].hdfsave();
                    B2_blocks_[ill * Nang + illp][m * Nchan2p + n].drop();
                }
            }
            
            // outer problem matrix : r1 -> inf, l2 bound
            B1_blocks_[ill * Nang + illp].resize(Nchan1 * Nchan1p);
            for (int m = 0; m < Nchan1; m++)
            for (int n = 0; n < Nchan1p; n++)
            {
                SymBandMatrix<Complex> subblock (Nspline_outer_x, order + 1);
                
                // channel-diagonal contribution
                if (ill == illp and m == n)
                {
                    subblock += (E_ + 1.0_z / (2.0_z * (l2 + m + 1.0_r) * (l2 + m + 1.0_r))) * Sx_outer
                             - 0.5_z * Dx_outer
                             - 0.5_z * (l1 * (l1 + 1.0_r)) * Mm2x_outer;
                }
                
                // channel-offdiagonal contribution
                for (int lambda = 1; lambda <= rad_->maxlambda(); lambda++) if (ang_->f(ill,illp,lambda) != 0.0_r)
                    subblock += Complex(Zp * ang_->f(ill,illp,lambda) * special::hydro_rho(l2 + m + 1, l2, l2p + n + 1, l2p, lambda)) * Mtr_mLm1_outer_x[lambda];
                
                // use the block
                B1_blocks_[ill * Nang + illp][m * Nchan1p + n].hdflink(format("blk-B1-%d-%d-%d-%d.ooc", ill, illp, m, n));
                B1_blocks_[ill * Nang + illp][m * Nchan1p + n] = std::move(subblock);
                if (cmd_->outofcore)
                {
                    B1_blocks_[ill * Nang + illp][m * Nchan1p + n].hdfsave();
                    B1_blocks_[ill * Nang + illp][m * Nchan1p + n].drop();
                }
            }
            
            // transition area r2 > r1, upper : psi_kl expressed in terms of F_nl for 'l' out of inner area
            for (int i = 0; i < Nspline_inner_x; i++)
            for (int j = 0; j < Nspline_inner_y; j++) // *
            for (int k = std::max(0, i - order); k <= std::min(i + order, Nspline_inner_x - 1); k++)
            for (int l = Nspline_inner_y; l <= j + order; l++)
            for (int n = 0; n < Nchan2p; n++)
            {
                std::size_t row = i * Nspline_inner_y + j;
                std::size_t col = A_size + (Nchan1p + n) * Nspline_outer_y + (l - Nspline_inner_y);
                
                Complex elem = 0; // A_ij,kl Xp_nk
                
                if (ill == illp)
                {
                    elem += E_ * rad_->S_full_x(i,k) * rad_->S_full_y(j,l)
                         - 0.5_r * rad_->D_full_x(i,k) * rad_->S_full_y(j,l)
                         - 0.5_r * rad_->S_full_x(i,k) * rad_->D_full_y(j,l)
                         - 0.5_r * (l1 * (l1 + 1.0_r)) * rad_->Mm2_full_x(i,k) * rad_->S_full_y(j,l)
                         - 0.5_r * (l2 * (l2 + 1.0_r)) * rad_->S_full_x(i,k) * rad_->Mm2_full_y(j,l)
                         + rad_->Mm1_tr_full_x(i,k) * rad_->S_full_y(j,l)
                         - inp_->Zp * rad_->S_full_x(i,k) * rad_->Mm1_tr_full_y(j,l);
                }
                
                Real r1 = rad_->bspline_full_x().t(std::min(i,k) + order + 1).real();
                Real r2 = rad_->bspline_full_y().t(std::min(j,l) + order + 1).real();
                
                for (int lambda = 0; lambda <= rad_->maxlambda(); lambda++) if (ang_->f(ill,illp,lambda) != 0)
                {
                    Real scale = special::pow_int(r1/r2, lambda) / r2;
                    elem += Zp * ang_->f(ill,illp,lambda) * scale * rad_->Mtr_L_full_x(lambda,i,k) * rad_->Mtr_mLm1_full_y(lambda,j,l);
                }
                
                Cu_blocks_[ill * Nang + illp].add(row, col, Xpx[l1p][n][k] * elem);
            }
            
            // transition area r1 > r2, upper : psi_kl expressed in terms of F_nk for 'k' out of inner area
            for (int i = 0; i < Nspline_inner_x; i++)
            for (int j = 0; j < Nspline_inner_y; j++) // *
            for (int k = Nspline_inner_x; k <= i + order; k++)
            for (int l = std::max(0, j - order); l <= std::min(j + order, Nspline_inner_y - 1); l++)
            for (int n = 0; n < Nchan1p; n++)
            {
                std::size_t row = i * Nspline_inner_x + j;
                std::size_t col = A_size + n * Nspline_outer_x + (k - Nspline_inner_x);
                
                Complex elem = 0; // A_ij,kl Xp_nl
                
                if (ill == illp)
                {
                    elem += E_ * rad_->S_full_x(i,k) * rad_->S_full_y(j,l)
                         - 0.5_r * rad_->D_full_x(i,k) * rad_->S_full_y(j,l)
                         - 0.5_r * rad_->S_full_x(i,k) * rad_->D_full_y(j,l)
                         - 0.5_r * (l1 * (l1 + 1.0_r)) * rad_->Mm2_full_x(i,k) * rad_->S_full_y(j,l)
                         - 0.5_r * (l2 * (l2 + 1.0_r)) * rad_->S_full_x(i,k) * rad_->Mm2_full_y(j,l)
                         + rad_->Mm1_tr_full_x(i,k) * rad_->S_full_y(j,l)
                         - Zp * rad_->S_full_x(i,k) * rad_->Mm1_tr_full_y(j,l);
                }
                
                Real r1 = rad_->bspline_full_x().t(std::min(i,k) + order + 1).real();
                Real r2 = rad_->bspline_full_y().t(std::min(j,l) + order + 1).real();
                
                for (int lambda = 0; lambda <= rad_->maxlambda(); lambda++) if (ang_->f(ill,illp,lambda) != 0)
                {
                    Real scale = special::pow_int(r2/r1, lambda) / r1;
                    elem += Zp * ang_->f(ill,illp,lambda) * scale * rad_->Mtr_mLm1_full_x(lambda,i,k) * rad_->Mtr_L_full_y(lambda,j,l);
                }
                
                Cu_blocks_[ill * Nang + illp].add(row, col, Xpy[l2p][n][l] * elem);
            }
            
            // transition area r2 > r1, lower : F_nl expressed in terms of psi_kl for 'l' out of outer area
            for (int m = 0; m < Nchan2; m++)
            for (int n = 0; n < Nchan2p; n++)
            for (int j = Nspline_inner_y; j < Nspline_full_y; j++) // *
            for (int k = 0; k < Nspline_inner_x; k++)
            for (int l = j - order; l < Nspline_inner_y; l++)
            {
                std::size_t row = A_size + (Nchan1 + m) * Nspline_outer_y + (j - Nspline_inner_y);
                std::size_t col = k * Nspline_inner_y + l;
                
                Complex elem = 0; // B_mj,nl Sp_nk
                
                if (ill == illp and m == n)
                {
                    elem += (E_ + 0.5_r / ((n + l1 + 1) * (n + l1 + 1))) * rad_->S_full_y(j,l)
                         - 0.5_r * rad_->D_full_y(j,l)
                         - 0.5_r * (l2 * (l2 + 1.0_r)) * rad_->Mm2_full_y(j,l);
                }
                
                Real r2 = rad_->bspline_full_y().t(std::min(j,l) + order + 1).real();
                
                for (int lambda = 1; lambda <= rad_->maxlambda(); lambda++) if (ang_->f(ill,illp,lambda) != 0)
                {
                    Real scale = special::pow_int(1/r2, lambda + 1);
                    elem += Zp * ang_->f(ill,illp,lambda) * scale * rad_->Mtr_mLm1_full_y(lambda,j,l) * special::hydro_rho(m + l1 + 1, l1, n + l1p + 1, l1p, lambda);
                }
                
                Cl_blocks_[ill * Nang + illp].add(row, col, Spx[l1p][n][k] * elem);
            }
            
            // transition area r1 > r2, lower : F_nk expressed in terms of psi_kl for 'k' out of outer area
            for (int m = 0; m < Nchan1; m++)
            for (int n = 0; n < Nchan1p; n++)
            for (int i = Nspline_inner_x; i < Nspline_full_x; i++) // *
            for (int k = i - order; k < Nspline_inner_x; k++)
            for (int l = 0; l < Nspline_inner_y; l++)
            {
                std::size_t row = A_size + m * Nspline_outer_x + (i - Nspline_inner_x);
                std::size_t col = k * Nspline_inner_y + l;
                
                Complex elem = 0; // B_mi,nk Sp_nl
                
                if (ill == illp and m == n)
                {
                    elem += (E_ + 0.5_r / ((n + l2 + 1) * (n + l2 + 1))) * rad_->S_full_x(i,k)
                         - 0.5_r * rad_->D_full_x(i,k)
                         - 0.5_r * (l1 * (l1 + 1.0_r)) * rad_->Mm2_full_x(i,k);
                }
                
                Real r1 = rad_->bspline_full_x().t(std::min(i,k) + order + 1).real();
                
                for (int lambda = 1; lambda <= rad_->maxlambda(); lambda++) if (ang_->f(ill,illp,lambda) != 0)
                {
                    Real scale = special::pow_int(1/r1, lambda + 1);
                    elem += Zp * ang_->f(ill,illp,lambda) * scale * rad_->Mtr_mLm1_full_x(lambda,i,k) * special::hydro_rho(m + l2 + 1, l2, n + l2p + 1, l2p, lambda);
                }
                
                Cl_blocks_[ill * Nang + illp].add(row, col, Spy[l2p][n][l] * elem);
            }
        }
    }
    
    std::cout << "done after " << t.nice_time() << std::endl;
    par_->wait();
}

void NoPreconditioner::rhs (BlockArray<Complex> & chi, int ie, int instate) const
{
    // shorthands
    int ni = std::get<0>(inp_->instates[instate]);
    int li = std::get<1>(inp_->instates[instate]);
    int mi = std::get<2>(inp_->instates[instate]);
    
    // shorthands
    int order = inp_->order;
    std::size_t Nspline_inner_x = rad_->bspline_inner_x().Nspline();
    std::size_t Nspline_inner_y = rad_->bspline_inner_y().Nspline();
    std::size_t Nspline_full_x  = rad_->bspline_full_x ().Nspline();
    std::size_t Nspline_full_y  = rad_->bspline_full_y ().Nspline();
    std::size_t Nspline_outer_x = Nspline_full_x - Nspline_inner_x;
    std::size_t Nspline_outer_y = Nspline_full_y - Nspline_inner_y;
    
    // impact momentum
    rArray ki = { std::sqrt(inp_->Etot[ie] + 1.0_r/(ni*ni)) };
    
    // calculate LU-decomposition of the overlap matrix
    CsrMatrix<LU_int_t,Complex> Sx_csr_full = rad_->S_full_x().tocoo<LU_int_t>().tocsr();
    CsrMatrix<LU_int_t,Complex> Sy_csr_full = rad_->S_full_y().tocoo<LU_int_t>().tocsr();
    std::shared_ptr<LUft> lu_Sx_full;
    std::shared_ptr<LUft> lu_Sy_full;
    lu_Sx_full.reset(LUft::Choose("lapack"));
    lu_Sy_full.reset(LUft::Choose("lapack"));
    lu_Sx_full->factorize(Sx_csr_full);
    lu_Sy_full->factorize(Sy_csr_full);
    
    // j-overlaps of shape [Nangmom × Nspline]
    cArray ji_overlaps_full_x = rad_->overlapj(rad_->bspline_full_x(), rad_->gaussleg_full_x(), inp_->maxell, ki, cmd_->fast_bessel);
    cArray ji_overlaps_full_y = rad_->overlapj(rad_->bspline_full_y(), rad_->gaussleg_full_y(), inp_->maxell, ki, cmd_->fast_bessel);
    if (not std::isfinite(ji_overlaps_full_x.norm()) or not std::isfinite(ji_overlaps_full_y.norm()))
        HexException("Unable to compute Riccati-Bessel function B-spline overlaps!");
    
    // j-expansions
    cArray ji_expansion_full_x = lu_Sx_full->solve(ji_overlaps_full_x, inp_->maxell + 1);
    cArray ji_expansion_full_y = lu_Sy_full->solve(ji_overlaps_full_y, inp_->maxell + 1);
    if (not std::isfinite(ji_expansion_full_x.norm()) or not std::isfinite(ji_expansion_full_y.norm()))
        HexException("Unable to expand Riccati-Bessel function in B-splines!");
    
    // (anti)symmetrization
    Real Sign = ((ang_->S() + ang_->Pi()) % 2 == 0) ? 1. : -1.;
    
    // for all segments constituting the RHS
    for (unsigned ill = 0; ill < ang_->states().size(); ill++) if (par_->isMyGroupWork(ill))
    {
        int l1 = ang_->states()[ill].first;
        int l2 = ang_->states()[ill].second;
        
        // get number of open channels in the outer region
        int Nchan1 = Nchan_[ill].first;     // r1 -> inf, l2 bound
        int Nchan2 = Nchan_[ill].second;    // r2 -> inf, l1 bound
        
        // setup storage
        cArray chi_block (Nspline_inner_x * Nspline_inner_y + Nchan1 * Nspline_outer_x + Nchan2 * Nspline_outer_y);
        
        // for all allowed angular momenta (by momentum composition) of the projectile
        for (int l = std::abs(li - ang_->L()); l <= li + ang_->L(); l++)
        {
            // skip wrong parity
            if ((ang_->L() + li + l) % 2 != ang_->Pi())
                continue;
            
            // compute energy- and angular momentum-dependent prefactor
            Complex prefactor = std::pow(1.0_i,l)
                              * std::sqrt(4.0_r * special::constant::pi * (2 * l + 1))
                              * (Real)special::ClebschGordan(li,mi, l,0, inp_->L,mi) / ki[0];
            
            // skip non-contributing terms
            if (prefactor == 0.0_r)
                continue;
            
            // calculate angular integrals
            rArray f1 (rad_->maxlambda() + 1), f2 (rad_->maxlambda() + 1);
            for (int lambda = 0; lambda <= rad_->maxlambda(); lambda++)
            {
                f1[lambda] = special::computef(lambda, l1, l2, li, l, inp_->L);
                f2[lambda] = special::computef(lambda, l1, l2, l, li, inp_->L);
                
                // abort if any of the coefficients is non-number (factorial overflow etc.)
                if (not std::isfinite(f1[lambda]))
                    HexException("Invalid result of computef(%d,%d,%d,%d,%d,%d)\n", lambda,l1,l2,li,l,inp_->L);
                if (not std::isfinite(f2[lambda]))
                    HexException("Invalid result of computef(%d,%d,%d,%d,%d,%d)\n", lambda,l1,l2,l,li,inp_->L);
            }
            
            // quadrature degree
            int points = order + li + l + 1;
            
            // precompute quadrature nodes and weights
            cArray xs ((rad_->bspline_full_x().Nreknot() - 1) * points), xws ((rad_->bspline_full_x().Nreknot() - 1) * points);
            cArray ys ((rad_->bspline_full_y().Nreknot() - 1) * points), yws ((rad_->bspline_full_y().Nreknot() - 1) * points);
            # pragma omp parallel for
            for (int ixknot = 0; ixknot < rad_->bspline_full_x().Nreknot() - 1; ixknot++)
                rad_->gaussleg_full_x().scaled_nodes_and_weights(points, rad_->bspline_full_x().t(ixknot), rad_->bspline_full_x().t(ixknot + 1), &xs[ixknot * points], &xws[ixknot * points]);
            # pragma omp parallel for
            for (int iyknot = 0; iyknot < rad_->bspline_full_y().Nreknot() - 1; iyknot++)
                rad_->gaussleg_full_y().scaled_nodes_and_weights(points, rad_->bspline_full_y().t(iyknot), rad_->bspline_full_y().t(iyknot + 1), &ys[iyknot * points], &yws[iyknot * points]);
            
            // precompute B-splines
            cArray B_x (rad_->bspline_full_x().Nspline() * (order + 1) * points);
            cArray B_y (rad_->bspline_full_y().Nspline() * (order + 1) * points);
            # pragma omp parallel for
            for (int ixspline = 0; ixspline < rad_->bspline_full_x().Nspline(); ixspline++)
            for (int ixknot = ixspline; ixknot <= ixspline + order and ixknot < rad_->bspline_full_x().Nreknot() - 1; ixknot++)
                rad_->bspline_full_x().B(ixspline, ixknot, points, &xs[ixknot * points], &B_x[(ixspline * (order + 1) + ixknot - ixspline) * points]);
            # pragma omp parallel for
            for (int iyspline = 0; iyspline < rad_->bspline_full_y().Nspline(); iyspline++)
            for (int iyknot = iyspline; iyknot <= iyspline + order and iyknot < rad_->bspline_full_y().Nreknot() - 1; iyknot++)
                rad_->bspline_full_y().B(iyspline, iyknot, points, &ys[iyknot * points], &B_y[(iyspline * (order + 1) + iyknot - iyspline) * points]);
            
            // precompute radial functions and Riccati-Bessel functions
            rArray Pi_x ((rad_->bspline_full_x().Nreknot() - 1) * points), ji_x ((rad_->bspline_full_x().Nreknot() - 1) * points);
            rArray Pi_y ((rad_->bspline_full_y().Nreknot() - 1) * points), ji_y ((rad_->bspline_full_y().Nreknot() - 1) * points);
            # pragma omp parallel for
            for (unsigned ix = 0; ix < xs.size(); ix++)
            {
                gsl_sf_result Pi;
                Pi_x[ix] = (gsl_sf_hydrogenicR_e(ni, li, 1., xs[ix].real(), &Pi) == GSL_EUNDRFLW ? 0. : xs[ix].real() * Pi.val);
                ji_x[ix] = special::ric_j(l, ki[0] * xs[ix].real());
            }
            # pragma omp parallel for
            for (unsigned iy = 0; iy < ys.size(); iy++)
            {
                gsl_sf_result Pi;
                Pi_y[iy] = (gsl_sf_hydrogenicR_e(ni, li, 1., ys[iy].real(), &Pi) == GSL_EUNDRFLW ? 0. : ys[iy].real() * Pi.val);
                ji_y[iy] = special::ric_j(l, ki[0] * ys[iy].real());
            }
            
            // precompute integral moments
            cArrays M_L_Px (rad_->maxlambda() + 1), M_mLm1_Px (rad_->maxlambda() + 1), M_L_jx (rad_->maxlambda() + 1), M_mLm1_jx (rad_->maxlambda() + 1);
            cArrays M_L_Py (rad_->maxlambda() + 1), M_mLm1_Py (rad_->maxlambda() + 1), M_L_jy (rad_->maxlambda() + 1), M_mLm1_jy (rad_->maxlambda() + 1);
            for (int lambda = 0; lambda <= rad_->maxlambda(); lambda++)
            {
                M_L_Px[lambda].resize(Nspline_full_x); M_mLm1_Px[lambda].resize(Nspline_full_x);
                M_L_jx[lambda].resize(Nspline_full_x); M_mLm1_jx[lambda].resize(Nspline_full_x);
                
                # pragma omp parallel for
                for (int ixspline = 0; ixspline < (int)Nspline_full_x; ixspline++)
                {
                    for (int ixknot = ixspline; ixknot < std::min(ixspline + order + 1, rad_->bspline_full_x().Nreknot() - 1); ixknot++)
                    if (rad_->bspline_full_x().t(ixknot).real() != rad_->bspline_full_x().t(ixknot + 1).real())
                    {
                        for (int ixpoint = 0; ixpoint < points; ixpoint++)
                        {
                            Real x = xs[ixknot * points + ixpoint].real();
                            Real t = rad_->bspline_full_x().t(ixspline + order + 1).real();
                            
                            Real L = special::pow_int(x / t, lambda);
                            Real mLm1 = special::pow_int(x / t, -lambda-1);
                            
                            Complex w = xws[ixknot * points + ixpoint];
                            Complex B = B_x[(ixspline * (order + 1) + ixknot - ixspline) * points + ixpoint];
                            
                            Real P = Pi_x[ixknot * points + ixpoint];
                            Real j = ji_x[ixknot * points + ixpoint];
                            
                            Real d = 1 /*damp(x, 0, distance)*/;
                            
                            M_L_Px[lambda][ixspline] += w * B * L * P * d;
                            M_L_jx[lambda][ixspline] += w * B * L * j * d;
                            M_mLm1_Px[lambda][ixspline] += w * B * mLm1 * P * d;
                            M_mLm1_jx[lambda][ixspline] += w * B * mLm1 * j * d;
                        }
                    }
                }
                
                M_L_Py[lambda].resize(Nspline_full_y); M_mLm1_Py[lambda].resize(Nspline_full_y);
                M_L_jy[lambda].resize(Nspline_full_y); M_mLm1_jy[lambda].resize(Nspline_full_y);
                
                # pragma omp parallel for
                for (int iyspline = 0; iyspline < (int)Nspline_full_y; iyspline++)
                {
                    for (int iyknot = iyspline; iyknot < std::min(iyspline + order + 1, rad_->bspline_full_y().Nreknot() - 1); iyknot++)
                    if (rad_->bspline_full_y().t(iyknot).real() != rad_->bspline_full_y().t(iyknot + 1).real())
                    {
                        for (int iypoint = 0; iypoint < points; iypoint++)
                        {
                            Real y = ys[iyknot * points + iypoint].real();
                            Real t = rad_->bspline_full_y().t(iyspline + order + 1).real();
                            
                            Real L = special::pow_int(y / t, lambda);
                            Real mLm1 = special::pow_int(y / t, -lambda-1);
                            
                            Complex w = yws[iyknot * points + iypoint];
                            Complex B = B_y[(iyspline * (order + 1) + iyknot - iyspline) * points + iypoint];
                            
                            Real P = Pi_y[iyknot * points + iypoint];
                            Real j = ji_y[iyknot * points + iypoint];
                            
                            Real d = 1 /*damp(y, 0, distance)*/;
                            
                            M_L_Py[lambda][iyspline] += w * B * L * P * d;
                            M_L_jy[lambda][iyspline] += w * B * L * j * d;
                            M_mLm1_Py[lambda][iyspline] += w * B * mLm1 * P * d;
                            M_mLm1_jy[lambda][iyspline] += w * B * mLm1 * j * d;
                        }
                    }
                }
            }
            
            // precompute hydrogen multipoles
            cArrays rho_l1 (Nchan2), rho_l2 (Nchan1);
            for (int ichan1 = 0; ichan1 < Nchan1; ichan1++)
            {
                rho_l2[ichan1].resize(rad_->maxlambda() + 1);
                for (int lambda = 1; lambda <= rad_->maxlambda(); lambda++)
                    rho_l2[ichan1][lambda] = special::hydro_rho(ichan1 + l2 + 1, l2, ni, li, lambda);
            }
            for (int ichan2 = 0; ichan2 < Nchan2; ichan2++)
            {
                rho_l1[ichan2].resize(rad_->maxlambda() + 1);
                for (int lambda = 1; lambda <= rad_->maxlambda(); lambda++)
                    rho_l1[ichan2][lambda] = special::hydro_rho(ichan2 + l1 + 1, l1, ni, li, lambda);
            }
            
            // for all B-spline pairs (elements of the right-hand side)
            # pragma omp parallel for
            for (std::size_t ispline = 0; ispline < Nspline_inner_x * Nspline_inner_y + Nchan1 * Nspline_outer_x + Nchan2 * Nspline_outer_y; ispline++)
            {
                // contributions to the element of the right-hand side
                Complex contrib_direct = 0, contrib_exchange = 0;
                
                // get B-spline indices
                if (ispline >= Nspline_inner_x * Nspline_inner_y and ispline < Nspline_inner_x * Nspline_inner_y + Nchan1 * Nspline_outer_x)
                {
                    // r1 > Ra, r2 < Ra
                    int ixspline = (ispline - Nspline_inner_x * Nspline_inner_y) % Nspline_outer_x + Nspline_inner_x;
                    int ichan1 = (ispline - Nspline_inner_x * Nspline_inner_y) / Nspline_outer_x;
                    
                    // calculate the exchange contribution
                    for (int lambda = 1; lambda <= rad_->maxlambda(); lambda++) if (f2[lambda] != 0)
                        contrib_exchange += f2[lambda] * rho_l2[ichan1][lambda] * M_mLm1_jx[lambda][ixspline] / special::pow_int(rad_->bspline_full_x().t(ixspline + order + 1).real(), lambda + 1);
                }
                else if (ispline >= Nspline_inner_x * Nspline_inner_y + Nchan1 * Nspline_outer_x)
                {
                    // r1 < Ra, r2 > Ra
                    int iyspline = (ispline - Nspline_inner_x * Nspline_inner_y - Nchan1 * Nspline_outer_x) % Nspline_outer_y + Nspline_inner_y;
                    int ichan2 = (ispline - Nspline_inner_x * Nspline_inner_x - Nchan1 * Nspline_outer_x) / Nspline_outer_y;
                    
                    // calculate the direct contribution
                    for (int lambda = 1; lambda <= rad_->maxlambda(); lambda++) if (f1[lambda] != 0)
                        contrib_direct += f1[lambda] * rho_l1[ichan2][lambda] * M_mLm1_jy[lambda][iyspline] / special::pow_int(rad_->bspline_full_y().t(iyspline + order + 1).real(), lambda + 1);
                }
                else /* if (ispline < Nspline_inner * Nspline_inner) */
                {
                    // r1 < Ra, r2 < Ra
                    int ixspline = ispline / Nspline_inner_y;
                    int iyspline = ispline % Nspline_inner_y;
                    
                    // non-overlapping B-splines ?
                    if (std::abs(ixspline - iyspline) > order)
                    {
                        // monopole contribution
                        if (ixspline > iyspline and f1[0] != 0)
                        {
                            contrib_direct   += f1[0] * (M_mLm1_Px[0][ixspline] * M_L_jy[0][iyspline] / rad_->bspline_full_x().t(ixspline + order + 1) - M_L_Px[0][ixspline] * M_mLm1_jy[0][iyspline] / rad_->bspline_full_y().t(iyspline + order + 1));
                            contrib_exchange += 0;
                        }
                        if (ixspline < iyspline and f2[0] != 0)
                        {
                            contrib_direct   += 0;
                            contrib_exchange += f2[0] * (M_L_jx[0][ixspline] * M_mLm1_Py[0][iyspline] / rad_->bspline_full_y().t(iyspline + order + 1) - M_mLm1_jx[0][ixspline] * M_L_Py[0][iyspline] / rad_->bspline_full_x().t(ixspline + order + 1));
                        }
                        
                        // multipole contributions
                        for (int lambda = 1; lambda <= rad_->maxlambda(); lambda++)
                        {
                            if (ixspline > iyspline)
                            {
                                Real scale = special::pow_int(rad_->bspline_full_y().t(iyspline + order + 1).real() / rad_->bspline_full_x().t(ixspline + order + 1).real(), lambda) / rad_->bspline_full_x().t(ixspline + order + 1).real();
                                
                                contrib_direct   += f1[lambda] * M_mLm1_Px[lambda][ixspline] * M_L_jy[lambda][iyspline] * scale;
                                contrib_exchange += f2[lambda] * M_mLm1_jx[lambda][ixspline] * M_L_Py[lambda][iyspline] * scale;
                            }
                            if (ixspline < iyspline)
                            {
                                Real scale = special::pow_int(rad_->bspline_full_x().t(ixspline + order + 1).real() / rad_->bspline_full_y().t(iyspline + order + 1).real(), lambda) / rad_->bspline_full_y().t(iyspline + order + 1).real();
                                
                                contrib_direct   += f1[lambda] * M_L_Px[lambda][ixspline] * M_mLm1_jy[lambda][iyspline] * scale;
                                contrib_exchange += f2[lambda] * M_L_jx[lambda][ixspline] * M_mLm1_Py[lambda][iyspline] * scale;
                            }
                        }
                    }
                    else
                    {
                        // for all knots
                        for (int ixknot = ixspline; ixknot <= ixspline + order and ixknot < rad_->bspline_full_x().Nreknot() - 1; ixknot++) if (rad_->bspline_full_x().t(ixknot).real() != rad_->bspline_full_x().t(ixknot + 1).real())
                        for (int iyknot = iyspline; iyknot <= iyspline + order and iyknot < rad_->bspline_full_y().Nreknot() - 1; iyknot++) if (rad_->bspline_full_y().t(iyknot).real() != rad_->bspline_full_y().t(iyknot + 1).real())
                        {
                            // off-diagonal contribution
                            if (ixknot != iyknot)
                            {
                                // for all quadrature points
                                for (int ix = 0; ix < points; ix++)
                                for (int iy = 0; iy < points; iy++)
                                {
                                    // radii
                                    Real rx = xs[ixknot * points + ix].real(), ry = ys[iyknot * points + iy].real(), rmin = std::min(rx,ry), rmax = std::max(rx,ry);
                                    
                                    // evaluated functions
                                    Complex Bx = B_x[(ixspline * (order + 1) + ixknot - ixspline) * points + ix];
                                    Complex By = B_y[(iyspline * (order + 1) + iyknot - iyspline) * points + iy];
                                    Real Pix = Pi_x[ixknot * points + ix];
                                    Real Piy = Pi_y[iyknot * points + iy];
                                    Real jix = ji_x[ixknot * points + ix];
                                    Real jiy = ji_y[iyknot * points + iy];
                                    Complex wx = xws[ixknot * points + ix];
                                    Complex wy = yws[iyknot * points + iy];
                                    
                                    // damp factor
                                    Real dampfactor = 1 /* damp(rx, ry, distance) */;
                                    
                                    // monopole contribution
                                    if (rx > ry and li == l1 and l == l2) contrib_direct   += Bx * By * (1.0_r/rx - 1.0_r/ry) * Pix * jiy * dampfactor * wx * wy;
                                    if (ry > rx and li == l2 and l == l1) contrib_exchange += Bx * By * (1.0_r/ry - 1.0_r/rx) * jix * Piy * dampfactor * wx * wy;
                                    
                                    // higher multipoles contribution
                                    for (int lambda = 1; lambda <= rad_->maxlambda(); lambda++)
                                    {
                                        Real multipole = special::pow_int(rmin/rmax, lambda) / rmax;
                                        if (f1[lambda] != 0) contrib_direct   += f1[lambda] * Bx * By * multipole * Pix * jiy * dampfactor * wx * wy;
                                        if (f2[lambda] != 0) contrib_exchange += f2[lambda] * Bx * By * multipole * jix * Piy * dampfactor * wx * wy;
                                    }
                                }
                            }
                            // diagonal contribution: needs to be integrated more carefully
                            else if (ixknot < rad_->bspline_full_x().Nreknot() - 1)
                            {
                                // for all quadrature points from the triangle x < y
                                for (int ix = 0; ix < points; ix++)
                                {
                                    cArray ys (points), yws (points), B_y (points);
                                    rad_->gaussleg_full_y().scaled_nodes_and_weights(points, xs[ixknot * points + ix], rad_->bspline_full_y().t(iyknot + 1), &ys[0], &yws[0]);
                                    rad_->bspline_full_y().B(iyspline, iyknot, points, &ys[0], &B_y[0]);
                                    
                                    for (int iy = 0; iy < points; iy++)
                                    {
                                        // radii
                                        Real rx = xs[ixknot * points + ix].real(), ry = ys[iy].real(), rmin = std::min(rx,ry), rmax = std::max(rx,ry);
                                        
                                        // evaluated functions
                                        Complex Bx = B_x[(ixspline * (order + 1) + ixknot - ixspline) * points + ix];
                                        Complex By = B_y[iy];
                                        Real Pix = Pi_x[ixknot * points + ix];
                                        gsl_sf_result piy;
                                        Real Piy = (gsl_sf_hydrogenicR_e(ni, li, 1., ry, &piy) == GSL_EUNDRFLW ? 0. : ry * piy.val);
                                        Real jix = ji_x[ixknot * points + ix];
                                        Real jiy = special::ric_j(l, ki[0] * ry);
                                        Complex wx = xws[ixknot * points + ix];
                                        Complex wy = yws[iy];
                                        
                                        // damp factor
                                        Real dampfactor = 1/*damp(rx, ry, distance)*/;
                                        
                                        // monopole contribution
                                        if (rx > ry and li == l1 and l == l2) contrib_direct   += Bx * By * (1.0_r/rx - 1.0_r/ry) * Pix * jiy * dampfactor * wx * wy;
                                        if (ry > rx and li == l2 and l == l1) contrib_exchange += Bx * By * (1.0_r/ry - 1.0_r/rx) * jix * Piy * dampfactor * wx * wy;
                                        
                                        // higher multipoles contribution
                                        for (int lambda = 1; lambda <= rad_->maxlambda(); lambda++)
                                        {
                                            Real multipole = special::pow_int(rmin/rmax, lambda) / rmax;
                                            if (f1[lambda] != 0) contrib_direct   += f1[lambda] * Bx * By * multipole * Pix * jiy * dampfactor * wx * wy;
                                            if (f2[lambda] != 0) contrib_exchange += f2[lambda] * Bx * By * multipole * jix * Piy * dampfactor * wx * wy;
                                        }
                                    }
                                }
                                
                                // for all quadrature points from the triangle x > y
                                for (int ix = 0; ix < points; ix++)
                                {
                                    cArray ys (points), yws (points), B_y (points);
                                    rad_->gaussleg_full_y().scaled_nodes_and_weights(points, rad_->bspline_full_y().t(iyknot), xs[ixknot * points + ix], &ys[0], &yws[0]);
                                    rad_->bspline_full_y().B(iyspline, iyknot, points, &ys[0], &B_y[0]);
                                    
                                    for (int iy = 0; iy < points; iy++)
                                    {
                                        // radii
                                        Real rx = xs[ixknot * points + ix].real(), ry = ys[iy].real(), rmin = std::min(rx,ry), rmax = std::max(rx,ry);
                                        
                                        // evaluated functions
                                        Complex Bx = B_x[(ixspline * (order + 1) + ixknot - ixspline) * points + ix];
                                        Complex By = B_y[iy];
                                        Real Pix = Pi_x[ixknot * points + ix];
                                        gsl_sf_result piy;
                                        Real Piy = (gsl_sf_hydrogenicR_e(ni, li, 1., ry, &piy) == GSL_EUNDRFLW ? 0. : ry * piy.val);
                                        Real jix = ji_x[ixknot * points + ix];
                                        Real jiy = special::ric_j(l, ki[0] * ry);
                                        Complex wx = xws[ixknot * points + ix];
                                        Complex wy = yws[iy];
                                        
                                        // damp factor
                                        Real dampfactor = 1/*damp(rx, ry, distance)*/;
                                        
                                        // monopole contribution
                                        if (rx > ry and li == l1 and l == l2) contrib_direct   += Bx * By * (1.0_r/rx - 1.0_r/ry) * Pix * jiy * dampfactor * wx * wy;
                                        if (ry > rx and li == l2 and l == l1) contrib_exchange += Bx * By * (1.0_r/ry - 1.0_r/rx) * jix * Piy * dampfactor * wx * wy;
                                        
                                        // higher multipoles contribution
                                        for (int lambda = 1; lambda <= rad_->maxlambda(); lambda++)
                                        {
                                            Real multipole = special::pow_int(rmin/rmax, lambda) / rmax;
                                            if (f1[lambda] != 0) contrib_direct   += f1[lambda] * Bx * By * multipole * Pix * jiy * dampfactor * wx * wy;
                                            if (f2[lambda] != 0) contrib_exchange += f2[lambda] * Bx * By * multipole * jix * Piy * dampfactor * wx * wy;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                }
                
                if (std::isnan(contrib_direct.real()))
                {
                    std::cout << "Nan!" << std::endl;
                    std::cout << "ispline = " << ispline << std::endl;
                }
                
                // update element of the right-hand side
                if (inp_->Zp > 0)
                {
                    chi_block[ispline] += -prefactor * contrib_direct;
                }
                else
                {
                    chi_block[ispline] += prefactor * (contrib_direct + Sign * contrib_exchange) / special::constant::sqrt_two;
                }
            }
        }
        
        // use the calculated block
        chi[ill] = chi_block;
        
        // optionally transfer to disk
        if (not chi.inmemory())
        {
            chi.hdfsave(ill);
            chi[ill].drop();
        }
    }
}

void NoPreconditioner::multiply (BlockArray<Complex> const & p, BlockArray<Complex> & q, MatrixSelection::Selection tri) const
{
    // shorthands
    unsigned Nspline_inner_x = rad_->bspline_inner_x().Nspline();
    unsigned Nspline_inner_y = rad_->bspline_inner_y().Nspline();
    unsigned Nspline_full_x  = rad_->bspline_full_x().Nspline();
    unsigned Nspline_full_y  = rad_->bspline_full_y().Nspline();
    unsigned Nspline_outer_x = Nspline_full_x - Nspline_inner_x;
    unsigned Nspline_outer_y = Nspline_full_y - Nspline_inner_y;
    unsigned Nang = ang_->states().size();
    
    // make sure no process is playing with the data
    par_->wait();
    
    // TODO : It is slightly more subtle to do this efficiently in out-of-core mode, so we are just
    //        loading the destination vectors here. But would it possible to rewrite the code to load
    //        only the necessary pieces of those vectors on the fly.
    
    // de-const-ed source vector reference
    BlockArray<Complex> & v = const_cast<BlockArray<Complex> &>(p);
    
    // in distributed case we first need to collect the whole source vector
    for (unsigned ill = 0; ill < Nang; ill++)
    {
        // calculate rank of the process that owns this source vector segment (use group master process)
        int owner = (ill % par_->Ngroup()) * par_->groupsize();
        
        // load the source segment from disk, if necessary
        if (par_->iproc() == owner and cmd_->outofcore)
            v.hdfload(ill);
        
        // broadcast the source segment from the owner process to all others
        par_->bcast(owner, v[ill]);
        
        // also load the destination segment
        if (par_->isMyGroupWork(ill) and cmd_->outofcore)
            q.hdfload(ill);
    }
    
    // for all angular block rows
    for (unsigned ill = 0; ill < Nang; ill++) if (par_->isMyGroupWork(ill))
    {
        std::memset(q[ill].data(), 0, q[ill].size() * sizeof(Complex));
        
        // for all angular blocks in a block row; only executed by one of the processes in a process group
        for (unsigned illp = 0; illp < Nang; illp++) if (par_->igroupproc() == (int)illp % par_->groupsize())
        {
            // determine which part of the block should be considered non-zero for a particlar selection
            MatrixSelection::Selection selection = tri;
            if (ill < illp) selection = (tri & MatrixSelection::StrictUpper ? MatrixSelection::Both : MatrixSelection::None);
            if (ill > illp) selection = (tri & MatrixSelection::StrictLower ? MatrixSelection::Both : MatrixSelection::None);
            
            // near-origin part multiplication
            if (cmd_->lightweight_full)
            {
                // only one-electron contribution; the rest is below
                calc_A_block(ill, illp, false).dot
                (
                    1.0_z, cArrayView(p[illp], 0, Nspline_inner_x * Nspline_inner_y),
                    1.0_z, cArrayView(q[ill],  0, Nspline_inner_x * Nspline_inner_y),
                    true,
                    selection
                );
            }
            else
            {
                // read matrix from disk
                if (cmd_->outofcore and cmd_->wholematrix)
                    const_cast<BlockSymBandMatrix<Complex> &>(A_blocks_[ill * Nang + illp]).hdfload();
                
                // full diagonal block multiplication
                A_blocks_[ill * Nang + illp].dot
                (
                    1.0_z, cArrayView(p[illp], 0, Nspline_inner_x * Nspline_inner_y),
                    1.0_z, cArrayView(q[ill],  0, Nspline_inner_x * Nspline_inner_y),
                    true,
                    selection
                );
                
                // release memory
                if (cmd_->outofcore and cmd_->wholematrix)
                    const_cast<BlockSymBandMatrix<Complex> &>(A_blocks_[ill * Nang + illp]).drop();
            }
            
            // channel expansion part multiplication
            if (not inp_->inner_only)
            {
                int Nchan1 = Nchan_[ill].first;     // # r1 -> inf; l2 bound
                int Nchan2 = Nchan_[ill].second;    // # r2 -> inf; l1 bound
                int Nchan1p = Nchan_[illp].first;   // # r1 -> inf; l2p bound
                int Nchan2p = Nchan_[illp].second;  // # r2 -> inf; l1p bound
                
                // r1 -> inf
                # pragma omp parallel for
                for (int m = 0; m < Nchan1; m++)
                for (int n = 0; n < Nchan1p; n++)
                {
                    // read matrix from disk
                    if (cmd_->outofcore)
                        const_cast<SymBandMatrix<Complex>&>(B1_blocks_[ill * Nang + illp][m * Nchan1p + n]).hdfload();
                    
                    // multiply
                    B1_blocks_[ill * Nang + illp][m * Nchan1p + n].dot
                    (
                        1.0_z, cArrayView(p[illp], Nspline_inner_x * Nspline_inner_y + n * Nspline_outer_x, Nspline_outer_x),
                        1.0_z, cArrayView(q[ill],  Nspline_inner_x * Nspline_inner_y + m * Nspline_outer_x, Nspline_outer_x)
                    );
                    
                    // release memory
                    if (cmd_->outofcore)
                        const_cast<SymBandMatrix<Complex>&>(B1_blocks_[ill * Nang + illp][m * Nchan1p + n]).drop();
                }
                
                // r2 -> inf
                # pragma omp parallel for
                for (int m = 0; m < Nchan2; m++)
                for (int n = 0; n < Nchan2p; n++)
                {
                    // read matrix from disk
                    if (cmd_->outofcore)
                        const_cast<SymBandMatrix<Complex>&>(B2_blocks_[ill * Nang + illp][m * Nchan2p + n]).hdfload();
                    
                    // multiply
                    B2_blocks_[ill * Nang + illp][m * Nchan2p + n].dot
                    (
                        1.0_z, cArrayView(p[illp], Nspline_inner_x * Nspline_inner_y + Nchan1p * Nspline_outer_x + n * Nspline_outer_y, Nspline_outer_y),
                        1.0_z, cArrayView(q[ill],  Nspline_inner_x * Nspline_inner_y + Nchan1  * Nspline_outer_x + m * Nspline_outer_y, Nspline_outer_y)
                    );
                    
                    // release memory
                    if (cmd_->outofcore) const_cast<SymBandMatrix<Complex>&>(B2_blocks_[ill * Nang + illp][m * Nchan2p + n]).drop();
                }
                
                // multiply by coupling matrices
                Cu_blocks_[ill * Nang + illp].dot(1.0_z, p[illp], 1.0_z, q[ill]);
                Cl_blocks_[ill * Nang + illp].dot(1.0_z, p[illp], 1.0_z, q[ill]);
            }
        }
    }
    
    // lightweight-full off-diagonal contribution
    if (cmd_->lightweight_full)
    {
        OMP_CREATE_LOCKS(Nang * Nspline_inner_x);
        
        int maxlambda = rad_->maxlambda();
        
        # pragma omp parallel for collapse (3) schedule (dynamic,1)
        for (int lambda = 0; lambda <= maxlambda; lambda++)
        for (unsigned i = 0; i < Nspline_inner_x; i++)
        for (unsigned d = 0; d <= (unsigned)inp_->order; d++)
        if (i + d < Nspline_inner_x)
        {
            unsigned k = i + d;
            
            SymBandMatrix<Complex> R = std::move(rad_->calc_R_tr_dia_block(lambda, i, k));
            
            for (unsigned ill = 0; ill < Nang; ill++) if (par_->isMyGroupWork(ill))
            for (unsigned illp = 0; illp < Nang; illp++) if (par_->igroupproc() == (int)illp % par_->groupsize())
            if (Real f = ang_->f(ill, illp, lambda))
            {
                // diagonal blocks
                if (d == 0)
                {
                    OMP_LOCK_LOCK(ill * Nspline_inner_x + i);
                    
                    R.dot
                    (
                        inp_->Zp * f, cArrayView(p[illp], i * Nspline_inner_y, Nspline_inner_y),
                        1.0_z, cArrayView(q[ill], i * Nspline_inner_y, Nspline_inner_y),
                        ill == illp ? tri : (ill < illp ? (tri & MatrixSelection::StrictUpper ? MatrixSelection::Both : MatrixSelection::None) :
                                                          (tri & MatrixSelection::StrictLower ? MatrixSelection::Both : MatrixSelection::None))
                    );
                    
                    OMP_UNLOCK_LOCK(ill * Nspline_inner_x + i);
                }
                
                // off-diagonal blocks
                else
                {
                    if (tri & MatrixSelection::StrictUpper)
                    {
                        OMP_LOCK_LOCK(ill * Nspline_inner_x + i);
                        
                        R.dot
                        (
                            inp_->Zp * f, cArrayView(p[illp], k * Nspline_inner_y, Nspline_inner_y),
                            1.0_z, cArrayView(q[ill], i * Nspline_inner_y, Nspline_inner_y)
                        );
                        
                        OMP_UNLOCK_LOCK(ill * Nspline_inner_x + i);
                    }
                    
                    if (tri & MatrixSelection::StrictLower)
                    {
                        OMP_LOCK_LOCK(ill * Nspline_inner_x + k);
                        
                        R.dot
                        (
                            inp_->Zp * f, cArrayView(p[illp], i * Nspline_inner_y, Nspline_inner_y),
                            1.0_z, cArrayView(q[ill], k * Nspline_inner_y, Nspline_inner_y)
                        );
                        
                        OMP_UNLOCK_LOCK(ill * Nspline_inner_x + k);
                    }
                }
            }
        }
        
        OMP_DELETE_LOCKS();
    }
    
    // release source vectors
    for (unsigned ill = 0; ill < Nang; ill++) if (cmd_->outofcore)
        v[ill].drop();
    
    // synchronize and release the result vectors
    for (unsigned ill = 0; ill < Nang; ill++) if (par_->isMyGroupWork(ill))
    {
        // synchronize (sum) across the group
        par_->syncsum_g(q[ill].data(), q[ill].size());
        
        // constrain
        if (const CGPreconditioner * cgprec = dynamic_cast<const CGPreconditioner*>(this))
            cgprec->CG_constrain(q[ill]);
        
        // release memory
        if (cmd_->outofcore)
        {
            if (par_->IamGroupMaster())
                q.hdfsave(ill);
            
            q[ill].drop();
        }
    }
}

void NoPreconditioner::precondition (const BlockArray< Complex >& r, BlockArray< Complex >& z) const
{
    z = r;
}

void NoPreconditioner::finish ()
{
    A_blocks_ .resize(0);
    B1_blocks_.resize(0);
    B2_blocks_.resize(0);
    Cu_blocks_.resize(0);
    Cl_blocks_.resize(0);
}

// --------------------------------------------------------------------------------- //

addClassToParentRunTimeSelectionTable(PreconditionerBase, NoPreconditioner)

// --------------------------------------------------------------------------------- //
