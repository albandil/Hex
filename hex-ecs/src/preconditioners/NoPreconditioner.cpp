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

#include "hex-arrays.h"
#include "hex-misc.h"

#include "gauss.h"
#include "parallel.h"
#include "preconditioners.h"
#include "radial.h"

#include "NoPreconditioner.h"

const std::string NoPreconditioner::prec_name = "none";
const std::string NoPreconditioner::prec_description = "\"Preconditioning\" by the identity matrix.";

void NoPreconditioner::setup ()
{
    rad_.setupOneElectronIntegrals(par_, cmd_);
    rad_.setupTwoElectronIntegrals(par_, cmd_);
}

void NoPreconditioner::update (Real E)
{
    // shorthands
    int order = inp_.order;
    int Nang = ang_.states().size();
    int Nspline_inner = rad_.bspline_inner().Nspline();
    int Nspline_full  = rad_.bspline_full().Nspline();
    int Nspline_outer = Nspline_full - Nspline_inner;
    std::size_t A_size = std::size_t(Nspline_inner) * std::size_t(Nspline_inner);
    
    // update energy
    E_ = E;
    
    // get maximal asymptotic principal quantum number
    max_n_ = (E_ >= 0 ? 0 : 1.0 / std::sqrt(-2 * E_));
    
    // update number of asymptotic channels
    Nchan_.clear();
    for (int ill = 0; ill < Nang; ill++)
    {
        int l1 = ang_.states()[ill].first;
        int l2 = ang_.states()[ill].second;
        
        // number of channels when r1 -> inf (i.e. second electron is bound)
        int Nchan1 = std::max(0, max_n_ - l2);
        
        // number of channels when r2 -> inf (i.e. first electron is bound)
        int Nchan2 = std::max(0, max_n_ - l1);
        
        Nchan_.push_back(std::make_pair(Nchan1, Nchan2));
    }
    
    // skip pre-calculation of the diagonal blocks in full lightweight mode
    if (cmd_.lightweight_full)
        return;
    
    std::cout << "\tPrecompute matrix blocks ... " << std::flush;
    Timer t;
    
    // LU-factorize the overlap matrix
    CsrMatrix<LU_int_t,Complex> csr_S = rad_.S_full().tocoo<LU_int_t>().tocsr();
    std::shared_ptr<LUft<LU_int_t,Complex>> lu_S = csr_S.factorize();
    
    // outer one-electron overlap matrix
    SymBandMatrix<Complex> S_outer (Nspline_outer, order + 1);
    S_outer.populate([&](int m, int n) { return rad_.S_full()(Nspline_inner + m, Nspline_inner + n); });
    
    // outer one-electron derivative matrix
    SymBandMatrix<Complex> D_outer (Nspline_outer, order + 1);
    D_outer.populate([&](int m, int n) { return rad_.D_full()(Nspline_inner + m, Nspline_inner + n); });
    
    // outer one-electron centrifugal moment matrix
    SymBandMatrix<Complex> Mm2_outer (Nspline_outer, order + 1);
    Mm2_outer.populate([&](int m, int n) { return rad_.Mm2_full()(Nspline_inner + m, Nspline_inner + n); });
    
    // outer one-electron multipole moment matrices
    std::vector<SymBandMatrix<Complex>> Mtr_mLm1_outer;
    for (int lambda = 0; lambda <= rad_.maxlambda(); lambda++)
    {
        Mtr_mLm1_outer.push_back(SymBandMatrix<Complex>(Nspline_outer, order + 1));
        Mtr_mLm1_outer.back().populate
        (
            [ & ] (int m, int n)
            {
                return rad_.Mtr_mLm1_full(lambda)(Nspline_inner + m, Nspline_inner + n)
                     * special::pow_int(rad_.bspline_full().t(Nspline_inner + std::min(m,n) + order + 1).real(), -lambda-1);
            }
        );
    }
    
    // setup blocks
    for (int ill = 0; ill < Nang; ill++)
    for (int illp = 0; illp < Nang; illp++)
    {
        // angular momenta
        int l1 = ang_.states()[ill].first;
        int l2 = ang_.states()[ill].second;
        int l1p = ang_.states()[illp].first;
        int l2p = ang_.states()[illp].second;
        
        // initialize diagonal block of the inner problem
        A_blocks_[ill * Nang + illp] = BlockSymBandMatrix<Complex>
        (
            Nspline_inner,              // block count
            inp_.order + 1,             // block structure half-bandwidth
            Nspline_inner,              // block size
            inp_.order + 1,             // block half-bandwidth
            !cmd_.outofcore,            // keep in memory?
            format("blk-A-%d-%d.ooc", ill, illp)  // scratch disk file name
        );
        
        // skip calculation if the disk file is already present
        /*if (cmd_.outofcore and cmd_.reuse_dia_blocks and dia_blocks_[ill].hdfcheck())
            continue;
        else if (cmd_.outofcore)
            dia_blocks_[ill].hdfinit();*/
        
        // for all sub-blocks
        # pragma omp parallel for if (A_blocks_[ill * Nang + illp].inmemory())
        for (int i = 0; i < Nspline_inner; i++)
        for (int d = 0; d <= order; d++)
        if (i + d < Nspline_inner)
        {
            int k = i + d;
            
            SymBandMatrix<Complex> subblock (Nspline_inner, order + 1);
            
            // one-electron part
            if (ill == illp)
            {
                subblock.populate
                (
                    [&](int j, int l)
                    {
                        return E * rad_.S_full()(i,k) * rad_.S_full()(j,l)
                               - (0.5_z * rad_.D_full()(i,k) - rad_.Mm1_tr_full()(i,k)) * rad_.S_full()(j,l)
                               - 0.5_r * l1 * (l1 + 1) * rad_.Mm2_full()(i,k) * rad_.S_full()(j,l)
                               - rad_.S_full()(i,k) * (0.5_z * rad_.D_full()(j,l) - rad_.Mm1_tr_full()(j,l))
                               - 0.5_r * l2 * (l2 + 1) * rad_.S_full()(i,k) * rad_.Mm2_full()(j,l);
                    }
                );
            }
            
            // two-electron part
            for (int lambda = 0; lambda <= rad_.maxlambda(); lambda++)
            {
                // check that the "f" coefficient is nonzero
                if (ang_.f(ill,illp,lambda) == 0.)
                    continue;
                
                // calculate two-electron term
                if (not cmd_.lightweight_radial_cache)
                {
                    // use precomputed block from scratch file or from memory
                    subblock.data() += -ang_.f(ill,illp,lambda) * rad_.R_tr_dia(lambda).getBlock(i * (order + 1) + d).slice(0, Nspline_inner * (order + 1));
                }
                else
                {
                    // compute the data anew
                    subblock.data() += -ang_.f(ill,illp,lambda) * rad_.calc_R_tr_dia_block(lambda, i, k).data().slice(0, Nspline_inner * (order + 1));
                }
            }
            
            // save block
            A_blocks_[ill * Nang + illp].setBlock(i * (order + 1) + d, subblock.data());
        }
        
        // get number of asymptotic bound channels
        int Nchan1 = Nchan_[ill].first;     // # r1 -> inf, l2 bound
        int Nchan2 = Nchan_[ill].second;    // # r2 -> inf, l1 bound
        int Nchan1p = Nchan_[illp].first;   // # r1 -> inf, l2p bound
        int Nchan2p = Nchan_[illp].second;  // # r2 -> inf, l1p bound
        
        // create inner-outer coupling blocks
        Cu_blocks_[ill * Nang + illp] = CooMatrix<LU_int_t,Complex>
        (
            A_size + (Nchan1 + Nchan2) * Nspline_outer,
            A_size + (Nchan1p + Nchan2p) * Nspline_outer
        );
        Cl_blocks_[ill * Nang + illp] = CooMatrix<LU_int_t,Complex>
        (
            A_size + (Nchan1 + Nchan2) * Nspline_outer,
            A_size + (Nchan1p + Nchan2p) * Nspline_outer
        );
        
        // calculate all missing hydrogen overlaps
        while (Sp.size() <= (unsigned)mmax(l1,l2,l1p,l2p)) Sp.push_back(cArrays());
        while (Xp.size() <= (unsigned)mmax(l1,l2,l1p,l2p)) Xp.push_back(cArrays());
        while (Sp[l1] .size() < (unsigned)Nchan2 ) Sp[l1] .push_back(cArray());
        while (Sp[l2] .size() < (unsigned)Nchan1 ) Sp[l2] .push_back(cArray());
        while (Sp[l1p].size() < (unsigned)Nchan2p) Sp[l1p].push_back(cArray());
        while (Sp[l2p].size() < (unsigned)Nchan1p) Sp[l2p].push_back(cArray());
        while (Xp[l1] .size() < (unsigned)Nchan2 ) Xp[l1] .push_back(cArray());
        while (Xp[l2] .size() < (unsigned)Nchan1 ) Xp[l2] .push_back(cArray());
        while (Xp[l1p].size() < (unsigned)Nchan2p) Xp[l1p].push_back(cArray());
        while (Xp[l2p].size() < (unsigned)Nchan1p) Xp[l2p].push_back(cArray());
        for (int n = 0; n < Nchan2; n++) if (Sp[l1][n].empty())
        {
            Sp[l1][n] = rad_.overlapP(rad_.bspline_full(), rad_.gaussleg_full(), l1 + n + 1, l1, weightEndDamp(rad_.bspline_full()));
            Xp[l1][n] = lu_S->solve(Sp[l1][n]);
        }
        for (int n = 0; n < Nchan1; n++) if (Sp[l2][n].empty())
        {
            Sp[l2][n] = rad_.overlapP(rad_.bspline_full(), rad_.gaussleg_full(), l2 + n + 1, l2, weightEndDamp(rad_.bspline_full()));
            Xp[l2][n] = lu_S->solve(Sp[l2][n]);
        }
        for (int n = 0; n < Nchan2p; n++) if (Sp[l1p][n].empty())
        {
            Sp[l1p][n] = rad_.overlapP(rad_.bspline_full(), rad_.gaussleg_full(), l1p + n + 1, l1p, weightEndDamp(rad_.bspline_full()));
            Xp[l1p][n] = lu_S->solve(Sp[l1p][n]);
        }
        for (int n = 0; n < Nchan1p; n++) if (Sp[l2p][n].empty())
        {
            Sp[l2p][n] = rad_.overlapP(rad_.bspline_full(), rad_.gaussleg_full(), l2p + n + 1, l2p, weightEndDamp(rad_.bspline_full()));
            Xp[l2p][n] = lu_S->solve(Sp[l2p][n]);
        }
        
        // setup stretched inner-outer problem
        if (not inp_.inner_only)
        {
            // outer problem matrix : r2 -> inf, l1 bound
            B2_blocks_[ill * Nang + illp].resize(Nchan2 * Nchan2p);
            for (int m = 0; m < Nchan2; m++)
            for (int n = 0; n < Nchan2p; n++)
            {
                SymBandMatrix<Complex> subblock (Nspline_outer, order + 1);
                
                // channel-diagonal contribution
                if (ill == illp and m == n)
                {
                    subblock += (E_ + 1.0_z / (2.0_z * (l1 + m + 1.0_r) * (l1 + m + 1.0_r))) * S_outer
                             - 0.5_z * D_outer
                             - 0.5_z * (l2 * (l2 + 1.0_r)) * Mm2_outer;
                }
                
                // channel-offdiagonal contribution
                for (int lambda = 1; lambda <= rad_.maxlambda(); lambda++) if (ang_.f(ill,illp,lambda) != 0.0_r)
                    subblock -= Complex(ang_.f(ill,illp,lambda) * special::hydro_rho(l1 + m + 1, l1, l1p + n + 1, l1p, lambda)) * Mtr_mLm1_outer[lambda];
                
                // use the block
                B2_blocks_[ill * Nang + illp][m * Nchan2p + n] = std::move(subblock);
            }
            
            // outer problem matrix : r1 -> inf, l2 bound
            B1_blocks_[ill * Nang + illp].resize(Nchan1 * Nchan1p);
            for (int m = 0; m < Nchan1; m++)
            for (int n = 0; n < Nchan1p; n++)
            {
                SymBandMatrix<Complex> subblock (Nspline_outer, order + 1);
                
                // channel-diagonal contribution
                if (ill == illp and m == n)
                {
                    subblock += (E_ + 1.0_z / (2.0_z * (l2 + m + 1.0_r) * (l2 + m + 1.0_r))) * S_outer
                             - 0.5_z * D_outer
                             - 0.5_z * (l1 * (l1 + 1.0_r)) * Mm2_outer;
                }
                
                // channel-offdiagonal contribution
                for (int lambda = 1; lambda <= rad_.maxlambda(); lambda++) if (ang_.f(ill,illp,lambda) != 0.0_r)
                    subblock -= Complex(ang_.f(ill,illp,lambda) * special::hydro_rho(l2 + m + 1, l2, l2p + n + 1, l2p, lambda)) * Mtr_mLm1_outer[lambda];
                
                // use the block
                B1_blocks_[ill * Nang + illp][m * Nchan1p + n] = std::move(subblock);
            }
            
            // transition area r2 > r1, upper : psi_kl expressed in terms of F_nl for 'l' out of inner area
            for (int i = 0; i < Nspline_inner; i++)
            for (int j = 0; j < Nspline_inner; j++) // *
            for (int k = std::max(0, i - order); k <= std::min(i + order, Nspline_inner - 1); k++)
            for (int l = Nspline_inner; l <= j + order; l++)
            for (int n = 0; n < Nchan2p; n++)
            {
                std::size_t row = i * Nspline_inner + j;
                std::size_t col = A_size + (Nchan1p + n) * Nspline_outer + (l - Nspline_inner);
                
                Complex elem = 0; // A_ij,kl Xp_nk
                
                if (ill == illp)
                {
                    elem += E_ * rad_.S_full()(i,k) * rad_.S_full()(j,l)
                         - 0.5_r * rad_.D_full()(i,k) * rad_.S_full()(j,l)
                         - 0.5_r * rad_.S_full()(i,k) * rad_.D_full()(j,l)
                         - 0.5_r * (l1 * (l1 + 1.0_r)) * rad_.Mm2_full()(i,k) * rad_.S_full()(j,l)
                         - 0.5_r * (l2 * (l2 + 1.0_r)) * rad_.S_full()(i,k) * rad_.Mm2_full()(j,l)
                         + rad_.Mm1_tr_full()(i,k) * rad_.S_full()(j,l)
                         + rad_.S_full()(i,k) * rad_.Mm1_tr_full()(j,l);
                }
                
                Real r1 = rad_.bspline_full().t(std::min(i,k) + order + 1).real();
                Real r2 = rad_.bspline_full().t(std::min(j,l) + order + 1).real();
                
                for (int lambda = 0; lambda <= rad_.maxlambda(); lambda++)
                {
                    Real scale = special::pow_int(r1/r2, lambda) / r2;
                    elem -= ang_.f(ill,illp,lambda) * scale * rad_.Mtr_L_full(lambda)(i,k) * rad_.Mtr_mLm1_full(lambda)(j,l);
                }
                
                Cu_blocks_[ill * Nang + illp].add(row, col, Xp[l1p][n][k] * elem);
            }
            
            // transition area r1 > r2, upper : psi_kl expressed in terms of F_nk for 'k' out of inner area
            for (int i = 0; i < Nspline_inner; i++)
            for (int j = 0; j < Nspline_inner; j++) // *
            for (int k = Nspline_inner; k <= i + order; k++)
            for (int l = std::max(0, j - order); l <= std::min(j + order, Nspline_inner - 1); l++)
            for (int n = 0; n < Nchan1p; n++)
            {
                std::size_t row = i * Nspline_inner + j;
                std::size_t col = A_size + n * Nspline_outer + (k - Nspline_inner);
                
                Complex elem = 0; // A_ij,kl Xp_nl
                
                if (ill == illp)
                {
                    elem += E_ * rad_.S_full()(i,k) * rad_.S_full()(j,l)
                         - 0.5_r * rad_.D_full()(i,k) * rad_.S_full()(j,l)
                         - 0.5_r * rad_.S_full()(i,k) * rad_.D_full()(j,l)
                         - 0.5_r * (l1 * (l1 + 1.0_r)) * rad_.Mm2_full()(i,k) * rad_.S_full()(j,l)
                         - 0.5_r * (l2 * (l2 + 1.0_r)) * rad_.S_full()(i,k) * rad_.Mm2_full()(j,l)
                         + rad_.Mm1_tr_full()(i,k) * rad_.S_full()(j,l)
                         + rad_.S_full()(i,k) * rad_.Mm1_tr_full()(j,l);
                }
                
                Real r1 = rad_.bspline_full().t(std::min(i,k) + order + 1).real();
                Real r2 = rad_.bspline_full().t(std::min(j,l) + order + 1).real();
                
                for (int lambda = 0; lambda <= rad_.maxlambda(); lambda++)
                {
                    Real scale = special::pow_int(r2/r1, lambda) / r1;
                    elem -= ang_.f(ill,illp,lambda) * scale * rad_.Mtr_mLm1_full(lambda)(i,k) * rad_.Mtr_L_full(lambda)(j,l);
                }
                
                Cu_blocks_[ill * Nang + illp].add(row, col, Xp[l2p][n][l] * elem);
            }
            
            // transition area r2 > r1, lower : F_nl expressed in terms of psi_kl for 'l' out of outer area
            for (int m = 0; m < Nchan2; m++)
            for (int n = 0; n < Nchan2p; n++)
            for (int j = Nspline_inner; j < Nspline_full; j++) // *
            for (int k = 0; k < Nspline_inner; k++)
            for (int l = j - order; l < Nspline_inner; l++)
            {
                std::size_t row = A_size + (Nchan1 + m) * Nspline_outer + (j - Nspline_inner);
                std::size_t col = k * Nspline_inner + l;
                
                Complex elem = 0; // B_mj,nl Sp_nk
                
                if (ill == illp and m == n)
                {
                    elem += (E_ + 0.5_r / ((n + l1 + 1) * (n + l1 + 1))) * rad_.S_full()(j,l)
                         - 0.5_r * rad_.D_full()(j,l)
                         - 0.5_r * (l2 * (l2 + 1.0_r)) * rad_.Mm2_full()(j,l);
                }
                
                Real r2 = rad_.bspline_full().t(std::min(j,l) + order + 1).real();
                
                for (int lambda = 1; lambda <= rad_.maxlambda(); lambda++)
                {
                    Real scale = special::pow_int(1/r2, lambda + 1);
                    elem -= ang_.f(ill,illp,lambda) * scale * rad_.Mtr_mLm1_full(lambda)(j,l) * special::hydro_rho(m + l1 + 1, l1, n + l1p + 1, l1p, lambda);
                }
                
                Cl_blocks_[ill * Nang + illp].add(row, col, Sp[l1p][n][k] * elem);
            }
            
            // transition area r1 > r2, lower : F_nk expressed in terms of psi_kl for 'k' out of outer area
            for (int m = 0; m < Nchan1; m++)
            for (int n = 0; n < Nchan1p; n++)
            for (int i = Nspline_inner; i < Nspline_full; i++) // *
            for (int k = i - order; k < Nspline_inner; k++)
            for (int l = 0; l < Nspline_inner; l++)
            {
                std::size_t row = A_size + m * Nspline_outer + (i - Nspline_inner);
                std::size_t col = k * Nspline_inner + l;
                
                Complex elem = 0; // B_mi,nk Sp_nl
                
                if (ill == illp and m == n)
                {
                    elem += (E_ + 0.5_r / ((n + l2 + 1) * (n + l2 + 1))) * rad_.S_full()(i,k)
                         - 0.5_r * rad_.D_full()(i,k)
                         - 0.5_r * (l1 * (l1 + 1.0_r)) * rad_.Mm2_full()(i,k);
                }
                
                Real r1 = rad_.bspline_full().t(std::min(i,k) + order + 1).real();
                
                for (int lambda = 1; lambda <= rad_.maxlambda(); lambda++)
                {
                    Real scale = special::pow_int(1/r1, lambda + 1);
                    elem -= ang_.f(ill,illp,lambda) * scale * rad_.Mtr_mLm1_full(lambda)(i,k) * special::hydro_rho(m + l2 + 1, l2, n + l2p + 1, l2p, lambda);
                }
                
                Cl_blocks_[ill * Nang + illp].add(row, col, Sp[l2p][n][l] * elem);
            }
        }
    }
    
    std::cout << "done after " << t.nice_time() << std::endl;
    par_.wait();
}

void NoPreconditioner::rhs (BlockArray<Complex> & chi, int ie, int instate) const
{
    // shorthands
    int ni = std::get<0>(inp_.instates[instate]);
    int li = std::get<1>(inp_.instates[instate]);
    int mi = std::get<2>(inp_.instates[instate]);
    
    // shorthands
    int order = rad_.bspline_inner().order();
    std::size_t Nspline_inner = rad_.bspline_inner().Nspline();
    std::size_t Nspline_full  = rad_.bspline_full ().Nspline();
    std::size_t Nspline_outer = Nspline_full - Nspline_inner;
    
    // impact momentum
    rArray ki = { std::sqrt(inp_.Etot[ie] + 1.0_r/(ni*ni)) };
    
    // calculate LU-decomposition of the overlap matrix
    CsrMatrix<LU_int_t,Complex> S_csr_full = rad_.S_full().tocoo<LU_int_t>().tocsr();
    std::shared_ptr<LUft<LU_int_t,Complex>> lu_S_full = S_csr_full.factorize();
    
    // j-overlaps of shape [Nangmom × Nspline]
    cArray ji_overlaps_full = rad_.overlapj(rad_.bspline_full(), rad_.gaussleg_full(), inp_.maxell, ki, weightEdgeDamp(rad_.bspline_full()));
    if (not std::isfinite(ji_overlaps_full.norm()))
        HexException("Unable to compute Riccati-Bessel function B-spline overlaps!");
    
    // j-expansions
    cArray ji_expansion_full = lu_S_full->solve(ji_overlaps_full, inp_.maxell + 1);
    if (not std::isfinite(ji_expansion_full.norm()))
        HexException("Unable to expand Riccati-Bessel function in B-splines!");
    
    // (anti)symmetrization
    Real Sign = ((ang_.S() + ang_.Pi()) % 2 == 0) ? 1. : -1.;
    
    // for all segments constituting the RHS
    for (unsigned ill = 0; ill < ang_.states().size(); ill++) if (par_.isMyGroupWork(ill))
    {
        int l1 = ang_.states()[ill].first;
        int l2 = ang_.states()[ill].second;
        
        // get number of open channels in the outer region
        int Nchan1 = Nchan_[ill].first;     // r1 -> inf, l2 bound
        int Nchan2 = Nchan_[ill].second;    // r2 -> inf, l1 bound
        
        // setup storage
        cArray chi_block (Nspline_inner * Nspline_inner + (Nchan1 + Nchan2) * Nspline_outer);
        
        // for all allowed angular momenta (by momentum composition) of the projectile
        for (int l = std::abs(li - ang_.L()); l <= li + ang_.L(); l++)
        {
            // skip wrong parity
            if ((ang_.L() + li + l) % 2 != ang_.Pi())
                continue;
            
            // compute energy- and angular momentum-dependent prefactor
            Complex prefactor = std::pow(1.0_i,l)
                              * std::sqrt(special::constant::two_pi * (2 * l + 1))
                              * (Real)special::ClebschGordan(li,mi, l,0, inp_.L,mi) / ki[0];
            
            // skip non-contributing terms
            if (prefactor == 0.0_r)
                continue;
            
            // calculate angular integrals
            rArray f1 (rad_.maxlambda() + 1), f2 (rad_.maxlambda() + 1);
            for (int lambda = 0; lambda <= rad_.maxlambda(); lambda++)
            {
                f1[lambda] = special::computef(lambda, l1, l2, li, l, inp_.L);
                f2[lambda] = special::computef(lambda, l1, l2, l, li, inp_.L);
                
                // abort if any of the coefficients is non-number (factorial overflow etc.)
                if (not std::isfinite(f1[lambda]))
                    HexException("Invalid result of computef(%d,%d,%d,%d,%d,%d)\n", lambda,l1,l2,li,l,inp_.L);
                if (not std::isfinite(f2[lambda]))
                    HexException("Invalid result of computef(%d,%d,%d,%d,%d,%d)\n", lambda,l1,l2,l,li,inp_.L);
            }
            
            // calculate the right-hand side
            if (cmd_.exact_rhs)
            {
                // quadrature degree
                int points = order + li + l + 1;
                
                // precompute quadrature nodes and weights
                cArray xs ((rad_.bspline_full().Nreknot() - 1) * points), xws ((rad_.bspline_full().Nreknot() - 1) * points);
                # pragma omp parallel for
                for (int ixknot = 0; ixknot < rad_.bspline_full().Nreknot() - 1; ixknot++)
                    rad_.gaussleg_full().scaled_nodes_and_weights(points, rad_.bspline_full().t(ixknot), rad_.bspline_full().t(ixknot + 1), &xs[ixknot * points], &xws[ixknot * points]);
                
                // precompute B-splines
                cArray B_x (rad_.bspline_full().Nspline() * (order + 1) * points);
                # pragma omp parallel for
                for (int ixspline = 0; ixspline < rad_.bspline_full().Nspline(); ixspline++)
                for (int ixknot = ixspline; ixknot <= ixspline + order and ixknot < rad_.bspline_full().Nreknot() - 1; ixknot++)
                    rad_.bspline_full().B(ixspline, ixknot, points, &xs[ixknot * points], &B_x[(ixspline * (order + 1) + ixknot - ixspline) * points]);
                
                // precompute radial functions and Riccati-Bessel functions
                rArray Pi_x ((rad_.bspline_full().Nreknot() - 1) * points);
                rArray ji_x ((rad_.bspline_full().Nreknot() - 1) * points);
                # pragma omp parallel for
                for (unsigned ix = 0; ix < xs.size(); ix++)
                {
                    gsl_sf_result Pi;
                    Pi_x[ix] = (gsl_sf_hydrogenicR_e(ni, li, 1., xs[ix].real(), &Pi) == GSL_EUNDRFLW ? 0. : xs[ix].real() * Pi.val);
                    ji_x[ix] = special::ric_j(l, ki[0] * xs[ix].real());
                }
                
                // damping distance
                Real distance = rad_.bspline_full().R0();
                
                // precompute integral moments
                cArrays M_L_P (rad_.maxlambda() + 1), M_mLm1_P (rad_.maxlambda() + 1);
                cArrays M_L_j (rad_.maxlambda() + 1), M_mLm1_j (rad_.maxlambda() + 1);
                for (int lambda = 0; lambda <= rad_.maxlambda(); lambda++)
                {
                    M_L_P[lambda].resize(Nspline_full); M_mLm1_P[lambda].resize(Nspline_full);
                    M_L_j[lambda].resize(Nspline_full); M_mLm1_j[lambda].resize(Nspline_full);
                    
                    # pragma omp parallel for
                    for (int ispline = 0; ispline < (int)Nspline_full; ispline++)
                    {
                        for (int iknot = ispline; iknot < std::min(ispline + order + 1, rad_.bspline_full().Nreknot() - 1); iknot++) if (rad_.bspline_full().t(iknot).real() != rad_.bspline_full().t(iknot + 1).real())
                        {
                            for (int ipoint = 0; ipoint < points; ipoint++)
                            {
                                Real x = xs[iknot * points + ipoint].real();
                                Real t = rad_.bspline_full().t(ispline + order + 1).real();
                                
                                Real L = special::pow_int(x / t, lambda);
                                Real mLm1 = special::pow_int(x / t, -lambda-1);
                                
                                Complex w = xws[iknot * points + ipoint];
                                Complex B = B_x[(ispline * (order + 1) + iknot - ispline) * points + ipoint];
                                
                                Real P = Pi_x[iknot * points + ipoint];
                                Real j = ji_x[iknot * points + ipoint];
                                
                                Real d = damp(x, 0, distance);
                                
                                M_L_P[lambda][ispline] += w * B * L * P * d;
                                M_L_j[lambda][ispline] += w * B * L * j * d;
                                M_mLm1_P[lambda][ispline] += w * B * mLm1 * P * d;
                                M_mLm1_j[lambda][ispline] += w * B * mLm1 * j * d;
                            }
                        }
                    }
                }
                
                // precompute hydrogen multipoles
                cArrays rho_l1 (Nchan2), rho_l2 (Nchan1);
                for (int ichan1 = 0; ichan1 < Nchan1; ichan1++)
                {
                    rho_l2[ichan1].resize(rad_.maxlambda() + 1);
                    for (int lambda = 1; lambda <= rad_.maxlambda(); lambda++)
                        rho_l2[ichan1][lambda] = special::hydro_rho(ichan1 + l2 + 1, l2, ni, li, lambda);
                }
                for (int ichan2 = 0; ichan2 < Nchan2; ichan2++)
                {
                    rho_l1[ichan2].resize(rad_.maxlambda() + 1);
                    for (int lambda = 1; lambda <= rad_.maxlambda(); lambda++)
                        rho_l1[ichan2][lambda] = special::hydro_rho(ichan2 + l1 + 1, l1, ni, li, lambda);
                }
                
                // for all B-spline pairs (elements of the right-hand side)
                # pragma omp parallel for schedule(dynamic,Nspline_inner)
                for (std::size_t ispline = 0; ispline < Nspline_inner * Nspline_inner + (Nchan1 + Nchan2) * Nspline_outer; ispline++)
                {
                    // contributions to the element of the right-hand side
                    Complex contrib_direct = 0, contrib_exchange = 0;
                    
                    // get B-spline indices
                    if (ispline >= Nspline_inner * Nspline_inner and ispline < Nspline_inner * Nspline_inner + Nchan1 * Nspline_outer)
                    {
                        // r1 > Ra, r2 < Ra
                        int ixspline = (ispline - Nspline_inner * Nspline_inner) % Nspline_outer + Nspline_inner;
                        int ichan1 = (ispline - Nspline_inner * Nspline_inner) / Nspline_outer;
                        
                        // calculate the exchange contribution
                        for (int lambda = 1; lambda <= rad_.maxlambda(); lambda++)
                            contrib_exchange += f2[lambda] * rho_l2[ichan1][lambda] * M_mLm1_j[lambda][ixspline] / special::pow_int(rad_.bspline_full().t(ixspline + order + 1).real(), lambda + 1);
                    }
                    else if (ispline >= Nspline_inner * Nspline_inner + Nchan1 * Nspline_outer)
                    {
                        // r1 < Ra, r2 > Ra
                        int iyspline = (ispline - Nspline_inner * Nspline_inner - Nchan1 * Nspline_outer) % Nspline_outer + Nspline_inner;
                        int ichan2 = (ispline - Nspline_inner * Nspline_inner - Nchan1 * Nspline_outer) / Nspline_outer;
                        
                        // calculate the direct contribution
                        for (int lambda = 1; lambda <= rad_.maxlambda(); lambda++)
                            contrib_direct += f1[lambda] * rho_l1[ichan2][lambda] * M_mLm1_j[lambda][iyspline] / special::pow_int(rad_.bspline_full().t(iyspline + order + 1).real(), lambda + 1);
                    }
                    else /* if (ispline < Nspline_inner * Nspline_inner) */
                    {
                        // r1 < Ra, r2 < Ra
                        int ixspline = ispline / Nspline_inner;
                        int iyspline = ispline % Nspline_inner;
                        
                        // non-overlapping B-splines ?
                        if (std::abs(ixspline - iyspline) > order)
                        {
                            // monopole contribution
                            if (ixspline > iyspline and f1[0] != 0)
                            {
                                contrib_direct   += f1[0] * (M_mLm1_P[0][ixspline] * M_L_j[0][iyspline] / rad_.bspline_full().t(ixspline + order + 1) - M_L_P[0][ixspline] * M_mLm1_j[0][iyspline] / rad_.bspline_full().t(iyspline + order + 1));
                                contrib_exchange += 0;
                            }
                            if (ixspline < iyspline and f2[0] != 0)
                            {
                                contrib_direct   += 0;
                                contrib_exchange += f2[0] * (M_L_j[0][ixspline] * M_mLm1_P[0][iyspline] / rad_.bspline_full().t(iyspline + order + 1) - M_mLm1_j[0][ixspline] * M_L_P[0][iyspline] / rad_.bspline_full().t(ixspline + order + 1));
                            }
                            
                            // multipole contributions
                            for (int lambda = 1; lambda <= rad_.maxlambda(); lambda++)
                            {
                                if (ixspline > iyspline)
                                {
                                    Real scale = special::pow_int(rad_.bspline_full().t(iyspline + order + 1).real() / rad_.bspline_full().t(ixspline + order + 1).real(), lambda) / rad_.bspline_full().t(ixspline + order + 1).real();
                                    
                                    contrib_direct   += f1[lambda] * M_mLm1_P[lambda][ixspline] * M_L_j[lambda][iyspline] * scale;
                                    contrib_exchange += f2[lambda] * M_mLm1_j[lambda][ixspline] * M_L_P[lambda][iyspline] * scale;
                                }
                                if (ixspline < iyspline)
                                {
                                    Real scale = special::pow_int(rad_.bspline_full().t(ixspline + order + 1).real() / rad_.bspline_full().t(iyspline + order + 1).real(), lambda) / rad_.bspline_full().t(iyspline + order + 1).real();
                                    
                                    contrib_direct   += f1[lambda] * M_L_P[lambda][ixspline] * M_mLm1_j[lambda][iyspline] * scale;
                                    contrib_exchange += f2[lambda] * M_L_j[lambda][ixspline] * M_mLm1_P[lambda][iyspline] * scale;
                                }
                            }
                        }
                        else
                        {
                            // for all knots
                            for (int ixknot = ixspline; ixknot <= ixspline + order and ixknot < rad_.bspline_full().Nreknot() - 1; ixknot++) if (rad_.bspline_full().t(ixknot).real() != rad_.bspline_full().t(ixknot + 1).real())
                            for (int iyknot = iyspline; iyknot <= iyspline + order and iyknot < rad_.bspline_full().Nreknot() - 1; iyknot++) if (rad_.bspline_full().t(iyknot).real() != rad_.bspline_full().t(iyknot + 1).real())
                            {
                                // off-diagonal contribution
                                if (ixknot != iyknot)
                                {
                                    // for all quadrature points
                                    for (int ix = 0; ix < points; ix++)
                                    for (int iy = 0; iy < points; iy++)
                                    {
                                        // radii
                                        Real rx = xs[ixknot * points + ix].real(), ry = xs[iyknot * points + iy].real(), rmin = std::min(rx,ry), rmax = std::max(rx,ry);
                                        
                                        // evaluated functions
                                        Complex Bx = B_x[(ixspline * (order + 1) + ixknot - ixspline) * points + ix];
                                        Complex By = B_x[(iyspline * (order + 1) + iyknot - iyspline) * points + iy];
                                        Real Pix = Pi_x[ixknot * points + ix];
                                        Real Piy = Pi_x[iyknot * points + iy];
                                        Real jix = ji_x[ixknot * points + ix];
                                        Real jiy = ji_x[iyknot * points + iy];
                                        Complex wx = xws[ixknot * points + ix];
                                        Complex wy = xws[iyknot * points + iy];
                                        
                                        // damp factor
                                        Real dampfactor = damp(rx, ry, distance);
                                        
                                        // monopole contribution
                                        if (rx > ry and li == l1 and l == l2) contrib_direct   += Bx * By * (1.0_r/rx - 1.0_r/ry) * Pix * jiy * dampfactor * wx * wy;
                                        if (ry > rx and li == l2 and l == l1) contrib_exchange += Bx * By * (1.0_r/ry - 1.0_r/rx) * jix * Piy * dampfactor * wx * wy;
                                        
                                        // higher multipoles contribution
                                        for (int lambda = 1; lambda <= rad_.maxlambda(); lambda++)
                                        {
                                            Real multipole = special::pow_int(rmin/rmax, lambda) / rmax;
                                            if (f1[lambda] != 0) contrib_direct   += f1[lambda] * Bx * By * multipole * Pix * jiy * dampfactor * wx * wy;
                                            if (f2[lambda] != 0) contrib_exchange += f2[lambda] * Bx * By * multipole * jix * Piy * dampfactor * wx * wy;
                                        }
                                    }
                                }
                                // diagonal contribution: needs to be integrated more carefully
                                else if (ixknot < rad_.bspline_full().Nreknot() - 1)
                                {
                                    // for all quadrature points from the triangle x < y
                                    for (int ix = 0; ix < points; ix++)
                                    {
                                        cArray ys (points), yws (points), B_y (points);
                                        rad_.gaussleg_full().scaled_nodes_and_weights(points, xs[ixknot * points + ix], rad_.bspline_full().t(iyknot + 1), &ys[0], &yws[0]);
                                        rad_.bspline_full().B(iyspline, iyknot, points, &ys[0], &B_y[0]);
                                        
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
                                            Real dampfactor = damp(rx, ry, distance);
                                            
                                            // monopole contribution
                                            if (rx > ry and li == l1 and l == l2) contrib_direct   += Bx * By * (1.0_r/rx - 1.0_r/ry) * Pix * jiy * dampfactor * wx * wy;
                                            if (ry > rx and li == l2 and l == l1) contrib_exchange += Bx * By * (1.0_r/ry - 1.0_r/rx) * jix * Piy * dampfactor * wx * wy;
                                            
                                            // higher multipoles contribution
                                            for (int lambda = 1; lambda <= rad_.maxlambda(); lambda++)
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
                                        rad_.gaussleg_full().scaled_nodes_and_weights(points, rad_.bspline_full().t(iyknot), xs[ixknot * points + ix], &ys[0], &yws[0]);
                                        rad_.bspline_full().B(iyspline, iyknot, points, &ys[0], &B_y[0]);
                                        
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
                                            Real dampfactor = damp(rx, ry, distance);
                                            
                                            // monopole contribution
                                            if (rx > ry and li == l1 and l == l2) contrib_direct   += Bx * By * (1.0_r/rx - 1.0_r/ry) * Pix * jiy * dampfactor * wx * wy;
                                            if (ry > rx and li == l2 and l == l1) contrib_exchange += Bx * By * (1.0_r/ry - 1.0_r/rx) * jix * Piy * dampfactor * wx * wy;
                                            
                                            // higher multipoles contribution
                                            for (int lambda = 1; lambda <= rad_.maxlambda(); lambda++)
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
                    
                    // update element of the right-hand side
                    chi_block[ispline] += prefactor * (contrib_direct + Sign * contrib_exchange);
                }
                
                chi[ill] = std::move(chi_block);
            }
            else
            {
                HexException("Please use --exact-rhs.");
            }
        }
        
        if (not chi.inmemory())
        {
            chi.hdfsave(ill);
            chi[ill].drop();
        }
    }
}

void NoPreconditioner::multiply (BlockArray<Complex> const & p, BlockArray<Complex> & q) const
{
    // shorthands
    unsigned Nspline_inner = rad_.bspline_inner().Nspline();
    unsigned Nspline_full  = rad_.bspline_full().Nspline();
    unsigned Nspline_outer = Nspline_full - Nspline_inner;
    unsigned Nang = ang_.states().size();
    
    for (unsigned ill = 0; ill < Nang; ill++)
    {
        std::memset(q[ill].data(), 0, q[ill].size() * sizeof(Complex));
        
        for (unsigned illp = 0; illp < Nang; illp++)
        {
            A_blocks_[ill * Nang + illp].dot
            (
                1.0_z, cArrayView(p[illp], 0, Nspline_inner * Nspline_inner),
                1.0_z, cArrayView(q[ill], 0, Nspline_inner * Nspline_inner)
            );
            
            if (not inp_.inner_only)
            {
                int Nchan1 = Nchan_[ill].first;     // # r1 -> inf; l2 bound
                int Nchan2 = Nchan_[ill].second;    // # r2 -> inf; l1 bound
                int Nchan1p = Nchan_[illp].first;   // # r1 -> inf; l2p bound
                int Nchan2p = Nchan_[illp].second;  // # r2 -> inf; l1p bound
                
                // r1 -> inf
                for (int m = 0; m < Nchan1; m++)
                for (int n = 0; n < Nchan1p; n++)
                {
                    B1_blocks_[ill * Nang + illp][m * Nchan1p + n].dot
                    (
                        1.0_z, cArrayView(p[illp], Nspline_inner * Nspline_inner + n * Nspline_outer, Nspline_outer),
                        1.0_z, cArrayView(q[ill], Nspline_inner * Nspline_inner + m * Nspline_outer, Nspline_outer)
                    );
                }
                
                // r2 -> inf
                for (int m = 0; m < Nchan2; m++)
                for (int n = 0; n < Nchan2p; n++)
                {
                    B2_blocks_[ill * Nang + illp][m * Nchan2p + n].dot
                    (
                        1.0_z, cArrayView(p[illp], Nspline_inner * Nspline_inner + (Nchan1p + n) * Nspline_outer, Nspline_outer),
                        1.0_z, cArrayView(q[ill], Nspline_inner * Nspline_inner + (Nchan1 + m) * Nspline_outer, Nspline_outer)
                    );
                }
                
                Cu_blocks_[ill * Nang + illp].dot(1.0_z, p[illp], 1.0_z, q[ill]);
                Cl_blocks_[ill * Nang + illp].dot(1.0_z, p[illp], 1.0_z, q[ill]);
            }
        }
    }
    
    // constrain the result
    if (const CGPreconditioner * cgprec = dynamic_cast<const CGPreconditioner*>(this))
    {
        for (std::size_t ill = 0; ill < Nang; ill++)
            cgprec->CG_constrain(q[ill]);
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
