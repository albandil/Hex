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

#include <iostream>

#include "../arrays.h"
#include "../gauss.h"
#include "../misc.h"
#include "../parallel.h"
#include "../preconditioners.h"
#include "../radial.h"

const std::string NoPreconditioner::name = "none";
const std::string NoPreconditioner::description = "\"Preconditioning\" by the identity matrix.";

void NoPreconditioner::setup ()
{
    // TODO : Determine which lambdas are needed by this process.
    // NOTE : At the moment each process holds in memory radial integrals
    //        for all lambdas, which needlessly raises memory requirements.
    Array<bool> lambdas (inp_.L + 2 * inp_.levels + 1, true);
    
    // compute one-electron radial integrals
    s_rad_.setupOneElectronIntegrals(cmd_);
    
    // compute two-eletron radial integrals
    s_rad_.setupTwoElectronIntegrals(par_, cmd_, lambdas);
}

void NoPreconditioner::update (double E)
{
    // update energy
    E_ = E;
    
    // skip pre-calculation of the diagonal blocks in full lightweight mode
    if (cmd_.lightweight_full)
        return;
    
    std::cout << "\tPrecompute diagonal blocks... " << std::flush;
    
    // some accelerator data for the calculation of the radial integrals
    cArrays Mtr_L(s_rad_.maxlambda() + 1), Mtr_mLm1(s_rad_.maxlambda() + 1);
    if (cmd_.lightweight_radial_cache)
    {
        for (int lambda = 0; lambda <= s_rad_.maxlambda(); lambda++)
            s_rad_.init_R_tr_dia_block(lambda, Mtr_L[lambda], Mtr_mLm1[lambda]);
    }
    
    // setup diagonal blocks
    # pragma omp parallel for if (cmd_.parallel_block)
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyWork(ill))
    {
        // angular momenta
        int l1 = l1_l2_[ill].first;
        int l2 = l1_l2_[ill].second;
        
        // initialize diagonal block
        dia_blocks_[ill] = BlockSymDiaMatrix
        (
            s_bspline_.Nspline(),       // block count (and size)
            s_rad_.S().nzpattern(),     // block structure
            s_rad_.S().diag(),          // non-zero diagonals
            not cmd_.outofcore,         // keep in memory?
            format("dblk-%d.ooc", ill)  // scratch disk file name
        );
        
        // skip calculation if the disk file is already present
        if (cmd_.outofcore and cmd_.reuse_dia_blocks and dia_blocks_[ill].hdfcheck())
            continue;
        
        // for all blocks
        for (unsigned iblock = 0; iblock < dia_blocks_[ill].structure().size(); iblock++)
        {
            // get block position
            int i = dia_blocks_[ill].structure()[iblock].first;
            int j = dia_blocks_[ill].structure()[iblock].second;
            
            // one-electron part
            cArray block = E * s_rad_.S()(i,j) * s_rad_.S().data();
            block -= (0.5 * s_rad_.D()(i,j) - s_rad_.Mm1_tr()(i,j)) * s_rad_.S().data();
            block -= 0.5 * l1 * (l1 + 1) * s_rad_.Mm2()(i,j) * s_rad_.S().data();
            block -= s_rad_.S()(i,j) * (0.5 * s_rad_.D().data() - s_rad_.Mm1_tr().data());
            block -= 0.5 * l2 * (l2 + 1) * s_rad_.S()(i,j) * s_rad_.Mm2().data();
            
            // two-electron part
            for (int lambda = 0; lambda <= s_rad_.maxlambda(); lambda++)
            {
                // calculate angular integral
                Complex f = special::computef(lambda,l1,l2,l1,l2,inp_.L);
                
                // check that the "f" coefficient is valid (no factorial overflow etc.)
                if (not Complex_finite(f))
                    Exception("Overflow in computation of f[%d](%d,%d,%d,%d).", inp_.L, l1, l2, l1, l2);
                
                // check that the "f" coefficient is nonzero
                if (f == 0.)
                    continue;
            
                // add two-electron term
                if (not cmd_.lightweight_radial_cache)
                {
                    // use precomputed block (from memory or from scratch file)
                    block += (-f) * s_rad_.R_tr_dia(lambda).getBlock(iblock);
                }
                else
                {
                    // compute the data anew
                    block += (-f) * s_rad_.calc_R_tr_dia_block(lambda, i, j, Mtr_L[lambda], Mtr_mLm1[lambda], true);
                }
            }
            
            // save block
            dia_blocks_[ill].setBlock(iblock,block);
        }
    }
    
    par_.wait();
    std::cout << "ok" << std::endl;
}

void NoPreconditioner::rhs (cArray & chi, int ie, int instate, int Spin) const
{
    // shorthands
    int li = std::get<1>(inp_.instates[instate]);
    int mi = std::get<2>(inp_.instates[instate]);
    
    // shorthands
    int Nspline = s_rad_.bspline().Nspline();
    
    // j-overlaps of shape [Nangmom × Nspline]
    cArray ji_overlaps = s_rad_.overlapj
    (
        inp_.maxell,
        inp_.ki.slice(ie, ie + 1), // use only one ki
        weightEdgeDamp(s_rad_.bspline())
    );
    ji_overlaps.hdfsave("ji_overlaps.hdf");
    if (not std::isfinite(ji_overlaps.norm()))
        Exception("Unable to compute Riccati-Bessel function B-spline overlaps!");
    
    // j-expansions
    cArray ji_expansion = s_rad_.S().tocoo().tocsr().solve(ji_overlaps, ji_overlaps.size() / Nspline);
    if (not std::isfinite(ji_expansion.norm()))
        Exception("Unable to expand Riccati-Bessel function in B-splines!");
    
    // compute P-overlaps and P-expansion
    cArray Pi_overlaps, Pi_expansion;
    Pi_overlaps = s_rad_.overlapP(inp_.ni, li, weightEndDamp(s_rad_.bspline()));
    Pi_expansion = s_rad_.S().tocoo().tocsr().solve(Pi_overlaps);
    if (not std::isfinite(Pi_expansion.norm()))
        Exception("Unable to expand hydrogen bound orbital in B-splines!");
    
    // for all segments constituting the RHS
    # pragma omp parallel for schedule (dynamic,1) if (cmd_.parallel_block)
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyWork(ill))
    {
        int l1 = l1_l2_[ill].first;
        int l2 = l1_l2_[ill].second;
        
        // setup storage
        cArray chi_block (Nspline * Nspline);
        
        // for all allowed angular momenta (by momentum composition) of the projectile
        for (int l = std::abs(li - inp_.L); l <= li + inp_.L; l++)
        {
            // skip wrong parity
            if ((inp_.L + li + l) % 2 != inp_.Pi)
                continue;
            
            // (anti)symmetrization
            double Sign = ((Spin + inp_.Pi) % 2 == 0) ? 1. : -1.;
            
            // compute energy- and angular momentum-dependent prefactor
            Complex prefactor = std::pow(Complex(0.,1.),l)
                              * std::sqrt(special::constant::two_pi * (2 * l + 1))
                              * special::ClebschGordan(li,mi, l,0, inp_.L,mi) / inp_.ki[ie];
            
            // skip non-contributing terms
            if (prefactor == 0.)
                continue;
            
            // pick the correct Bessel function expansion
            cArrayView Ji_expansion (ji_expansion, l * Nspline, Nspline);
            
            // compute outer products of B-spline expansions
            cArray Pj1 = outer_product(Pi_expansion, Ji_expansion);
            cArray Pj2 = outer_product(Ji_expansion, Pi_expansion);
            
            // skip angular forbidden right hand sides
            for (int lambda = 0; lambda <= s_rad_.maxlambda(); lambda++)
            {
                // calculate angular integrals
                double f1 = special::computef(lambda, l1, l2, li, l, inp_.L);
                double f2 = special::computef(lambda, l1, l2, l, li, inp_.L);
                
                // abort if any of the coefficients is non-number (factorial overflow etc.)
                if (not std::isfinite(f1))
                    Exception("Invalid result of computef(%d,%d,%d,%d,%d,%d)\n", lambda,l1,l2,li,l,inp_.L);
                if (not std::isfinite(f2))
                    Exception("Invalid result of computef(%d,%d,%d,%d,%d,%d)\n", lambda,l1,l2,l,li,inp_.L);
                
                // add multipole terms (direct/exchange)
                if (not cmd_.lightweight_radial_cache)
                {
                    if (f1 != 0.) chi_block += (       prefactor * f1) * s_rad_.R_tr_dia(lambda).dot(Pj1, cmd_.parallel_dot);
                    if (f2 != 0.) chi_block += (Sign * prefactor * f2) * s_rad_.R_tr_dia(lambda).dot(Pj2, cmd_.parallel_dot);
                }
                else
                {
                    if (f1 != 0.) chi_block += (       prefactor * f1) * s_rad_.apply_R_matrix(lambda, cArrays(1, Pj1))[0];
                    if (f2 != 0.) chi_block += (Sign * prefactor * f2) * s_rad_.apply_R_matrix(lambda, cArrays(1, Pj2))[0];
                }
            }
            
            // add monopole terms (direct/exchange)
            if (li == l1 and l == l2)
                chi_block += (-prefactor       ) * outer_product(s_rad_.S().dot(Pi_expansion), s_rad_.Mm1_tr().dot(Ji_expansion));
            if (li == l2 and l == l1)
                chi_block += (-prefactor * Sign) * outer_product(s_rad_.Mm1_tr().dot(Ji_expansion), s_rad_.S().dot(Pi_expansion));
            
            // update the right-hand side
            # pragma omp critical
            {
                // resize if necessary
                if (chi.size() < (ill / par_.Nproc() + 1) * Nspline * Nspline)
                    chi.resize((ill / par_.Nproc() + 1) * Nspline * Nspline);
                
                // copy data
                cArrayView(chi, ill / par_.Nproc() * Nspline * Nspline, Nspline * Nspline) = chi_block;
            }
        }
    }
}

void NoPreconditioner::multiply (const cArrayView p, cArrayView q) const
{
    // shorthands
    int Nspline = s_rad_.bspline().Nspline();
    int Nang = l1_l2_.size();
    int Nchunk = Nspline * Nspline;
    
    // clear output array
    std::memset(q.data(), 0, q.size() * sizeof(Complex));

    if (not cmd_.lightweight_radial_cache)
    {
        //
        // Simple multiplication by the super-matrix.
        //
        
        // multiply "p" by the diagonal blocks
        # pragma omp parallel for schedule (dynamic,1) if (cmd_.parallel_block)
        for (int ill = 0;  ill < Nang;  ill++)
        if (par_.isMyWork(ill))
        {
            // copy-from segment of "p"
            cArrayView p_block (p, (ill / par_.Nproc()) * Nchunk, Nchunk);
            
            // copy-to segment of "q"
            cArrayView q_block (q, (ill / par_.Nproc()) * Nchunk, Nchunk);
            
            // multiply
            q_block = dia_blocks_[ill].dot(p_block, cmd_.parallel_dot);
        }
        
        // multiply "p" by the off-diagonal blocks
        for (int ill = 0;  ill < Nang;  ill++)
        {
            // product of line of blocks with the source vector (-> one segment of destination vector)
            cArray product (Nchunk);
            
            # pragma omp parallel for schedule (dynamic,1) collapse (2) if (cmd_.parallel_block)
            for (int illp = 0; illp < Nang; illp++)
            for (int lambda = 0; lambda <= s_rad_.maxlambda(); lambda++)
            {
                // skip diagonal
                if (ill == illp)
                    continue;
                
                // skip segments not owned by this process (will be computed by other processes)
                if (not par_.isMyWork(illp))
                    continue;
                
                // row multi-index
                int l1 = l1_l2_[ill].first;
                int l2 = l1_l2_[ill].second;
                
                // column multi-index
                int l1p = l1_l2_[illp].first;
                int l2p = l1_l2_[illp].second;
                
                // calculate angular integral
                double f = special::computef(lambda, l1, l2, l1p, l2p, inp_.L);
                if (not std::isfinite(f))
                    Exception("Invalid result of computef(%d,%d,%d,%d,%d,%d).", lambda, l1, l2, l1p, l2p, inp_.L);
                
                // check non-zero
                if (f == 0.)
                    continue;
                
                // source segment of "p"
                cArrayView p_block (p, (illp / par_.Nproc()) * Nchunk, Nchunk);
                
                // calculate product
                cArray p0 = std::move( (-f) * s_rad_.R_tr_dia(lambda).dot(p_block, cmd_.parallel_dot) );
                
                // update collected product
                # pragma omp critical
                product += p0;
            }
            
            // sum all contributions to this destination segment on its owner node
            par_.sum(product.data(), Nchunk, ill % par_.Nproc());
            
            // finally, owner will update its segment
            if (par_.isMyWork(ill))
                cArrayView (q, (ill / par_.Nproc()) * Nchunk, Nchunk) += product;
        }
    }
    else
    {
        //
        // Lightweight mode multiplication
        //
        
        // get one-electron matrix structure
        std::vector<std::pair<int,int>> structure = s_rad_.S().nzpattern();
        
        // multiply "p" by the diagonal super-blocks
        # pragma omp parallel for schedule (dynamic,1) if (cmd_.parallel_block)
        for (int ill = 0;  ill < Nang;  ill++)
        if (par_.isMyWork(ill))
        {
            // copy-from segment of "p"
            cArrayView p_block (p, (ill / par_.Nproc()) * Nchunk, Nchunk);
            
            // copy-to segment of "q"
            cArrayView q_block (q, (ill / par_.Nproc()) * Nchunk, Nchunk);
            
            // get block angular momemnta
            int l1 = l1_l2_[ill].first;
            int l2 = l1_l2_[ill].second;
            
            // multiply 'p_block' by the diagonal block (except for the two-electron term)
            q_block  = kron_dot(Complex(E_) * s_rad_.S_d(), s_rad_.S_d(), p_block);
            q_block -= kron_dot(Complex(0.5) * s_rad_.D_d() - s_rad_.Mm1_tr_d() + Complex(0.5 * l1 * (l1 + 1.)) * s_rad_.Mm2_d(), s_rad_.S_d(), p_block);
            q_block -= kron_dot(s_rad_.S_d(), Complex(0.5) * s_rad_.D_d() - s_rad_.Mm1_tr_d() + Complex(0.5 * l2 * (l2 + 1.)) * s_rad_.Mm2_d(), p_block);
        }
        
        // auxiliary buffers
        cArray buffer (Nspline * Nspline);
        cArrays Mtr_L(s_rad_.maxlambda() + 1), Mtr_mLm1(s_rad_.maxlambda() + 1);
        for (int lambda = 0; lambda <= s_rad_.maxlambda(); lambda++)
            s_rad_.init_R_tr_dia_block(lambda, Mtr_L[lambda], Mtr_mLm1[lambda]);
        
        // multiply by offdiagonal blocks
        for (int illp = 0; illp < Nang; illp++)
        {
            // owner of this block will copy it to the buffer
            if (par_.isMyWork(illp))
                std::memcpy(&buffer[0], p.data() + (illp / par_.Nproc()) * Nspline * Nspline, Nspline * Nspline * sizeof(Complex));
            
            // owner will broadcast the source segment 'p(illp)' to all non-owners
            par_.bcast(illp % par_.Nproc(), buffer);
            
            // everyone will update its result segments; for all angular momentum transfers
            for (int lambda = 0; lambda <= s_rad_.maxlambda(); lambda++)
            {
                // precalculate all angular factors 'f(lambda,ill,illp)'
                rArray fs (Nang);
                for (int ill = 0; ill < Nang; ill++)
                {
                    fs[ill] = special::computef(lambda, l1_l2_[ill].first, l1_l2_[ill].second, l1_l2_[illp].first, l1_l2_[illp].second, inp_.L);
                    if (not std::isfinite(fs[ill]))
                        Exception("Failed to evaluate the angular integral f[%d](%d,%d,%d,%d;%d).", lambda, l1_l2_[ill].first, l1_l2_[ill].second, l1_l2_[illp].first, l1_l2_[illp].second, inp_.L);
                }
                
                // for all sub-blocks of the radial matrix
                for (std::pair<int,int> bpos : structure)
                {
                    // calculate the radial sub-block
                    cArray R_block_ik = s_rad_.calc_R_tr_dia_block(lambda, bpos.first, bpos.second, Mtr_L[lambda], Mtr_mLm1[lambda], cmd_.parallel_dot);
                    
                    // assemble all symmetries (one or two positions)
                    std::vector<std::pair<int,int>> syms = { bpos };
                    if (bpos.first != bpos.second)
                        syms.push_back(std::make_pair(bpos.second, bpos.first));
                    
                    // multiply the appropriate segments of the source vector (which is contained in 'buffer')
                    for (std::pair<int,int> sym : syms)
                    {
                        cArray product = SymDiaMatrix::sym_dia_dot(Nspline, s_rad_.S().diag(), R_block_ik.data(), buffer.data() + sym.second * Nspline);
                        
                        // for all destination super-segments
                        for (int ill = 0; ill < Nang; ill++)
                        {
                            // if eligible, update the sub-segment in the super-segment 'q(ill)'
                            if (fs[ill] != 0 and par_.isMyWork(ill))
                                cArrayView(q, (ill / par_.Nproc()) * Nspline * Nspline + sym.first * Nspline, Nspline) -= fs[ill] * product;
                        }
                        
                    } // for sub-block positions
                    
                } // for sub-blocks
                
            } // for lambda
            
        } // for illp
        
    } // if not cmd_.lightweight_radial_cache
}
