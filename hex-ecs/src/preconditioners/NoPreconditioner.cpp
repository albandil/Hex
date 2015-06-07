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

#ifdef _OPENMP
    #include <omp.h>
    #define OMP_prepare omp_lock_t writelock; omp_init_lock(&writelock)
    #define OMP_exclusive_in omp_set_lock(&writelock)
    #define OMP_exclusive_out omp_unset_lock(&writelock)
    #define OMP_clean omp_destroy_lock(&writelock)
#else
    #define OMP_prepare
    #define OMP_exclusive_in
    #define OMP_exclusive_out
    #define OMP_clean
#endif

#include "../arrays.h"
#include "../gauss.h"
#include "../misc.h"
#include "../parallel.h"
#include "../preconditioners.h"
#include "../radial.h"

const std::string NoPreconditioner::prec_name = "none";
const std::string NoPreconditioner::prec_description = "\"Preconditioning\" by the identity matrix.";

void NoPreconditioner::setup ()
{
    // TODO : Determine which lambdas are needed by this process.
    // NOTE : At the moment each process holds in memory radial integrals
    //        for all lambdas, which needlessly raises memory requirements.
    Array<bool> lambdas (inp_.L + 2 * inp_.levels + 1, true);
    
    // compute one-electron radial integrals
    s_rad_.setupOneElectronIntegrals(par_, cmd_);
    
    // compute two-eletron radial integrals
    s_rad_.setupTwoElectronIntegrals(par_, cmd_, lambdas);
}

void NoPreconditioner::update (double E)
{
    OMP_prepare;
    
    // shorthands
    unsigned Nspline = s_bspline_.Nspline();
    unsigned order = s_bspline_.order();
    
    // update energy
    E_ = E;
    
    // skip pre-calculation of the diagonal blocks in full lightweight mode
    if (cmd_.lightweight_full)
        return;
    
    std::cout << "\tPrecompute diagonal blocks... " << std::flush;
    
    // setup diagonal blocks
    # pragma omp parallel for if (cmd_.parallel_block)
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyWork(ill))
    {
        // angular momenta
        int l1 = l1_l2_[ill].first;
        int l2 = l1_l2_[ill].second;
        
        // initialize diagonal block
        dia_blocks_[ill] = BlockSymBandMatrix<Complex>
        (
            s_bspline_.Nspline(),       // block count (and size)
            s_bspline_.order() + 1,     // half-bandwidth
            !cmd_.outofcore,            // keep in memory?
            format("dblk-%d.ooc", ill)  // scratch disk file name
        );
        
        // skip calculation if the disk file is already present
        if (cmd_.outofcore and cmd_.reuse_dia_blocks and dia_blocks_[ill].hdfcheck())
            continue;
        else if (cmd_.outofcore)
            dia_blocks_[ill].hdfinit();
        
        // for all blocks
        # pragma omp parallel for schedule (dynamic,1)
        for (unsigned i = 0; i < Nspline; i++)
        for (unsigned d = 0; d <= order; d++)
        if (i + d < Nspline)
        {
            unsigned j = i + d;
            
            // one-electron part
            Complex half (0.5,0.0);
            SymBandMatrix<Complex> block = E * s_rad_.S()(i,j) * s_rad_.S();
            block -= (half * s_rad_.D()(i,j) - s_rad_.Mm1_tr()(i,j)) * s_rad_.S();
            block -= 0.5 * l1 * (l1 + 1) * s_rad_.Mm2()(i,j) * s_rad_.S();
            block -= s_rad_.S()(i,j) * (half * s_rad_.D() - s_rad_.Mm1_tr());
            block -= 0.5 * l2 * (l2 + 1) * s_rad_.S()(i,j) * s_rad_.Mm2();
            
            // two-electron part
            for (int lambda = 0; lambda <= s_rad_.maxlambda(); lambda++)
            {
                // calculate angular integral
                Complex f = special::computef(lambda,l1,l2,l1,l2,inp_.L);
                
                // check that the "f" coefficient is valid (no factorial overflow etc.)
                if (not Complex_finite(f))
                    HexException("Overflow in computation of f[%d](%d,%d,%d,%d).", inp_.L, l1, l2, l1, l2);
                
                // check that the "f" coefficient is nonzero
                if (f == 0.)
                    continue;
                
                // calculate two-electron term
                if (not cmd_.lightweight_radial_cache)
                {
                    // use precomputed block ... 
                    if (not cmd_.cache_all_radint) // ... from scratch file
                    {
                        OMP_exclusive_in;
                        block.data() += (-f) * s_rad_.R_tr_dia(lambda).getBlock(i * (order + 1) + d);
                        OMP_exclusive_out;
                    }
                    else // ... from memory
                    {
                        block.data() += (-f) * s_rad_.R_tr_dia(lambda).getBlock(i * (order + 1) + d);
                    }
                }
                else
                {
                    // compute the data anew
                    block += (-f) * s_rad_.calc_R_tr_dia_block(lambda, i, j);
                }
            }
            
            // save block
            OMP_exclusive_in;
            dia_blocks_[ill].setBlock(i * (order + 1) + d, block.data());
            OMP_exclusive_out;
        }
    }
    
    par_.wait();
    std::cout << "ok" << std::endl;
    
    OMP_clean;
}

void NoPreconditioner::rhs (BlockArray<Complex> & chi, int ie, int instate, int Spin) const
{
    OMP_prepare;
    
    // shorthands
    int ni = std::get<0>(inp_.instates[instate]);
    int li = std::get<1>(inp_.instates[instate]);
    int mi = std::get<2>(inp_.instates[instate]);
    
    // shorthands
    int Nspline = s_rad_.bspline().Nspline();
    
    // impact momentum
    rArray ki = { std::sqrt(inp_.Etot[ie] + 1./(ni*ni)) };
    
    // j-overlaps of shape [Nangmom × Nspline]
    cArray ji_overlaps = s_rad_.overlapj(inp_.maxell, ki, weightEdgeDamp(s_rad_.bspline()));
    ji_overlaps.hdfsave("ji_overlaps.hdf");
    if (not std::isfinite(ji_overlaps.norm()))
        HexException("Unable to compute Riccati-Bessel function B-spline overlaps!");
    
    // j-expansions
    cArray ji_expansion = s_rad_.S().tocoo<LU_int_t>().tocsr().solve(ji_overlaps, ji_overlaps.size() / Nspline);
    if (not std::isfinite(ji_expansion.norm()))
        HexException("Unable to expand Riccati-Bessel function in B-splines!");
    
    // compute P-overlaps and P-expansion
    cArray Pi_overlaps, Pi_expansion;
    Pi_overlaps = s_rad_.overlapP(ni, li, weightEndDamp(s_rad_.bspline()));
    Pi_expansion = s_rad_.S().tocoo<LU_int_t>().tocsr().solve(Pi_overlaps);
    if (not std::isfinite(Pi_expansion.norm()))
        HexException("Unable to expand hydrogen bound orbital in B-splines!");
    
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
                              * special::ClebschGordan(li,mi, l,0, inp_.L,mi) / ki[0];
            
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
                    HexException("Invalid result of computef(%d,%d,%d,%d,%d,%d)\n", lambda,l1,l2,li,l,inp_.L);
                if (not std::isfinite(f2))
                    HexException("Invalid result of computef(%d,%d,%d,%d,%d,%d)\n", lambda,l1,l2,l,li,inp_.L);
                
                // add multipole terms (direct/exchange)
                if (not cmd_.lightweight_radial_cache)
                {
                    if (f1 != 0.) chi_block += (       prefactor * f1) * s_rad_.R_tr_dia(lambda).dot(Pj1, true);
                    if (f2 != 0.) chi_block += (Sign * prefactor * f2) * s_rad_.R_tr_dia(lambda).dot(Pj2, true);
                }
                else
                {
                    if (f1 != 0.) chi_block += (       prefactor * f1) * s_rad_.apply_R_matrix(lambda, Pj1);
                    if (f2 != 0.) chi_block += (Sign * prefactor * f2) * s_rad_.apply_R_matrix(lambda, Pj2);
                }
            }
            
            // add monopole terms (direct/exchange)
            if (li == l1 and l == l2)
                chi_block += (-prefactor       ) * outer_product(s_rad_.S().dot(Pi_expansion), s_rad_.Mm1_tr().dot(Ji_expansion));
            if (li == l2 and l == l1)
                chi_block += (-prefactor * Sign) * outer_product(s_rad_.Mm1_tr().dot(Ji_expansion), s_rad_.S().dot(Pi_expansion));
            
            // update the right-hand side
            OMP_exclusive_in;
            chi[ill] = chi_block;
            if (not chi.inmemory())
            {
                chi.hdfsave(ill);
                chi[ill].drop();
            }
            OMP_exclusive_out;
        }
    }
    
    OMP_clean;
}

void NoPreconditioner::multiply (BlockArray<Complex> const & p, BlockArray<Complex> & q) const
{
    OMP_prepare;
    
    // shorthands
    int Nspline = s_rad_.bspline().Nspline();
    int order = s_rad_.bspline().order();
    int Nang = l1_l2_.size();
    int Nchunk = Nspline * Nspline;
    
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
            // load data from scratch disk
            if (cmd_.outofcore)
            {
                const_cast<BlockArray<Complex>&>(p).hdfload(ill);
                
                q.hdfload(ill);
                
                if (cmd_.wholematrix)
                    const_cast<BlockSymBandMatrix<Complex>&>(dia_blocks_[ill]).hdfload();
            }
            
            // multiply
            q[ill] = dia_blocks_[ill].dot(p[ill], true);
            
            // unload data
            if (cmd_.outofcore)
            {
                const_cast<BlockArray<Complex>&>(p)[ill].drop();
                
                q.hdfsave(ill);
                q[ill].drop();
                
                if (cmd_.wholematrix)
                    const_cast<BlockSymBandMatrix<Complex>&>(dia_blocks_[ill]).drop();
            }
        }
        
        // multiply "p" by the off-diagonal blocks - single proces
        if (par_.Nproc() == 1)
        for (int lambda = 0; lambda <= s_rad_.maxlambda(); lambda++)
        {
            // load data from scratch disk
            if (not cmd_.cache_own_radint and cmd_.wholematrix)
                const_cast<BlockSymBandMatrix<Complex>&>(s_rad_.R_tr_dia(lambda)).hdfload();
            
            // update all blocks with this multipole potential matrix
            for (int ill = 0;  ill < Nang;  ill++)
            {
                if (cmd_.outofcore)
                    q.hdfload(ill);
                
                for (int illp = 0; illp < Nang; illp++)
                {
                    // skip diagonal
                    if (ill == illp)
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
                        HexException("Invalid result of computef(%d,%d,%d,%d,%d,%d).", lambda, l1, l2, l1p, l2p, inp_.L);
                    
                    // check non-zero
                    if (f == 0.)
                        continue;
                    
                    // load data
                    if (cmd_.outofcore)
                        const_cast<BlockArray<Complex>&>(p).hdfload(illp);
                    
                    // calculate product
                    q[ill] += (-f) * s_rad_.R_tr_dia(lambda).dot(p[illp], true);
                    
                    // unload data
                    if (cmd_.outofcore)
                        const_cast<BlockArray<Complex>&>(p)[illp].drop();
                }
                
                if (cmd_.outofcore)
                {
                    q.hdfsave(ill);
                    q[ill].drop();
                }
            }
            
            // release radial integrals
            if (not cmd_.cache_own_radint)
                const_cast<BlockSymBandMatrix<Complex>&>(s_rad_.R_tr_dia(lambda)).drop();
        }
        
        // multiply "p" by the off-diagonal blocks multiprocess
        if (par_.Nproc() > 1)
        for (int ill = 0;  ill < Nang;  ill++)
        {
            // product of line of blocks with the source vector (-> one segment of destination vector)
            cArray product (Nchunk);
            
            // load data
            if (cmd_.outofcore and par_.isMyWork(ill))
                q.hdfload(ill);
            
            // maximal multipole
            int maxlambda = s_rad_.maxlambda();
            
            # pragma omp parallel for schedule (dynamic,1) if (cmd_.parallel_block)
            for (int illp = 0; illp < Nang; illp++)
            if (par_.isMyWork(illp))
            {
                if (cmd_.outofcore)
                    const_cast<BlockArray<Complex>&>(p).hdfload(illp);
                
                for (int lambda = 0; lambda <= maxlambda; lambda++)
                {
                    // skip diagonal
                    if (ill == illp)
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
                        HexException("Invalid result of computef(%d,%d,%d,%d,%d,%d).", lambda, l1, l2, l1p, l2p, inp_.L);
                    
                    // check non-zero
                    if (f == 0.)
                        continue;
                    
                    // calculate product
                    cArray p0 = std::move( (-f) * s_rad_.R_tr_dia(lambda).dot(p[illp], true) );
                    
                    // update collected product
                    OMP_exclusive_in;
                    product += p0;
                    OMP_exclusive_out;
                }
                
                if (cmd_.outofcore)
                    const_cast<BlockArray<Complex>&>(p)[illp].drop();
            }
            
            // sum all contributions to this destination segment on its owner node
            par_.sum(product.data(), Nchunk, ill % par_.Nproc());
            
            // finally, owner will update its segment (and move back to disk, if OOC)
            if (par_.isMyWork(ill))
            {
                q[ill] += product;
                
                if (cmd_.outofcore)
                {
                    q.hdfsave(ill);
                    q[ill].drop();
                }
            }
        }
    }
    else
    {
        //
        // Lightweight mode multiplication
        //
        
        // multiply "p" by the diagonal super-blocks
        # pragma omp parallel for schedule (dynamic,1) if (cmd_.parallel_block)
        for (int ill = 0;  ill < Nang;  ill++)
        if (par_.isMyWork(ill))
        {
            if (cmd_.outofcore)
            {
                const_cast<BlockArray<Complex>&>(p).hdfload(ill);
                q.hdfload(ill);
            }
            
            // get block angular momemnta
            int l1 = l1_l2_[ill].first;
            int l2 = l1_l2_[ill].second;
            
            // multiply 'p_block' by the diagonal block (except for the two-electron term)
            q[ill]  = kron_dot(Complex(E_) * s_rad_.S(), s_rad_.S(), p[ill]);
            q[ill] -= kron_dot(Complex(0.5) * s_rad_.D() - s_rad_.Mm1_tr() + Complex(0.5*l1*(l1+1)) * s_rad_.Mm2(), s_rad_.S(), p[ill]);
            q[ill] -= kron_dot(s_rad_.S(), Complex(0.5) * s_rad_.D() - s_rad_.Mm1_tr() + Complex(0.5*l2*(l2+1)) * s_rad_.Mm2(), p[ill]);
            
            if (cmd_.outofcore)
            {
                const_cast<BlockArray<Complex>&>(p)[ill].drop();
                q.hdfsave(ill);
                q[ill].drop();
            }
        }
        
        // auxiliary buffers
        cArray buffer (Nang * Nspline);
        
        // precalculate angular integrals
        double fs[Nang][Nang][s_rad_.maxlambda() + 1];
        for (int ill  = 0; ill  < Nang; ill ++)
        for (int illp = 0; illp < Nang; illp++)
        for (int lambda = 0; lambda <= s_rad_.maxlambda(); lambda++)
        {
            fs[ill][illp][lambda] = special::computef(lambda, l1_l2_[ill].first, l1_l2_[ill].second, l1_l2_[illp].first, l1_l2_[illp].second, inp_.L);
            if (not std::isfinite(fs[ill][illp][lambda]))
                HexException("Failed to evaluate the angular integral f[%d](%d,%d,%d,%d;%d).", lambda, l1_l2_[ill].first, l1_l2_[ill].second, l1_l2_[illp].first, l1_l2_[illp].second, inp_.L);
        }
        
#ifdef _OPENMP
        // I/O access locks
        std::vector<omp_lock_t> locks(Nang);
        for (omp_lock_t & lock : locks)
            omp_init_lock(&lock);
#endif
        
        // for all source vector sub-segments
        for (int k = 0; k < Nspline; k++)
        {
            // copy owned sub-segments to the buffer
            for (int illp = 0; illp < Nang; illp++) if (par_.isMyWork(illp))
            {
                std::memcpy
                (
                    &buffer[0] + illp * Nspline,
                    p.segment(illp, k * Nspline, Nspline).ptr(),
                    Nspline * sizeof(Complex)
                );
            }
            
            // synchronize source sub-segments buffer across processes
            par_.sync(&buffer[0], Nspline, Nang);
            
            // auxiliary variables
            int min_i = std::max(0, k - order);
            int max_i = std::min(k + order, Nspline - 1);
            int maxlambda = s_rad_.maxlambda();
            
            // for all destination sub-blocks
            #pragma omp parallel for collapse (2)
            for (int i = min_i; i <= max_i; i++)
            {
                // for all potential multipoles
                for (int lambda = 0; lambda <= maxlambda; lambda++)
                {
                    // calculate the radial sub-block
                    SymBandMatrix<Complex> R_block_ik = s_rad_.calc_R_tr_dia_block(lambda, i, k);
                    
                    // apply all superblocks
                    for (int ill  = 0; ill  < Nang; ill ++) if (par_.isMyWork(ill))
                    {
                        // collected products of superblocks in this row
                        cArray product (Nspline);
                        
                        // for all superblocks in this row
                        for (int illp = 0; illp < Nang; illp++) if (fs[ill][illp][lambda] != 0.)
                        {
                            // multiply sub-segment by the R block
                            product += fs[ill][illp][lambda] * R_block_ik.dot(cArrayView(buffer, illp * Nspline, Nspline));
                        }
                        
                        // atomic update of the owned sub-segment
#ifdef _OPENMP
                        omp_set_lock(&locks[ill]);
#endif
                        q.setSegment(ill, i * Nspline, Nspline, q.segment(ill, i * Nspline, Nspline)() - product);
#ifdef _OPENMP
                        omp_unset_lock(&locks[ill]);
#endif
                    }
                }
            }
        }
        
#ifdef _OPENMP
        for (omp_lock_t & lock : locks)
            omp_destroy_lock(&lock);
#endif
    } // if not cmd_.lightweight_radial_cache
    
    OMP_clean;
}
