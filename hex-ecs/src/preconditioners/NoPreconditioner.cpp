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
    rad_.setupOneElectronIntegrals(par_, cmd_);
    rad_.setupTwoElectronIntegrals(par_, cmd_);
}

void NoPreconditioner::update (double E)
{
    OMP_prepare;
    
    // shorthands
    unsigned order = inp_.order;
    unsigned Nspline_atom = bspline_atom_.Nspline();
    unsigned Nspline_proj = bspline_proj_.Nspline();
    
    // update energy
    E_ = E;
    
    // skip pre-calculation of the diagonal blocks in full lightweight mode
    if (cmd_.lightweight_full)
        return;
    
    std::cout << "\tPrecompute diagonal blocks... " << std::flush;
    
    // setup diagonal blocks
    # pragma omp parallel for if (cmd_.parallel_multiply)
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyGroupWork(ill))
    {
        // angular momenta
        int l1 = l1_l2_[ill].first;
        int l2 = l1_l2_[ill].second;
        
        // initialize diagonal block
        dia_blocks_[ill] = BlockSymBandMatrix<Complex>
        (
            Nspline_atom,               // block count
            inp_.order + 1,             // block structure half-bandwidth
            Nspline_proj,               // block size
            inp_.order + 1,             // block half-bandwidth
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
        for (unsigned i = 0; i < Nspline_atom; i++)
        for (unsigned d = 0; d <= order; d++)
        if (i + d < Nspline_atom)
        {
            unsigned j = i + d;
            
            // one-electron part
            Complex half (0.5,0.0);
            SymBandMatrix<Complex> block = E * rad_.S_atom()(i,j) * rad_.S_proj();
            block -= (half * rad_.D_atom()(i,j) - rad_.Mm1_tr_atom()(i,j)) * rad_.S_proj();
            block -= 0.5 * l1 * (l1 + 1) * rad_.Mm2_atom()(i,j) * rad_.S_proj();
            block -= rad_.S_atom()(i,j) * (half * rad_.D_proj() - rad_.Mm1_tr_proj());
            block -= 0.5 * l2 * (l2 + 1) * rad_.S_atom()(i,j) * rad_.Mm2_proj();
            
            // two-electron part
            for (int lambda = 0; lambda <= rad_.maxlambda(); lambda++)
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
                        block.data() += (-f) * rad_.R_tr_dia(lambda).getBlock(i * (order + 1) + d);
                        OMP_exclusive_out;
                    }
                    else // ... from memory
                    {
                        block.data() += (-f) * rad_.R_tr_dia(lambda).getBlock(i * (order + 1) + d);
                    }
                }
                else
                {
                    // compute the data anew
                    block += (-f) * rad_.calc_R_tr_dia_block(lambda, i, j);
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

void NoPreconditioner::rhs (BlockArray<Complex> & chi, int ie, int instate, int Spin, Bspline const & bspline_proj_full) const
{
    // shorthands
    int ni = std::get<0>(inp_.instates[instate]);
    int li = std::get<1>(inp_.instates[instate]);
    int mi = std::get<2>(inp_.instates[instate]);
    
    // shorthands
    int Nspline_atom = rad_.bspline_atom().Nspline();
    int Nspline_proj = rad_.bspline_proj().Nspline();
    
    // impact momentum
    rArray ki = { std::sqrt(inp_.Etot[ie] + 1./(ni*ni)) };
    
    // radial information for full projectil B-spline basis (used only to expand Riccati-Bessel function)
    RadialIntegrals radf (rad_.bspline_atom(), bspline_proj_full, 0);
    radf.verbose(false);
    radf.setupOneElectronIntegrals(par_, cmd_);
    
    // calculate LU-decomposition of the overlap matrix
    CsrMatrix<LU_int_t,Complex> S_csr_atom = rad_.S_atom().tocoo<LU_int_t>().tocsr();
    CsrMatrix<LU_int_t,Complex> S_csr_proj = radf.S_proj().tocoo<LU_int_t>().tocsr();
    std::shared_ptr<LUft<LU_int_t,Complex>> lu_S_atom = S_csr_atom.factorize();
    std::shared_ptr<LUft<LU_int_t,Complex>> lu_S_proj = S_csr_proj.factorize();
    
    // j-overlaps of shape [Nangmom × Nspline]
    cArray ji_overlaps_atom = rad_.overlapj(rad_.bspline_atom(), rad_.gaussleg_atom(), inp_.maxell, ki, weightEdgeDamp(rad_.bspline_atom()));
    cArray ji_overlaps_proj = radf.overlapj(radf.bspline_proj(), radf.gaussleg_proj(), inp_.maxell, ki, weightEdgeDamp(radf.bspline_proj()));
    if (not std::isfinite(ji_overlaps_atom.norm()) or not std::isfinite(ji_overlaps_proj.norm()))
        HexException("Unable to compute Riccati-Bessel function B-spline overlaps!");
    
    // j-expansions
    cArray ji_expansion_atom = lu_S_atom->solve(ji_overlaps_atom, inp_.maxell + 1);
    cArray ji_expansion_proj = lu_S_proj->solve(ji_overlaps_proj, inp_.maxell + 1);
    if (not std::isfinite(ji_expansion_atom.norm()) or not std::isfinite(ji_expansion_proj.norm()))
        HexException("Unable to expand Riccati-Bessel function in B-splines!");
    
    // compute P-overlaps and P-expansion
    cArray Pi_overlaps_atom = rad_.overlapP(rad_.bspline_atom(), rad_.gaussleg_atom(), ni, li, weightEndDamp(rad_.bspline_atom()));
    cArray Pi_overlaps_proj = radf.overlapP(radf.bspline_proj(), radf.gaussleg_proj(), ni, li, weightEndDamp(radf.bspline_proj()));
    cArray Pi_expansion_atom = lu_S_atom->solve(Pi_overlaps_atom);
    cArray Pi_expansion_proj = lu_S_proj->solve(Pi_overlaps_proj);
    if (not std::isfinite(Pi_expansion_atom.norm()) or not std::isfinite(Pi_expansion_proj.norm()))
        HexException("Unable to expand hydrogen bound orbital in B-splines!");
    
    // truncate the projectile expansions for non-origin panels
    if (rad_.bspline_proj().Nspline() != radf.bspline_proj().Nspline())
    {
        // truncate all Riccati-Bessel function expansion (for various angular momenta)
        cArrays ji_expansion_proj_trunc;
        for (int i = 0; i <= inp_.maxell; i++)
            ji_expansion_proj_trunc.push_back(cArrayView(ji_expansion_proj, (i + 1) * radf.bspline_proj().Nspline() - Nspline_proj, Nspline_proj));
        ji_expansion_proj = join(ji_expansion_proj_trunc);
        
        // truncate the hydrogen radial orbital expansion
        Pi_expansion_proj = Pi_expansion_proj.slice(Pi_expansion_proj.size() - Nspline_proj, Pi_expansion_proj.size());
    }
    
    // for all segments constituting the RHS
    # pragma omp parallel for schedule (dynamic,1) if (cmd_.parallel_multiply)
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyGroupWork(ill))
    {
        int l1 = l1_l2_[ill].first;
        int l2 = l1_l2_[ill].second;
        
        // setup storage
        cArray chi_block (Nspline_atom * Nspline_proj);
        
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
            cArrayView Ji_expansion_atom (ji_expansion_atom, l * Nspline_atom, Nspline_atom);
            cArrayView Ji_expansion_proj (ji_expansion_proj, l * Nspline_proj, Nspline_proj);
            
            // compute outer products of B-spline expansions
            cArray Pj1 = outer_product(Pi_expansion_atom, Ji_expansion_proj);
            cArray Pj2 = outer_product(Ji_expansion_atom, Pi_expansion_proj);
            
            // skip angular forbidden right hand sides
            for (int lambda = 0; lambda <= rad_.maxlambda(); lambda++)
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
                    if (f1 != 0.) chi_block += (       prefactor * f1) * rad_.R_tr_dia(lambda).dot(Pj1, true);
                    if (f2 != 0.) chi_block += (Sign * prefactor * f2) * rad_.R_tr_dia(lambda).dot(Pj2, true);
                }
                else
                {
                    if (f1 != 0.) chi_block += (       prefactor * f1) * rad_.apply_R_matrix(lambda, Pj1);
                    if (f2 != 0.) chi_block += (Sign * prefactor * f2) * rad_.apply_R_matrix(lambda, Pj2);
                }
            }
            
            // add monopole terms (direct/exchange)
            if (li == l1 and l == l2)
                chi_block += (-prefactor       ) * outer_product(rad_.S_atom().dot(Pi_expansion_atom), rad_.Mm1_tr_proj().dot(Ji_expansion_proj));
            if (li == l2 and l == l1)
                chi_block += (-prefactor * Sign) * outer_product(rad_.Mm1_tr_atom().dot(Ji_expansion_atom), rad_.S_proj().dot(Pi_expansion_proj));
        }
        
        // update the right-hand side
        chi[ill] = chi_block;
        if (not chi.inmemory())
        {
            chi.hdfsave(ill);
            chi[ill].drop();
        }
    }
}

void NoPreconditioner::multiply (BlockArray<Complex> const & p, BlockArray<Complex> & q) const
{
    OMP_prepare;
    
    // shorthands
    int order = rad_.bspline_atom().order();
    int Nspline_atom = rad_.bspline_atom().Nspline();
    int Nspline_proj = rad_.bspline_proj().Nspline();
    int Nang = l1_l2_.size();
    int Nchunk = Nspline_atom * Nspline_proj;
    
    // precalculate angular integrals
    rArray ffs (Nang * Nang * (rad_.maxlambda() + 1));
    # define fs(i,j,k) (ffs[k + (j + i * Nang) * (rad_.maxlambda() + 1)])
    for (int ill  = 0; ill  < Nang; ill ++)
    for (int illp = 0; illp < Nang; illp++)
    for (int lambda = 0; lambda <= rad_.maxlambda(); lambda++)
    {
        fs(ill,illp,lambda) = special::computef(lambda, l1_l2_[ill].first, l1_l2_[ill].second, l1_l2_[illp].first, l1_l2_[illp].second, inp_.L);
        if (not std::isfinite(fs(ill,illp,lambda)))
            HexException("Failed to evaluate the angular integral f[%d](%d,%d,%d,%d;%d).", lambda, l1_l2_[ill].first, l1_l2_[ill].second, l1_l2_[illp].first, l1_l2_[illp].second, inp_.L);
    }
    
    if (not cmd_.lightweight_radial_cache)
    {
        //
        // Simple multiplication by the super-matrix.
        //
        
        // multiply "p" by the diagonal blocks
        # pragma omp parallel for schedule (dynamic,1) if (cmd_.parallel_multiply)
        for (int ill = 0;  ill < Nang;  ill++)
        if (par_.isMyGroupWork(ill))
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
        for (int lambda = 0; lambda <= rad_.maxlambda(); lambda++)
        {
            // load data from scratch disk
            if (not cmd_.cache_own_radint and cmd_.wholematrix)
                const_cast<BlockSymBandMatrix<Complex>&>(rad_.R_tr_dia(lambda)).hdfload();
            
            // update all blocks with this multipole potential matrix
            for (int ill = 0;  ill < Nang;  ill++)
            {
                if (cmd_.outofcore)
                    q.hdfload(ill);
                
                # pragma omp parallel for schedule (dynamic,1) if (cmd_.parallel_multiply)
                for (int illp = 0; illp < Nang; illp++)
                {
                    // skip diagonal
                    if (ill == illp)
                        continue;
                    
                    // check non-zero
                    if (fs(ill,illp,lambda) == 0.)
                        continue;
                    
                    // load data
                    if (cmd_.outofcore)
                        const_cast<BlockArray<Complex>&>(p).hdfload(illp);
                    
                    // calculate product
                    q[ill] += (-fs(ill,illp,lambda)) * rad_.R_tr_dia(lambda).dot(p[illp], true);
                    
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
                const_cast<BlockSymBandMatrix<Complex>&>(rad_.R_tr_dia(lambda)).drop();
        }
        
        // multiply "p" by the off-diagonal blocks (multi-process)
        if (par_.Nproc() > 1)
        for (int ill = 0;  ill < Nang;  ill++)
        {
            // product of line of blocks with the source vector (-> one segment of destination vector)
            cArray product (Nchunk);
            
            // load data
            if (cmd_.outofcore and par_.isMyGroupWork(ill))
                q.hdfload(ill);
            
            // maximal multipole
            int maxlambda = rad_.maxlambda();
            
            // for all source segments that this group owns
            # pragma omp parallel for schedule (dynamic,1) if (cmd_.parallel_multiply)
            for (int illp = 0; illp < Nang; illp++) if (par_.isMyGroupWork(illp))
            {
                // load the segment, if on disk
                if (cmd_.outofcore)
                    const_cast<BlockArray<Complex>&>(p).hdfload(illp);
                
                // for all multipoles (each group's process has a different work)
                for (int lambda = 0; lambda <= maxlambda; lambda++) if (lambda % par_.groupsize() == par_.igroupproc())
                {
                    // skip diagonal
                    if (ill == illp)
                        continue;
                    
                    // check non-zero
                    if (fs(ill,illp,lambda) == 0.)
                        continue;
                    
                    // calculate product
                    cArray p0 = std::move( (-fs(ill,illp,lambda)) * rad_.R_tr_dia(lambda).dot(p[illp], true) );
                    
                    // update collected product
                    OMP_exclusive_in;
                    product += p0;
                    OMP_exclusive_out;
                }
                
                if (cmd_.outofcore)
                    const_cast<BlockArray<Complex>&>(p)[illp].drop();
            }
            
            // sum all group processes contributions on that group master process
            par_.sum_g(product.data(), product.size(), 0);
            
            // sum all master segments to the master of the group that owns the destination segment
            par_.mastersum(product.data(), product.size(), ill % par_.Ngroup());
            
            // redistribute the summed segment over the whole owning group
            par_.bcast_g(ill % par_.Ngroup(), 0, &product[0], product.size());
            
            // finally, owner will update its segment (and move back to disk, if OOC)
            if (par_.isMyGroupWork(ill))
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
        # pragma omp parallel for schedule (dynamic,1) if (cmd_.parallel_multiply)
        for (int ill = 0;  ill < Nang;  ill++) if (par_.isMyGroupWork(ill))
        {
            std::cout << "before 1el: " << p[ill].norm() << std::endl;
            if (cmd_.outofcore)
            {
                const_cast<BlockArray<Complex>&>(p).hdfload(ill);
                q.hdfload(ill);
            }
            
            // get block angular momemnta
            int l1 = l1_l2_[ill].first;
            int l2 = l1_l2_[ill].second;
            
            // multiply 'p_block' by the diagonal block (except for the two-electron term)
            q[ill]  = kron_dot(Complex(E_) * rad_.S_atom(), rad_.S_proj(), p[ill]);
            q[ill] -= kron_dot(Complex(0.5) * rad_.D_atom() - rad_.Mm1_tr_atom() + Complex(0.5*l1*(l1+1)) * rad_.Mm2_atom(), rad_.S_proj(), p[ill]);
            q[ill] -= kron_dot(rad_.S_atom(), Complex(0.5) * rad_.D_proj() - rad_.Mm1_tr_proj() + Complex(0.5*l2*(l2+1)) * rad_.Mm2_proj(), p[ill]);
            
            if (cmd_.outofcore)
            {
                const_cast<BlockArray<Complex>&>(p)[ill].drop();
                q.hdfsave(ill);
                q[ill].drop();
            }
            
            std::cout << "after 1el: " << q[ill].norm() << std::endl;
        }
        
        // auxiliary buffers
        cArray buffer (Nang * Nspline_proj);
        
#ifdef _OPENMP
        // I/O access locks
        std::vector<omp_lock_t> locks(Nang);
        for (omp_lock_t & lock : locks)
            omp_init_lock(&lock);
#endif
        
        // for all source vector sub-segments
        for (int k = 0; k < Nspline_atom; k++)
        {
            // copy owned sub-segments to the buffer
            for (int illp = 0; illp < Nang; illp++) if (par_.isMyGroupWork(illp) and par_.IamGroupMaster())
            {
                std::memcpy
                (
                    &buffer[0] + illp * Nspline_proj,
                    p.segment(illp, k * Nspline_proj, Nspline_proj).ptr(),
                    Nspline_proj * sizeof(Complex)
                );
            }
            
            // synchronize source sub-segments buffer across groups' masters
            par_.sync_m(&buffer[0], Nspline_proj, Nang);
            
            // broadcast buffer to all members of the group
            par_.bcast_g(par_.igroup(), 0, &buffer[0], Nang * Nspline_proj);
            
            // auxiliary variables
            int min_i = std::max(0, k - order);
            int max_i = std::min(k + order, Nspline_atom - 1);
            int maxlambda = rad_.maxlambda();
            
            // for all destination sub-blocks
            #pragma omp parallel for collapse (2)
            for (int i = min_i; i <= max_i; i++)
            {
                // for all potential multipoles
                for (int lambda = 0; lambda <= maxlambda; lambda++)
                {
                    // calculate the radial sub-block
                    SymBandMatrix<Complex> R_block_ik = rad_.calc_R_tr_dia_block(lambda, i, k);
                    
                    // apply all superblocks
                    for (int ill = 0; ill  < Nang; ill ++) if (par_.isMyGroupWork(ill))
                    {
                        // collected products of superblocks in this row
                        cArray product (Nspline_proj);
                        
                        // for all superblocks in this row
                        for (int illp = 0; illp < Nang; illp++) if (fs(ill,illp,lambda) != 0.)
                        {
                            // multiply sub-segment by the R block
                            product += fs(ill,illp,lambda) * R_block_ik.dot(cArrayView(buffer, illp * Nspline_proj, Nspline_proj));
                        }
                        
                        // atomic update of the owned sub-segment
#ifdef _OPENMP
                        omp_set_lock(&locks[ill]);
#endif
                        q.setSegment(ill, i * Nspline_proj, Nspline_proj, q.segment(ill, i * Nspline_proj, Nspline_proj)() - product);
#ifdef _OPENMP
                        omp_unset_lock(&locks[ill]);
#endif
                    }
                }
            }
        }
        
        for (int ill = 0;  ill < Nang;  ill++)
            std::cout << "after 2el: " << q[ill].norm() << std::endl;
        
#ifdef _OPENMP
        for (omp_lock_t & lock : locks)
            omp_destroy_lock(&lock);
#endif
    } // if not cmd_.lightweight_radial_cache
    
    OMP_clean;
}

void NoPreconditioner::finish ()
{
    dia_blocks_.resize(0);
}
