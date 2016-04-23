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

#include "hex-arrays.h"
#include "hex-misc.h"

#include "gauss.h"
#include "parallel.h"
#include "preconditioners.h"
#include "radial.h"

const std::string NoPreconditioner::prec_name = "none";
const std::string NoPreconditioner::prec_description = "\"Preconditioning\" by the identity matrix.";

void NoPreconditioner::setup ()
{
    rad_.setupOneElectronIntegrals(par_, cmd_);
    rad_.setupTwoElectronIntegrals(par_, cmd_);
    rad_.setupRadialEigenstates(par_, cmd_);
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
    for (unsigned ill = 0; ill < ang_.states().size(); ill++) if (par_.isMyGroupWork(ill))
    {
        // angular momenta
        int l1 = ang_.states()[ill].first;
        int l2 = ang_.states()[ill].second;
        
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
                // check that the "f" coefficient is nonzero
                if (ang_.f(ill,ill,lambda) == 0.)
                    continue;
                
                // calculate two-electron term
                if (not cmd_.lightweight_radial_cache)
                {
                    // use precomputed block ... 
                    if (not cmd_.cache_all_radint) // ... from scratch file
                    {
                        OMP_exclusive_in;
                        block.data() += (-ang_.f(ill,ill,lambda)) * rad_.R_tr_dia(lambda).getBlock(i * (order + 1) + d);
                        OMP_exclusive_out;
                    }
                    else // ... from memory
                    {
                        block.data() += (-ang_.f(ill,ill,lambda)) * rad_.R_tr_dia(lambda).getBlock(i * (order + 1) + d);
                    }
                }
                else
                {
                    // compute the data anew
                    block += Complex(-ang_.f(ill,ill,lambda)) * rad_.calc_R_tr_dia_block(lambda, i, j);
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

void NoPreconditioner::rhs (BlockArray<Complex> & chi, int ie, int instate) const
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
    
    // radial information for full projectile B-spline basis (used only to expand Riccati-Bessel function)
    RadialIntegrals radf (rad_.bspline_atom(), rad_.bspline_proj_full(), rad_.bspline_proj_full(), 0, 0);
    radf.verbose(false);
    radf.setupOneElectronIntegrals(par_, cmd_);
    
    // calculate LU-decomposition of the overlap matrix
    CsrMatrix<LU_int_t,Complex> S_csr_atom = rad_.S_atom().tocoo<LU_int_t>().tocsr();
    CsrMatrix<LU_int_t,Complex> S_csr_proj = radf.S_proj().tocoo<LU_int_t>().tocsr();
    std::shared_ptr<LUft<LU_int_t,Complex>> lu_S_atom = S_csr_atom.factorize();
    std::shared_ptr<LUft<LU_int_t,Complex>> lu_S_proj = S_csr_proj.factorize();
    
    // weight functions
    double Rp = cmd_.polarization;
    double R0 = bspline_atom_.R0();
    auto pol_pot = [Rp] (double r) { return -std::expm1(-special::pow_int(r/Rp,4)) / (r*r); };
    auto w_edge = [R0] (Complex z) { return damp(z.real(), R0); };
    auto w_polar = [R0,&w_edge,&pol_pot] (Complex z) { double r = z.real(); return w_edge(z) * pol_pot(r); };
    
    // j-overlaps of shape [Nangmom × Nspline]
    cArray ji_overlaps_atom = rad_.overlapj(rad_.bspline_atom(), rad_.gaussleg_atom(), inp_.maxell, ki, w_edge);
    cArray ji_overlaps_proj = radf.overlapj(radf.bspline_proj(), radf.gaussleg_proj(), inp_.maxell, ki, w_edge);
    cArray ji_overlaps_polar = rad_.overlapj(rad_.bspline_proj(), rad_.gaussleg_proj(), inp_.maxell, ki, w_polar);
    if (not std::isfinite(ji_overlaps_atom.norm()) or not std::isfinite(ji_overlaps_proj.norm()))
        HexException("Unable to compute Riccati-Bessel function B-spline overlaps!");
    
    // j-expansions
    cArray ji_expansion_atom = lu_S_atom->solve(ji_overlaps_atom, inp_.maxell + 1);
    cArray ji_expansion_proj = lu_S_proj->solve(ji_overlaps_proj, inp_.maxell + 1);
    cArray ji_expansion_polar = lu_S_atom->solve(ji_overlaps_polar, inp_.maxell + 1);
    if (not std::isfinite(ji_expansion_atom.norm()) or not std::isfinite(ji_expansion_proj.norm()) or not std::isfinite(ji_expansion_polar.norm()))
        HexException("Unable to expand Riccati-Bessel function in B-splines!");
    
    rArray grid = linspace(0., bspline_atom_.Rmax(), 10000);
    for (int ell = 0; ell <= inp_.maxell; ell++)
    {
        write_array
        (
            grid,
            bspline_atom_.zip
            (
                ji_expansion_polar.slice(ell * Nspline_atom, (ell + 1) * Nspline_atom),
                grid
            ),
            format("ji-%d.dat", ell)
        );
    }
    
    // compute P-overlaps and P-expansion
    cArray Pi_expansion_atom = rad_.getstate(ni,li,rad_.eigenstates(li));
    cArray Pi_expansion_proj = rad_.getstate(ni,li,rad_.eigenstates(li)); // FIXME : Proj
    double Eni = rad_.eigenenergy(ni, li).real();
    
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
    
    // (anti)symmetrization w.r.t. spin multiplet
    double Sign = special::pow_int(-1, ang_.S());
    
    // for all segments constituting the RHS
    for (unsigned ill = 0; ill < ang_.states().size(); ill++) if (par_.isMyGroupWork(ill))
    {
        int l1 = ang_.states()[ill].first;
        int l2 = ang_.states()[ill].second;
        
        cArray chi_block (Nspline_atom * Nspline_proj);
        
        for (unsigned illp = 0; illp < ang_.states().size(); illp++)
        {
            int l1p = ang_.states()[illp].first;
            int l2p = ang_.states()[illp].second;
            
            // calculate other angular integrals
            rArray f1 (rad_.maxlambda() + 1), f2 (rad_.maxlambda() + 1);
            for (int lambda = 0; lambda <= rad_.maxlambda(); lambda++)
            {
                f1[lambda] = special::computef(lambda, l1, l2, l1p, l2p, inp_.L);
                f2[lambda] = special::computef(lambda, l1, l2, l2p, l1p, inp_.L);
                
                // abort if any of the coefficients is non-number (factorial overflow etc.)
                if (not std::isfinite(f1[lambda]))
                    HexException("Invalid result of computef(%d,%d,%d,%d,%d,%d)\n", lambda,l1,l2,l1p,l2p,inp_.L);
                if (not std::isfinite(f2[lambda]))
                    HexException("Invalid result of computef(%d,%d,%d,%d,%d,%d)\n", lambda,l1,l2,l2p,l1p,inp_.L);
            }
            
            // for all eigenstates
            for (int i = 0; i < Nspline_atom; i++)
            {
                // get "pricipal quantum number"
                int n = l1p + i + 1;
                
                // for all multipoles
                for (int lambda = 0; lambda <= ang_.maxlambda(); lambda++)
                {
                    // hydrogen orbital
                    if (n == ni)
                    {
                        if (l1p == li and f1[lambda] != 0)
                        {
                            Complex prefactor = 4 * special::constant::pi / ki[0] * special::pow_int(1.0_i, l2p) * std::sqrt((2*l2p + 1) / (4 * special::constant::pi)) / special::constant::sqrt_two;
                            cArray Ji_expansion = prefactor * ji_expansion_atom.slice(l2p * Nspline_atom, (l2p + 1) * Nspline_atom);
                            cArray Pj1 = outer_product(Pi_expansion_atom, Ji_expansion);
                            
                            double CG = special::ClebschGordan(li, mi, l2p, 0, ang_.L(), mi);
                            
                            if (cmd_.lightweight_radial_cache)
                                rad_.apply_R_matrix(lambda, CG * f1[lambda], Pj1, 1., chi_block);
                            else
                                rad_.R_tr_dia(lambda).dot(CG * f1[lambda], Pj1, 1., chi_block, true);
                            
                            if (lambda == 0)
                                kron_dot(1., chi_block, -CG, Pj1, rad_.S_atom(), rad_.Mm1_tr_proj());
                            
                            if (lambda == 1 and cmd_.polarization > 0 and f1[lambda] != 0)
                                kron_dot(1., chi_block, -CG * f1[lambda], Pj1, rad_.M1_atom(), rad_.Xi());
                        }
                        if (l2p == li and f2[lambda] != 0)
                        {
                            Complex prefactor = 4 * special::constant::pi / ki[0] * special::pow_int(1.0_i, l1p) * std::sqrt((2*l1p + 1) / (4 * special::constant::pi)) / special::constant::sqrt_two;
                            cArray Ji_expansion = prefactor * ji_expansion_atom.slice(l1p * Nspline_atom, (l1p + 1) * Nspline_atom);
                            cArray Pj2 = outer_product(Ji_expansion, Pi_expansion_atom);
                            
                            double CG = special::ClebschGordan(l1p, 0, li, mi, ang_.L(), mi);
                            
                            if (cmd_.lightweight_radial_cache)
                                rad_.apply_R_matrix(lambda, Sign * CG * f2[lambda], Pj2, 1., chi_block);
                            else
                                rad_.R_tr_dia(lambda).dot(Sign * CG * f2[lambda], Pj2, 1., chi_block, true);
                            
                            if (lambda == 0)
                                kron_dot(1., chi_block, -CG * Sign, Pj2, rad_.Mm1_tr_proj(), rad_.S_atom());
                            
                            if (lambda == 1 and cmd_.polarization > 0 and f2[lambda] != 0)
                                kron_dot(1., chi_block, -CG * f2[lambda] * Sign, Pj2, rad_.Xi(), rad_.M1_atom());
                        }
                    }
                    
                    // hydrogen orbital polarization
                    else if (cmd_.polarization > 0)
                    {
                        cArray BState1 = rad_.getstate(n, l1p, rad_.eigenstates(l1p));
                        Complex r_elem1 = (Pi_expansion_atom.conj() | rad_.M1_atom() | BState1);
                        double En1 = rad_.eigenenergy(n, l1p).real();
                        BState1 *=  r_elem1 / (En1 - Eni);
                        for (int ell : std::vector<int>{ l2p - 1, l2p + 1 }) if (ell >= 0 and f1[lambda] != 0)
                        {
                            Complex prefactor = std::sqrt(2 * special::constant::pi * (2*ell + 1)) / ki[0] * special::pow_int(1.0_i, ell) * special::computef(lambda, l1p, l2p, li, ell, ang_.L());
                            if (prefactor == 0.)
                                continue;
                            
                            cArray Ji_expansion = ji_expansion_polar.slice(ell * Nspline_atom, (ell + 1) * Nspline_atom);
                            cArray Pj1 = outer_product(-BState1, prefactor * Ji_expansion);
                            
                            double CG = special::ClebschGordan(li, mi, ell, 0, ang_.L(), mi);
                            
                            if (cmd_.lightweight_radial_cache)
                                rad_.apply_R_matrix(lambda, CG * f1[lambda], Pj1, 1., chi_block);
                            else
                                rad_.R_tr_dia(lambda).dot(CG * f1[lambda], Pj1, 1., chi_block, true);
                            
                            if (lambda == 0)
                                kron_dot(1., chi_block, -CG, Pj1, rad_.S_atom(), rad_.Mm1_tr_proj());
                            
                            if (lambda == 1 and f1[lambda] != 0)
                                kron_dot(1., chi_block, -CG * f1[lambda], Pj1, rad_.M1_atom(), rad_.Xi());
                        }
                        
                        cArray BState2 = rad_.getstate(n, l2p, rad_.eigenstates(l2p));
                        Complex r_elem2 = (Pi_expansion_atom.conj() | rad_.M1_atom() | BState2);
                        double En2 = rad_.eigenenergy(n, l2p).real();
                        BState2 *=  r_elem2 / (En2 - Eni);
                        for (int ell : std::vector<int>{ l1p - 1, l1p + 1 }) if (ell >= 0 and f2[lambda] != 0)
                        {
                            Complex prefactor = std::sqrt(2 * special::constant::pi * (2*ell + 1)) / ki[0] * special::pow_int(1.0_i, ell) * special::computef(lambda, l1p, l2p, ell, li, ang_.L());
                            if (prefactor == 0.)
                                continue;
                            
                            cArray Ji_expansion = ji_expansion_polar.slice(ell * Nspline_atom, (ell + 1) * Nspline_atom);
                            cArray Pj2 = outer_product(prefactor * Ji_expansion, -BState2);
                            
                            double CG = special::ClebschGordan(ell, 0, li, mi, ang_.L(), mi);
                            
                            if (cmd_.lightweight_radial_cache)
                                rad_.apply_R_matrix(lambda, CG * f2[lambda] * Sign, Pj2, 1., chi_block);
                            else
                                rad_.R_tr_dia(lambda).dot(CG * f2[lambda] * Sign, Pj2, 1., chi_block, true);
                            
                            if (lambda == 0)
                                kron_dot(1., chi_block, -CG * Sign, Pj2, rad_.Mm1_tr_proj(), rad_.S_atom());
                            
                            if (lambda == 1 and f2[lambda] != 0)
                                kron_dot(1., chi_block, -CG * f2[lambda] * Sign, Pj2, rad_.Xi(), rad_.M1_atom());
                        }
                    }
                }
            }
        }
        
        // update the right-hand side
        chi[ill] = chi_block;
        chi[ill].hdfsave(format("chi-%d-%d.hdf", ang_.S(), ill));
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
    unsigned order = rad_.bspline_atom().order();
    unsigned Nspline_atom = rad_.bspline_atom().Nspline();
    unsigned Nspline_proj = rad_.bspline_proj().Nspline();
    unsigned Nang = ang_.states().size();
    unsigned Nchunk = Nspline_atom * Nspline_proj;
    
    cArray chunk (Nchunk);
    
    if (not cmd_.lightweight_radial_cache)
    {
        //
        // Simple multiplication by the super-matrix.
        //
        
        // multiply "p" by the diagonal blocks
        for (unsigned ill = 0;  ill < Nang;  ill++)
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
            dia_blocks_[ill].dot(1., p[ill], 0., q[ill], true);
            
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
        for (int lambda = 1; lambda <= rad_.maxlambda(); lambda++)
        {
            // load off-diagonal block from scratch disk
            if (not cmd_.cache_own_radint and cmd_.wholematrix)
                const_cast<BlockSymBandMatrix<Complex>&>(rad_.R_tr_dia(lambda)).hdfload();
            
            // for all source vector segments
            for (unsigned illp = 0;  illp < Nang;  illp++)
            {
                // load source segment, if needed
                if (cmd_.outofcore)
                    const_cast<BlockArray<Complex>&>(p).hdfload(illp);
                
                // multiply by the off-diagonal block
                rad_.R_tr_dia(lambda).dot(1., p[illp], 0., chunk, true);
                
                // update all destination vector segments
                for (unsigned ill = 0; ill < Nang; ill++)
                {
                    // skip diagonal and non-contributing transfers
                    if (ill == illp or ang_.f(ill,illp,lambda) == 0.)
                        continue;
                    
                    // load destination segment, if needed
                    if (cmd_.outofcore)
                        q.hdfload(ill);
                    
                    // update destination segment
                    q[ill] -= ang_.f(ill,illp,lambda) * chunk;
                    
                    // save and unload destination segment, if needed
                    if (cmd_.outofcore)
                        q.hdfsave(ill), q[ill].drop();
                }
                
                // uload the source segment
                if (cmd_.outofcore)
                    const_cast<BlockArray<Complex>&>(p)[illp].drop();
            }
            
            // release radial integrals
            if (not cmd_.cache_own_radint)
                const_cast<BlockSymBandMatrix<Complex>&>(rad_.R_tr_dia(lambda)).drop();
        }
        
        // multiply "p" by the off-diagonal blocks (multi-process)
        if (par_.Nproc() > 1)
        for (int lambda = 1; lambda <= rad_.maxlambda(); lambda++) 
        {
            // Every group of MPI processes cares for several angular momentum segments of the vectors 'p' and 'q'.
            // There are at most N segments per group, where N = ceil(A / G). Here A is the total number of blocks
            // and G is the number of groups.
            unsigned N = (ang_.states().size() + par_.Ngroup() - 1) / par_.Ngroup();
            
            // There are M processes per group, so we can process all segments in P = ceil(N / M) passes.
            unsigned P = (N + par_.groupsize() - 1) / par_.groupsize();
            
            // Working array.
            cArray work (Nchunk);
            
            // for all source segments
            for (unsigned pass = 0; pass < P; pass++)
            {
                // get angular momentum state index
                unsigned illp = (pass * par_.groupsize() + par_.igroupproc()) * par_.Ngroup() + par_.igroup();
                
                // load the source segment, if needed
                if (illp < ang_.states().size() and cmd_.outofcore and par_.isMyGroupWork(illp))
                    const_cast<BlockArray<Complex>&>(p).hdfload(illp);
                
                // calculate product, if applicable
                if (illp < ang_.states().size())
                    rad_.R_tr_dia(lambda).dot(1., p[illp], 0., work, true);
                
                // for all destination segments
                for (unsigned ill = 0; ill < Nang; ill++)
                {
                    // calculate this process' contribution
                    if (illp >= ang_.states().size() or ill == illp)
                        chunk.fill(0.);
                    else
                        chunk = -ang_.f(ill,illp,lambda) * work;
                    
                    // sum results of all group's processes to the group master proccess
                    par_.sum_g(chunk.data(), chunk.size(), 0);
                    
                    // sum all masters' segments to the master of the group that owns the destination segment
                    par_.mastersum(chunk.data(), chunk.size(), ill % par_.Ngroup());
                    
                    // redistribute the summed segment over the whole owning group
                    par_.bcast_g(ill % par_.Ngroup(), 0, &chunk[0], chunk.size());
                    
                    // finally, owning group will update the segment
                    if (par_.isMyGroupWork(ill))
                    {
                        // load destination segment, if needed
                        if (cmd_.outofcore)
                            q.hdfload(ill);
                        
                        // update the destination vector segment
                        q[ill] += chunk;
                        
                        // write to disk, if needed
                        if (cmd_.outofcore)
                        {
                            // only the group master process will write the chunk in case of shared scratch
                            if (not cmd_.shared_scratch or par_.IamGroupMaster())
                                q.hdfsave(ill);
                            
                            // unload from memory
                            q[ill].drop();
                        }
                    }
                }
                
                // unload the source segment
                if (illp < ang_.states().size() and cmd_.outofcore and par_.isMyGroupWork(illp))
                    const_cast<BlockArray<Complex>&>(p)[illp].drop();
            }
        }
    }
    else
    {
        //
        // Lightweight mode multiplication
        //
        
        // multiply "p" by the diagonal super-blocks
        for (unsigned ill = 0;  ill < Nang;  ill++) if (par_.isMyGroupWork(ill))
        {
            if (cmd_.outofcore)
                const_cast<BlockArray<Complex>&>(p).hdfload(ill), q.hdfload(ill);
            
            // get block angular momemnta
            int l1 = ang_.states()[ill].first;
            int l2 = ang_.states()[ill].second;
            
            // multiply 'p_block' by the diagonal block (except for the two-electron term)
            kron_dot(0., q[ill],  1., p[ill], Complex(E_) * rad_.S_atom(), rad_.S_proj());
            kron_dot(1., q[ill], -1., p[ill], Complex(0.5) * rad_.D_atom() - rad_.Mm1_tr_atom() + Complex(0.5*l1*(l1+1)) * rad_.Mm2_atom(), rad_.S_proj());
            kron_dot(1., q[ill], -1., p[ill], rad_.S_atom(), Complex(0.5) * rad_.D_proj() - rad_.Mm1_tr_proj() + Complex(0.5*l2*(l2+1)) * rad_.Mm2_proj());
            
            if (cmd_.outofcore)
                const_cast<BlockArray<Complex>&>(p)[ill].drop(), q.hdfsave(ill), q[ill].drop();
        }
        
        // auxiliary buffers
        cArray buffer (Nang * Nspline_proj), products (Nang * Nspline_proj), updates (Nang * Nspline_proj);
        
        // for all source vector sub-segments
        for (unsigned k = 0; k < Nspline_atom; k++)
        {
            // copy owned sub-segments to the buffer
            for (unsigned illp = 0; illp < Nang; illp++) if (par_.isMyGroupWork(illp) and par_.IamGroupMaster())
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
            
            // for all destination sub-segments
            for (unsigned i = (k > order ? k - order : 0); i <= k + order and i < Nspline_atom; i++)
            {
                // for all potential multipoles
                for (int lambda = 0; lambda <= rad_.maxlambda(); lambda++)
                {
                    // calculate the radial sub-block
                    SymBandMatrix<Complex> R_block_ik = rad_.calc_R_tr_dia_block(lambda, i, k);
                    
                    // calculate all sub-products of R * p
                    # pragma omp parallel for
                    for (unsigned illp = 0; illp < Nang; illp++)
                        cArrayView(products, illp * Nspline_proj, Nspline_proj) = R_block_ik.dot(cArrayView(buffer, illp * Nspline_proj, Nspline_proj));
                    
                    // clear updates
                    std::memset(&updates[0], 0, updates.size() * sizeof(Complex));
                    
                    // sum with angular integrals -> f * R * p
                    # pragma omp parallel for schedule (dynamic, 1) // !! no collapse !!
                    for (unsigned ill = 0; ill < Nang; ill++)
                    for (unsigned illp = 0; illp < Nang; illp++)
                    if (ang_.f(ill,illp,lambda) != 0.)
                        blas::xpby(cArrayView(updates, ill * Nspline_proj, Nspline_proj), ang_.f(ill,illp,lambda), cArrayView(products, illp * Nspline_proj, Nspline_proj)); // (not threaded)
                    
                    // update the owned sub-segment
                    # pragma omp parallel for
                    for (unsigned ill = 0; ill < Nang; ill++)
                    if (par_.isMyGroupWork(ill))
                    {
                        // get i-th section of ill-th segment (load it from disk, if needed)
                        TmpNumberArray<Complex> q_section = q.segment(ill, i * Nspline_proj, Nspline_proj);
                        
                        // update section
                        cArray section = std::move(q_section.view() - cArrayView(updates, ill * Nspline_proj, Nspline_proj));
                        
                        // store section (save to disk, if needed)
                        q.setSegment(ill, i * Nspline_proj, Nspline_proj, section);
                    }
                }
            } // for all destination sub-segments
        } // for all source vector sub-segments
    } // if not cmd_.lightweight_radial_cache
}

void NoPreconditioner::precondition (const BlockArray< Complex >& r, BlockArray< Complex >& z) const
{
    z = r;
}

void NoPreconditioner::finish ()
{
    dia_blocks_.resize(0);
}
