//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2018, Jakub Benda, Charles University in Prague                    //
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
#include <cstdio>

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"
#include "hex-itersolve.h"
#include "hex-misc.h"
#include "hex-vtkfile.h"

// --------------------------------------------------------------------------------- //

#include "CGPreconditioner.h"
#include "HybPreconditioner.h"

// --------------------------------------------------------------------------------- //

std::string CGPreconditioner::description () const
{
    return "Block inversion using plain conjugate gradients. Use --tolerance option to set the termination tolerance.";
}

void CGPreconditioner::update (Real E)
{
    // update parent
    NoPreconditioner::update(E);
    
    // in arrowhead mode calculate LU decompositions of the diagonal B blocks
    if (cmd_->arrowhead)
    {
        luB1_.resize(ang_->states().size());
        luB2_.resize(ang_->states().size());
        
        // for all angular states
        for (unsigned ill = 0; ill < ang_->states().size(); ill++)
        {
            // get number of asymptotic channels
            int Nchan1 = Nchan_[ill].first;
            int Nchan2 = Nchan_[ill].second;
            
            // reserve space for LU decompositions
            luB1_[ill].resize(Nchan1);
            luB2_[ill].resize(Nchan2);
            
            // for all asymptotic channels of the first electron
            for (int ichan1 = 0; ichan1 < Nchan1; ichan1++)
            {
                luB1_[ill][ichan1].reset(LUft::Choose("lapack"));
                luB1_[ill][ichan1]->factorize(B1_blocks_[ill][ichan1 * Nchan1 + ichan1].tocoo<LU_int_t>().tocsr());
            }
            
            // for all asymptotic channels of the second electron
            for (int ichan2 = 0; ichan2 < Nchan2; ichan2++)
            {
                luB2_[ill][ichan2].reset(LUft::Choose("lapack"));
                luB2_[ill][ichan2]->factorize(B2_blocks_[ill][ichan2 * Nchan2 + ichan2].tocoo<LU_int_t>().tocsr());
            }
        }
    }
}

int CGPreconditioner::solve_channels (int ill, const cArrayView r, cArrayView z) const
{
    std::cout << std::endl << "   CGPreconditioner::solve_channels" << std::endl;
    std::cout << "    - |r| = " << r.norm() << std::endl;
    std::cout << "    - |z| = " << z.norm() << std::endl;
    
    int Nchan1 = Nchan_[ill].first;
    int Nchan2 = Nchan_[ill].second;
    
    std::cout << "    - Nchan1 = " << Nchan1 << std::endl;
    std::cout << "    - Nchan2 = " << Nchan2 << std::endl;
    
    std::size_t Nspline_full_x  = rad_panel_->bspline_x().Nspline();
    std::size_t Nspline_inner_x = rad_panel_->bspline_x().hash() == rad_full_->bspline().hash() ? rad_inner_->bspline().Nspline() : rad_panel_->bspline_x().Nspline();
    std::size_t Nspline_outer_x = Nspline_full_x - Nspline_inner_x;
    
    std::size_t Nspline_full_y  = rad_panel_->bspline_y().Nspline();
    std::size_t Nspline_inner_y = rad_panel_->bspline_y().hash() == rad_full_->bspline().hash() ? rad_inner_->bspline().Nspline() : rad_panel_->bspline_y().Nspline();
    std::size_t Nspline_outer_y = Nspline_full_y - Nspline_inner_y;
    
    // solve using the CG solver
    ConjugateGradients < Complex, cArray, cArrayView > CG;
    CG.reset();
    CG.verbose              = true; // <- temporary
    CG.scalar_product       = [ ](const cArrayView a, const cArrayView b) { return ( a | b ); };
    CG.compute_norm         = [ ](const cArrayView a) { return a.norm(); };
    CG.axby                 = [ ](Complex a, cArrayView x, Complex b, const cArrayView y) { x = a * x + b * y; };
    CG.new_array            = [ ](std::size_t n, std::string name) { return cArray(n); };
    CG.apply_preconditioner = [&](const cArrayView a, cArrayView b)
    {
        for (int ichan1 = 0; ichan1 < Nchan1; ichan1++)
        {
            cArrayView src (a, ichan1 * Nspline_outer_x, Nspline_outer_x);
            cArrayView dst (b, ichan1 * Nspline_outer_x, Nspline_outer_x);
            luB1_[ill][ichan1]->solve(src, dst, 1);
        }
        for (int ichan2 = 0; ichan2 < Nchan2; ichan2++)
        {
            cArrayView src (a, Nchan1 * Nspline_outer_x + ichan2 * Nspline_outer_y, Nspline_outer_y);
            cArrayView dst (b, Nchan1 * Nspline_outer_x + ichan2 * Nspline_outer_y, Nspline_outer_y);
            luB2_[ill][ichan2]->solve(src, dst, 1);
        }
        
        std::cout << "      prec " << a.norm() << " " << b.norm() << std::endl;
    };
    CG.matrix_multiply      = [&](const cArrayView a, cArrayView b)
    {
        for (int ichan1 = 0; ichan1 < Nchan1; ichan1 ++)
        for (int ichan1p= 0; ichan1p< Nchan1; ichan1p++)
        {
            cArrayView src (a, ichan1p * Nspline_outer_x, Nspline_outer_x);
            cArrayView dst (b, ichan1  * Nspline_outer_x, Nspline_outer_x);
            B1_blocks_[ill][ichan1 * Nchan1 + ichan1p].dot(1., src, 1., dst);
        }
        for (int ichan2 = 0; ichan2 < Nchan2; ichan2 ++)
        for (int ichan2p= 0; ichan2p< Nchan2; ichan2p++)
        {
            cArrayView src (a, Nchan1 * Nspline_outer_x + ichan2p * Nspline_outer_y, Nspline_outer_y);
            cArrayView dst (b, Nchan1 * Nspline_outer_x + ichan2  * Nspline_outer_y, Nspline_outer_y);
            B1_blocks_[ill][ichan2 * Nchan2 + ichan2p].dot(1., src, 1., dst);
        }
        
        std::cout << "      mmul " << a.norm() << " " << b.norm() << std::endl;
    };
    
    int n = CG.solve(r, z, cmd_->prec_itertol, 0, 100);
    
    std::cout << "    - |z| = " << z.norm() << std::endl;
    std::cout << "    - n = " << n << std::endl;
    
    return n;
}

int CGPreconditioner::solve_block (int ill, const cArrayView r, cArrayView z) const
{
    std::cout << std::endl << "CGPreconditioner::solve_block" << std::endl;
    
    // shorthands
    int Nspline_inner_x = rad_panel_->bspline_x().hash() == rad_full_->bspline().hash() ? rad_inner_->bspline().Nspline() : rad_panel_->bspline_x().Nspline();
    int Nspline_inner_y = rad_panel_->bspline_y().hash() == rad_full_->bspline().hash() ? rad_inner_->bspline().Nspline() : rad_panel_->bspline_y().Nspline();
    
    // prepare the block-preconditioner for run
    this->CG_init(ill);
    
    // maximal number of nested iterations, never more than 10 thousand (also prevents the int to overflow)
    int max_iterations = cmd_->max_sub_iter > 0 ? cmd_->max_sub_iter : std::min<std::size_t>(10000, std::size_t(Nspline_inner_x) * std::size_t(Nspline_inner_y));
    
    // adjust max iterations for ILU-preconditioned blocks
    if (HybCGPreconditioner const * hp = dynamic_cast<HybCGPreconditioner const*>(this))
    {
        if (hp->ilu_needed(ill) and cmd_->ilu_max_iter > 0)
            max_iterations = cmd_->ilu_max_iter;
    }
    else if (dynamic_cast<ILUCGPreconditioner const*>(this) != nullptr)
    {
        if (cmd_->ilu_max_iter > 0)
            max_iterations = cmd_->ilu_max_iter;
    }
    
    // views to the source and solution; FIXME : combination arrowhead + multiple RHSs not implemented
    cArrayView rview (r);
    cArrayView zview (z);
    
    // views of the inner region part of sources and solution
    cArrayView inner_sources   (r, 0, Nspline_inner_x * Nspline_inner_y);
    cArrayView inner_solutions (z, 0, Nspline_inner_x * Nspline_inner_y);
    
    // views of the outer region part of sources and solution
    cArrayView channel_sources   (r, Nspline_inner_x * Nspline_inner_y, r.size() - Nspline_inner_x * Nspline_inner_y);
    cArrayView channel_solutions (z, Nspline_inner_x * Nspline_inner_y, r.size() - Nspline_inner_x * Nspline_inner_y);
    
    // auxiliary arrays used by arrowhead decomposition
    cArray modified_inner_sources, modified_channel_sources;
    
    std::cout << " - #r = " << r.size() << std::endl;
    std::cout << " - |r| = " << r.norm() << std::endl;
    std::cout << " - #inner_sources = " << inner_sources.size() << std::endl;
    std::cout << " - |inner_sources| = " << inner_sources.norm() << std::endl;
    std::cout << " - #channel_sources = " << channel_sources.size() << std::endl;
    std::cout << " - |channel_sources| = " << channel_sources.norm() << std::endl;
    
    // modify vectors for arrowhead decomposition
    if (cmd_->arrowhead)
    {
        // In the case of arrowhead decomposition we solve only for the inner wave function
        // iteratively. The outer (channel) wave functions are solved once and for all before
        // (to correct the source) and after.
        
        
        // solve detached outer problem
        solve_channels(ill, channel_sources, channel_solutions);
        std::cout << " - |channel_solutions| = " << channel_solutions.norm() << std::endl;
        
        // correct right-hand side of the inner problem
        modified_inner_sources = r;
        Cu_blocks_[ill].dot(-1., z, 1., modified_inner_sources);
        
        std::cout << " - |modified_inner_sources| = " << modified_inner_sources.norm() << std::endl;
        
        // redirect views to inner data only (which is what preconditioners expect in arrowhead mode)
        rview.reset(Nspline_inner_x * Nspline_inner_y, modified_inner_sources.data());
        zview.reset(Nspline_inner_x * Nspline_inner_y, z.data());
    }
    
    // solve using the CG solver
    ConjugateGradients < Complex, cArray, cArrayView > CG;
    CG.reset();
    CG.verbose              = cmd_->sub_prec_verbose;
    CG.apply_preconditioner = [&](const cArrayView a, cArrayView b)
                              {
                                  par_->wait_g();
                                  Timer timer;
                                  this->CG_prec(ill, a, b);
                                  us_prec_ += timer.microseconds();
                              };
    CG.matrix_multiply      = [&](const cArrayView a, cArrayView b)
                              {
                                  par_->wait_g();
                                  Timer timer;
                                  this->CG_mmul(ill, a, b);
                                  us_mmul_ += timer.microseconds();
                              };
    CG.scalar_product       = [&](const cArrayView a, const cArrayView b)
                              {
                                  par_->wait_g();
                                  Timer timer;
                                  Complex prod = this->CG_scalar_product(a, b);
                                  us_spro_ += timer.microseconds();
                                  return prod;
                              };
    CG.compute_norm         = [&](const cArrayView a)
                              {
                                  par_->wait_g();
                                  Timer timer;
                                  Real nrm = this->CG_compute_norm(a);
                                  us_norm_ += timer.microseconds();
                                  return nrm;
                              };
    CG.axby                 = [&](Complex a, cArrayView x, Complex b, const cArrayView y)
                              {
                                  par_->wait_g();
                                  Timer timer;
                                  this->CG_axby_operation(a, x, b, y);
                                  us_axby_ += timer.microseconds();
                              };
    CG.new_array            = [&](std::size_t n, std::string name)
                              {
                                  return cArray(n);
                              };
    CG.constrain            = [&](cArrayView r)
                              {
                                  this->CG_constrain(r);
                              };
    int n = CG.solve(rview, zview, cmd_->prec_itertol, 0, max_iterations);
    
    if (n >= max_iterations)
    {
        if (cmd_->fail_on_sub_iter)
            HexException("Error: Maximal number of iterations (%d) reached in the sub-preconditioner.", max_iterations);
        else
            std::cout << "\tWarning: Maximal number of iterations (" << max_iterations << ") reached in the sub-preconditioner." << std::endl;
    }
    
    std::cout << " - |inner_solutions| = " << zview.norm() << std::endl;
    
    // finalize arrowhead solution
    if (cmd_->arrowhead)
    {
        
        // correct right-hand sides of the outer problems using the inner solution
        cArray modified_channel_sources = channel_sources;
        Cl_blocks_[ill].dot(-1., zview, 1., modified_channel_sources);
        
        // solve channels
        solve_channels(ill, modified_channel_sources, channel_solutions);
        
        std::cout << " - final |channel_solutions| = " << channel_solutions.norm() << std::endl;
    }
    
    std::cout << " - final |z| = " << z.norm() << std::endl;
    
    // release block-preconditioner block-specific data
    this->CG_exit(ill);
    
    return n;
}

void CGPreconditioner::precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const
{
    // clear timing information
    us_axby_ = us_mmul_ = us_norm_ = us_prec_ = us_spro_ = 0;
    
    // apply SSOR TODO : fix for multi-rhs
    if (cmd_->ssor > 0)
    {
        // working arrays
        cArrays y (r.size()), x (r.size());
        
        // forward SOR
        for (int ill = 0; ill < (int)ang_->states().size(); ill++)
        {
            // start with right-hand side
            y[ill] = r[ill];
            
            // subtract lower block diagonals for ill-th block row
            for (int illp = 0; illp < ill; illp++)
            for (int lambda = 0; lambda <= rad_full_->maxlambda(); lambda++)
            if (ang_->f(ill, illp, lambda) != 0)
                rad_full_->R_tr_dia(lambda).dot(-ang_->f(ill, illp, lambda), y[illp], 1., y[ill], true);
            
            // use (preconditioned) conjugate gradients to invert a diagonal block
            x[ill].resize(y[ill].size());
            n_[ill] = solve_block(ill, cmd_->ssor * y[ill], x[ill]);
        }
        
        // normalize
        for (int ill = 0; ill < (int)ang_->states().size(); ill++)
        {
            this->CG_mmul(ill, (2.0_r - cmd_->ssor) / cmd_->ssor * x[ill], y[ill]);
        }
        
        // backward SOR TODO : fix for multi-rhs
        for (int ill = (int)ang_->states().size() - 1; ill >= 0; ill--)
        {
            // subtract upper block diagonals for ill-th block row
            for (int illp = ill + 1; illp < (int)ang_->states().size(); illp++)
            for (int lambda = 0; lambda <= rad_full_->maxlambda(); lambda++)
            if (ang_->f(ill, illp, lambda) != 0)
                rad_full_->R_tr_dia(lambda).dot(-ang_->f(ill, illp, lambda), y[illp], 1., y[ill], true);
            
            // use (preconditioned) conjugate gradients to invert a diagonal block
            n_[ill] += solve_block(ill, cmd_->ssor * y[ill], z[ill]);
        }
    }
    else
    {
#ifndef DISABLE_PARALLEL_PRECONDITION
        // NOTE : If the BLAS is multi-threaded, this will result in nested parallelism.
        //        Some combinations of system kernel and OpenMP implementation lead to system
        //        freeze. It can be avoided by using a serial BLAS.
        # pragma omp parallel for schedule (dynamic, 1) if (cmd_->parallel_precondition && cmd_->groupsize == 1)
#endif
        for (int ill = 0; ill < (int)ang_->states().size(); ill++) if (par_->isMyGroupWork(ill))
        {
            // load segment, if necessary
            if (cmd_->outofcore)
            {
                const_cast<BlockArray<Complex>&>(r).hdfload(ill);
                z.hdfload(ill);
            }
            
            // invert diagonal block
            n_[ill] = solve_block(ill, r[ill], z[ill]);
            
            // unload segment
            if (cmd_->outofcore)
            {
                const_cast<BlockArray<Complex>&>(r)[ill].drop();
                
                if (not cmd_->shared_scratch or par_->IamGroupMaster())
                {
                    z.hdfsave(ill);
                    z.drop(ill);
                }
            }
        }
    }
    
    // broadcast inner preconditioner iterations
    par_->sync_m(n_.data(), 1, ang_->states().size());
    par_->bcast_g(par_->igroup(), 0, n_.data(), ang_->states().size());
    
    // inner preconditioner info (max and avg number of iterations)
    std::cout << " | ";
    std::cout << std::setw(5) << (*std::min_element(n_.begin(), n_.end()));
    std::cout << std::setw(5) << (*std::max_element(n_.begin(), n_.end()));
    std::cout << std::setw(5) << format("%g", std::accumulate(n_.begin(), n_.end(), 0) / float(n_.size()));
    
    // preconditioner timing
    std::size_t us_total = us_axby_ + us_mmul_ + us_norm_ + us_prec_ + us_spro_;
    if (us_total != 0)
    std::cout << " [prec: " << format("%2d", int(us_prec_ * 100. / us_total)) << "%"
              << ", mmul: " << format("%2d", int(us_mmul_ * 100. / us_total)) << "%"
              << ", axby: " << format("%2d", int(us_axby_ * 100. / us_total)) << "%"
              << ", norm: " << format("%2d", int(us_norm_ * 100. / us_total)) << "%"
              << ", spro: " << format("%2d", int(us_spro_ * 100. / us_total)) << "%"
              << "]";
}

void CGPreconditioner::CG_init (int iblock) const
{
    if (cmd_->lightweight_simple)
    {
        const_cast<BlockSymBandMatrix<Complex> &>(A_blocks_[iblock * ang_->states().size() + iblock])
            = calc_A_block(iblock, iblock, true);
    }
}

void CGPreconditioner::CG_mmul (int iblock, const cArrayView p, cArrayView q) const
{
    std::size_t Nspline_full_x  = rad_panel_->bspline_x().Nspline();
    std::size_t Nspline_inner_x = rad_panel_->bspline_x().hash() == rad_full_->bspline().hash() ? rad_inner_->bspline().Nspline() : rad_panel_->bspline_x().Nspline();
    std::size_t Nspline_outer_x = Nspline_full_x - Nspline_inner_x;
    
    std::size_t Nspline_full_y  = rad_panel_->bspline_y().Nspline();
    std::size_t Nspline_inner_y = rad_panel_->bspline_y().hash() == rad_full_->bspline().hash() ? rad_inner_->bspline().Nspline() : rad_panel_->bspline_y().Nspline();
    std::size_t Nspline_outer_y = Nspline_full_y - Nspline_inner_y;
    
    std::size_t Nini = p.size() / block_rank_[iblock];
    std::size_t Nang = ang_->states().size();
    std::size_t iang = iblock * Nang + iblock;
    
    // FIXME : multi-rhs in arrowhead mode
    if (cmd_->arrowhead)
        Nini = 1;
    
    q.fill(0.0_z);
    
    for (std::size_t ini = 0; ini < Nini; ini++)
    {
        std::size_t offset = block_rank_[iblock] * ini;
        
        // inner-region subset of the vectors
        cArrayView p_inner (p, offset, Nspline_inner_x * Nspline_inner_y);
        cArrayView q_inner (q, offset, Nspline_inner_x * Nspline_inner_y);
        
        if (cmd_->lightweight_full and not cmd_->lightweight_simple)
        {
            // get block angular momemnta
            int l1 = ang_->states()[iblock].first;
            int l2 = ang_->states()[iblock].second;
            
            // one-electron matrices
            SymBandMatrix<Complex> const & Sx = rad_panel_->S_x();
            SymBandMatrix<Complex> const & Sy = rad_panel_->S_y();
            SymBandMatrix<Complex> Hx = 0.5_z * rad_panel_->D_x() + (0.5_z * (l1 * (l1 + 1.0_r))) * rad_panel_->Mm2_x() + Complex(inp_->Za *   -1.0_r) * rad_panel_->Mm1_x();
            SymBandMatrix<Complex> Hy = 0.5_z * rad_panel_->D_y() + (0.5_z * (l2 * (l2 + 1.0_r))) * rad_panel_->Mm2_y() + Complex(inp_->Za * inp_->Zp) * rad_panel_->Mm1_y();
            
            // multiply 'p' by the diagonal block
            // - except for the two-electron term
            // - restrict to inner region
            kron_dot(0., q_inner, E_, p_inner, Sx, Sy, Nspline_inner_x, Nspline_inner_y);
            kron_dot(1., q_inner, -1, p_inner, Hx, Sy, Nspline_inner_x, Nspline_inner_y);
            kron_dot(1., q_inner, -1, p_inner, Sx, Hy, Nspline_inner_x, Nspline_inner_y);
            
            // multiply 'p' by the two-electron integrals
            for (int lambda = 0; lambda <= rad_inner_->maxlambda(); lambda++)
            {
                // calculate angular integral
                Real f = special::computef(lambda, l1, l2, l1, l2, inp_->L);
                if (not std::isfinite(f))
                    HexException("Invalid result of computef(%d,%d,%d,%d,%d,%d).", lambda, l1, l2, l1, l2, inp_->L);
                
                // multiply
                if (f != 0.)
                {
                    rad_panel_->apply_R_matrix(lambda, inp_->Zp * f, p_inner, 1.0, q_inner, Nspline_inner_x, Nspline_inner_y);
                }
            }
        }
        else
        {
            if (not cmd_->lightweight_simple and cmd_->outofcore and cmd_->wholematrix) const_cast<BlockSymBandMatrix<Complex> &>(A_blocks_[iang]).hdfload();
            
            A_blocks_[iang].dot
            (
                1.0_r, p_inner,
                1.0_r, q_inner,
                !cmd_->parallel_precondition
            );
            
            if (not cmd_->lightweight_simple and cmd_->outofcore and cmd_->wholematrix) const_cast<BlockSymBandMatrix<Complex> &>(A_blocks_[iang]).drop();
        }
        
        if (not inp_->inner_only and not cmd_->arrowhead)
        {
            std::size_t Nchan1 = Nchan_[iblock].first;
            std::size_t Nchan2 = Nchan_[iblock].second;
            
            # pragma omp parallel for if (!cmd_->parallel_precondition)
            for (std::size_t m = 0; m < Nchan1; m++)
            for (std::size_t n = 0; n < Nchan1; n++)
            {
                if (cmd_->outofcore) const_cast<SymBandMatrix<Complex>&>(B1_blocks_[iang][m * Nchan1 + n]).hdfload();
                B1_blocks_[iang][m * Nchan1 + n].dot
                (
                    1.0_r, cArrayView(p, offset + Nspline_inner_x * Nspline_inner_y + n * Nspline_outer_x, Nspline_outer_x),
                    1.0_r, cArrayView(q, offset + Nspline_inner_x * Nspline_inner_y + m * Nspline_outer_x, Nspline_outer_x)
                );
                if (cmd_->outofcore) const_cast<SymBandMatrix<Complex>&>(B1_blocks_[iang][m * Nchan1 + n]).drop();
            }
            
            # pragma omp parallel for if (!cmd_->parallel_precondition)
            for (std::size_t m = 0; m < Nchan2; m++)
            for (std::size_t n = 0; n < Nchan2; n++)
            {
                if (cmd_->outofcore) const_cast<SymBandMatrix<Complex>&>(B2_blocks_[iang][m * Nchan2 + n]).hdfload();
                B2_blocks_[iang][m * Nchan2 + n].dot
                (
                    1.0_r, cArrayView(p, offset + Nspline_inner_x * Nspline_inner_y + Nchan1 * Nspline_outer_x + n * Nspline_outer_y, Nspline_outer_y),
                    1.0_r, cArrayView(q, offset + Nspline_inner_x * Nspline_inner_y + Nchan1 * Nspline_outer_x + m * Nspline_outer_y, Nspline_outer_y)
                );
                if (cmd_->outofcore) const_cast<SymBandMatrix<Complex>&>(B2_blocks_[iang][m * Nchan2 + n]).drop();
            }
            
            Cu_blocks_[iang].dot(1.0_r, cArrayView(p, offset, block_rank_[iblock]), 1.0_r, cArrayView(q, offset, block_rank_[iblock]));
            Cl_blocks_[iang].dot(1.0_r, cArrayView(p, offset, block_rank_[iblock]), 1.0_r, cArrayView(q, offset, block_rank_[iblock]));
        }
        
        std::cout << "(mmul) |p| = " << p_inner.norm() << std::endl;
        std::cout << "(mmul) |q| = " << q_inner.norm() << std::endl;
        
        if (cmd_->arrowhead) // FIXME : multi-rhs in arrohead mode
        {
            cArray tmp1 (block_rank_[iblock]);
            cArray tmp2 (block_rank_[iblock]);
            
            cArrayView tmp1_inner (tmp1, 0, Nspline_inner_x * Nspline_inner_y);
            cArrayView tmp2_inner (tmp2, 0, Nspline_inner_x * Nspline_inner_y);
            
            cArrayView tmp1_outer (tmp1, Nspline_inner_x * Nspline_inner_y, tmp1.size() - Nspline_inner_x * Nspline_inner_y);
            cArrayView tmp2_outer (tmp2, Nspline_inner_x * Nspline_inner_y, tmp2.size() - Nspline_inner_x * Nspline_inner_y);
            
            std::cout << "(mmul) |Cl| = " << Cl_blocks_[iblock].v().norm() << std::endl;
            std::cout << "(mmul) |Cu| = " << Cu_blocks_[iblock].v().norm() << std::endl;
            
            if (p.norm() != 0)
            {
                VTKRectGridFile vtk;
                
                rArray grid = linspace(0., rad_inner().bspline().R2(), 500);
                cArray eval = rad_inner().bspline().zip(p, grid, grid);
                
                vtk.setGridX(grid);
                vtk.setGridY(grid);
                vtk.setGridZ({ 0. });
                
                vtk.appendScalarAttribute("RePsi", realpart(eval));
                vtk.appendScalarAttribute("ImPsi", imagpart(eval));
                
                vtk.writePoints("p.vtk");
                
                std::exit(0);
            }
            
            tmp1_inner = p;
            Cl_blocks_[iblock].dot(1., tmp1, 0., tmp2);
            std::cout << "(mmul) | Cl p | = " << tmp2.norm() << std::endl;
            
            solve_channels(iblock, tmp2_outer, tmp1_outer);
            std::cout << "(mmul) | B-1 Cl p | = " << tmp1_outer.norm() << std::endl;
            
            Cu_blocks_[iblock].dot(-1., tmp1, 1., tmp2);
            q = tmp2_inner;
            
            std::cout << "(mmul) -> |q| = " << q.norm() << std::endl;
        }
    }
}

void CGPreconditioner::CG_prec (int iblock, const cArrayView r, cArrayView z) const
{
    z = r;
}

void CGPreconditioner::CG_exit (int iblock) const
{
    if (cmd_->lightweight_simple)
    {
        const_cast<BlockSymBandMatrix<Complex> &>(A_blocks_[iblock * ang_->states().size() + iblock])
            .drop();
    }
}

void CGPreconditioner::finish ()
{
    n_.fill(-1);
    luB1_.clear();
    luB2_.clear();
    NoPreconditioner::finish();
}

Real CGPreconditioner::CG_compute_norm (const cArrayView a) const
{
    // compute norm (part will be computed by every process in group, result will be synchronized)
        
    // calculate number of elements every group's process will sum
    std::size_t N = (a.size() + par_->groupsize() - 1) / par_->groupsize();
    
    // calculate part of the norm
    Real norm2 = 0;
    for (std::size_t i = par_->igroupproc() * N; i < (par_->igroupproc() + 1) * N and i < a.size(); i++)
        norm2 += sqrabs(a[i]);
    
    // sum across group and return result
    par_->syncsum_g(&norm2, 1);
    return std::sqrt(norm2);
}

Complex CGPreconditioner::CG_scalar_product (const cArrayView a, const cArrayView b) const
{
    // compute scalar product (part will be computed by every process in group, result will be synchronized)
    
    assert(a.size() == b.size());
    
    // calculate number of elements every group's process will sum
    std::size_t N = (a.size() + par_->groupsize() - 1) / par_->groupsize();
    
    // calculate part of the norm
    Complex prod = 0;
    for (std::size_t i = par_->igroupproc() * N; i < (par_->igroupproc() + 1) * N and i < a.size(); i++)
        prod += a[i] * b[i];
    
    // sum across group and return result
    par_->syncsum_g(&prod, 1);
    return prod;
}

void CGPreconditioner::CG_axby_operation (Complex a, cArrayView x, Complex b, const cArrayView y) const
{
    // a*x + b*y operation
    
    assert(x.size() == y.size());
    
    // number of elements each node should process (rounded up)
    std::size_t N = (x.size() + par_->groupsize() - 1) / par_->groupsize();
    
    // beginning and end of the segment to be processed by this group member
    std::size_t begin = par_->igroupproc() * N;
    std::size_t end = std::min(x.size(), (par_->igroupproc() + 1) * N);
    
    // calculate the linear combination
    for (std::size_t i = begin; i < end; i++)
        x[i] = a * x[i] + b * y[i];
    
    // synchronize the segments
    for (int inode = 0; inode < par_->groupsize(); inode++)
    {
        begin = inode * N;
        end = std::min(x.size(), (inode + 1) * N);
        
        // broadcast the data to all members of the group
        par_->bcast_g
        (
            par_->igroup(),         // broadcast within this node's group
            inode,                  // owners rank within the group
            &x[0] + begin,          // data pointer
            end - begin             // number of elements to synchronize
        );
    }
}

void CGPreconditioner::CG_constrain (cArrayView r) const
{
    // leave the resudual as it is
}

// --------------------------------------------------------------------------------- //

addClassToParentRunTimeSelectionTable(PreconditionerBase, CGPreconditioner)

// --------------------------------------------------------------------------------- //
