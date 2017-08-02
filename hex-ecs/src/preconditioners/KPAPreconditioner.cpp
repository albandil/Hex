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
#include "hex-blas.h"
#include "hex-hydrogen.h"
#include "hex-misc.h"

// --------------------------------------------------------------------------------- //

#include "gauss.h"
#include "radial.h"

// --------------------------------------------------------------------------------- //

#include "KPAPreconditioner.h"

// --------------------------------------------------------------------------------- //

std::string KPACGPreconditioner::description () const
{
    return "Block inversion using conjugate gradients preconditioned by Kronecker product approximation.";
}

void KPACGPreconditioner::setup ()
{
    NoPreconditioner::setup();
    
    refcount_atom_.resize(Hl_[0].size());
    refcount_proj_.resize(Hl_[1].size());
    
    refcount_atom_.fill(0);
    refcount_proj_.fill(0);
    
    if (verbose_)
        std::cout << "Set up KPA preconditioner" << std::endl << std::endl;
    
    // get maximal number of threads that will run the preconditioning routines concurrently
    unsigned n = 1;
#ifdef _OPENMP
    if (cmd_->parallel_precondition)
        n = omp_get_max_threads();
#endif
    
    // allocate workspaces
    workspace_.resize(n);
}

void KPACGPreconditioner::CG_init (int iblock) const
{
    // initialize parent
    CGPreconditioner::CG_init(iblock);
    
    // get block angular momenta
    int l1 = ang_->states()[iblock].first;
    int l2 = ang_->states()[iblock].second;
    
    // initialize self
    lock_kpa_access();
    {
        // update reference count
        refcount_atom_[l1]++;
        refcount_proj_[l2]++;
        
        // load preconditioners from disk if needed
        if (refcount_atom_[l1] == 1 and not Hl_[0][l1].hdfload())
            HexException("Failed to read preconditioner matrix for l = %d.", l1);
        if (refcount_proj_[l2] == 1 and not Hl_[1][l2].hdfload())
            HexException("Failed to read preconditioner matrix for l = %d.", l2);
    }
    unlock_kpa_access();
    
    // dimension of the matrices
    std::size_t Nspline_inner_x = rad_inner_->bspline_x().Nspline();
    std::size_t Nspline_inner_y = rad_inner_->bspline_y().Nspline();
    std::size_t size = Nspline_inner_x * Nspline_inner_y;
    
    // allocate thread workspace
    unsigned ithread = 0;
#ifdef _OPENMP
    if (cmd_->parallel_precondition)
        ithread = omp_get_thread_num();
#endif
    workspace_[ithread].resize(size);
}

void KPACGPreconditioner::CG_prec (int iblock, const cArrayView r, cArrayView z) const
{
    // get angular momenta of this block
    int l1 = ang_->states()[iblock].first;
    int l2 = ang_->states()[iblock].second;
    
    // get number of right-hand sides (initial states)
    unsigned Nini = r.size() / block_rank_[iblock];
    
    // get the right radial integrals (interior of the full domain, or full panel)
    RadialIntegrals const * rint;
    if (rad_panel_->bspline_x().hash() == rad_full_->bspline_x().hash() and
        rad_panel_->bspline_y().hash() == rad_full_->bspline_y().hash())
        rint = rad_inner_;
    else
        rint = rad_panel_;
    
    // dimension of the matrices
    std::size_t Nspline_inner_x = rint->bspline_x().Nspline();
    std::size_t Nspline_inner_y = rint->bspline_y().Nspline();
    
    // get workspace
    int ithread = 0;
#ifdef _OPENMP
    if (cmd_->parallel_precondition)
        ithread = omp_get_thread_num();
#endif
    
    // for all initial states
    for (unsigned ini = 0; ini < Nini; ini++)
    {
        // solution offset
        std::size_t offset = ini * block_rank_[iblock];
        
        // encapsulated memory chunks
        ColMatrixView<Complex> U (Nspline_inner_x, Nspline_inner_y, workspace_[ithread]);
        RowMatrixView<Complex> R (Nspline_inner_x, Nspline_inner_y, cArrayView(r, offset, block_rank_[iblock]));
        RowMatrixView<Complex> Z (Nspline_inner_x, Nspline_inner_y, cArrayView(z, offset, block_rank_[iblock]));
        
        // first Kronecker product
        {
            // U = (AV)B
            RowMatrixView<Complex> A (Nspline_inner_x, Nspline_inner_x, Hl_[0][l1].invCl_invsqrtS.data());
            ColMatrixView<Complex> B (Nspline_inner_y, Nspline_inner_y, Hl_[1][l2].invCl_invsqrtS.data());
            
            blas::gemm(1., A, R, 0., Z); // N³ operations
            blas::gemm(1., Z, B, 0., U); // N³ operations
        }
        
        // Divide elements by the diagonal; N² operations
        # pragma omp parallel for
        for (std::size_t i = 0; i < Nspline_inner_x; i++) 
        for (std::size_t j = 0; j < Nspline_inner_y; j++)
            Z(i,j) = U(i,j) / (E_ - Hl_[0][l1].Dl[i] - Hl_[1][l2].Dl[j]);
        
        // second Kronecker product
        {
            // W = (AU)B
            RowMatrixView<Complex> A (Nspline_inner_x, Nspline_inner_x, Hl_[0][l1].invsqrtS_Cl.data());
            ColMatrixView<Complex> B (Nspline_inner_y, Nspline_inner_y, Hl_[1][l2].invsqrtS_Cl.data());
            
            blas::gemm(1., A, Z, 0., U);    // N³ operations
            blas::gemm(1., U, B, 0., Z);    // N³ operations
        }
    }
}

void KPACGPreconditioner::CG_exit (int iblock) const
{
    // exit self
    lock_kpa_access();
    {
        // get block angular momenta
        int l1 = ang_->states()[iblock].first;
        int l2 = ang_->states()[iblock].second;
        
        // update reference count (allow drop only in out-of-core)
        if (refcount_atom_[l1] > 1 or not cmd_->outofcore)
            refcount_atom_[l1]--;
        if (refcount_proj_[l2] > 1 or not cmd_->outofcore)
            refcount_proj_[l2]--;
        
        // release memory
        if (refcount_atom_[l1] == 0)
            Hl_[0][l1].drop();
        if (refcount_proj_[l2] == 0)
            Hl_[1][l2].drop();
    }
    unlock_kpa_access();
    
    // exit parent
    CGPreconditioner::CG_exit(iblock);
}

void KPACGPreconditioner::finish ()
{
    CGPreconditioner::finish();
}

void KPACGPreconditioner::lock_kpa_access () const
{
#ifdef _OPENMP
    omp_set_lock(&lck_);
#endif
}

void KPACGPreconditioner::unlock_kpa_access () const
{
#ifdef _OPENMP
    omp_unset_lock(&lck_);
#endif    
}

// --------------------------------------------------------------------------------- //

addClassToParentRunTimeSelectionTable(PreconditionerBase, KPACGPreconditioner)

// --------------------------------------------------------------------------------- //
