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
    
    std::cout << "Set up KPA preconditioner" << std::endl << std::endl;
    
    // get maximal number of threads that will run the preconditioning routines concurrently
    unsigned n = 1;
#ifdef _OPENMP
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
    
    // calculate workspace size
    std::size_t size = rad_->bspline_inner().Nspline() * rad_->bspline_inner().Nspline();
    
    // allocate thread workspace
    unsigned ithread = 0;
#ifdef _OPENMP
    ithread = omp_get_thread_num();
#endif
    workspace_[ithread].resize(size);
}

void KPACGPreconditioner::CG_mmul (int iblock, const cArrayView p, cArrayView q) const
{
    // multiply by components
    if (cmd_->kpa_simple_rad or cmd_->lightweight_full)
    {
        // get block angular momemnta
        int l1 = ang_->states()[iblock].first;
        int l2 = ang_->states()[iblock].second;
        
        // multiply 'p' by the diagonal block (except for the two-electron term)
        kron_dot(0., q,  1., p, Complex(E_) * rad_->S_inner(), rad_->S_inner());
        kron_dot(1., q, -1., p, Complex(0.5) * rad_->D_inner() - rad_->Mm1_tr_inner() + Complex(0.5*(l1+1)*l1) * rad_->Mm2_inner(), rad_->S_inner());
        kron_dot(1., q, -1., p, rad_->S_inner(), Complex(0.5) * rad_->D_inner() + Complex(inp_->Zp) * rad_->Mm1_tr_inner() + Complex(0.5*(l2+1)*l2) * rad_->Mm2_inner());
        
        // multiply 'p' by the two-electron integrals
        for (int lambda = 0; lambda <= rad_->maxlambda(); lambda++)
        {
            // calculate angular integral
            Real f = special::computef(lambda, l1, l2, l1, l2, inp_->L);
            if (not std::isfinite(f))
                HexException("Invalid result of computef(%d,%d,%d,%d,%d,%d).", lambda, l1, l2, l1, l2, inp_->L);
            
            // multiply
            if (f != 0.)
                rad_->apply_R_matrix(lambda, inp_->Zp * f, p, 1., q, cmd_->kpa_simple_rad);
        }
    }
    
    // or let the parent do the multiplication
    else
    {
        CGPreconditioner::CG_mmul(iblock, p, q);
    }
}

void KPACGPreconditioner::CG_prec (int iblock, const cArrayView r, cArrayView z) const
{
    // get angular momenta of this block
    int l1 = ang_->states()[iblock].first;
    int l2 = ang_->states()[iblock].second;
    
    // dimension of the matrices
    std::size_t Nspline_inner = rad_->bspline_inner().Nspline();
    
    // get workspace
    int ithread = 0;
#ifdef _OPENMP
    ithread = omp_get_thread_num();
#endif
    
    // encapsulated memory chunks
    ColMatrix<Complex> U (Nspline_inner, Nspline_inner, workspace_[ithread]);
    RowMatrixView<Complex> R (Nspline_inner, Nspline_inner, r);
    RowMatrixView<Complex> Z (Nspline_inner, Nspline_inner, z);
    
    // first Kronecker product
    {
        // U = (AV)B
        RowMatrixView<Complex> A (Nspline_inner, Nspline_inner, Hl_[0][l1].invCl_invsqrtS.data());
        ColMatrixView<Complex> B (Nspline_inner, Nspline_inner, Hl_[1][l2].invCl_invsqrtS.data());
        
        blas::gemm(1., A, R, 0., Z); // N³ operations
        blas::gemm(1., Z, B, 0., U); // N³ operations
    }
    
    // Divide elements by the diagonal; N² operations
    # pragma omp parallel for
    for (std::size_t i = 0; i < Nspline_inner; i++) 
    for (std::size_t j = 0; j < Nspline_inner; j++)
        Z(i,j) = U(i,j) / (E_ - Hl_[0][l1].Dl[i] - Hl_[1][l2].Dl[j]);
    
    // second Kronecker product
    {
        // W = (AU)B
        RowMatrixView<Complex> A (Nspline_inner, Nspline_inner, Hl_[0][l1].invsqrtS_Cl.data());
        ColMatrixView<Complex> B (Nspline_inner, Nspline_inner, Hl_[1][l2].invsqrtS_Cl.data());
        
        blas::gemm(1., A, Z, 0., U);    // N³ operations
        blas::gemm(1., U, B, 0., Z);    // N³ operations
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

void KPACGPreconditioner::CG_constrain (cArrayView r) const
{
    // nothing
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
