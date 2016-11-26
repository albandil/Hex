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

void KPACGPreconditioner::sData::hdflink (const char* file)
{
    filename = file;
}

bool KPACGPreconditioner::sData::hdfcheck (const char* file) const
{
    // open HDF file for reading
    HDFFile hdf ((file == nullptr ? filename.c_str() : file), HDFFile::readonly);
    return hdf.valid();
}

bool KPACGPreconditioner::sData::hdfload (const char* file)
{
    // open HDF file for reading
    HDFFile hdf ((file == nullptr ? filename.c_str() : file), HDFFile::readonly);
    if (not hdf.valid())
        return false;
    
    // read size
    unsigned size;
    if (not hdf.read("n", &size, 1))
        return false;
    
    // read matrices
    invCl_invsqrtS = RowMatrix<Complex>(size,size);
    if (not hdf.read("invCl_invsqrtS", invCl_invsqrtS.data().data(), size * size))
        return false;
    invsqrtS_Cl = RowMatrix<Complex>(size,size);
    if (not hdf.read("invsqrtS_Cl", invsqrtS_Cl.data().data(), size * size))
        return false;
    
    // read eigenvalues
    Dl.resize(size);
    if (not hdf.read("Dl", Dl.data(), size))
    {
        std::cout << "Failed to read Dl from " << hdf.name() << ": " << hdf.error() << std::endl;
        return false;
    }
    
    return true;
}

bool KPACGPreconditioner::sData::hdfsave (const char* file) const
{
    // open HDF file for writing
    HDFFile hdf ((file == nullptr ? filename.c_str() : file), HDFFile::overwrite);
    if (not hdf.valid())
        return false;
    
    // write size
    unsigned size = Dl.size();
    if (not hdf.write("n", &size, 1))
        return false;
    
    // write matrices
    if (not hdf.write("invCl_invsqrtS", invCl_invsqrtS.data().data(), size * size))
        return false;
    if (not hdf.write("invsqrtS_Cl", invsqrtS_Cl.data().data(), size * size))
        return false;
    
    // write eigenvalues
    if (not hdf.write("Dl", Dl.data(), size))
        return false;
    
    return true;
}

void KPACGPreconditioner::sData::drop ()
{
    invCl_invsqrtS.drop();
    invsqrtS_Cl.drop();
    Dl.clear();
}

void KPACGPreconditioner::prepare
(
    std::vector<Data> & prec,
    std::size_t Nspline,
    SymBandMatrix<Complex> const & mS,
    SymBandMatrix<Complex> const & mD,
    SymBandMatrix<Complex> const & mMm1_tr,
    SymBandMatrix<Complex> const & mMm2,
    Array<bool> done,
    std::set<int> comp_l,
    std::set<int> needed_l
)
{
    Timer timer;
    cArray D;
    ColMatrix<Complex> S = mS.torow().T(), CR, invCR, invsqrtS(Nspline,Nspline);
    
    // diagonalize overlap matrix
    if (not all(done))
    {
        std::cout << "\t\t- overlap matrix factorization" << std::endl;
        
        S.diagonalize(D, nullptr, &CR);
        CR.invert(invCR);
        
        // Now S = CR * (D * CR⁻¹)
        std::cout << "\t\t\t- time: " << timer.nice_time() << std::endl;
        for (std::size_t i = 0; i < Nspline * Nspline; i++)
            invCR.data()[i] *= D[i % Nspline];
        
        // S = S - CR * invCR
        blas::gemm(-1., CR, invCR, 1., S);
        std::cout << "\t\t\t- residual: " << S.data().norm() << std::endl;
        
        // compute √S⁻¹
        for (std::size_t i = 0; i < Nspline * Nspline; i++)
            invCR.data()[i] /= std::pow(D.data()[i % Nspline], 1.5);
        blas::gemm(1., CR, invCR, 0., invsqrtS);
    }
    
    // diagonalize one-electron hamiltonians for all angular momenta
    for (int l : comp_l)
    {
        // skip loaded
        if (done[l])
            continue;
        
        // reset timer
        std::cout << "\t\t- one-electron Hamiltonian factorization (l = " << l << ")" << std::endl;
        timer.reset();
        
        // compose the symmetrical one-electron hamiltonian
        ColMatrix<Complex> tHl = (Complex(0.5) * mD - mMm1_tr + Complex(0.5*l*(l+1)) * mMm2).torow().T();
        
        // symmetrically transform by inverse square root of the overlap matrix, tHl <- invsqrtS * tHl * invsqrtS
        blas::gemm(1., invsqrtS, tHl, 0., S);
        blas::gemm(1., S, invsqrtS, 0., tHl);
        
        // diagonalize the transformed matrix
        tHl.diagonalize(D, nullptr, &CR);
        CR.invert(invCR);
        
        // analyze bound states
        int maxn = 0;
        for (unsigned i = 0; i < Nspline; i++)
        {
            Complex E = D[i];
            if (E.real() < 0)
            {
                int n = std::floor(0.5 + 1.0 / std::sqrt(-2.0 * E.real()));
                if (n > 0)
                {
                    Real E0 = -0.5 / (n * n);
                    if (std::abs(E0 - E) < 1e-3 * std::abs(E0))
                        maxn = std::max(maxn, n);
                }
            }
        }
        
        // store the preconditioner data
        prec[l].Dl = D;
        prec[l].invsqrtS_Cl = RowMatrix<Complex>(Nspline,Nspline);
        prec[l].invCl_invsqrtS = RowMatrix<Complex>(Nspline,Nspline);
        blas::gemm(1., invsqrtS, CR, 0., prec[l].invsqrtS_Cl);
        blas::gemm(1., invCR, invsqrtS, 0., prec[l].invCl_invsqrtS);
        prec[l].hdfsave();
        prec[l].drop();
        
        // Now Hl = ClR * D * ClR⁻¹
        std::cout << "\t\t\t- time: " << timer.nice_time() << std::endl;
        for (std::size_t i = 0; i < Nspline * Nspline; i++)
            invCR.data()[i] *= D[i % Nspline];
        
        // Hl <- Hl - CR * invCR
        blas::gemm(-1., CR, invCR, 1., tHl);
        std::cout << "\t\t\t- residual: " << tHl.data().norm() << std::endl;
        std::cout << "\t\t\t- bound states with energy within 0.1 % from exact value: " << l + 1 << " <= n <= " << maxn << std::endl;
    }
    
    // wait for completition of diagonalization on other nodes
    par_->wait();

    // load all preconditioner matrices needed by this MPI node
//     for (int l : needed_l)
//         if (not cmd_.outofcore and not prec[l].hdfload())
//             HexException("Failed to read preconditioner matrix for l = %d.", l);
}

void KPACGPreconditioner::setup ()
{
    NoPreconditioner::setup();
    
    std::cout << "Set up KPA preconditioner" << std::endl;
    
    //
    // Preparations (presence checking, loading, ...).
    //
    
        // compose list of angular momenta NEEDED by this MPI node
        std::set<int> needed_l;
        for (int l = 0; l <= inp_->maxell; l++)
        {
            // check if this angular momentum is needed by some of the blocks owned by this process
            bool need_this_l = false;
            for (unsigned ill = 0; ill < ang_->states().size(); ill++)
            {
                if (par_->isMyWork(ill) and (ang_->states()[ill].first == l or ang_->states()[ill].second == l))
                    need_this_l = true;
            }
            
            // append the 'l' if needed
            if (need_this_l)
                needed_l.insert(l);
        }
        
        // compose list of angular momenta that this MPI node will COMPUTE
        std::set<int> comp_l;
        for (int l = 0; l <= inp_->maxell; l++)
        {
            // compute l-th preconditioner matrices if either
            // a) nodes share scratch and even distribution assigns l-th preconditioner matrices to this node, or
            if (cmd_->shared_scratch and par_->isMyWork(l))
                comp_l.insert(l);
            // b) nodes do not share scratch and this node will need l-th preconditioner matrices
            if (not cmd_->shared_scratch and needed_l.find(l) != needed_l.end())
                comp_l.insert(l);
        }
        
    //
    // Check presence of the preconditioner data
    //
        
        // status array indicating necessity to calculate the preconditioner matrix for given 'l'
        Array<bool> done_atom (inp_->maxell + 1, true), done_proj (inp_->maxell + 1, true);
        
        // "to compute matrices": link them to scratch disk files and check presence
        for (int l : comp_l)
        {
            prec_atom_[l].hdflink(format("kpa-%d-%.4lx.hdf",l,rad_->bspline_inner_x().hash()).c_str());
            prec_proj_[l].hdflink(format("kpa-%d-%.4lx.hdf",l,rad_->bspline_inner_y().hash()).c_str());
            done_atom[l] = prec_atom_[l].hdfcheck();
            done_atom[l] = prec_proj_[l].hdfcheck();
        }
        
        // "needed matrices": link them to scratch disk files and check presence, load if present
        for (int l : needed_l)
        {
            prec_atom_[l].hdflink(format("kpa-%d-%.4lx.hdf",l,rad_->bspline_inner_x().hash()).c_str());
            prec_proj_[l].hdflink(format("kpa-%d-%.4lx.hdf",l,rad_->bspline_inner_y().hash()).c_str());
            done_atom[l] = prec_atom_[l].hdfcheck();
            done_proj[l] = prec_proj_[l].hdfcheck();
            
            if (done_atom[l])
            {
                std::cout << "\t- atomic preconditioner data for l = " << l
                           << " present in \"" << prec_atom_[l].filename << "\"" << std::endl;
            }
            
            if (done_proj[l])
            {
                std::cout << "\t- projectile preconditioner data for l = " << l
                           << " present in \"" << prec_atom_[l].filename << "\"" << std::endl;
            }
        }
    
    //
    // Calculation of the preconditioner for atomic basis.
    //
        
        if (not all(done_atom))
        {
            std::cout << std::endl << "\tPrepare preconditioner matrices for atomic grid" << std::endl;
            prepare
            (
                prec_atom_, rad_->bspline_inner_x().Nspline(),
                rad_->S_inner_x(), rad_->D_inner_x(), rad_->Mm1_tr_inner_x(), rad_->Mm2_inner_x(),
                done_atom, comp_l, needed_l
            );
            std::cout << std::endl << "\tPrepare preconditioner matrices for projectile grid" << std::endl;
            prepare
            (
                prec_proj_, rad_->bspline_inner_y().Nspline(),
                rad_->S_inner_y(), rad_->D_inner_y(), Complex(-inp_->Zp) * rad_->Mm1_tr_inner_y(), rad_->Mm2_inner_y(),
                done_proj, comp_l, needed_l
            );
        }
    
    std::cout << std::endl;
    
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
        if (refcount_atom_[l1] == 1 and not prec_atom_[l1].hdfload())
            HexException("Failed to read preconditioner matrix for l = %d.", l1);
        if (refcount_proj_[l2] == 1 and not prec_proj_[l2].hdfload())
            HexException("Failed to read preconditioner matrix for l = %d.", l2);
    }
    unlock_kpa_access();
    
    // calculate workspace size
    std::size_t size = rad_->bspline_inner_x().Nspline() * rad_->bspline_inner_y().Nspline();
    
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
        kron_dot(0., q,  1., p, Complex(E_) * rad_->S_inner_x(), rad_->S_inner_y());
        kron_dot(1., q, -1., p, Complex(0.5) * rad_->D_inner_x() - rad_->Mm1_tr_inner_x() + Complex(0.5*(l1+1)*l1) * rad_->Mm2_inner_x(), rad_->S_inner_y());
        kron_dot(1., q, -1., p, rad_->S_inner_x(), Complex(0.5) * rad_->D_inner_y() + Complex(inp_->Zp) * rad_->Mm1_tr_inner_y() + Complex(0.5*(l2+1)*l2) * rad_->Mm2_inner_y());
        
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
    std::size_t Nspline_inner_x = rad_->bspline_inner_x().Nspline();
    std::size_t Nspline_inner_y = rad_->bspline_inner_y().Nspline();
    
    // get workspace
    int ithread = 0;
#ifdef _OPENMP
    ithread = omp_get_thread_num();
#endif
    
    // encapsulated memory chunks
    ColMatrix<Complex> U (Nspline_inner_x, Nspline_inner_y, workspace_[ithread]);
    RowMatrixView<Complex> R (Nspline_inner_x, Nspline_inner_y, r);
    RowMatrixView<Complex> Z (Nspline_inner_x, Nspline_inner_y, z);
    
    // first Kronecker product
    {
        // U = (AV)B
        RowMatrixView<Complex> A (Nspline_inner_x, Nspline_inner_x, prec_atom_[l1].invCl_invsqrtS.data());
        ColMatrixView<Complex> B (Nspline_inner_y, Nspline_inner_y, prec_proj_[l2].invCl_invsqrtS.data());
        
        blas::gemm(1., A, R, 0., Z); // N³ operations
        blas::gemm(1., Z, B, 0., U); // N³ operations
    }
    
    // Divide elements by the diagonal; N² operations
    # pragma omp parallel for
    for (std::size_t i = 0; i < Nspline_inner_x; i++) 
    for (std::size_t j = 0; j < Nspline_inner_y; j++)
        Z(i,j) = U(i,j) / (E_ - prec_atom_[l1].Dl[i] - prec_proj_[l2].Dl[j]);
    
    // second Kronecker product
    {
        // W = (AU)B
        RowMatrixView<Complex> A (Nspline_inner_x, Nspline_inner_x, prec_atom_[l1].invsqrtS_Cl.data());
        ColMatrixView<Complex> B (Nspline_inner_y, Nspline_inner_y, prec_proj_[l2].invsqrtS_Cl.data());
        
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
            prec_atom_[l1].drop();
        if (refcount_proj_[l2] == 0)
            prec_proj_ [l2].drop();
    }
    unlock_kpa_access();
    
    // exit parent
    CGPreconditioner::CG_exit(iblock);
}

void KPACGPreconditioner::finish ()
{
    prec_atom_.clear();
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
