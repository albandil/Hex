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
#include "hex-blas.h"
#include "hex-hydrogen.h"
#include "hex-misc.h"

#include "gauss.h"
#include "preconditioners.h"
#include "radial.h"

#ifndef NO_LAPACK

const std::string KPACGPreconditioner::prec_name = "KPA";
const std::string KPACGPreconditioner::prec_description = 
    "Block inversion using conjugate gradients preconditioned by Kronecker product approximation.";

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
    // get reference to √S⁻¹
    ColMatrix<Complex> const & invsqrtS = rad_.inveigenstates(ang_.maxell() + 1);
    
    // for all angular momenta
    for (int l : comp_l)
    {
        // skip loaded
        if (done[l])
            continue;
        
        std::cout << "\t- compose preconditioner for l = " << l << std::endl;
        
        // store the preconditioner data
        prec[l].Dl = rad_.eigenenergies(l);
        prec[l].invsqrtS_Cl = RowMatrix<Complex>(Nspline,Nspline);
        prec[l].invCl_invsqrtS = RowMatrix<Complex>(Nspline,Nspline);
        blas::gemm(1., invsqrtS, rad_.eigenstates(l), 0., prec[l].invsqrtS_Cl);
        blas::gemm(1., rad_.inveigenstates(l), invsqrtS, 0., prec[l].invCl_invsqrtS);
        prec[l].hdfsave();
        prec[l].drop();
    }
    
    // wait for completition of diagonalization on other nodes
    par_.wait();

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
        for (int l = 0; l <= inp_.maxell; l++)
        {
            // check if this angular momentum is needed by some of the blocks owned by this process
            bool need_this_l = false;
            for (unsigned ill = 0; ill < ang_.states().size(); ill++)
            {
                if (par_.isMyWork(ill) and (ang_.states()[ill].first == l or ang_.states()[ill].second == l))
                    need_this_l = true;
            }
            
            // append the 'l' if needed
            if (need_this_l)
                needed_l.insert(l);
        }
        
        // compose list of angular momenta that this MPI node will COMPUTE
        std::set<int> comp_l;
        for (int l = 0; l <= inp_.maxell; l++)
        {
            // compute l-th preconditioner matrices if either
            // a) nodes share scratch and even distribution assigns l-th preconditioner matrices to this node, or
            if (cmd_.shared_scratch and par_.isMyWork(l))
                comp_l.insert(l);
            // b) nodes do not share scratch and this node will need l-th preconditioner matrices
            if (not cmd_.shared_scratch and needed_l.find(l) != needed_l.end())
                comp_l.insert(l);
        }
        
    //
    // Check presence of the atomic electron preconditioner
    //
        
        // status array indicating necessity to calculate the preconditioner matrix for given 'l'
        Array<bool> done_atom (inp_.maxell + 1, true);
        
        // "to compute matrices": link them to scratch disk files and check presence
        for (int l : comp_l)
        {
            prec_atom_[l].hdflink(format("kpa-%d-%.4lx.hdf",l,rad_.bspline_atom().hash()).c_str());
            done_atom[l] = prec_atom_[l].hdfcheck();
        }
        
        // "needed matrices": link them to scratch disk files and check presence, load if present
        for (int l : needed_l)
        {
            prec_atom_[l].hdflink(format("kpa-%d-%.4lx.hdf",l,rad_.bspline_atom().hash()).c_str());
            done_atom[l] = prec_atom_[l].hdfcheck();
            
            if (done_atom[l])
            {
                std::cout << "\t- atomic preconditioner data for l = " << l
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
                prec_atom_, bspline_atom_.Nspline(),
                rad_.S_atom(), rad_.D_atom(), rad_.Mm1_tr_atom(), rad_.Mm2_atom(),
                done_atom, comp_l, needed_l
            );
        }
        
        std::cout << std::endl;
        
    //
    // Check presence of the projectile electron preconditioner
    //
        
        // status array indicating necessity to calculate the preconditioner matrix for given 'l'
        Array<bool> done_proj (inp_.maxell + 1, true);
        
        // "to compute matrices": link them to scratch disk files and check presence
        for (int l : comp_l)
        {
            prec_proj_[l].hdflink(format("kpa-%d-%.4lx.hdf",l,rad_.bspline_proj().hash()).c_str());
            done_proj[l] = prec_proj_[l].hdfcheck();
        }
        
        // "needed matrices": link them to scratch disk files and check presence, load if present
        for (int l : needed_l)
        {
            prec_proj_[l].hdflink(format("kpa-%d-%.4lx.hdf",l,rad_.bspline_atom().hash()).c_str());
            done_proj[l] = prec_proj_[l].hdfcheck();
            
            if (done_proj[l])
            {
                std::cout << "\t- projectile preconditioner data for l = " << l
                           << " present in \"" << prec_proj_[l].filename << "\"" << std::endl;
            }
        }
    
    //
    // Calculation of the preconditioner for atomic basis.
    //
    
        if (not all(done_proj))
        {
            std::cout << std::endl << "\tPrepare preconditioner matrices for projectile grid" << std::endl;
            prepare
            (
                prec_proj_, bspline_proj_.Nspline(),
                rad_.S_proj(), rad_.D_proj(), rad_.Mm1_tr_proj(), rad_.Mm2_proj(),
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

void KPACGPreconditioner::update (double E)
{
    // update parent
    CGPreconditioner::update(E);
    
    // determine knot where all open bound-state channels decrease below the drop tolerance
    if (E >= 0 or cmd_.kpa_drop < 0)
    {
        // use whole matrices
        maxknot_ = -1;
    }
    else if (cmd_.kpa_drop > 0)
    {
        // get knot from user-supplied distance
        maxknot_  = rad_.bspline_atom().knot(cmd_.kpa_drop);
        
        std::cout << "\tKPA: dropping splines beyond " << cmd_.kpa_drop << " a.u.; B-spline index "
                  << maxknot_ << " and up (of " << rad_.bspline_atom().Nspline() << ")" << std::endl;
    }
    else
    {
        // get highest open bound state channel
        int n = std::floor(std::sqrt(-1.0 / (2.0 * E)));
        
        // get smallest non-zero knot separation
        double h = special::constant::Inf;
        for (int i = 0; i < rad_.bspline_atom().Nspline(); i++)
            if (rad_.bspline_atom().t(i) != rad_.bspline_atom().t(i + 1))
                h = std::min(h, std::abs(rad_.bspline_atom().t(i) - rad_.bspline_atom().t(i + 1)));
        
        // hunt & bisect for drop tolerance, start at classical turning point
        double r1, r2, r3;
        for (r1 = 2 * n, r3 = 2 * r1;
             std::abs(Hydrogen::P(n, 0, r3)) > cmd_.droptol;
             r1 = r3, r3 *= 2);
        for (r2 = 0.5 * (r1 + r3);
             std::abs(r1 - r3) > h;
             (std::abs(Hydrogen::P(n, 0, r2)) > cmd_.droptol ? (r1 = r2) : (r3 = r2)), r2 = 0.5 * (r1 + r3));
        
        // get knot
        maxknot_  = rad_.bspline_atom().knot(r2);
        
        std::cout << "\tKPA: dropping splines beyond " << r2 << " a.u. (" << n << "s's " << cmd_.droptol << " threshold); B-spline index "
                  << maxknot_ << " and up (of " << rad_.bspline_atom().Nspline() << ")" << std::endl;
    }
}

void KPACGPreconditioner::rhs (BlockArray<Complex>& chi, int ienergy, int instate) const
{
    // let parent bake the right-hand side
    CGPreconditioner::rhs(chi, ienergy, instate);
    
    // constrain all blocks from the very beginning
    for (unsigned iblock = 0; iblock < ang_.states().size(); iblock++)
    {
        if (cmd_.outofcore)
            chi[iblock].hdfload();
        
        this->CG_constrain(chi[iblock]);
        
        if (cmd_.outofcore)
            chi[iblock].hdfsave(), chi[iblock].drop();
    }
}

void KPACGPreconditioner::CG_init (int iblock) const
{
    // initialize parent
    CGPreconditioner::CG_init(iblock);
    
    // initialize self
    {
        // get block angular momenta
        int l1 = ang_.states()[iblock].first;
        int l2 = ang_.states()[iblock].second;
        
        // load preconditioner from disk
        if (not prec_atom_[l1].hdfload())
            HexException("Failed to read preconditioner matrix for l = %d.", l1);
        if (not prec_proj_[l2].hdfload())
            HexException("Failed to read preconditioner matrix for l = %d.", l2);
    }
    
    // calculate workspace size
    std::size_t size = bspline_atom_.Nspline() * bspline_proj_.Nspline();
    if (maxknot_ > 0)
        size *= 2;
    
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
    if (cmd_.kpa_simple_rad or cmd_.lightweight_full)
    {
        // get block angular momemnta
        int l1 = ang_.states()[iblock].first;
        int l2 = ang_.states()[iblock].second;
        
        // multiply 'p' by the diagonal block (except for the two-electron term)
        kron_dot(0., q,  1., p, Complex(E_) * rad_.S_atom(), rad_.S_proj());
        kron_dot(1., q, -1., p, Complex(0.5) * rad_.D_atom() - rad_.Mm1_tr_atom() + Complex(0.5*(l1+1)*l1) * rad_.Mm2_atom(), rad_.S_proj());
        kron_dot(1., q, -1., p, rad_.S_atom(), Complex(0.5) * rad_.D_proj() - rad_.Mm1_tr_proj() + Complex(0.5*(l2+1)*l2) * rad_.Mm2_proj());
        
        // multiply 'p' by the two-electron integrals
        for (int lambda = 0; lambda <= rad_.maxlambda(); lambda++)
        {
            // calculate angular integral
            double f = special::computef(lambda, l1, l2, l1, l2, inp_.L);
            if (not std::isfinite(f))
                HexException("Invalid result of computef(%d,%d,%d,%d,%d,%d).", lambda, l1, l2, l1, l2, inp_.L);
            
            // multiply
            if (f != 0.)
                rad_.apply_R_matrix(lambda, -f, p, 1., q, cmd_.kpa_simple_rad);
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
    int l1 = ang_.states()[iblock].first;
    int l2 = ang_.states()[iblock].second;
    
    // dimension of the matrices
    std::size_t Nspline_atom = bspline_atom_.Nspline();
    std::size_t Nspline_proj = bspline_proj_.Nspline();
    
    // get workspace
    int ithread = 0;
#ifdef _OPENMP
    ithread = omp_get_thread_num();
#endif
    Complex * work = workspace_[ithread].data();
    
    // encapsulated memory regions
    cArrayView U_data (Nspline_atom * Nspline_proj, work);
    ColMatrixView<Complex> U (Nspline_atom, Nspline_proj, U_data);
    RowMatrixView<Complex> R (Nspline_atom, Nspline_proj, r);
    RowMatrixView<Complex> Z (Nspline_atom, Nspline_proj, z);
    
    // multiply by the first Kronecker product
    if (maxknot_ < 0)
    {
        // U = (AV)B
        RowMatrixView<Complex> A (Nspline_atom, Nspline_atom, prec_atom_[l1].invCl_invsqrtS.data());
        ColMatrixView<Complex> B (Nspline_proj, Nspline_proj, prec_proj_[l2].invCl_invsqrtS.data());
        
        blas::gemm(1., A, R, 0., Z); // N³ operations
        blas::gemm(1., Z, B, 0., U); // N³ operations
    }
    else
    {
        //  ┏━━━━━━━━━━━━┓     ┏━━━━━━━━━━━━┓  ┏━┱──────────┐  ┏━━━━━━━━━━━━┓     
        //  ┃            ┃     ┃            ┃  ┃ ┃          │  ┡━━━━━━━━━━━━┩     
        //  ┃            ┃     ┃            ┃  ┃ ┃          │  │            │     First multiplication, then second !
        //  ┃     U₁     ┃  =  ┃     A₁     ┃  ┃ ┃   V      │  │     B₁'    │     A₁ ... prec_atom_[l1].invCl_invsqrtS
        //  ┃            ┃     ┃            ┃  ┃ ┃          │  │            │     B₁ ... prec_proj_[l2].invCl_invsqrtS [transposed]
        //  ┃            ┃     ┃            ┃  ┃ ┃          │  │            │
        //  ┗━━━━━━━━━━━━┛     ┗━━━━━━━━━━━━┛  ┗━┹──────────┘  └────────────┘
        //
        //  ┏━━━━━━━━━━━━┓     ┏━┱──────────┐  ┌─┲━━━━━━━━━━┓  ┌────────────┐     
        //  ┃            ┃     ┃ ┃          │  │ ┗━━━━━━━━━━┩  ┢━━━━━━━━━━━━┪     
        //  ┃            ┃     ┃ ┃          │  │            │  ┃            ┃     Second multiplocation, then first !
        //  ┃     U₁     ┃ +=  ┃ ┃   A₁     │  │     V      │  ┃     B₁'    ┃     A₁ ... prec_atom_[l1].invCl_invsqrtS
        //  ┃            ┃     ┃ ┃          │  │            │  ┃            ┃     B₁ ... prec_proj_[l2].invCl_invsqrtS [transposed]
        //  ┃            ┃     ┃ ┃          │  │            │  ┃            ┃
        //  ┗━━━━━━━━━━━━┛     ┗━┹──────────┘  └────────────┘  ┗━━━━━━━━━━━━┛
        //
        //
        
        // create views of the sub-matrices
        RowMatrixView<Complex> A_full (Nspline_atom, Nspline_atom, prec_atom_[l1].invCl_invsqrtS.data(), Nspline_atom);
        RowMatrixView<Complex> A_left (Nspline_atom, maxknot_, prec_atom_[l1].invCl_invsqrtS.data(), Nspline_atom);
        ColMatrixView<Complex> B_top  (maxknot_, Nspline_proj, prec_proj_[l2].invCl_invsqrtS.data(), Nspline_proj);
        ColMatrixView<Complex> B_bot  (Nspline_proj - maxknot_, Nspline_proj, cArrayView(prec_proj_[l2].invCl_invsqrtS.data(), maxknot_, (Nspline_proj - maxknot_) * Nspline_proj), Nspline_proj);
        RowMatrixView<Complex> V_left (Nspline_atom, maxknot_, r, Nspline_proj);
        RowMatrixView<Complex> V_top  (maxknot_, Nspline_proj - maxknot_, cArrayView(r, maxknot_, maxknot_ * Nspline_proj), Nspline_proj);
        
        // create view of working memory
        cArrayView Xa_data (Nspline_atom * maxknot_, work + Nspline_atom * Nspline_proj);
        cArrayView Xp_data (Nspline_proj * maxknot_, work + Nspline_atom * Nspline_proj);
        ColMatrixView<Complex> AV (Nspline_atom, maxknot_, Xa_data);
        ColMatrixView<Complex> VB (maxknot_, Nspline_proj, Xp_data);
        
        // calculate the matrix products
        blas::gemm(1., A_full, V_left, 0., AV); // N²m operations
        blas::gemm(1., AV, B_top, 0., U);       // N²m operations
        blas::gemm(1., V_top, B_bot, 0., VB);   // mN(N-m) operations
        blas::gemm(1., A_left, VB, 1., U);      // N²m operations
    }
    
    //  ┏━━━━━━━━━━━━┓     ┏━━━━━━━━━━━━┓     ┏━━━━━━━━━━━━┓     
    //  ┃            ┃     ┃            ┃     ┃            ┃     
    //  ┃            ┃     ┃            ┃     ┃            ┃     Hadamard product.
    //  ┃     U₂     ┃  =  ┃     D      ┃  o  ┃     U₁     ┃     Scaling by the diagonal.
    //  ┃            ┃     ┃            ┃     ┃            ┃ 
    //  ┃            ┃     ┃            ┃     ┃            ┃ 
    //  ┗━━━━━━━━━━━━┛     ┗━━━━━━━━━━━━┛     ┗━━━━━━━━━━━━┛ 
    
    // Divide elements by the diagonal; N² operations
    # pragma omp parallel for
    for (std::size_t i = 0; i < Nspline_atom; i++) 
    for (std::size_t j = 0; j < Nspline_proj; j++)
    {
        U(i,j) /= (E_ - prec_atom_[l1].Dl[i] - prec_proj_[l2].Dl[j]);
        Z(i,j) = U(i,j); // <-- this is used in full product
    }
    
    // multiply by the second Kronecker product
    if (maxknot_ < 0)
    {
        // W = (AU)B
        RowMatrixView<Complex> A (Nspline_atom, Nspline_atom, prec_atom_[l1].invsqrtS_Cl.data());
        ColMatrixView<Complex> B (Nspline_proj, Nspline_proj, prec_proj_[l2].invsqrtS_Cl.data());
        
        blas::gemm(1., A, Z, 0., U);    // N³ operations
        blas::gemm(1., U, B, 0., Z);    // N³ operations
    }
    else
    {
        //  ┏━┱──────────┐     ┏━━━━━━━━━━━━┓   ┏━━━━━━━━━━━━┓  ┏━┱──────────┐
        //  ┃ ┃          │     ┃            ┃   ┃            ┃  ┃ ┃          │
        //  ┃ ┃          │     ┃            ┃   ┃            ┃  ┃ ┃          │    Second multiplication, then first !
        //  ┃ ┃   W      │  =  ┃     A₂     ┃   ┃     U₂     ┃  ┃ ┃   B₂'    │    A₂ ... prec_atom_[l1].invsqrtS_Cl
        //  ┃ ┃          │     ┃            ┃   ┃            ┃  ┃ ┃          │    B₂ ... prec_atom_[l2].invsqrtS_Cl [transposed]
        //  ┃ ┃          │     ┃            ┃   ┃            ┃  ┃ ┃          │
        //  ┗━┹──────────┘     ┗━━━━━━━━━━━━┛   ┗━━━━━━━━━━━━┛  ┗━┹──────────┘
        
        //  ┌─┲━━━━━━━━━━┓     ┏━━━━━━━━━━━━┓   ┏━━━━━━━━━━━━┓  ┌─┲━━━━━━━━━━┓
        //  │ ┗━━━━━━━━━━┩     ┡━━━━━━━━━━━━┩   ┃            ┃  │ ┃          ┃
        //  │            │     │            │   ┃            ┃  │ ┃          ┃    First multiplication, then second !
        //  │     W      │  =  │     A₂     │   ┃     U₂     ┃  │ ┃   B₂'    ┃    A₂ ... prec_atom_[l1].invsqrtS_Cl
        //  │            │     │            │   ┃            ┃  │ ┃          ┃    B₂ ... prec_atom_[l2].invsqrtS_Cl [transposed]
        //  │            │     │            │   ┃            ┃  │ ┃          ┃
        //  └────────────┘     └────────────┘   ┗━━━━━━━━━━━━┛  └─┺━━━━━━━━━━┛
        
        // create views of the sub-matrices
        RowMatrixView<Complex> A_full (Nspline_atom, Nspline_atom, prec_atom_[l1].invsqrtS_Cl.data(), Nspline_atom);
        RowMatrixView<Complex> A_top  (maxknot_, Nspline_atom, prec_atom_[l1].invsqrtS_Cl.data(), Nspline_atom);
        ColMatrixView<Complex> B_left (Nspline_proj, maxknot_, prec_proj_[l2].invsqrtS_Cl.data(), Nspline_proj);
        ColMatrixView<Complex> B_right(Nspline_proj, Nspline_proj - maxknot_, cArrayView(prec_proj_[l2].invsqrtS_Cl.data(), maxknot_ * Nspline_proj, (Nspline_proj - maxknot_) * Nspline_proj), Nspline_proj);
        RowMatrixView<Complex> W_left (Nspline_atom, maxknot_, z, Nspline_proj);
        RowMatrixView<Complex> W_top  (maxknot_, Nspline_proj - maxknot_, cArrayView(z, maxknot_, maxknot_ * Nspline_proj), Nspline_proj);
        
        // create view of working memory
        cArrayView Xa_data (Nspline_atom * maxknot_, work + Nspline_atom * Nspline_proj);
        cArrayView Xp_data (Nspline_proj * maxknot_, work + Nspline_atom * Nspline_proj);
        ColMatrixView<Complex> UB (Nspline_atom, maxknot_, Xa_data);
        ColMatrixView<Complex> AU (maxknot_, Nspline_proj, Xp_data);
        
        // calculate the matrix products
        blas::gemm(1., U, B_left, 0., UB);      // N²m operations
        blas::gemm(1., A_full, UB, 0., W_left); // N²m operations
        blas::gemm(1., A_top, U, 0., AU);       // N²m operations
        blas::gemm(1., AU, B_right, 0., W_top); // m(N-m)N operations
    }
    
    // erase thresholded part of the output
    if (maxknot_ > 0)
    {
        # pragma omp parallel for
        for (std::size_t i = maxknot_; i < Nspline_atom; i++) 
        for (std::size_t j = maxknot_; j < Nspline_proj; j++)
            Z(i,j) = 0.;
    }
}

void KPACGPreconditioner::CG_exit (int iblock) const
{
    // exit self
    {
        // get block angular momenta
        int l1 = ang_.states()[iblock].first;
        int l2 = ang_.states()[iblock].second;
        
        // release memory
        prec_atom_[l1].drop();
        prec_proj_[l2].drop();
    }
    
    // exit parent
    CGPreconditioner::CG_exit(iblock);
}

void KPACGPreconditioner::finish ()
{
    prec_atom_.clear();
    prec_proj_.clear();
    CGPreconditioner::finish();
}

void KPACGPreconditioner::CG_constrain (cArrayView r) const
{
    if (maxknot_ > 0)
    {
        std::size_t Nspline_atom = bspline_atom_.Nspline();
        std::size_t Nspline_proj = bspline_proj_.Nspline();
        
        # pragma omp parallel for
        for (std::size_t i = maxknot_; i < Nspline_atom; i++) 
        for (std::size_t j = maxknot_; j < Nspline_proj; j++)
            r[i * Nspline_proj + j] = 0.;
    }
}

#endif
