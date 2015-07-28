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
#include "../preconditioners.h"
#include "../radial.h"

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
    int Nspline,
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
    ColMatrix<Complex> S = mS.torow().T(), CR, invCR, invsqrtS;
    
    std::ofstream ofs ("S2.png");
    mS.torow().plot_abs(ofs);
    write_array(S.T().data(), "S2.txt");
    ofs.close();
    
    // diagonalize overlap matrix
    if (not all(done))
    {
        std::cout << "\t\t- overlap matrix factorization" << std::endl;
        
        S.diagonalize(D, nullptr, &CR);
        CR.invert(invCR);
        
        // Now S = CR * (D * CR⁻¹)
        std::cout << "\t\t\ttime: " << timer.nice_time() << std::endl;
        for (std::size_t i = 0; i < (std::size_t)Nspline * (std::size_t)Nspline; i++)
            invCR.data()[i] *= D[i % Nspline];
        std::cout << "\t\t\tresidual: " << (S - CR * invCR).data().norm() << std::endl;
        S = ColMatrix<Complex>();
        
        // compute √S⁻¹
        for (std::size_t i = 0; i < (std::size_t)Nspline * (std::size_t)Nspline; i++)
            invCR.data()[i] /= std::pow(D.data()[i % Nspline], 1.5);
        invsqrtS = std::move(CR * invCR);
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
        
        // symmetrically transform by inverse square root of the overlap matrix
        tHl = std::move(invsqrtS * tHl * invsqrtS);
        
        // diagonalize the transformed matrix
        tHl.diagonalize(D, nullptr, &CR);
        CR.invert(invCR);
        
        // store the preconditioner data
        prec[l].Dl = D;
        prec[l].invsqrtS_Cl = std::move(RowMatrix<Complex>(invsqrtS * CR));
        prec[l].invCl_invsqrtS = std::move(RowMatrix<Complex>(invCR * invsqrtS));
        prec[l].hdfsave();
        prec[l].drop();
        
        // Now Hl = ClR * D * ClR⁻¹
        std::cout << "\t\t\t- time: " << timer.nice_time() << std::endl;
        for (std::size_t i = 0; i < (std::size_t)Nspline * (std::size_t)Nspline; i++)
            invCR.data()[i] *= D[i % Nspline];
        std::cout << "\t\t\t- residual: " << (tHl - CR * invCR).data().norm() << std::endl;
    }
    
    // wait for completition of diagonalization on other nodes
    par_.wait();

    // load all preconditioner matrices needed by this MPI node
    for (int l : needed_l)
        if (not cmd_.outofcore and not prec[l].hdfload())
            HexException("Failed to read preconditioner matrix for l = %d.", l);
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
            for (unsigned ill = 0; ill < l1_l2_.size(); ill++)
            {
                if (par_.isMyWork(ill) and (l1_l2_[ill].first == l or l1_l2_[ill].second == l))
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
        
        // status array indicating necessity to calculate the preconditioner matrix for given 'l'
        Array<bool> done (inp_.maxell + 1, true);
        
        // "to compute matrices": link them to scratch disk files and check presence
        for (int l : comp_l)
        {
            prec_atom_[l].hdflink(format("kpa-%d-%.4lx.hdf",l,rad_.bspline_atom().hash()).c_str());
            prec_proj_[l].hdflink(format("kpa-%d-%.4lx.hdf",l,rad_.bspline_proj().hash()).c_str());
            
            done[l] = (prec_atom_[l].hdfcheck() and prec_proj_[l].hdfcheck());
        }
        
        // "needed matrices": link them to scratch disk files and check presence, load if present
        for (int l : needed_l)
        {
            prec_atom_[l].hdflink(format("kpa-%d-%.4lx.hdf",l,rad_.bspline_atom().hash()).c_str());
            prec_proj_[l].hdflink(format("kpa-%d-%.4lx.hdf",l,rad_.bspline_proj().hash()).c_str());
            
            done[l] = (prec_atom_[l].hdfcheck() and prec_proj_[l].hdfcheck());
            
            if (cmd_.outofcore and done[l])
            {
                std::cout << "\t- preconditioner data for l = " << l
                          << " present in \"" << prec_atom_[l].filename
                          << "\" and \"" << prec_proj_[l].filename
                          << "\"" << std::endl;
            }
            if (not cmd_.outofcore and (done[l] = (prec_atom_[l].hdfload() and prec_proj_[l].hdfload())))
            {
                std::cout << "\t- preconditioner data for l = " << l
                          << " loaded from \"" << prec_atom_[l].filename
                          << "\" and \"" << prec_proj_[l].filename
                          << "\"" << std::endl;
            }
        }
        
        // if all preconditioners have been loaded, exit this routine
        if (all(done))
            std::cout << std::endl;
    
    //
    // Calculation of the preconditioner for atomic basis.
    //
        
        // prepare preconditioner for atomic basis
        if (not all(done))
            std::cout << std::endl << "\tPrepare preconditioner matrices for atomic grid" << std::endl;
        prepare
        (
            prec_atom_, bspline_atom_.Nspline(),
            rad_.S_atom(), rad_.D_atom(), rad_.Mm1_tr_atom(), rad_.Mm2_atom(),
            done, comp_l, needed_l
        );
        
        // prepare preconditioner for projectile basis
        if (not all(done))
            std::cout << std::endl << "\tPrepare preconditioner matrices for projectile grid" << std::endl;
        prepare
        (
            prec_proj_, bspline_proj_.Nspline(),
            rad_.S_proj(), rad_.D_proj(), rad_.Mm1_tr_proj(), rad_.Mm2_proj(),
            done, comp_l, needed_l
        );
        
        std::cout << std::endl;
}

void KPACGPreconditioner::CG_init (int iblock) const
{
    // initialize parent
    CGPreconditioner::CG_init(iblock);
    
    // initialize self
    if (cmd_.outofcore)
    {
        // get block angular momenta
        int l1 = l1_l2_[iblock].first;
        int l2 = l1_l2_[iblock].second;
        
        // load preconditioner from disk
        if (not prec_atom_[l1].hdfload())
            HexException("Failed to read preconditioner matrix for l = %d.", l1);
        if (not prec_proj_[l2].hdfload())
            HexException("Failed to read preconditioner matrix for l = %d.", l2);
    }
}

void KPACGPreconditioner::CG_mmul (int iblock, const cArrayView p, cArrayView q) const
{
    // multiply by components
    if (cmd_.kpa_simple_rad or cmd_.lightweight_full)
    {
        // get block angular momemnta
        int l1 = l1_l2_[iblock].first;
        int l2 = l1_l2_[iblock].second;
        
        // multiply 'p' by the diagonal block (except for the two-electron term)
        q  = kron_dot(Complex(E_) * rad_.S_atom(), rad_.S_proj(), p);
        q -= kron_dot(Complex(0.5) * rad_.D_atom() - rad_.Mm1_tr_atom() + Complex(0.5*(l1+1)*l1) * rad_.Mm2_atom(), rad_.S_proj(), p);
        q -= kron_dot(rad_.S_atom(), Complex(0.5) * rad_.D_proj() - rad_.Mm1_tr_proj() + Complex(0.5*(l2+1)*l2) * rad_.Mm2_proj(), p);
        
        // multiply 'p' by the two-electron integrals
        for (int lambda = 0; lambda <= rad_.maxlambda(); lambda++)
        {
            // calculate angular integral
            double f = special::computef(lambda, l1, l2, l1, l2, inp_.L);
            if (not std::isfinite(f))
                HexException("Invalid result of computef(%d,%d,%d,%d,%d,%d).", lambda, l1, l2, l1, l2, inp_.L);
            
            // multiply
            if (f != 0.)
                q -= rad_.apply_R_matrix(lambda, f * p, cmd_.kpa_simple_rad);
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
    int l1 = l1_l2_[iblock].first;
    int l2 = l1_l2_[iblock].second;
    
    // dimension of the matrices
    std::size_t Nspline_atom = bspline_atom_.Nspline();
    std::size_t Nspline_proj = bspline_proj_.Nspline();
    
    // multiply by the first Kronecker product
    dense_kron_dot
    (
        Nspline_atom, Nspline_atom, prec_atom_[l1].invCl_invsqrtS.data().data(),
        Nspline_proj, Nspline_proj, prec_proj_[l2].invCl_invsqrtS.data().data(),
        r.data(), z.data()
    );
    
    // divide by the diagonal
    # pragma omp parallel for
    for (std::size_t i = 0; i < Nspline_atom; i++) 
    for (std::size_t j = 0; j < Nspline_proj; j++)
        z[i * Nspline_proj + j] /= E_ - prec_atom_[l1].Dl[i] - prec_proj_[l2].Dl[j];
    
    // multiply by the second Kronecker product
    dense_kron_dot
    (
        Nspline_atom, Nspline_atom, prec_atom_[l1].invsqrtS_Cl.data().data(),
        Nspline_proj, Nspline_proj, prec_proj_[l2].invsqrtS_Cl.data().data(),
        z.data(), z.data()
    );
}

void KPACGPreconditioner::CG_exit (int iblock) const
{
    // exit self
    if (cmd_.outofcore)
    {
        // get block angular momenta
        int l1 = l1_l2_[iblock].first;
        int l2 = l1_l2_[iblock].second;
        
        // release memory
        prec_atom_[l1].drop();
        prec_proj_[l2].drop();
    }
    
    // exit parent
    CGPreconditioner::CG_exit(iblock);
}

#endif
