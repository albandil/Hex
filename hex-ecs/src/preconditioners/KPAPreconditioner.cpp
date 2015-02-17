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
#include <set>

#include "../arrays.h"
#include "../gauss.h"
#include "../misc.h"
#include "../preconditioners.h"
#include "../radial.h"

#ifndef NO_LAPACK

const std::string KPACGPreconditioner::name = "KPA";
const std::string KPACGPreconditioner::description = "Block inversion using conjugate gradients preconditioned by Kronecker product approximation.";

void KPACGPreconditioner::sData::hdflink (const char* file)
{
    filename = file;
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
        return false;
    
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

void KPACGPreconditioner::setup ()
{
    NoPreconditioner::setup();
    
    std::cout << "Set up KPA preconditioner" << std::endl;
    
    // compose list of angular momenta needed by this processor
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
    
    // link preconditioner matrices to disk files and try to load precomputed
    Array<bool> done (inp_.maxell + 1, true);
    for (int l : needed_l)
    {
        prec_[l].hdflink(format("kpa-%d-%d.hdf",l,par_.iproc()).c_str());
        done[l] = prec_[l].hdfload();
        if (done[l])
            std::cout << "\t- preconditioner data for l = " << l << " loaded from \"" << prec_[l].filename << "\"" << std::endl;
    }
    
    // if all preconditioners have been loaded, exit this routine
    if (all(done))
    {
        std::cout << std::endl;
        return;
    }
    
    std::cout << "\t- Overlap matrix factorization" << std::endl;
    
    Timer timer;
    
    // diagonalize overlap matrix
    ColMatrix<Complex> S = s_rad_.S().torow().T();
    ColMatrix<Complex> CL, CR; cArray DS;
    std::tie(DS,CL,CR) = S.diagonalize();
    ColMatrix<Complex> invCR = CR.invert();
    
    // convert eigenvalues to diagonal matrix
    ColMatrix<Complex> DSmat(DS.size()), invDSmat(DS.size());
    ColMatrix<Complex> DSsqrtmat(DS.size()), invDSsqrtmat(DS.size());
    for (unsigned i = 0; i < DS.size(); i++)
    {
        DSmat(i,i) = DS[i];
        DSsqrtmat(i,i) = std::sqrt(DS[i]);
        invDSmat(i,i) = 1.0 / DS[i];
        invDSsqrtmat(i,i) = 1.0 / std::sqrt(DS[i]);
    }
    
    // Now S = CR * DSmat * CR⁻¹
    std::cout << "\t\ttime: " << timer.nice_time() << std::endl;
    std::cout << "\t\tresidual: " << cArray((RowMatrix<Complex>(S) - RowMatrix<Complex>(CR) * DSmat * invCR).data()).norm() << std::endl;
    
    // compute √S and √S⁻¹
    RowMatrix<Complex> sqrtS = RowMatrix<Complex>(CR) * DSsqrtmat * invCR;
    RowMatrix<Complex> invsqrtS = RowMatrix<Complex>(CR) * invDSsqrtmat * invCR;
    
    // necessary Kronecker product
    SymDiaMatrix half_D_minus_Mm1_tr = 0.5 * s_rad_.D() - s_rad_.Mm1_tr();
    
    // diagonalize one-electron hamiltonians for all angular momenta
    for (int l : needed_l)
    {
        // skip loaded
        if (done[l])
            continue;
        
        // reset timer
        std::cout << "\t- One-electron Hamiltonian factorization (l = " << l << ")" << std::endl;
        timer.reset();
        
        // compose the one-electron hamiltonian
        ColMatrix<Complex> Hl ( (half_D_minus_Mm1_tr + (0.5*l*(l+1)) * rad().Mm2()).torow() );
        
        // symmetrically transform by inverse square root of the overlap matrix
        RowMatrix<Complex> tHl = invsqrtS * Hl * ColMatrix<Complex>(invsqrtS);
        
        // diagonalize the transformed matrix
        ColMatrix<Complex> ClL, ClR;
        std::tie(prec_[l].Dl,ClL,ClR) = ColMatrix<Complex>(tHl).diagonalize();
        ColMatrix<Complex> invClR = ClR.invert();
        
        // store the data
        prec_[l].invsqrtS_Cl = invsqrtS * ClR;
        prec_[l].invCl_invsqrtS = RowMatrix<Complex>(invClR) * ColMatrix<Complex>(invsqrtS);
        
        // covert Dl to matrix form and print verification
        ColMatrix<Complex> Dlmat(prec_[l].Dl.size()), invDlmat(prec_[l].Dl.size());
        for (unsigned i = 0; i < prec_[l].Dl.size(); i++)
        {
            Dlmat(i,i) = prec_[l].Dl[i];
            invDlmat(i,i) = 1.0 / prec_[l].Dl[i];
        }
        
        // Now Hl = ClR * Dlmat * ClR⁻¹
        std::cout << "\t\t- time: " << timer.nice_time() << std::endl;
        std::cout << "\t\t- residual: " << cArray((tHl - RowMatrix<Complex>(ClR) * Dlmat * invClR).data()).norm() << std::endl;
        
        // Save to disk for possible future reuse.
        prec_[l].hdfsave();
    }
    
    std::cout << std::endl;
}

void KPACGPreconditioner::CG_mmul (int iblock, const cArrayView p, cArrayView q) const
{
    // let the parent do it if lightweight mode is off
    if (not cmd_.lightweight_full)
    {
        CGPreconditioner::CG_mmul(iblock, p, q);
    }
    else
    {
        // get block angular momemnta
        int l1 = l1_l2_[iblock].first;
        int l2 = l1_l2_[iblock].second;
        
        // multiply 'p' by the diagonal block (except for the two-electron term)
        q  = kron_dot(Complex(E_) * s_rad_.S_d(), s_rad_.S_d(), p);
        q -= kron_dot(Complex(0.5) * s_rad_.D_d() - s_rad_.Mm1_tr_d() + Complex(0.5 * (l1 + 1.) * l1) * s_rad_.Mm2_d(), s_rad_.S_d(), p);
        q -= kron_dot(s_rad_.S_d(), Complex(0.5) * s_rad_.D_d() - s_rad_.Mm1_tr_d() + Complex(0.5 * (l2 + 1.) * l2) * s_rad_.Mm2_d(), p);
        
        // multiply 'p' by the two-electron integrals
        for (int lambda = 0; lambda <= s_rad_.maxlambda(); lambda++)
        {
            // calculate angular integral
            double f = special::computef(lambda, l1, l2, l1, l2, inp_.L);
            if (not std::isfinite(f))
                HexException("Invalid result of computef(%d,%d,%d,%d,%d,%d).", lambda, l1, l2, l1, l2, inp_.L);
            
            // multiply
            if (f != 0.)
                q -= s_rad_.apply_R_matrix(lambda, cArrays(1, f * p))[0];
        }
    }
}

void KPACGPreconditioner::CG_prec (int iblock, const cArrayView r, cArrayView z) const
{
    try
    {
        // get angular momenta of this block
        int l1 = l1_l2_[iblock].first;
        int l2 = l1_l2_[iblock].second;
        
        // dimension of the matrices
        int Nspline = s_bspline_.Nspline();
        
        // multiply by the first Kronecker product
        z = kron_dot(prec_[l1].invCl_invsqrtS, prec_[l2].invCl_invsqrtS, r);
        
        // divide by the diagonal
        # pragma omp parallel for collapse (2) if (cmd_.parallel_dot)
        for (int i = 0; i < Nspline; i++) 
        for (int j = 0; j < Nspline; j++)
            z[i * Nspline + j] /= E_ - prec_[l1].Dl[i] - prec_[l2].Dl[j];
        
        // multiply by the second Kronecker product
        z = kron_dot(prec_[l1].invsqrtS_Cl, prec_[l2].invsqrtS_Cl, z);
    }
    catch (std::exception & e)
    {
        HexException("Standard exception in KPA preconditioner: %s.", e.what());
    }
    catch (...)
    {
        HexException("Unknown exception in KPA preconditioner.");
    }
}

#endif
