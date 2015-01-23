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

const std::string KPACGPreconditioner::name = "KPA";
const std::string KPACGPreconditioner::description = "Block inversion using conjugate gradients preconditioned by Kronecker product approximation.";

void KPACGPreconditioner::setup ()
{
    NoPreconditioner::setup();
    
    std::cout << "Set up KPA preconditioner" << std::endl;
    std::cout << "\t- Overlap matrix factorization" << std::endl;
    
    Timer timer;
    
    // resize arrays
    invsqrtS_Cl_.resize(inp_.maxell + 1);
    invCl_invsqrtS_.resize(inp_.maxell + 1);
    Dl_.resize(inp_.maxell + 1);
    
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
    for (int l = 0; l <= inp_.maxell; l++)
    {
        // check if this angular momentum is needed by some of the blocks owned by this process
        bool need_l = false;
        for (unsigned ill = 0; ill < l1_l2_.size(); ill++)
        {
            if (par_.isMyWork(ill) and (l1_l2_[ill].first == l or l1_l2_[ill].second == l))
                need_l = true;
        }
        
        // skip the angular momentum if no owned diagonal block needs it
        if (not need_l)
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
        std::tie(Dl_[l],ClL,ClR) = ColMatrix<Complex>(tHl).diagonalize();
        ColMatrix<Complex> invClR = ClR.invert();
        
        // store the data
        invsqrtS_Cl_[l] = invsqrtS * ClR;
        invCl_invsqrtS_[l] = RowMatrix<Complex>(invClR) * ColMatrix<Complex>(invsqrtS);
        
        // covert Dl to matrix form and print verification
        ColMatrix<Complex> Dlmat(Dl_[l].size()), invDlmat(Dl_[l].size());
        for (unsigned i = 0; i < Dl_[l].size(); i++)
        {
            Dlmat(i,i) = Dl_[l][i];
            invDlmat(i,i) = 1.0 / Dl_[l][i];
        }
        
        // Now Hl = ClR * Dlmat * ClR⁻¹
        std::cout << "\t\t- time: " << timer.nice_time() << std::endl;
        std::cout << "\t\t- residual: " << cArray((tHl - RowMatrix<Complex>(ClR) * Dlmat * invClR).data()).norm() << std::endl;
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
                Exception("Invalid result of computef(%d,%d,%d,%d,%d,%d).", lambda, l1, l2, l1, l2, inp_.L);
            
            // multiply
            if (f != 0.)
                q -= s_rad_.apply_R_matrix(lambda, cArrays(1, f * p))[0];
        }
    }
}

void KPACGPreconditioner::CG_prec (int iblock, const cArrayView r, cArrayView z) const
{
    // get angular momenta of this block
    int l1 = l1_l2_[iblock].first;
    int l2 = l1_l2_[iblock].second;
    
    // dimension of the matrices
    int Nspline = s_bspline_.Nspline();
    
    // multiply by the first Kronecker product
    z = kron_dot(invCl_invsqrtS_[l1], invCl_invsqrtS_[l2], r);
    
    // divide by the diagonal
    # pragma omp parallel for collapse (2) if (cmd_.parallel_dot)
    for (int i = 0; i < Nspline; i++) 
    for (int j = 0; j < Nspline; j++)
        z[i * Nspline + j] /= E_ - Dl_[l1][i] - Dl_[l2][j];
    
    // multiply by the second Kronecker product
    z = kron_dot(invsqrtS_Cl_[l1], invsqrtS_Cl_[l2], z);
}

#endif
