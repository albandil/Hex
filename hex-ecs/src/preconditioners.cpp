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

#ifdef _OPENMP
    #include <omp.h>
#endif

#include "arrays.h"
#include "gauss.h"
#include "itersolve.h"
#include "misc.h"
#include "preconditioners.h"
#include "radial.h"

const std::string NoPreconditioner::name = "none";
const std::string NoPreconditioner::description = "\"Preconditioning\" by the identity matrix.";

void NoPreconditioner::setup ()
{
    // TODO : Determine which lambdas are needed by this process.
    // NOTE : At the moment each process holds in memory radial integrals
    //        for all lambdas, which needlessly raises memory requirements.
    Array<bool> lambdas (inp_.L + 2 * inp_.levels + 1, true);
    
    // compute one-electron radial integrals
    s_rad_.setupOneElectronIntegrals();
    
    // compute two-eletron radial integrals
    s_rad_.setupTwoElectronIntegrals(par_, cmd_, lambdas);
}

void NoPreconditioner::update (double E)
{
    // update energy
    E_ = E;
    
    // skip pre-calculation of the diagonal blocks in lightweight mode
    if (cmd_.lightweight)
        return;
    
    std::cout << "\tPrecompute diagonal blocks... " << std::flush;
    
    // setup diagonal blocks
    # pragma omp parallel for if (cmd_.parallel_block)
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyWork(ill))
    {
        // angular momenta
        int l1 = l1_l2_[ill].first;
        int l2 = l1_l2_[ill].second;
        
        // initialize diagonal block
        dia_blocks_[ill] = BlockSymDiaMatrix
        (
            s_bspline_.Nspline(),       // block count (and size)
            s_rad_.S().nzpattern(),     // block structure
            s_rad_.S().diag(),          // non-zero diagonals
            not cmd_.outofcore,         // keep in memory?
            format("dblk-%d.ooc", ill)  // scratch disk file name
        );
        
        // one-electron parts; for all blocks
        for (unsigned iblock = 0; iblock < dia_blocks_[ill].structure().size(); iblock++)
        {
            // get block position
            int i = dia_blocks_[ill].structure()[iblock].first;
            int j = dia_blocks_[ill].structure()[iblock].second;
            
            cArray block = E * s_rad_.S()(i,j) * s_rad_.S().data();
            block -= (0.5 * s_rad_.D()(i,j) - s_rad_.Mm1_tr()(i,j)) * s_rad_.S().data();
            block -= 0.5 * l1 * (l1 + 1) * s_rad_.Mm2()(i,j) * s_rad_.S().data();
            block -= s_rad_.S()(i,j) * (0.5 * s_rad_.D().data() - s_rad_.Mm1_tr().data());
            block -= 0.5 * l2 * (l2 + 1) * s_rad_.S()(i,j) * s_rad_.Mm2().data();
            
            // two-electron part
            for (int lambda = 0; lambda <= s_rad_.maxlambda(); lambda++)
            {
                // calculate angular integral
                Complex f = special::computef(lambda,l1,l2,l1,l2,inp_.L);
                
                // check that the "f" coefficient is valid (no factorial overflow etc.)
                if (not Complex_finite(f))
                    Exception("Overflow in computation of f[%d](%d,%d,%d,%d).", inp_.L, l1, l2, l1, l2);
                
                // check that the "f" coefficient is nonzero
                if (f == 0.)
                    continue;
            
                // add two-electron term
                block += (-f) * s_rad_.R_tr_dia(lambda).getBlock(iblock);
            }
            
            // save block
            dia_blocks_[ill].setBlock(iblock,block);
        }
    }
    
    par_.wait();
    std::cout << "ok" << std::endl;
}

void NoPreconditioner::rhs (cArrayView chi, int ie, int instate, int Spin) const
{
    // shorthands
    int li = std::get<1>(inp_.instates[instate]);
    int mi = std::get<2>(inp_.instates[instate]);
    
    // shorthands
    int Nspline = s_rad_.bspline().Nspline();
    
    // j-overlaps of shape [Nangmom × Nspline]
    cArray ji_overlaps = s_rad_.overlapj
    (
        inp_.maxell,
        inp_.ki.slice(ie, ie + 1), // use only one ki
        weightEdgeDamp(s_rad_.bspline())
    );
    ji_overlaps.hdfsave("ji_overlaps.hdf");
    if (not std::isfinite(ji_overlaps.norm()))
        Exception("Unable to compute Riccati-Bessel function B-spline overlaps!");
    
    // j-expansions
    cArray ji_expansion = s_rad_.S().tocoo().tocsr().solve(ji_overlaps, ji_overlaps.size() / Nspline);
    if (not std::isfinite(ji_expansion.norm()))
        Exception("Unable to expand Riccati-Bessel function in B-splines!");
    
    // compute P-overlaps and P-expansion
    cArray Pi_overlaps, Pi_expansion;
    Pi_overlaps = s_rad_.overlapP(inp_.ni, li, weightEndDamp(s_rad_.bspline()));
    Pi_expansion = s_rad_.S().tocoo().tocsr().solve(Pi_overlaps);
    if (not std::isfinite(Pi_expansion.norm()))
        Exception("Unable to expand hydrogen bound orbital in B-splines!");
    
    //
    // Add two-electron part.
    //
    
    // for all segments constituting the RHS
    # pragma omp parallel for if (cmd_.parallel_block)
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++)
    {
        int l1 = l1_l2_[ill].first;
        int l2 = l1_l2_[ill].second;
        
        // setup storage
        cArrayView chi_block (chi, ill * Nspline * Nspline, Nspline * Nspline);
        chi_block.fill(0);
        
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
                              * special::ClebschGordan(li,mi, l,0, inp_.L,mi) / inp_.ki[ie];
            
            // skip non-contributing terms
            if (prefactor == 0.)
                continue;
            
            // pick the correct Bessel function expansion
            cArrayView Ji_expansion (ji_expansion, l * Nspline, Nspline);
            
            // compute outer products of B-spline expansions
            cArray Pj1 = outer_product(Pi_expansion, Ji_expansion);
            cArray Pj2 = outer_product(Ji_expansion, Pi_expansion);
            
            // skip angular forbidden right hand sides
            for (int lambda = 0; lambda <= s_rad_.maxlambda(); lambda++) if (par_.isMyWork(lambda))
            {
                // calculate angular integrals
                double f1 = special::computef(lambda, l1, l2, li, l, inp_.L);
                double f2 = special::computef(lambda, l1, l2, l, li, inp_.L);
                
                // abort if any of the coefficients is non-number (factorial overflow etc.)
                if (not std::isfinite(f1))
                    Exception("Invalid result of computef(%d,%d,%d,%d,%d,%d)\n", lambda,l1,l2,li,l,inp_.L);
                if (not std::isfinite(f2))
                    Exception("Invalid result of computef(%d,%d,%d,%d,%d,%d)\n", lambda,l1,l2,l,li,inp_.L);
                
                // add multipole terms (direct/exchange)
                if (not cmd_.lightweight)
                {
                    if (f1 != 0.) chi_block += (       prefactor * f1) * s_rad_.R_tr_dia(lambda).dot(Pj1, cmd_.parallel_dot);
                    if (f2 != 0.) chi_block += (Sign * prefactor * f2) * s_rad_.R_tr_dia(lambda).dot(Pj2, cmd_.parallel_dot);
                }
                else
                {
                    if (f1 != 0.) chi_block += (       prefactor * f1) * s_rad_.apply_R_matrix(lambda, cArrays(1, Pj1))[0];
                    if (f2 != 0.) chi_block += (Sign * prefactor * f2) * s_rad_.apply_R_matrix(lambda, cArrays(1, Pj2))[0];
                }
            }
        }
    }
    
    // sum all process' contributions
    par_.syncsum(chi.data(), chi.size());
    
    //
    // Add one-electron part.
    //
    
    # pragma omp parallel for
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++)
    {
        int l1 = l1_l2_[ill].first;
        int l2 = l1_l2_[ill].second;
        
        // setup storage
        cArrayView chi_block (chi, ill * Nspline * Nspline, Nspline * Nspline);
        
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
                              * special::ClebschGordan(li,mi, l,0, inp_.L,mi) / inp_.ki[ie];
            
            // skip non-contributing terms
            if (prefactor == 0.)
                continue;
            
            // pick the correct Bessel function expansion
            cArrayView Ji_expansion (ji_expansion, l * Nspline, Nspline);
            
            // compute outer products of B-spline expansions
            cArray Pj1 = outer_product(Pi_expansion, Ji_expansion);
            cArray Pj2 = outer_product(Ji_expansion, Pi_expansion);
            
            // add monopole terms (direct/exchange)
            if (li == l1 and l == l2)
                chi_block += (-prefactor       ) * outer_product(s_rad_.S().dot(Pi_expansion), s_rad_.Mm1_tr().dot(Ji_expansion));
            if (li == l2 and l == l1)
                chi_block += (-prefactor * Sign) * outer_product(s_rad_.Mm1_tr().dot(Ji_expansion), s_rad_.S().dot(Pi_expansion));
        }
    }
}

void NoPreconditioner::multiply (const cArrayView p, cArrayView q) const
{
    // shorthands
    int Nspline = s_rad_.bspline().Nspline();
    int Nang = l1_l2_.size();
    int Nchunk = Nspline * Nspline;
    
    // clear output array
    std::memset(q.data(), 0, q.size() * sizeof(Complex));

    if (not cmd_.lightweight)
    {
        //
        // Simple multiplication by the super-matrix.
        //
        
        // multiply "p" by the diagonal blocks
        # pragma omp parallel for schedule (dynamic,1) if (cmd_.parallel_block)
        for (int ill = 0;  ill < Nang;  ill++)
        if (par_.isMyWork(ill))
        {
            // copy-from segment of "p"
            cArrayView p_block (p, ill * Nchunk, Nchunk);
            
            // copy-to segment of "q"
            cArrayView q_block (q, ill * Nchunk, Nchunk);
            
            // multiply
            q_block = dia_blocks_[ill].dot(p_block, cmd_.parallel_dot);
        }
        
        // multiply "p" by the off-diagonal blocks
        # pragma omp parallel for schedule (dynamic,1) if (cmd_.parallel_block)
        for (int lambda = 0; lambda <= s_rad_.maxlambda(); lambda++)
        if (par_.isMyWork(lambda))
        {
            // for all blocks of the super-matrix
            for (int ill = 0;  ill < Nang;  ill++)
            for (int illp = 0; illp < Nang; illp++)
            {
                // skip diagonal
                if (ill == illp)
                    continue;
                
                // row multi-index
                int l1 = l1_l2_[ill].first;
                int l2 = l1_l2_[ill].second;
                
                // column multi-index
                int l1p = l1_l2_[illp].first;
                int l2p = l1_l2_[illp].second;
                
                // calculate angular integral
                double f = special::computef(lambda, l1, l2, l1p, l2p, inp_.L);
                if (not std::isfinite(f))
                    Exception("Invalid result of computef(%d,%d,%d,%d,%d,%d).", lambda, l1, l2, l1p, l2p, inp_.L);
                
                // check non-zero
                if (f == 0.)
                    continue;
                
                // source segment of "p"
                cArrayView p_block (p, illp * Nchunk, Nchunk);
                
                // product
                cArray product = std::move( (-f) * s_rad_.R_tr_dia(lambda).dot(p_block, cmd_.parallel_dot) );
                
                // update array
                # pragma omp critical
                cArrayView (q, ill * Nchunk, Nchunk) += product;
            }
        }
    }
    else
    {
        //
        // Lightweight mode multiplication
        //
        
        // multiply "p" by the diagonal blocks
        # pragma omp parallel for schedule (dynamic,1) if (cmd_.parallel_block)
        for (int ill = 0;  ill < Nang;  ill++)
        if (par_.isMyWork(ill))
        {
            // copy-from segment of "p"
            cArrayView p_block (p, ill * Nchunk, Nchunk);
            
            // copy-to segment of "q"
            cArrayView q_block (q, ill * Nchunk, Nchunk);
            
            // get block angular momemnta
            int l1 = l1_l2_[ill].first;
            int l2 = l1_l2_[ill].second;
            
            // multiply 'p_block' by the diagonal block (except for the two-electron term)
            q_block  = kron_dot(Complex(E_) * s_rad_.S_d(), s_rad_.S_d(), p_block);
            q_block -= kron_dot(Complex(0.5) * s_rad_.D_d() - s_rad_.Mm1_tr_d() + Complex(0.5 * l1 * (l1 + 1.)) * s_rad_.Mm2_d(), s_rad_.S_d(), p_block);
            q_block -= kron_dot(s_rad_.S_d(), Complex(0.5) * s_rad_.D_d() - s_rad_.Mm1_tr_d() + Complex(0.5 * l2 * (l2 + 1.)) * s_rad_.Mm2_d(), p_block);
        }
        
        // multiply "p" by the off-diagonal blocks
        # pragma omp parallel for schedule (dynamic,1) if (cmd_.parallel_block)
        for (int lambda = 0; lambda <= s_rad_.maxlambda(); lambda++)
        if (par_.isMyWork(lambda))
        {
            // vectors to be multiplied
            cArrays src;
            
            // destination segments
            iArray dst;
            
            // assemble all source vectors to be multiplied
            for (int ill = 0;  ill < Nang;  ill++)
            for (int illp = 0; illp < Nang; illp++)
            {
                // row and column multi-indices
                int l1 = l1_l2_[ill].first, l2 = l1_l2_[ill].second;
                int l1p = l1_l2_[illp].first, l2p = l1_l2_[illp].second;
                
                // calculate angular integral
                double f = special::computef(lambda, l1, l2, l1p, l2p, inp_.L);
                if (not std::isfinite(f))
                    Exception("Invalid result of computef(%d,%d,%d,%d,%d,%d).", lambda, l1, l2, l1p, l2p, inp_.L);
                
                // append source segment of "p" and the destination position
                if (f != 0.)
                {
                    src.push_back((-f) * cArrayView(p, illp * Nchunk, Nchunk));
                    dst.push_back(ill);
                }
            }
            
            // multiply all source segments by R^lambda
            src = std::move( s_rad_.apply_R_matrix(lambda, src) );
            
            // update destination array
            # pragma omp critical
            for (unsigned i = 0; i < src.size(); i++)
                cArrayView(q, dst[i] * Nchunk, Nchunk) += src[i];
        }
    }
    
//     std::cout << "multiply : q.norm() = " << std::sqrt(sum(sqrabs(q))) << std::endl;
    
    // synchronize across processes by summing individual contributions
    par_.syncsum(q.data(), Nchunk * l1_l2_.size());
}

void NoPreconditioner::precondition (const cArrayView r, cArrayView z) const
{
    z = r;
}

const std::string CGPreconditioner::name = "cg";
const std::string CGPreconditioner::description = "Block inversion using plain conjugate gradients. Use --tolerance option to set the termination tolerance.";

void CGPreconditioner::precondition (const cArrayView r, cArrayView z) const
{
    // shorthands
    int Nspline = s_rad_.bspline().Nspline();
    
    // iterations
    iArray n (l1_l2_.size());
    
    # pragma omp parallel for schedule (dynamic, 1) if (cmd_.parallel_block)
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyWork(ill))
    {
        // create segment views
        cArrayView rview (r, ill * Nspline * Nspline, Nspline * Nspline);
        cArrayView zview (z, ill * Nspline * Nspline, Nspline * Nspline);
        
        // wrappers around the callbacks
        auto inner_mmul = [&](const cArrayView a, cArrayView b) { this->CG_mmul(ill, a, b); };
        auto inner_prec = [&](const cArrayView a, cArrayView b) { this->CG_prec(ill, a, b); };
        
        // solve using the CG solver
        n[ill] = cg_callbacks < cArray, cArrayView >
        (
            rview,                  // rhs
            zview,                  // solution
            cmd_.prec_itertol,      // preconditioner tolerance
            0,                      // min. iterations
            Nspline * Nspline,      // max. iteration
            inner_prec,             // preconditioner
            inner_mmul,             // matrix multiplication
            false                   // verbose output
        );
    }
    
    // broadcast inner preconditioner iterations
    par_.sync(n.data(), 1, l1_l2_.size());
    
    // inner preconditioner info (max and avg number of iterations)
    std::cout << " | ";
    std::cout << std::setw(5) << (*std::min_element(n.begin(), n.end()));
    std::cout << std::setw(5) << (*std::max_element(n.begin(), n.end()));
    std::cout << std::setw(5) << format("%g", std::accumulate(n.begin(), n.end(), 0) / float(n.size()));
    
    // synchronize data across processes
    par_.sync(z.data(), Nspline * Nspline, l1_l2_.size());
}

void CGPreconditioner::CG_mmul (int iblock, const cArrayView p, cArrayView q) const
{
    q = dia_blocks_[iblock].dot(p, cmd_.parallel_dot);
}

void CGPreconditioner::CG_prec (int iblock, const cArrayView r, cArrayView z) const
{
    z = r;
}

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
    if (not cmd_.lightweight)
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

const std::string ILUCGPreconditioner::name = "ILU";
const std::string ILUCGPreconditioner::description = "Block inversion using conjugate gradients preconditioned by Incomplete LU. The drop tolerance can be given as the --droptol parameter.";

void ILUCGPreconditioner::update (double E)
{
    if (E != E_)
    {
        // release outdated LU factorizations
        for (auto & lu : lu_)
        {
            lu.drop();
            lu.unlink();
        }
        
        // release outdated CSR diagonal blocks
        for (auto & csr : csr_blocks_)
        {
            csr.drop();
            csr.unlink();
        }
    }
    
    // update parent
    CGPreconditioner::update(E);
}

void ILUCGPreconditioner::CG_prec (int iblock, const cArrayView r, cArrayView z) const
{
    // load data from linked disk files
    if (cmd_.outofcore)
    {
        csr_blocks_[iblock].hdfload();
        # pragma omp critical
        lu_[iblock].silent_load();
    }
    
    // check that the factorization is loaded
    if (lu_[iblock].size() == 0)
    {
        // create CSR block
        csr_blocks_[iblock] = dia_blocks_[iblock].tocoo().tocsr();
        
        // start timer
        Timer timer;
        
        // factorize the block and store it
        lu_[iblock].transfer(csr_blocks_[iblock].factorize(droptol_));
        
        // print time and memory info for this block (one thread at a time)
        # pragma omp critical
        std::cout << std::endl << std::setw(37) << format
        (
            "\tLU #%d (%d,%d) in %d:%02d (%d MiB)",
            iblock, l1_l2_[iblock].first, l1_l2_[iblock].second,    // block identification (id, ℓ₁, ℓ₂)
            timer.seconds() / 60, timer.seconds() % 60,             // factorization time
            lu_[iblock].size() / 1048576                            // final memory size
        );
        
        // save the diagonal block
        if (cmd_.outofcore)
        {
            csr_blocks_[iblock].hdflink(format("csr-%d.ooc", iblock));
            csr_blocks_[iblock].hdfsave();
        }
    }
    
    // precondition by LU
    z = lu_[iblock].solve(r);
    
    // release memory
    if (cmd_.outofcore)
    {
        // link to a disk file and save (if not already done)
        if (lu_[iblock].name().size() == 0)
        {
            lu_[iblock].link(format("lu-%d.ooc", iblock));
            # pragma omp critical
            lu_[iblock].save();
        }
        
        // release memory objects
        lu_[iblock].drop();
        csr_blocks_[iblock].drop();
    }
}
