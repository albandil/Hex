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

#include "hex-hydrogen.h"

#include "coupled-solver.h"

CoupledSolver::CoupledSolver
(
    CommandLine & cmd,
    InputFile const & inp,
    Parallel const & par,
    AngularBasis const & ang,
    std::vector<Bspline> const & bspline,
    std::vector<Bspline> const & bspline_full
) : cmd_(cmd), inp_(inp), par_(par), ang_(ang), bspline_(bspline),
    rad_(bspline.back(), bspline.back(), bspline_full.back(), ang.maxlambda() + 1)
{
    // nothing to do
}

void CoupledSolver::choose_preconditioner (int ipanel)
{
    
}

void CoupledSolver::setup_preconditioner ()
{
    if (not (cmd_.itinerary & CommandLine::StgRadial))
    {
        std::cout << "Skipped computation of radial integrals." << std::endl;
        return;
    }
    
    // setup the preconditioner (compute radial integrals etc.)
//     prec_->setup();
    rad_.setupOneElectronIntegrals(par_, cmd_);
    rad_.setupTwoElectronIntegrals(par_, cmd_);
}

void CoupledSolver::solve ()
{
    if (not (cmd_.itinerary & CommandLine::StgSolve))
    {
        std::cout << "Skipped solution of the equation." << std::endl;
        return;
    }
    
    // shorthands
    int Nspline = bspline_[0].Nspline();
    
    // print some information on the Hamiltonian
    std::cout << "Hamiltonian size: " << Nspline * Nspline * ang_.states().size() << std::endl;
    
    double E = special::constant::Nan;
    int iterations_done = 0, computations_done = 0;
    
    for (unsigned ie = 0; ie < inp_.Etot.size(); ie++)
    {
        // print progress information
        std::cout << "\nSolving the system for Etot[" << ie << "] = " << inp_.Etot[ie] << " ("
                  << int(std::trunc(ie * 100. / inp_.Etot.size() + 0.5)) << " % finished" << std::endl;
        
        // we may have already computed all solutions for this energy... is it so?
        std::vector<std::pair<int,int>> work;
        for (unsigned instate = 0; instate < inp_.instates.size(); instate++)
        for (unsigned Spin = 0; Spin <= 1; Spin++)
        {
            // decode initial state
            int ni = std::get<0>(inp_.instates[instate]);
            int li = std::get<1>(inp_.instates[instate]);
            int mi = std::get<2>(inp_.instates[instate]);
            
            // skip energy-forbidden states
            if (inp_.Etot[ie] <= -1./(ni*ni))
            {
                std::cout << "\tSkip initial state " << Hydrogen::stateName(ni,li,mi) << " (S = " << Spin
                          << ") : not allowed by total E." << std::endl;;
                continue;
            }
            
            // check if the right hand side will be zero for this instate
            bool allowed = false;
            for (int l = std::abs(li - inp_.L); l <= li + inp_.L; l++)
            {
                // does this combination conserve parity?
                if ((inp_.L + li + l) % 2 != inp_.Pi)
                    continue;
                
                // does this combination have valid 'mi' for this partial wave?
                if (special::ClebschGordan(li,mi,l,0,inp_.L,mi) != 0)
                    allowed = true;
            }
            
            // skip angular forbidden states
            if (not allowed)
            {
                std::cout << "\tSkip initial state " << Hydrogen::stateName(ni,li,mi) << " (S = " << Spin
                          << ") : not allowed by total L, Pi and nL." << std::endl;
                continue;
            }
            
            // check if there is some precomputed solution on the disk
            SolutionIO reader (inp_.L, Spin, inp_.Pi, ni, li, mi, inp_.Etot[ie], ang_.states());
            std::size_t size = reader.check();
            
            // solution has the expected size
            if (size == (std::size_t)Nspline * (std::size_t)Nspline)
            {
                std::cout << "\tSolution for initial state " << Hydrogen::stateName(ni,li,mi) << " (S = " << Spin << ") found." << std::endl;
                continue;
            }
            
            // add work
            work.push_back(std::make_pair(instate,Spin));
        }
        
        // skip this energy if nothing to compute
        if (work.empty())
        {
            std::cout << "\tAll solutions for Etot[" << ie << "] = " << inp_.Etot[ie] << " loaded." << std::endl;
            continue;
        }
        
        // update the preconditioner, if this is the first energy to compute or it changed from previous iteration
        E = 0.5 * inp_.Etot[ie];
        
        // for all initial states
        for (auto workitem : work)
        {
            // decode initial state
            int instate = std::get<0>(workitem);
            ang_.S() = std::get<1>(workitem);
            int ni = std::get<0>(inp_.instates[instate]);
            int li = std::get<1>(inp_.instates[instate]);
            int mi = std::get<2>(inp_.instates[instate]);
            
            // create right hand side
            cArray chi (ang_.states().size() * Nspline * Nspline);
            std::cout << "\tCreate right-hand side for initial state " << Hydrogen::stateName(ni,li,mi) << " and total spin S = " << ang_.S() << " ... " << std::flush;
            rhs(chi, ie, instate);
            std::cout << "ok" << std::endl;
            
            // compute and check norm of the right hand side vector
            double chi_norm = chi.norm();
            if (chi_norm == 0.)
            {
                // this should not happen, hopefully we already checked
                std::cout << "\t! Right-hand-side is zero, check L, Pi and nL." << std::endl;
                continue;
            }
            if (not std::isfinite(chi_norm))
            {
                // this is a numerical problem, probably in evaluation of special functions (P, j)
                std::cout << "\t! Right hand side has invalid norm (" << chi_norm << ")." << std::endl;
                continue;
            }
            
            // prepare solution vector
            cArray psi (chi.size());
            
            std::cout << "\tSolve for initial state " << Hydrogen::stateName(ni,li,mi) << " and total spin S = " << ang_.S() << "." << std::endl;
            
            // calculate matrix elements
            Array<NumberArray<LU_int_t>> column_indices (ang_.states().size() * Nspline * Nspline);
            Array<NumberArray<Complex>> elements (ang_.states().size() * Nspline * Nspline);
            # pragma omp parallel for collapse (3)
            for (std::size_t ill = 0; ill < ang_.states().size(); ill++)
            for (std::size_t i = 0; i < Nspline; i++)
            for (std::size_t j = 0; j < Nspline; j++)
            {
                // calculate row index
                std::size_t row_index = (ill * Nspline + i) * Nspline + j;
                
                // angular momenta
                unsigned l1 = ang_.states()[ill].first;
                unsigned l2 = ang_.states()[ill].second;
                
                for (std::size_t illp = 0; illp < ang_.states().size(); illp++)
                for (std::size_t k = 0; k < Nspline; k++)
                for (std::size_t l = 0; l < Nspline; l++)
                {
                    // calculate column index
                    std::size_t col_index = (illp * Nspline + k) * Nspline + l;
                    
                    // matrix element
                    Complex a = 0;
                    
                    // diagonal contribution
                    if (ill == illp)
                    {
                        a = E * rad_.S_atom()(i,k) * rad_.S_proj()(j,l);
                        a -= 0.5 * rad_.D_atom()(i,k) * rad_.S_proj()(j,l);
                        a -= 0.5 * rad_.S_atom()(i,k) * rad_.D_proj()(j,l);
                        a -= 0.5 * l1 * (l1 + 1.) * rad_.Mm2_atom()(i,k) * rad_.S_proj()(j,l);
                        a -= 0.5 * l2 * (l2 + 1.) * rad_.S_atom()(i,k) * rad_.Mm2_proj()(j,l);
                        a += rad_.Mm1_tr_atom()(i,k) * rad_.S_proj()(j,l);
                        a += rad_.S_atom()(i,k) * rad_.Mm1_tr_proj()(j,l);
                    }
                    
                    // off-diagonal contribution
                    for (unsigned lambda = 0; lambda <= ang_.maxlambda(); lambda++)
                    {
                        double f = ang_.f(ill, illp, lambda);
                        if (f != 0.)
                            a -= f * rad_.computeR(lambda, i, j, k, l);
                    }
                    
                    // insert the new element
                    if (a != 0.)
                    {
                        column_indices[row_index].push_back(col_index);
                        elements[row_index].push_back(a);
                    }
                }
            }
            
            // collect row pointers
            NumberArray<LU_int_t> row_ptrs = { 0 };
            for (std::size_t i = 1; i <= column_indices.size(); i++)
                row_ptrs.push_back(row_ptrs.back() + column_indices[i - 1].size());
            
            // create the matrix in CSR format
            CsrMatrix<LU_int_t,Complex> A
            (
                chi.size(), chi.size(),
                row_ptrs,
                join(column_indices),
                join(elements)
            );
            
            // factorize matrix
            std::shared_ptr<LUft<LU_int_t,Complex>> LU = A.factorize(1e-15);
            
            // solve the system
            LU->solve(chi, psi, 1);
            
            std::cout << "\tSolution found." << std::endl;
            
            // save solution to disk (if valid)
            SolutionIO reader (ang_.L(), ang_.S(), ang_.Pi(), ni, li, mi, inp_.Etot[ie], ang_.states());
            if (std::isfinite(psi.norm()))
            {
                for (unsigned ill = 0; ill < ang_.states().size(); ill++)
                {
                    if (par_.isMyGroupWork(ill) and par_.IamGroupMaster())
                    {
                        // save origin panel
//                         if (ipanel_ == 0)
//                         {
                            // write the updated solution to disk
                            if (not reader.save_segment(psi.slice(ill * Nspline * Nspline, (ill + 1) * Nspline * Nspline), ill))
                                HexException("Failed to save solution to disk - the data are lost!");
//                         }
                        
                        // connect further panels to the previous solution
//                         else
//                         {
//                             HexException("Connection not implemented for coupled solver.");
//                         }
                    }
                }
            }
            
            // reset some one-solution command line flags
            cmd_.reuse_dia_blocks = false;
            
        } // end of For Spin, instate
        
    } // end of For ie = 0, ..., inp.Ei.size() - 1
    
    // wait for completition of all processes before next step
    par_.wait();
    
    std::cout << std::endl << "All solutions computed." << std::endl;
    if (computations_done > 0)
        std::cout << "\t(typically " << iterations_done / computations_done << " CG iterations per solution)" << std::endl;
}

void CoupledSolver::finish ()
{
//     prec_->finish();
}

void CoupledSolver::rhs (cArray& rhs, int ie, int instate) const
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
    
    // radial information for full projectil B-spline basis (used only to expand Riccati-Bessel function)
    RadialIntegrals radf (rad_.bspline_atom(), rad_.bspline_proj_full(), rad_.bspline_proj_full(), 0);
    radf.verbose(false);
    radf.setupOneElectronIntegrals(par_, cmd_);
    
    // calculate LU-decomposition of the overlap matrix
    CsrMatrix<LU_int_t,Complex> S_csr_atom = rad_.S_atom().tocoo<LU_int_t>().tocsr();
    CsrMatrix<LU_int_t,Complex> S_csr_proj = radf.S_proj().tocoo<LU_int_t>().tocsr();
    std::shared_ptr<LUft<LU_int_t,Complex>> lu_S_atom = S_csr_atom.factorize();
    std::shared_ptr<LUft<LU_int_t,Complex>> lu_S_proj = S_csr_proj.factorize();
    
    // j-overlaps of shape [Nangmom × Nspline]
    cArray ji_overlaps_atom = rad_.overlapj(rad_.bspline_atom(), rad_.gaussleg_atom(), inp_.maxell, ki, weightEdgeDamp(rad_.bspline_atom()));
    cArray ji_overlaps_proj = radf.overlapj(radf.bspline_proj(), radf.gaussleg_proj(), inp_.maxell, ki, weightEdgeDamp(radf.bspline_proj()));
    if (not std::isfinite(ji_overlaps_atom.norm()) or not std::isfinite(ji_overlaps_proj.norm()))
        HexException("Unable to compute Riccati-Bessel function B-spline overlaps!");
    
    // j-expansions
    cArray ji_expansion_atom = lu_S_atom->solve(ji_overlaps_atom, inp_.maxell + 1);
    cArray ji_expansion_proj = lu_S_proj->solve(ji_overlaps_proj, inp_.maxell + 1);
    if (not std::isfinite(ji_expansion_atom.norm()) or not std::isfinite(ji_expansion_proj.norm()))
        HexException("Unable to expand Riccati-Bessel function in B-splines!");
    
    // compute P-overlaps and P-expansion
    cArray Pi_overlaps_atom = rad_.overlapP(rad_.bspline_atom(), rad_.gaussleg_atom(), ni, li, weightEndDamp(rad_.bspline_atom()));
    cArray Pi_overlaps_proj = radf.overlapP(radf.bspline_proj(), radf.gaussleg_proj(), ni, li, weightEndDamp(radf.bspline_proj()));
    cArray Pi_expansion_atom = lu_S_atom->solve(Pi_overlaps_atom);
    cArray Pi_expansion_proj = lu_S_proj->solve(Pi_overlaps_proj);
    if (not std::isfinite(Pi_expansion_atom.norm()) or not std::isfinite(Pi_expansion_proj.norm()))
        HexException("Unable to expand hydrogen bound orbital in B-splines!");
    
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
    
    // for all segments constituting the RHS
    # pragma omp parallel for schedule (dynamic,1) if (cmd_.parallel_multiply)
    for (unsigned ill = 0; ill < ang_.states().size(); ill++) if (par_.isMyGroupWork(ill))
    {
        int l1 = ang_.states()[ill].first;
        int l2 = ang_.states()[ill].second;
        
        // setup storage
        cArray chi_block (Nspline_atom * Nspline_proj);
        
        // for all allowed angular momenta (by momentum composition) of the projectile
        for (int l = std::abs(li - ang_.L()); l <= li + ang_.L(); l++)
        {
            // skip wrong parity
            if ((ang_.L() + li + l) % 2 != ang_.Pi())
                continue;
            
            // (anti)symmetrization
            double Sign = ((ang_.S() + ang_.Pi()) % 2 == 0) ? 1. : -1.;
            
            // compute energy- and angular momentum-dependent prefactor
            Complex prefactor = std::pow(Complex(0.,1.),l)
                              * std::sqrt(special::constant::two_pi * (2 * l + 1))
                              * special::ClebschGordan(li,mi, l,0, inp_.L,mi) / ki[0];
            
            // skip non-contributing terms
            if (prefactor == 0.)
                continue;
            
            // pick the correct Bessel function expansion
            cArrayView Ji_expansion_atom (ji_expansion_atom, l * Nspline_atom, Nspline_atom);
            cArrayView Ji_expansion_proj (ji_expansion_proj, l * Nspline_proj, Nspline_proj);
            
            // compute outer products of B-spline expansions
            cArray Pj1 = outer_product(Pi_expansion_atom, Ji_expansion_proj);
            cArray Pj2 = outer_product(Ji_expansion_atom, Pi_expansion_proj);
            
            // skip angular forbidden right hand sides
            for (int lambda = 0; lambda <= rad_.maxlambda(); lambda++)
            {
                // calculate angular integrals
                double f1 = special::computef(lambda, l1, l2, li, l, inp_.L);
                double f2 = special::computef(lambda, l1, l2, l, li, inp_.L);
                
                // abort if any of the coefficients is non-number (factorial overflow etc.)
                if (not std::isfinite(f1))
                    HexException("Invalid result of computef(%d,%d,%d,%d,%d,%d)\n", lambda,l1,l2,li,l,inp_.L);
                if (not std::isfinite(f2))
                    HexException("Invalid result of computef(%d,%d,%d,%d,%d,%d)\n", lambda,l1,l2,l,li,inp_.L);
                
                // add multipole terms (direct/exchange)
                if (not cmd_.lightweight_radial_cache)
                {
                    if (f1 != 0.) rad_.R_tr_dia(lambda).dot(       prefactor * f1, Pj1, 1., chi_block, true);
                    if (f2 != 0.) rad_.R_tr_dia(lambda).dot(Sign * prefactor * f2, Pj2, 1., chi_block, true);
                }
                else
                {
                    if (f1 != 0.) rad_.apply_R_matrix(lambda,        prefactor * f1, Pj1, 1., chi_block);
                    if (f2 != 0.) rad_.apply_R_matrix(lambda, Sign * prefactor * f2, Pj2, 1., chi_block);
                }
            }
            
            // add monopole terms (direct/exchange)
            if (li == l1 and l == l2)
                chi_block += (-prefactor       ) * outer_product(rad_.S_atom().dot(Pi_expansion_atom), rad_.Mm1_tr_proj().dot(Ji_expansion_proj));
            if (li == l2 and l == l1)
                chi_block += (-prefactor * Sign) * outer_product(rad_.Mm1_tr_atom().dot(Ji_expansion_atom), rad_.S_proj().dot(Pi_expansion_proj));
        }
        
        // update the right-hand side
        cArrayView(rhs, ill * Nspline_atom * Nspline_proj, Nspline_atom * Nspline_proj) = chi_block;
    }
}
