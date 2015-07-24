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

#include "hydrogen.h"
#include "solver.h"

Solver::Solver
(
    CommandLine & cmd,
    InputFile const & inp,
    Parallel const & par,
    std::vector<std::pair<int,int>> const & angs,
    std::vector<Bspline> const & bspline,
    std::vector<int> const & Nspline_fixed
) : cmd_(cmd), inp_(inp), par_(par), angs_(angs), bspline_(bspline), Nspline_fixed_(Nspline_fixed),
    prec_(nullptr), ipanel_(0)
{
    
}

void Solver::choose_preconditioner (int ipanel)
{
    // save panel choice
    ipanel_ = ipanel;
    
    // create the preconditioner
    prec_ = Preconditioners::choose(par_, inp_, angs_, bspline_[0], bspline_[ipanel_], cmd_);
    
    // check success
    if (prec_ == nullptr)
        HexException("Preconditioner %d not implemented.", cmd_.preconditioner);
}

void Solver::setup_preconditioner ()
{
    if (not (cmd_.itinerary & CommandLine::StgRadial))
    {
        std::cout << "Skipped computation of radial integrals." << std::endl;
        return;
    }
    
    // setup the preconditioner (compute radial integrals etc.)
    prec_->setup();
}

void Solver::solve ()
{
    if (not (cmd_.itinerary & CommandLine::StgSolve))
    {
        std::cout << "Skipped solution of the equation." << std::endl;
        return;
    }
    
    // shorthands
    int Nspline_atom = bspline_[0].Nspline();
    int Nspline_proj = bspline_[ipanel_].Nspline();
    
    // print some information on the Hamiltonian
    std::cout << "Hamiltonian size: " << Nspline_atom * Nspline_proj * angs_.size() << std::endl;
    
    // wrap member functions to lambda-functions for use in the CG solver
    auto apply_preconditioner = [&](BlockArray<Complex> const & r, BlockArray<Complex> & z) -> void { this->apply_preconditioner_(r,z); };
    auto matrix_multiply = [&](BlockArray<Complex> const & p, BlockArray<Complex> & q) -> void { this->matrix_multiply_(p,q); };
    auto scalar_product = [&](BlockArray<Complex> const & x, BlockArray<Complex> const & y) -> Complex { return this->scalar_product_(x,y); };
    auto compute_norm = [&](BlockArray<Complex> const & r) -> double { return this->compute_norm_(r); };
    auto axby_operation = [&](Complex a, BlockArray<Complex> & x, Complex b, BlockArray<Complex> const & y) -> void { this->axby_operation_(a,x,b,y); };
    auto new_array = [&](std::size_t N, std::string name) -> BlockArray<Complex> { return this->new_array_(N,name); };
    auto constrain = [&](BlockArray<Complex> & r) -> void { constrain_(r); };
    
    double E = special::constant::Nan;
    int iterations_done = 0, computations_done = 0;
    
    for (unsigned ie = 0; ie < inp_.Etot.size(); ie++)
    {
        // print progress information
        std::cout << "\nSolving the system for Etot[" << ie << "] = " << inp_.Etot[ie] << " ("
                  << int(std::trunc(ie * 100. / inp_.Etot.size() + 0.5)) << " % finished, typically "
                  << (computations_done == 0 ? 0 : iterations_done / computations_done)
                  << " CG iterations per energy)" << std::endl;
        
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
            if (inp_.Etot[ie] < -1./(ni*ni))
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
            SolutionIO reader (inp_.L, Spin, inp_.Pi, ni, li, mi, inp_.Etot[ie], angs_);
            std::size_t size = reader.check();
            
            // solution has the expected size
            if (size == (std::size_t)Nspline_atom * (std::size_t)Nspline_proj)
                continue;
            
            // solution is smaller
            if (0 < size and size < (std::size_t)Nspline_atom * (std::size_t)Nspline_proj)
            {
                // is the size equal to previous panel basis size?
                if (ipanel_ > 0 and size == (std::size_t)Nspline_atom * (std::size_t)bspline_[ipanel_-1].Nspline())
                {
                    // use that older solution as a constraint for the current solution (see further down)
                }
                else
                {
                    // inform user that there is already some incompatible solution file
                    std::cout << "Warning: Solution for initial state (" << ni << "," << li << "," << mi << "), S = " << Spin << ", Etot = " << inp_.Etot[ie]
                            << " found, but has a smaller block size (" << size << " < " << Nspline_atom * Nspline_proj << ") and will be recomputed." << std::endl;
                }
            }
            
            // solution is larger
            if (size > (std::size_t)Nspline_atom * (std::size_t)Nspline_proj)
            {
                // is the size equal to some higher panel basis size?
                bool higher_panel_match = false;
                for (int ipanel = ipanel_ + 1; ipanel < cmd_.panels; ipanel++)
                {
                    if (size == (std::size_t)Nspline_atom * (std::size_t)bspline_[ipanel].Nspline())
                        higher_panel_match = true;
                }
                
                // consider this solution all right if sizes match for some further panel
                if (higher_panel_match)
                    continue;
                
                // inform user that there is already some incompatible solution file
                std::cout << "Warning: Solution for initial state (" << ni << "," << li << "," << mi << "), S = " << Spin << ", Etot = " << inp_.Etot[ie]
                          << " found, but has a larger block size (" << size << " > " << Nspline_atom * Nspline_proj << ") and will be recomputed." << std::endl;
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
        if (not (E == 0.5 * inp_.Etot[ie]))
            prec_->update(E = 0.5 * inp_.Etot[ie]);
        
        // for all initial states
        for (auto workitem : work)
        {
            // decode initial state
            int instate = std::get<0>(workitem);
            int Spin = std::get<1>(workitem);
            int ni = std::get<0>(inp_.instates[instate]);
            int li = std::get<1>(inp_.instates[instate]);
            int mi = std::get<2>(inp_.instates[instate]);
            
            // create right hand side
            BlockArray<Complex> chi (angs_.size(), !cmd_.outofcore, "cg-b");
            if (not cmd_.cont)
            {
                std::cout << "\tCreate right-hand side for initial state " << Hydrogen::stateName(ni,li,mi) << " and total spin S = " << Spin << " ... " << std::flush;
                prec_->rhs(chi, ie, instate, Spin);
                std::cout << "ok" << std::endl;
            }
            
            // compute and check norm of the right hand side vector
            double chi_norm = compute_norm(chi);
            if (chi_norm == 0.)
            {
                // this should not happen, hopefully we already checked
                std::cout << "\t! Right-hand-side is zero, check L, Pi and nL." << std::endl;
                computations_done++;
                continue;
            }
            if (not std::isfinite(chi_norm))
            {
                // this is a numerical problem, probably in evaluation of special functions (P, j)
                std::cout << "\t! Right hand side has invalid norm (" << chi_norm << ")." << std::endl;
                computations_done++;
                continue;
            }
            
            // load / reset solver state
            if (cmd_.cont) CG_.recover(); else CG_.reset();
            
            // prepare solution vector
            BlockArray<Complex> psi (std::move(new_array(angs_.size(),"cg-x")));
            
            // update constraint, if any
            SolutionIO loader (inp_.L, Spin, inp_.Pi, ni, li, mi, inp_.Etot[ie], angs_);
            if (ipanel_ > 0 and loader.check() == (std::size_t)Nspline_atom * (std::size_t)bspline_[ipanel_-1].Nspline())
            {
                // previous solution
                BlockArray<Complex> prev (angs_.size());
                
                // previous projectile basis size
                std::size_t Nspline_proj_prev = bspline_[ipanel_-1].Nspline();
                
                // map data from solution of the previous parcel to the new solution vector
                for (int ill = 0; ill < (int)angs_.size(); ill++)
                {
                    // read block from disk
                    loader.load(prev, ill);
                    
                    // copy all segments
                    for (int i = 0; i < Nspline_atom; i++)
                    {
                        psi.setSegment
                        (
                            ill,                // block index
                            i * Nspline_proj,   // segment placement position
                            Nspline_proj_prev,  // number of elements to copy
                            prev.segment
                            (
                                ill,                    // block index
                                i * Nspline_proj_prev,  // segment read position
                                Nspline_proj_prev       // number of elements to copy
                            )()                         // <-- dereference from TmpNumberArray to ArrayView
                        );
                    }
                    
                    // unload block
                    prev[ill].drop();
                }
            }
            
            // launch the linear system solver
            unsigned max_iter = (inp_.maxell + 1) * Nspline_atom;
            std::cout << "\tStart linear solver with tolerance " << cmd_.itertol << " for initial state " << Hydrogen::stateName(ni,li,mi) << " and total spin S = " << Spin << "." << std::endl;
            std::cout << "\t   i | time        | residual        | min  max  avg  block precond. iter." << std::endl;
            unsigned iterations = CG_.solve
            (
                chi,                    // right-hand side
                psi,                    // on input, the initial guess, on return, the solution
                cmd_.itertol,           // requested precision, |A·x - b|² < ε·|b|²
                0,                      // minimal iteration count
                max_iter,               // maximal iteration count
                apply_preconditioner,   // preconditioner callback
                matrix_multiply,        // matrix multiplication callback
                true,                   // verbose output
                compute_norm,           // how to evaluate norm of an array
                scalar_product,         // how to calculate scalar product of two arrays
                axby_operation,         // ax + by
                new_array,              // recipe for creation of a new array
                constrain               // linear constraint on the solution
            );
            
            if (iterations >= max_iter)
                std::cout << "\tConvergence too slow... The saved solution will be probably non-converged." << std::endl;
            else
                std::cout << "\tSolution converged." << std::endl;
            
            // update progress
            iterations_done += iterations;
            computations_done++;
            
            // save solution to disk (if valid)
            SolutionIO reader (inp_.L, Spin, inp_.Pi, ni, li, mi, inp_.Etot[ie], angs_);
            if (std::isfinite(compute_norm(psi)))
            {
                for (unsigned ill = 0; ill < angs_.size(); ill++)
                {
                    if (par_.isMyGroupWork(ill) and par_.IamGroupMaster() and not reader.save(psi, ill))
                        HexException("Failed to save solution to disk - the data are lost!");
                }
            }
            
            // reset some one-solution command line flags
            cmd_.reuse_dia_blocks = false;
            cmd_.cont = false;
            
        } // end of For Spin, instate
        
    } // end of For ie = 0, ..., inp.Ei.size() - 1
    
    // wait for completition of all processes before next step
    par_.wait();
    
    std::cout << std::endl << "All solutions computed." << std::endl;
    if (computations_done > 0)
        std::cout << "\t(typically " << iterations_done / computations_done << " CG iterations per solution)" << std::endl;
}

void Solver::apply_preconditioner_ (BlockArray<Complex> const & r, BlockArray<Complex> & z) const
{
    // save state of the solver for possible restart
    if (cmd_.outofcore)
        CG_.dump();
    
    // MPI-distributed preconditioning
    prec_->precondition(r, z);
}

void Solver::matrix_multiply_ (BlockArray<Complex> const & p, BlockArray<Complex> & q) const
{
    // MPI-distributed multiplication
    prec_->multiply(p, q);
}

Complex Solver::scalar_product_ (BlockArray<Complex> const & x, BlockArray<Complex> const & y) const
{
    // compute node-local scalar product
    Complex prod = 0;
    
    // for all segments
    for (std::size_t i = 0; i < x.size(); i++) if (par_.isMyGroupWork(i))
    {
        if (not x.inmemory()) const_cast<BlockArray<Complex>&>(x).hdfload(i);
        if (not y.inmemory()) const_cast<BlockArray<Complex>&>(y).hdfload(i);
        
        prod += (x[i]|y[i]);
        
        if (not x.inmemory()) const_cast<BlockArray<Complex>&>(x)[i].drop();
        if (not y.inmemory()) const_cast<BlockArray<Complex>&>(y)[i].drop();
    }
    
    // colect products from other nodes
    par_.syncsum(&prod, 1);
    
    // return global scalar product
    return prod / double(cmd_.groupsize);
}

double Solver::compute_norm_ (BlockArray<Complex> const & r) const
{
    // compute node-local norm of 'r'
    double rnorm2 = 0;
    for (std::size_t i = 0; i < r.size(); i++) if (par_.isMyGroupWork(i))
    {
        if (not r.inmemory())
            const_cast<BlockArray<Complex>&>(r).hdfload(i);
        
        rnorm2 += r[i].sqrnorm();
        
        if (not r.inmemory())
            const_cast<BlockArray<Complex>&>(r)[i].drop();
    }
    
    // collect norms from other nodes
    par_.syncsum(&rnorm2, 1);
    
    // return global norm
    return std::sqrt(rnorm2 / par_.groupsize());
}

void Solver::axby_operation_ (Complex a, BlockArray<Complex> & x, Complex b, BlockArray<Complex> const & y) const
{
    // only references blocks that are local to this MPI node
    for (std::size_t i = 0; i < x.size(); i++) if (par_.isMyGroupWork(i))
    {
        if (not x.inmemory()) x.hdfload(i);
        if (not y.inmemory()) const_cast<BlockArray<Complex>&>(y).hdfload(i);
        
        for (std::size_t j = 0; j < x[i].size(); j++)
            x[i][j] = a * x[i][j] + b * y[i][j];
        
        if (not x.inmemory()) { x.hdfsave(i); x[i].drop(); }
        if (not y.inmemory()) const_cast<BlockArray<Complex>&>(y)[i].drop();
    }
}

BlockArray<Complex> Solver::new_array_ (std::size_t N, std::string name) const
{
    // create a new block array and initialize blocks local to this MPI node
    BlockArray<Complex> array (N, !cmd_.outofcore, name);
    
    // initialize all blocks ('resize' automatically zeroes added elements)
    for (std::size_t i = 0; i < N; i++) if (par_.isMyGroupWork(i))
    {
        // allocate memory
        array[i].resize(bspline_[0].Nspline() * bspline_[ipanel_].Nspline());
        
        // take care of out-of-core mode
        if (not array.inmemory())
        {
            // save to disk if we are not continuing an older calculation
            if (not cmd_.cont)
                array.hdfsave(i);
            
            // release memory
            array[i].drop();
        }
    }
    
    return array;
}

void Solver::constrain_ (BlockArray<Complex> & r) const
{
/*
    // do not constrain the first panel
    if (ipanel_ == 0)
        return;
    
    // get number of elements fixed from the previous panel
    std::size_t Nspline_fixed = Nspline_fixed_[ipanel_ - 1];
    
    // for all blocks
    for (int ill = 0; ill < (int)angs_.size(); ill++)
    {
        // load the block from disk, if needed
        if (not r.inmemory())
            r.hdfload(ill);
        
        // erase leading part of the residual vector
        cArrayView(r[ill], 0, Nspline_fixed).fill(0);
        
        // release the memory
        if (not r.inmemory())
            r[ill].drop();
    }
*/
}
