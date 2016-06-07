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

#include "hex-hydrogen.h"

#include "solver.h"

Solver::Solver
(
    CommandLine & cmd,
    InputFile const & inp,
    Parallel const & par,
    AngularBasis const & ang,
    Bspline const & bspline_inner,
    Bspline const & bspline_full
) : cmd_(cmd), inp_(inp), par_(par), ang_(ang),
    bspline_inner_(bspline_inner),
    bspline_full_ (bspline_full),
    prec_(nullptr)
{
    // nothing to do
}

void Solver::choose_preconditioner ()
{
    // create the preconditioner
    prec_ = Preconditioners::choose(par_, inp_, ang_, bspline_inner_, bspline_full_, cmd_);
    
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
    std::size_t Nspline_inner = bspline_inner_.Nspline();
    std::size_t Nspline_full  = bspline_full_ .Nspline();
    std::size_t Nspline_outer = Nspline_full - Nspline_inner;
    
    // print Hamiltonian size as a number with thousands separator (apostroph used)
    class MyNumPunct : public std::numpunct<char>
    {
        protected:
            virtual char do_thousands_sep() const { return '\''; }
            virtual std::string do_grouping() const { return "\03"; }
    };
    std::cout.imbue(std::locale(std::locale::classic(), new MyNumPunct));
    std::cout << "Inner region hamiltonian size: " << Nspline_inner * Nspline_inner * ang_.states().size() << std::endl;
    std::cout.imbue(std::locale::classic());
    
    // wrap member functions to lambda-functions for use in the CG solver
    auto apply_preconditioner = [&](BlockArray<Complex> const & r, BlockArray<Complex> & z) { this->apply_preconditioner_(r,z); };
    auto matrix_multiply      = [&](BlockArray<Complex> const & p, BlockArray<Complex> & q) { this->matrix_multiply_(p,q); };
    auto scalar_product       = [&](BlockArray<Complex> const & x, BlockArray<Complex> const & y) { return this->scalar_product_(x,y); };
    auto compute_norm         = [&](BlockArray<Complex> const & r) { return this->compute_norm_(r); };
    auto axby_operation       = [&](Complex a, BlockArray<Complex> & x, Complex b, BlockArray<Complex> const & y) { this->axby_operation_(a,x,b,y); };
    auto new_array            = [&](std::size_t N, std::string name) { return this->new_array_(N,name); };
    auto process_solution     = [&](unsigned iteration, BlockArray<Complex> const & x) { return this->process_solution_(iteration,x); };
    
    E_ = special::constant::Nan;
    int iterations_done = 0, computations_done = 0;
    
    for (unsigned ie = 0; ie < inp_.Etot.size(); ie++)
    {
        // print progress information
        std::cout << "\nSolving the system for Etot[" << ie << "] = " << inp_.Etot[ie] << " ("
                  << int(std::trunc(ie * 100. / inp_.Etot.size() + 0.5)) << " % finished, typically "
                  << (computations_done == 0 ? 0 : iterations_done / computations_done)
                  << " CG iterations per energy)" << std::endl;
        
        // check applicability of the projectile basis extension
        if (not inp_.inner_only and inp_.Etot[ie] >= 0)
            HexException("Projectile basis extension cannot be used for energies above ionization.");
        
        // get maximal asymptotic principal quantum number
        int max_n = (inp_.Etot[ie] >= 0 ? 0 : 1.0 / std::sqrt(-inp_.Etot[ie]));
        
        // get asymptotical bound states for each of the angular momentum pairs
        for (unsigned ill = 0; ill < ang_.states().size(); ill++)
        {
            int l1 = ang_.states()[ill].first;
            int l2 = ang_.states()[ill].second;
            
            // store all principal quantum numbers for given l₁ or l₂ and the total energy
            bstates_.push_back
            (
                std::make_pair
                (
                    l1 + 1 > max_n ? iArray{} : linspace(l1 + 1, max_n, max_n - l1),
                    l2 + 1 > max_n ? iArray{} : linspace(l2 + 1, max_n, max_n - l2)
                )
            );
        }
        
        // print system information
        std::cout.imbue(std::locale(std::locale::classic(), new MyNumPunct));
        std::cout << "\tFull hamiltonian / solution size: " << std::accumulate
        (
            bstates_.begin(), bstates_.end(),
            Nspline_inner * Nspline_inner * ang_.states().size(),
            [&](std::size_t n, std::pair<iArray,iArray> const & p) { return n + (p.first.size() + p.second.size()) * Nspline_outer; }
        ) << std::endl;
        std::cout.imbue(std::locale::classic());
        
        // we may have already computed all solutions for this energy... is it so?
        std::vector<std::pair<int,int>> work;
        for (unsigned instate = 0; instate < inp_.instates.size(); instate++)
        for (unsigned Spin : inp_.Spin)
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
            for (unsigned ill = 0; ill < ang_.states().size(); ill++) if (ang_.states()[ill].first == li)
            {
                // get partia wave
                int l = ang_.states()[ill].second;
                
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
            /*std::size_t size = reader.check();*/
            
            // TODO : solution has the expected size
            /*if (size == ? and not cmd_.refine_solution)
            {
                std::cout << "\tSolution for initial state " << Hydrogen::stateName(ni,li,mi) << " (S = " << Spin << ") found." << std::endl;
                continue;
            }*/
            
            // TODO : solution is smaller
            /*if (0 < size and size < (std::size_t)Nspline_atom * (std::size_t)bspline_full_[ipanel_].Nspline())
            {
                // is the size equal to previous panel basis size?
                if (ipanel_ > 0 and size == (std::size_t)Nspline_atom * (std::size_t)bspline_full_[ipanel_-1].Nspline())
                {
                    // use that older solution in construction of the new one (later)
                }
                else if (ipanel_ > 0)
                {
                    // inform user that there is already some incompatible solution file
                    std::cout << "Warning: Solution for initial state (" << ni << "," << li << "," << mi << "), S = " << Spin << ", Etot = " << inp_.Etot[ie]
                              << " found, but has a smaller block size (" << size << " < " << (std::size_t)Nspline_atom * bspline_full_[ipanel_-1].Nspline() << ") and will be recomputed." << std::endl;
                }
                else
                {
                    // the same for first panel
                    std::cout << "Warning: Solution for initial state (" << ni << "," << li << "," << mi << "), S = " << Spin << ", Etot = " << inp_.Etot[ie]
                              << " found, but has a smaller block size (" << size << " < " << (std::size_t)Nspline_atom * (std::size_t)bspline_full_[ipanel_].Nspline() << ") and will be recomputed." << std::endl;
                }
            }*/
            
            // TODO : solution is larger
            /*if (size > (std::size_t)Nspline_atom * (std::size_t)bspline_full_[ipanel_].Nspline())
            {
                // is the size equal to some higher panel basis size?
                bool higher_panel_match = false;
                for (int ipanel = ipanel_ + 1; ipanel < cmd_.panels; ipanel++)
                {
                    if (size == (std::size_t)Nspline_atom * (std::size_t)bspline_full_[ipanel].Nspline())
                        higher_panel_match = true;
                }
                
                // consider this solution all right if sizes match for some further panel
                if (higher_panel_match)
                {
                    std::cout << "\tSolution for initial state " << Hydrogen::stateName(ni,li,mi) << " (S = " << Spin << ") present and propagated." << std::endl;
                    continue;
                }
                
                // inform user that there is already some incompatible solution file
                std::cout << "Warning: Solution for initial state (" << ni << "," << li << "," << mi << "), S = " << Spin << ", Etot = " << inp_.Etot[ie]
                          << " found, but has a larger block size (" << size << " > " << (std::size_t)Nspline_atom * Nspline_proj << ") and will be recomputed." << std::endl;
            }*/
            
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
        if (not (E_ == 0.5 * inp_.Etot[ie]))
            prec_->update(E_ = 0.5 * inp_.Etot[ie]);
        
        // for all initial states
        for (auto workitem : work)
        {
            // decode initial state
            int instate = std::get<0>(workitem);
            ang_.S() = std::get<1>(workitem);
            ni_ = std::get<0>(inp_.instates[instate]);
            li_ = std::get<1>(inp_.instates[instate]);
            mi_ = std::get<2>(inp_.instates[instate]);
            
            // create right hand side
            BlockArray<Complex> chi (ang_.states().size(), !cmd_.outofcore, "cg-b");
            if (not cmd_.cont)
            {
                Timer t;
                std::cout << "\tCreate right-hand side for initial state " << Hydrogen::stateName(ni_,li_,mi_) << " and total spin S = " << ang_.S() << " ... " << std::flush;
                
                // use the preconditioner setup routine
                prec_->rhs(chi, ie, instate);
                
                std::cout << "done after " << t.nice_time() << std::endl;
            }
            
            // compute and check norm of the right hand side vector
            Real chi_norm = compute_norm(chi);
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
            if (cmd_.cont /*or cmd_.refine_solution*/) CG_.recover(); else CG_.reset();
            
            // prepare solution vector
            BlockArray<Complex> psi (std::move(new_array(ang_.states().size(),"cg-x")));
            
            // load initial guess
            if (not cmd_.cont and ie > 0 and cmd_.carry_initial_guess)
            {
                SolutionIO prev_sol_reader (ang_.L(), ang_.S(), ang_.Pi(), ni_, li_, mi_, inp_.Etot[ie-1], ang_.states());
                prev_sol_reader.load(psi);
            }
            if (not cmd_.cont and cmd_.refine_solution)
            {
                SolutionIO prev_sol_reader (ang_.L(), ang_.S(), ang_.Pi(), ni_, li_, mi_, inp_.Etot[ie], ang_.states());
                prev_sol_reader.load(psi);
            }
            
            // launch the linear system solver
            Timer t;
            unsigned max_iter = (inp_.maxell + 1) * (std::size_t)Nspline_inner;
            std::cout << "\tStart linear solver with tolerance " << cmd_.itertol << " for initial state " << Hydrogen::stateName(ni_,li_,mi_) << " and total spin S = " << ang_.S() << "." << std::endl;
            std::cout << "\t   i | time        | residual        | min  max  avg  block precond. iter." << std::endl;
            CG_.apply_preconditioner = apply_preconditioner;
            CG_.matrix_multiply      = matrix_multiply;
            CG_.verbose              = true;
            CG_.compute_norm         = compute_norm;
            CG_.scalar_product       = scalar_product;
            CG_.axby                 = axby_operation;
            CG_.new_array            = new_array;
            CG_.process_solution     = process_solution;
            unsigned iterations = CG_.solve(chi, psi, cmd_.itertol, 0, max_iter);
            
            if (iterations >= max_iter)
                std::cout << "\tConvergence too slow... The saved solution will be probably non-converged." << std::endl;
            else
                std::cout << "\tSolution converged after " << t.nice_time() << "." << std::endl;
            
            // update progress
            iterations_done += iterations;
            computations_done++;
            
            // save solution to disk (if valid)
            SolutionIO reader (ang_.L(), ang_.S(), ang_.Pi(), ni_, li_, mi_, 2 * E_, ang_.states());
            if (std::isfinite(compute_norm(psi)))
            {
                for (unsigned ill = 0; ill < ang_.states().size(); ill++)
                {
                    if (par_.isMyGroupWork(ill) and par_.IamGroupMaster())
                    {
                        if (not psi.inmemory())
                            psi.hdfload(ill);
                        
                        // write the updated solution to disk
                        if (not reader.save(psi, ill))
                            HexException("Failed to save solution to disk - the data are lost!");
                        
                        if (not psi.inmemory())
                            psi[ill].drop();
                    }
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
    return prod / Real(cmd_.groupsize);
}

Real Solver::compute_norm_ (BlockArray<Complex> const & r) const
{
    // compute node-local norm of 'r'
    Real rnorm2 = 0;
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
    // just a quick check: N must be equal to the number of blocks
    if (N != ang_.states().size())
        HexException("Runtime error: Wrong number of angular blocks: %d != %d.", N, ang_.states().size());
    
    // create a new block array and initialize blocks local to this MPI node
    BlockArray<Complex> array (N, !cmd_.outofcore, name);
    
    // sizes
    std::size_t Nspline_inner = bspline_inner_.Nspline();
    std::size_t Nspline_full = bspline_full_.Nspline();
    std::size_t Nspline_outer = Nspline_full - Nspline_inner;
    
    // initialize all blocks ('resize' automatically zeroes added elements)
    for (std::size_t i = 0; i < N; i++) if (par_.isMyGroupWork(i))
    {
        // get correct size of the block depending on the calculation mode (with/without projectile basis extension)
        std::size_t size = Nspline_inner * Nspline_inner;
        if (not inp_.inner_only)
            size += (bstates_[i].first.size() + bstates_[i].second.size()) * Nspline_outer;
        
        // allocate memory
        array[i].resize(size);
        
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

void Solver::process_solution_ (unsigned iteration, BlockArray<Complex> const & x) const
{
    if (cmd_.write_intermediate_solutions)
    {
        SolutionIO writer (ang_.L(), ang_.S(), ang_.Pi(), ni_, li_, mi_, 2 * E_, ang_.states(), format("tmp-%d", iteration));
        writer.save(x);
    }
}

void Solver::finish ()
{
    prec_->finish();
}
