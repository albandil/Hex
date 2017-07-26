//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2017, Jakub Benda, Charles University in Prague                    //
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
#include "hex-misc.h"

// --------------------------------------------------------------------------------- //

#include "amplitudes.h"
#include "solver.h"

// --------------------------------------------------------------------------------- //

Solver::Solver
(
    CommandLine        & cmd,
    InputFile    const & inp,
    Parallel     const & par,
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
    prec_ = PreconditionerBase::Choose
    (
        cmd_.preconditioner,
        cmd_, inp_,  par_, ang_,
        bspline_inner_, // inner region basis
        bspline_full_,  // full basis
        bspline_full_,  // x axis full basis
        bspline_full_   // y axis full basis
    );
    
    // check success
    if (prec_ == nullptr)
        HexException("Preconditioner %s not implemented.", cmd_.preconditioner.c_str());
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
        
        // get asymptotical bound states for each of the angular momentum pairs
        channels_.resize(ang_.states().size());
        bstates_.clear();
        for (unsigned ill = 0; ill < ang_.states().size(); ill++)
        {
            int l1 = ang_.states()[ill].first;
            int l2 = ang_.states()[ill].second;
            
            // get number of bound states for each particle at this energy
            channels_[ill] = prec_->bstates(std::max(inp_.Etot[ie], inp_.channel_max_E), l1, l2);
            
            // the number of scattering channels
            std::swap(channels_[ill].first, channels_[ill].second);
            
            // store all principal quantum numbers for given l1 or l2
            bstates_.push_back
            (
                std::make_pair
                (
                    linspace<int>(l1 + 1, l1 + channels_[ill].second, channels_[ill].second),
                    linspace<int>(l2 + 1, l2 + channels_[ill].first,  channels_[ill].first)
                )
            );
        }
        
        // calculate size of the hamiltonian
        std::size_t Hsize = Nspline_inner * Nspline_inner * ang_.states().size();
        for (std::pair<iArray,iArray> const & p : bstates_)
            Hsize += (p.first.size() + p.second.size()) * Nspline_outer;
        
        // print system information
        std::cout.imbue(std::locale(std::locale::classic(), new MyNumPunct));
        std::cout << "\tFull hamiltonian / solution size: " << Hsize << std::endl;
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
            SolutionIO reader (inp_.L, Spin, inp_.Pi, ni, li, mi, inp_.Etot[ie], ang_.states(), channels_);
            std::size_t size = 0;
            reader.check(SolutionIO::All, size);
            
            // sum partial sizes in distributed case
            if (not cmd_.shared_scratch)
            {
                par_.mastersum(&size, 1, 0);
                par_.bcast(0, &size, 1);
            }
            
            // solution has the expected size
            if (size == Hsize and not cmd_.refine_solution)
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
            SolutionIO reader (ang_.L(), ang_.S(), ang_.Pi(), ni_, li_, mi_, 2 * E_, ang_.states(), channels_);
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
    // compute scalar product
    Complex prod = 0;
    
    // make sure no process is playing with the data
    par_.wait();
    
    // for all segments owned by current process' group
    for (std::size_t i = 0; i < x.size(); i++) if (par_.isMyGroupWork(i))
    {
        // load segments, if needed
        if (not x.inmemory()) const_cast<BlockArray<Complex>&>(x).hdfload(i);
        if (not y.inmemory()) const_cast<BlockArray<Complex>&>(y).hdfload(i);
        
        // choose a chunk for this member of the group
        std::size_t jmin = x[i].size() * (par_.igroupproc()    ) / par_.groupsize();
        std::size_t jmax = x[i].size() * (par_.igroupproc() + 1) / par_.groupsize();
        
        // update the scalar product
        for (std::size_t j = jmin; j < x[i].size() and j < jmax; j++)
            prod += x[i][j] * y[i][j];
        
        // release memory
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
    // compute norm
    Real rnorm2 = 0;
    
    // make sure no process is playing with the data
    par_.wait();
    
    // compute node-local norm of 'r'
    for (std::size_t i = 0; i < r.size(); i++) if (par_.isMyGroupWork(i))
    {
        // load segment, if needed
        if (not r.inmemory()) const_cast<BlockArray<Complex>&>(r).hdfload(i);
        
        // choose a chunk for this member of the group
        std::size_t jmin = r[i].size() * (par_.igroupproc()    ) / par_.groupsize();
        std::size_t jmax = r[i].size() * (par_.igroupproc() + 1) / par_.groupsize();
        
        // update the square norm
        for (std::size_t j = jmin; j < r[i].size() and j < jmax; j++)
            rnorm2 += sqrabs(r[i][j]);
        
        // release memory
        if (not r.inmemory()) const_cast<BlockArray<Complex>&>(r)[i].drop();
    }
    
    // collect norms from other nodes
    par_.syncsum(&rnorm2, 1);
    
    // return global norm
    return std::sqrt(rnorm2 / par_.groupsize());
}

void Solver::axby_operation_ (Complex a, BlockArray<Complex> & x, Complex b, BlockArray<Complex> const & y) const
{
    // make sure no process is playing with the data
    par_.wait();
    
    // only references blocks that are local to this MPI node
    for (std::size_t i = 0; i < x.size(); i++) if (par_.isMyGroupWork(i))
    {
        // load segments, if needed
        if (not x.inmemory()) x.hdfload(i);
        if (not y.inmemory()) const_cast<BlockArray<Complex>&>(y).hdfload(i);
        
        // choose a chunk for this member of the group
        std::size_t jmin = x[i].size() * (par_.igroupproc()    ) / par_.groupsize();
        std::size_t jmax = x[i].size() * (par_.igroupproc() + 1) / par_.groupsize();
        
        // update the linear combination
        for (std::size_t j = jmin; j < jmax and j < x[i].size(); j++)
            x[i][j] = a * x[i][j] + b * y[i][j];
        
        // synchronize across the group
        for (int inode = 0; inode < par_.groupsize(); inode++)
        {
            std::size_t jmin = x[i].size() * (inode    ) / par_.groupsize();
            std::size_t jmax = x[i].size() * (inode + 1) / par_.groupsize();
            
            jmax = std::min(x[i].size(), jmax);
            
            par_.bcast_g(par_.igroup(), inode, x[i].data() + jmin, jmax - jmin);
        }
        
        // release memory (groupmaster optionally saves the data to disk)
        if (not x.inmemory()) { if (par_.IamGroupMaster()) x.hdfsave(i); x[i].drop(); }
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
    std::string dir = format("iter-%d", iteration);
    
    if (cmd_.write_intermediate_solutions or cmd_.runtime_postprocess)
    {
        create_directory(dir);
        
        // write the solution
        SolutionIO writer (ang_.L(), ang_.S(), ang_.Pi(), ni_, li_, mi_, 2 * E_, ang_.states(), channels_, dir + "/psi");
        writer.save(x);
    }
    
    if (cmd_.runtime_postprocess)
    {
        create_directory(dir);
        
        // extract amplitudes
        Amplitudes ampl (bspline_inner_, bspline_full_, inp_, par_, cmd_, ang_.states());
        ampl.verbose(false);
        ampl.extract(dir);
        ampl.writeSQL_files(dir);
        ampl.writeICS_files(dir);
    }
}

void Solver::finish ()
{
    prec_->finish();
    delete prec_;
    prec_ = nullptr;
}
