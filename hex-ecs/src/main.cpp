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

#include <cstdio>
#include <cmath>
#include <csignal>
#include <vector>
#include <iostream>
#include <string>
#include <tuple>

#include <gsl/gsl_errno.h>

#ifdef _OPENMP
    #include <omp.h>
#endif

#include "amplitudes.h"
#include "arrays.h"
#include "bspline.h"
#include "hydrogen.h"
#include "io.h"
#include "itersolve.h"
#include "misc.h"
#include "parallel.h"
#include "preconditioners.h"
#include "radial.h"
#include "version.h"

ConjugateGradients < cBlockArray, cBlockArray& > CG;

void signal_handler (int param)
{
    // save solver state
    CG.dump();
    
    // terminate the program
    std::exit(EXIT_FAILURE);
}

int main (int argc, char* argv[])
{
// ------------------------------------------------------------------------- //
//                                                                           //
// Preparations                                                              //
//                                                                           //
// ------------------------------------------------------------------------- //
    
    //
    // Program initialization ---------------------------------------------- //
    //
    
    // display logo
    std::cout << logo(" ") << std::endl;
    std::cout << "=== Exterior complex scaling in B-splines ===" << std::endl << std::endl;
    
    // echo command line
    std::cout << "Command line used" << std::endl;
    std::cout << "\t";
    for (int iarg = 0; iarg < argc; iarg++)
        std::cout << argv[iarg] << " ";
    std::cout << std::endl << std::endl;
    
    // turn off GSL exceptions
    gsl_set_error_handler_off();
    
    // error handler
    std::signal(SIGINT, &signal_handler);
    
    // disable buffering of the standard output (-> immediate logging)
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    
    // get input from command line
    CommandLine cmd (argc, argv);
    
    // check some exclusive options
    if (cmd.parallel_block and cmd.lightweight_radial_cache)
        HexException("The options --parallel-block and --lightweight-radial-cache/--lightweight-full can't be used together because of different multiplication scheme.");
    
#ifdef _OPENMP
    // set OpenMP parallel nesting (avoid oversubscription)
    // - disable for concurrent diagonal block preconditioning
    // - enable for sequential diagonal block preconditioning
    omp_set_nested(!cmd.parallel_block);
    # pragma omp parallel
    # pragma omp master
    {
        std::cout << "OpenMP environment" << std::endl;
        std::cout << "\tthreads: " << omp_get_num_threads() << std::endl;
        std::cout << "\tnesting: " << (omp_get_nested() ? "on" : "off") << std::endl;
    }
#endif
    
    // setup MPI
    Parallel par (&argc, &argv, cmd.parallel);
    
    // check input file
    if (not cmd.inputfile.is_open())
    {
        cmd.inputfile.open("ecs.inp");
        if (not cmd.inputfile.good())
        {
            std::cout << "Cannot open the file \"ecs.inp\"." << std::endl;
            std::cout << std::endl;
            std::cout << "Either (1) provide input settings in the file \"ecs.inp\", " << std::endl;
            std::cout << "    or (2) give another name using the '--input' command line option." << std::endl;
            std::cout << std::endl;
            std::exit(EXIT_FAILURE);
        };
    }
    
    // get input from input file
    InputFile inp (cmd.inputfile);
    
    // is there something to compute?
    if (inp.Etot.empty() or inp.instates.empty() or inp.outstates.empty())
    {
        std::cout << "Nothing to compute." << std::endl;
        std::exit(EXIT_SUCCESS);
    }
    
    // 
    // Setup B-spline environment ------------------------------------------ //
    //
    
    double R0 = inp.rknots.back();       // end of real grid
    double Rmax = inp.cknots.back();     // end of complex grid
    
    // check ordering
    if (R0 >= Rmax)
        HexException("ERROR: Rmax = %g (end of grid) must be greater than R0 = %g (end of real grid)!", Rmax, R0);
    
    // create B-splines
    Bspline bspline (inp.order, inp.rknots, inp.ecstheta, inp.cknots);
    
    // shortcuts
    int Nspline = bspline.Nspline();
    
    // info
    std::cout << "B-spline count: " << Nspline << std::endl << std::endl;
    
    //
    // Setup angular data -------------------------------------------------- //
    //
    
    std::cout << "Setting up the coupled angular states..." << std::endl;
    
    // coupled angular momentum pairs
    std::vector<std::pair<int,int>> coupled_states;
    
    // for given L, Π and levels list all available (ℓ₁ℓ₂) pairs
    for (int ell = 0; ell <= inp.levels; ell++)
    {
        std::cout << "\t-> [" << ell << "] ";
        
        // get sum of the angular momenta for this angular level
        int sum = 2 * ell + inp.L + inp.Pi;
        
        // for all angular momentum pairs that do compose L
        for (int l1 = ell; l1 <= sum - ell; l1++)
        {
            std::cout << "(" << l1 << "," << sum - l1 << ") ";
            coupled_states.push_back(std::make_pair(l1, sum - l1));
        }
        std::cout << std::endl;
    }
    
    std::cout << "\t-> The matrix of the set contains " << coupled_states.size()
              << " diagonal blocks." << std::endl;
    
    // skip if there is nothing to compute
    if (coupled_states.empty())
        std::exit(EXIT_SUCCESS);
    
    std::cout << std::endl;
    
    //
    // Zip solution file into VTK geometry if told so
    //
    
    if (cmd.zipfile.size() != 0 and par.IamMaster())
    {
        zip_solution (cmd, bspline, coupled_states);
        std::cout << std::endl << "Done." << std::endl << std::endl;
        return 0;
    }
    
    //
    // Create and the preconditioner.
    //
    
    PreconditionerBase * prec = Preconditioners::choose (par, inp, coupled_states, bspline, cmd);
    if (prec == nullptr)
        HexException("Preconditioner %d not implemented.", cmd.preconditioner);
    
// ------------------------------------------------------------------------- //
//                                                                           //
// StgRadial                                                                 //
//                                                                           //
// ------------------------------------------------------------------------- //
    
if (cmd.itinerary & CommandLine::StgRadial)
{
    // setup the preconditioner (compute radial integrals etc.)
    prec->setup();
}
else
{
    std::cout << "Skipped computation of radial integrals." << std::endl;
}

// ------------------------------------------------------------------------- //
//                                                                           //
// StgSolve                                                                  //
//                                                                           //
// ------------------------------------------------------------------------- //
    
if (cmd.itinerary & CommandLine::StgSolve)
{
    // CG preconditioner callback
    auto apply_preconditioner = [ & ](BlockArray<Complex> const & r, BlockArray<Complex> & z) -> void
    {
        // MPI-distributed preconditioning
        prec->precondition(r, z);
    };
    
    // CG matrix multiplication callback
    auto matrix_multiply = [ & ](BlockArray<Complex> const & p, BlockArray<Complex> & q) -> void
    {
        // MPI-distributed multiplication
        prec->multiply(p, q);
    };
    
    // CG scalar product function callback
    auto scalar_product = [ & ](BlockArray<Complex> const & x, BlockArray<Complex> const & y) -> Complex
    {
        // compute node-local scalar product
        Complex prod = 0;
        
        // for all segments
        for (std::size_t i = 0; i < x.size(); i++) if (par.isMyWork(i))
        {
            if (not x.inmemory()) const_cast<BlockArray<Complex>&>(x).hdfload(i);
            if (not y.inmemory()) const_cast<BlockArray<Complex>&>(y).hdfload(i);
            
            prod += (x[i]|y[i]);
            
            if (not x.inmemory()) const_cast<BlockArray<Complex>&>(x)[i].drop();
            if (not y.inmemory()) const_cast<BlockArray<Complex>&>(y)[i].drop();
        }
        
        // colect products from other nodes
        par.syncsum(&prod, 1);
        
        // return global scalar product
        return prod;
    };
    
    // CG norm function that broadcasts master's result to all nodes
    auto compute_norm = [ & ](BlockArray<Complex> const & r) -> double
    {
        // compute node-local norm of 'r'
        double rnorm2 = 0;
        for (std::size_t i = 0; i < r.size(); i++) if (par.isMyWork(i))
        {
            if (not r.inmemory())
                const_cast<BlockArray<Complex>&>(r).hdfload(i);
            
            rnorm2 += r[i].sqrnorm();
            
            if (not r.inmemory())
                const_cast<BlockArray<Complex>&>(r)[i].drop();
        }
        
        // collect norms from other nodes
        par.syncsum(&rnorm2, 1);
        
        // return global norm
        return std::sqrt(rnorm2);
    };
    
    // CG linear combination
    auto axby_operation = [ & ](Complex a, BlockArray<Complex> & x, Complex b, BlockArray<Complex> const & y) -> void
    {
        // only references blocks that are local to this MPI node
        for (std::size_t i = 0; i < x.size(); i++) if (par.isMyWork(i))
        {
            if (not x.inmemory()) x.hdfload(i);
            if (not y.inmemory()) const_cast<BlockArray<Complex>&>(y).hdfload(i);
            
            for (std::size_t j = 0; j < x[i].size(); j++)
                x[i][j] = a * x[i][j] + b * y[i][j];
            
            if (not x.inmemory()) { x.hdfsave(i); x[i].drop(); }
            if (not y.inmemory()) const_cast<BlockArray<Complex>&>(y)[i].drop();
        }
    };
    
    // CG new array
    auto new_array = [ & ](std::size_t N, std::string name) -> BlockArray<Complex>
    {
        // create a new block array and initialize blocks local to this MPI node
        BlockArray<Complex> array (N, !cmd.outofcore, name);
        
        // initialize all blocks ('resize' automatically zeroes added elements)
        for (std::size_t i = 0; i < N; i++) if (par.isMyWork(i))
        {
            array[i].resize(Nspline * Nspline);
            
            if (not array.inmemory())
            {
                if (not cmd.cont)
                    array.hdfsave(i);
                array[i].drop();
            }
        }
        
        return array;
    };
    
    //
    // Solve the equations
    //
    
    std::cout << "Hamiltonian size: " << Nspline * Nspline * coupled_states.size() << std::endl;
    double E = special::constant::Nan;
    int iterations_done = 0, computations_done = 0;
    for (unsigned ie = 0; ie < inp.Etot.size(); ie++)
    {
        // print progress information
        std::cout << "\nSolving the system for Etot[" << ie << "] = " << inp.Etot[ie] << " ("
                  << int(trunc(ie * 100. / inp.Etot.size() + 0.5)) << " % finished, typically "
                  << (computations_done == 0 ? 0 : iterations_done / computations_done)
                  << " CG iterations per energy)" << std::endl;
        
        // we may have already computed all solutions for this energy... is it so?
        std::vector<std::pair<int,int>> work;
        for (unsigned instate = 0; instate < inp.instates.size(); instate++)
        for (unsigned Spin = 0; Spin <= 1; Spin++)
        {
            // decode initial state
            int ni = std::get<0>(inp.instates[instate]);
            int li = std::get<1>(inp.instates[instate]);
            int mi = std::get<2>(inp.instates[instate]);
            
            // skip energy-forbidden states
            if (inp.Etot[ie] < -1./(ni*ni))
            {
                std::cout << "\tSkip initial state " << Hydrogen::stateName(ni,li,mi) << " (S = " << Spin
                          << ") : not allowed by total E." << std::endl;;
                continue;
            }
            
            // check if the right hand side will be zero for this instate
            bool allowed = false;
            for (int l = std::abs(li - inp.L); l <= li + inp.L; l++)
            {
                // does this combination conserve parity?
                if ((inp.L + li + l) % 2 != inp.Pi)
                    continue;
                
                // does this combination have valid 'mi' for this partial wave?
                if (special::ClebschGordan(li,mi,l,0,inp.L,mi) != 0)
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
            SolutionIO reader (inp.L, Spin, inp.Pi, ni, li, mi, inp.Etot[ie], coupled_states, Nspline);
            if (reader.check())
                continue;
            
            // add work
            work.push_back(std::make_pair(instate,Spin));
        }
        
        // skip this energy if nothing to compute
        if (work.empty())
        {
            std::cout << "\tAll solutions for Etot[" << ie << "] = " << inp.Etot[ie] << " loaded." << std::endl;
            continue;
        }
        
        // update the preconditioner, if this is the first energy to compute or it changed from previous iteration
        if (not (E == 0.5 * inp.Etot[ie]))
            prec->update(E = 0.5 * inp.Etot[ie]);
        
        // for all initial states
        for (auto workitem : work)
        {
            // decode initial state
            int instate = std::get<0>(workitem);
            int Spin = std::get<1>(workitem);
            int ni = std::get<0>(inp.instates[instate]);
            int li = std::get<1>(inp.instates[instate]);
            int mi = std::get<2>(inp.instates[instate]);
            
            // create right hand side
            BlockArray<Complex> chi (coupled_states.size(), !cmd.outofcore, "cg-b");
            if (not cmd.cont)
            {
                std::cout << "\tCreate right-hand side for initial state " << Hydrogen::stateName(ni,li,mi) << " and total spin S = " << Spin << " ... " << std::flush;
                prec->rhs(chi, ie, instate, Spin);
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
            
            // load solver state
            if (cmd.cont)
                CG.recover();
            
            // prepare solution vector
            BlockArray<Complex> psi (std::move(new_array(coupled_states.size(),"cg-x")));
            
            // launch the linear system solver
            unsigned max_iter = (inp.maxell + 1) * Nspline;
            std::cout << "\tStart linear solver with tolerance " << cmd.itertol << " for initial state " << Hydrogen::stateName(ni,li,mi) << " and total spin S = " << Spin << "." << std::endl;
            std::cout << "\t   i | time        | residual        | min  max  avg  block precond. iter." << std::endl;
            unsigned iterations = CG.solve
            (
                chi,                    // right-hand side
                psi,                    // on input, the initial guess, on return, the solution
                cmd.itertol,            // requested precision, |A·x - b|² < ε·|b|²
                0,                      // minimal iteration count
                max_iter,               // maximal iteration count
                apply_preconditioner,   // preconditioner callback
                matrix_multiply,        // matrix multiplication callback
                true,                   // verbose output
                compute_norm,           // how to evaluate norm of an array
                scalar_product,         // how to calculate scalar product of two arrays
                axby_operation,         // ax + by
                new_array               // recipe for creation of a new array
            );
            
            if (iterations >= max_iter)
                std::cout << "\tConvergence too slow... The saved solution will be probably non-converged." << std::endl;
            else
                std::cout << "\tSolution converged." << std::endl;
            
            // update progress
            iterations_done += iterations;
            computations_done++;
            
            // save solution to disk (if valid)
            SolutionIO reader (inp.L, Spin, inp.Pi, ni, li, mi, inp.Etot[ie], coupled_states, Nspline);
            if (std::isfinite(compute_norm(psi)))
            {
                for (unsigned ill = 0; ill < coupled_states.size(); ill++)
                    if (par.isMyWork(ill))
                        reader.save(psi, ill);
            }
            
            // reset some one-solution command line flags
            cmd.reuse_dia_blocks = false;
            cmd.cont = false;
            
        } // end of For Spin, instate
        
    } // end of For ie = 0, ..., inp.Ei.size() - 1
    
    std::cout << std::endl << "All solutions computed." << std::endl;
    if (computations_done > 0)
        std::cout << "\t(typically " << iterations_done / computations_done << " CG iterations per solution)" << std::endl;
}
else // i.e. no StgSolve
{
    std::cout << "Skipped solution of the equation." << std::endl;
}

// ------------------------------------------------------------------------- //
//                                                                           //
// StgExtract                                                                //
//                                                                           //
// ------------------------------------------------------------------------- //

if (cmd.itinerary & CommandLine::StgExtract)
{
    // extract amplitudes
    Amplitudes ampl (bspline, inp, par, cmd, coupled_states);
    ampl.extract();
    
    // write T-matrices to a text file as SQL statements 
    ampl.writeSQL_files();
    
    // write integral cross sections to a text file
    ampl.writeICS_files();
}
else
{
    std::cout << "Skipped extraction of amplitudes." << std::endl;
}

// ------------------------------------------------------------------------- //
//                                                                           //
// End                                                                       //
//                                                                           //
// ------------------------------------------------------------------------- //

    std::cout << std::endl << "Done." << std::endl << std::endl;
    
    return 0;
}
