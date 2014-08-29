/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2014                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cstdio>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <tuple>

#include <gsl/gsl_errno.h>

#include "amplitudes.h"
#include "arrays.h"
#include "bspline.h"
#include "complex.h"
#include "io.h"
#include "itersolve.h"
#include "misc.h"
#include "parallel.h"
#include "preconditioners.h"
#include "radial.h"
#include "version.h"

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
    std::cout << "Command line used:" << std::endl;
    std::cout << "\t";
    for (int iarg = 0; iarg < argc; iarg++)
        std::cout << argv[iarg] << " ";
    std::cout << std::endl << std::endl;
    
    // turn off GSL and HDF exceptions
    gsl_set_error_handler_off();
    H5::Exception::dontPrint();
    
    // disable buffering of the standard output (-> immediate logging)
    setvbuf(stdout, nullptr, _IONBF, 0);
    
    // get input from command line
    CommandLine cmd (argc, argv);
    
    // setup MPI
    Parallel par (&argc, &argv, cmd.parallel);
    
    // check input file
    if (not cmd.inputfile.is_open())
    {
        cmd.inputfile.open("ecs.inp");
        if (not cmd.inputfile.good())
            throw exception("Input error: Cannot open the file \"ecs.inp\".");
    }
    
    // get input from input file
    InputFile inp (cmd.inputfile);
    
    // is there something to compute?
    if (inp.Ei.empty() or inp.instates.empty() or inp.outstates.empty())
    {
        std::cout << "Nothing to compute." << std::endl;
        exit(0);
    }
    
    // 
    // Setup B-spline environment ------------------------------------------ //
    //
    
    double R0 = inp.rknots.back();       // end of real grid
    double Rmax = inp.cknots.back();     // end of complex grid
    
    // check ordering
    if (R0 >= Rmax)
        throw exception ("ERROR: Rmax = %g (end of grid) must be greater than R0 = %g (end of real grid)!", Rmax, R0);
    
    // create B-splines
    Bspline bspline (inp.order, inp.rknots, inp.ecstheta, inp.cknots);
    
    // shortcuts
    int Nspline = bspline.Nspline();
    
    // info
    std::cout << "B-spline solver count: " << Nspline << "\n\n";
    
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
        std::cout << "\n";
    }
    
    std::cout << "\t-> The matrix of the set contains " << coupled_states.size()
              << " diagonal blocks." << std::endl;
    
    // skip if there is nothing to compute
    if (coupled_states.empty())
        std::exit(0);
    
    std::cout << "\n";
    
    //
    // Zip solution file into VTK geometry if told so
    //
    
    if (cmd.zipfile.size() != 0 and par.IamMaster())
    {
        zip_solution (cmd, bspline, coupled_states);
        std::cout << std::endl << "Done." << std::endl << std::endl;
        return 0;
    }
    
// ------------------------------------------------------------------------- //
//                                                                           //
// StgSolve                                                                  //
//                                                                           //
// ------------------------------------------------------------------------- //
    
if (cmd.itinerary & CommandLine::StgSolve)
{
    // create and initialize the preconditioner
    PreconditionerBase * prec = Preconditioners::choose (par, inp, coupled_states, bspline, cmd);
    if (prec == nullptr)
        throw exception ("Preconditioner %d not implemented.", cmd.preconditioner);
    prec->setup();
    
    // CG preconditioner callback
    auto apply_preconditioner = [ & ](cArray const & r, cArray & z) -> void { prec->precondition(r, z); };
    
    // CG matrix multiplication callback
    auto matrix_multiply = [ & ](cArray const & p, cArray & q) -> void { prec->multiply(p, q); };
    
    //
    // Solve the equations
    //
    
    std::cout << "Hamiltonian size: " << Nspline * Nspline * coupled_states.size() << "\n";
    double prevE = special::constant::Nan, E = special::constant::Nan;
    int iterations_done = 0, computations_done = 0;
    for (unsigned ie = 0; ie < inp.Ei.size(); ie++)
    {
        // print progress information
        std::cout << "\nSolving the system for Ei[" << ie << "] = " << inp.Ei[ie] << " ("
                  << int(trunc(ie * 100. / inp.Ei.size() + 0.5)) << " % finished, typically "
                  << (computations_done == 0 ? 0 : iterations_done / computations_done)
                  << " CG iterations per energy)" << std::endl;
        
        cArray current_solution, previous_solution;
        
        // we may have already computed all solutions for this energy... is it so?
        bool all_done = true;
        for (auto instate : inp.instates)
        for (unsigned Spin = 0; Spin <= 1; Spin++)
        {
            int li = std::get<1>(instate);
            int mi = std::get<2>(instate);
            
            // skip angular forbidden states
            bool allowed = false;
            for (int l = abs(li - inp.L); l <= li + inp.L; l++)
                allowed = allowed or special::ClebschGordan(li,mi,l,0,inp.L,mi);
            if (not allowed)
            {
                std::cout << "\tInitial state li=" << li << ", mi=" << mi << " will be skipped (not allowed by total angular variables).\n";
                continue;
            }
            
            // check if there is some precomputed solution on the disk
            SolutionIO reader (inp.L, Spin, inp.Pi, inp.ni, li, mi, inp.Ei[ie]);
            if (not reader.load(current_solution))
                all_done = false;
        }
        if (all_done)
        {
            std::cout << "\tAll solutions for Ei[" << ie << "] = " << inp.Ei[ie] << " loaded." << std::endl;
            continue;
        }
        
        // for all initial states
        for (unsigned instate = 0; instate < inp.instates.size(); instate++)
        for (unsigned Spin = 0; Spin <= 1; Spin++)
        {
            int li = std::get<1>(inp.instates[instate]);
            int mi = std::get<2>(inp.instates[instate]);
            
            // check if the right hand side will be zero for this instate
            bool allowed = false;
            for (int l = abs(li - inp.L); l <= li + inp.L; l++)
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
                continue;
            
            // we may have already computed solution for this state and energy... is it so?
            SolutionIO reader (inp.L, Spin, inp.Pi, inp.ni, li, mi, inp.Ei[ie]);
            if (reader.load(current_solution))
                continue;
            
            // get total energy of the system
            prevE = E; E = 0.5 * (inp.Ei[ie] - 1./(inp.ni * inp.ni));
            
            // update the preconditioner, if this is the first energy to compute or it changed from previous iteration
            if (not (E == prevE)) // does not use 'if (E != prevE)' to work with initial Nan values
                prec->update(E);
            
            // create right hand side
            std::cout << "\tCreate RHS for li = " << li << ", mi = " << mi << ", S = " << Spin << "\n";
            cArray chi (coupled_states.size() * Nspline * Nspline);
            prec->rhs(chi, ie, instate, Spin);
            
            // compute and check norm of the right hand side vector
            double chi_norm = chi.norm();
            if (chi_norm == 0.)
            {
                // this should not happen, hopefully we already checked
                std::cout << "\t! Right-hand-side is zero, check L, Π and nL." << std::endl;
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
            
            // custom conjugate gradients callback-based solver
            current_solution = cArray(chi.size());
            unsigned max_iter = (inp.maxell + 1) * Nspline;
            std::cout << "\tStart CG callback with tolerance " << cmd.itertol << std::endl;
            std::cout << "\t   i | time        | residual        | max & avg block precond. iter." << std::endl;
            unsigned iterations = cg_callbacks<cArray,cArrayView>
            (
                chi,                    // right-hand side
                current_solution,       // on input, the initial guess, on return, the solution
                cmd.itertol,            // requested precision, |A·x - b|² < ε·|b|²
                0,                      // minimal iteration count
                max_iter,               // maximal iteration count
                apply_preconditioner,   // preconditioner callback
                matrix_multiply,        // matrix multiplication callback
                true                    // verbose output
            );
            
            if (iterations >= max_iter)
                std::cout << "\tConvergence too slow... The saved solution will be probably non-converged." << std::endl;
            else
                std::cout << "\tEnd CG callback" << std::endl;
            
            // update progress
            iterations_done += iterations;
            computations_done++;
            
            // save solution to disk (if valid)
            if (std::isfinite(current_solution.norm()))
                reader.save(current_solution);
            
        } // end of For Spin, instate
        
    } // end of For ie = 0, ..., inp.Ei.size() - 1
    
    std::cout << "All solutions computed." << std::endl;
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
    Amplitudes ampl (bspline, inp, par, coupled_states);
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
