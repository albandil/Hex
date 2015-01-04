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
#include <vector>
#include <iostream>
#include <string>
#include <tuple>

#include <gsl/gsl_errno.h>

#include "amplitudes.h"
#include "angular.h"
#include "arrays.h"
#include "bspline.h"
#include "complex.h"
#include "hydrogen.h"
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
    std::cout << "=== Relativistic exterior complex scaling in B-splines ===" << std::endl << std::endl;
    
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
        cmd.inputfile.open("recs.inp");
        if (not cmd.inputfile.good())
            throw exception("Input error: Cannot open the file \"recs.inp\".");
    }
    
    // get input from input file
    InputFile inp (cmd.inputfile);
    
    // is there something to compute?
    if (inp.Ei.empty() or inp.instates.empty() or inp.outstates.empty())
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
    
    std::cout << "Setting up the angular expansion basis..." << std::endl;
    
    // coupled angular momentum pairs
    AngularBasis ang;
    
    // for given J list all available (L,S) pairs
    AngularState state;
    for (state.S = 0; state.S <= 1; state.S++)
    for (state.L = std::abs(inp.J - state.S); state.L <= inp.J + state.S; state.L++)
    {
        std::cout << "    (L,S) = (" << state.L << "," << state.S << ")" << std::endl;
        
        // for given L, Π and levels list all available (ℓ₁ℓ₂) pairs
        for (int ell = 0; ell <= inp.levels; ell++)
        {
            std::cout << "        -> [" << ell << "] ";
            
            // get sum of the angular momenta for this angular level
            int sum = 2 * ell + state.L;
            
            // for all angular momentum pairs that do compose L
            for (state.l1 = ell; state.l1 <= sum - ell; state.l1++)
            {
                state.l2 = sum - state.l1;
                std::cout << "(" << state.l1 << "," << state.l2 << ") ";
                ang.push_back(state);
            }
            std::cout << std::endl;
        }
    }
    
    std::cout << "    The matrix of the set contains " << ang.size()
              << " diagonal blocks." << std::endl;
    
    // skip if there is nothing to compute
    if (ang.empty())
        std::exit(EXIT_SUCCESS);
    
    std::cout << std::endl;
    
    //
    // Zip solution file into VTK geometry if told so
    //
    
    if (cmd.zipfile.size() != 0 and par.IamMaster())
    {
        zip_solution (cmd, bspline, ang);
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
    PreconditionerBase * prec = Preconditioners::choose(par, inp, ang, bspline, cmd);
    if (prec == nullptr)
        throw exception ("Preconditioner %d not implemented.", cmd.preconditioner);
    prec->setup();
    
    // CG preconditioner callback
    auto apply_preconditioner = [ & ](cArray const & r, cArray & z) -> void
    {
        prec->precondition(r, z);
    };
    
    // CG matrix multiplication callback
    auto matrix_multiply = [ & ](cArray const & p, cArray & q) -> void
    {
        prec->multiply(p, q);
    };
    
    //
    // Solve the equations
    //
    
    std::cout << "Hamiltonian size: " << Nspline * Nspline * ang.size() << "\n";
    int iterations_done = 0, computations_done = 0;
    double prevE = special::constant::Nan, E = special::constant::Nan;
    for (unsigned ie = 0; ie < inp.Ei.size(); ie++)
    {
        // print progress information
        std::cout << "\nSolving the system for Ei[" << ie << "] = " << inp.Ei[ie] << " ("
                  << int(trunc(ie * 100. / inp.Ei.size() + 0.5)) << " % finished, typically "
                  << (computations_done == 0 ? 0 : iterations_done / computations_done)
                  << " CG iterations per energy)" << std::endl;
        
        cArray current_solution, previous_solution;
        
        // for all initial states
        for (unsigned instate = 0; instate < inp.instates.size(); instate++)
        {
            // initial quantum numbers
            int ni = std::get<0>(inp.instates[instate]);
            int li = std::get<1>(inp.instates[instate]);
            int two_ji = std::get<2>(inp.instates[instate]); double ji = 0.5 * two_ji;
            int two_mi = std::get<3>(inp.instates[instate]); double mi = 0.5 * two_mi;
            
            // we may have already computed solution for this state and energy... is it so?
            SolutionIO reader (inp.J, inp.M, ni, li, two_ji, two_mi, inp.Ei[ie]);
            if (reader.load(current_solution))
            {
                std::cout << "\tSolution already exists in file \"" << reader.name() << "\"" << std::endl;
                computations_done++;
                continue;
            }
            
            // get total energy of the system
            prevE = E;
            E = 0.5 * (inp.Ei[ie] - 1./(ni * ni));
            
            // add spin-orbital interaction energy contribution to non-S-states
            if (li > 0)
            {
                E += 0.5 * special::constant::alpha_sqr * (ji*(ji+1) - li*(li+1) - 0.75) // = ½(Zα)² L·S
                     / (ni*ni*ni * li * (li+0.5) * (li+1)); // = ⟨ni,li|r⁻³|ni,li⟩
            }
            
            // update the preconditioner, if this is the first energy to compute or it changed from previous iteration
            if (not (E == prevE)) // does not use 'if (E != prevE)' to work with initial Nan values
                prec->update(E);
            
            // create right hand side
            std::cout << "\tCreate RHS for ni = " << ni << ", li = " << li << ", ji = " << ji << ", mi = " << mi << ", σi = " << inp.M - mi << std::endl;
            cArray chi (ang.size() * Nspline * Nspline);
            prec->rhs(chi, ie, instate);
            if (chi.norm() == 0.)
            {
                std::cout << "\t! Right-hand-side is zero (probably due to incompatible angular settings)." << std::endl;
                computations_done++;
                continue;
            }
            
            // custom conjugate gradients callback-based solver
            current_solution = cArray(chi.size());
            unsigned max_iter = (inp.maxell + 1) * Nspline;
            std::cout << "\tStart CG callback with tolerance " << cmd.itertol << std::endl;
            std::cout << "\t   i | time        | residual        | min max avg block precond. iter." << std::endl;
            unsigned iterations = cg_callbacks<cArray,cArrayView>
            (
                chi,                    // right-hand side
                current_solution,       // on input, the initial guess, on return, the solution
                cmd.itertol,            // requested precision, |A·x - b|² < ε·|b|²
                0,                      // minimal iteration count
                max_iter,               // maximal iteration count
                apply_preconditioner,   // preconditioner callback
                matrix_multiply         // matrix multiplication callback
            );
            std::cout << "\tEnd CG callback\n";
            
            // update progress
            iterations_done += iterations;
            computations_done++;
            
            // save solution to disk
            reader.save(current_solution);
            
        } // end of For Spin, instate
        
    } // end of For ie = 0, ..., inp.Ei.size() - 1
    
    std::cout << "All solutions computed." << std::endl;
    std::cout << "\t(typically " << iterations_done / computations_done << " CG iterations per solution)\n";
}
else
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
    Amplitudes ampl (bspline, inp, par, ang);
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
