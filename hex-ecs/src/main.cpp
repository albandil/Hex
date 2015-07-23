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

#ifdef _OPENMP
    #include <omp.h>
#endif

#include "amplitudes.h"
#include "arrays.h"
#include "bspline.h"
#include "io.h"
#include "misc.h"
#include "parallel.h"
#include "radial.h"
#include "solver.h"
#include "version.h"

int main (int argc, char* argv[])
{
    //
    // Program initialization
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
        
        // disable buffering of the standard output (-> immediate logging)
        std::setvbuf(stdout, nullptr, _IONBF, 0);
        
        // get input from command line
        CommandLine cmd (argc, argv);
    
    //
    // Setup parallel environments
    //
    
        // check some exclusive options
        if (cmd.parallel_block and cmd.lightweight_radial_cache)
            HexException("The options --parallel-block and --lightweight-radial-cache/--lightweight-full can't be used together because of different multiplication scheme.");
        if (cmd.factorizer == LUFT_SUPERLU_DIST and not cmd.parallel)
            HexException("You need to run the program using MPI launcher and with --mpi option to use the distributed SuperLU.");
        
        // setup MPI
        Parallel par (&argc, &argv, cmd.parallel, cmd.groupsize);
        
        // print information
        if (par.active())
        {
            std::cout << "MPI environment" << std::endl;
            std::cout << "\tthis process ID:  " << 1 + par.iproc() << " / " << par.Nproc() << std::endl;
            std::cout << "\tbelongs to group: " << 1 + par.igroup() << " / " << par.Ngroup() << std::endl;
            std::cout << "\tID within group:  " << 1 + par.igroupproc() << " / " << par.groupsize() << std::endl;
            std::cout << std::endl;
        }
    
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
    
    //
    // Read input file
    //
    
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
                return EXIT_FAILURE;
            };
        }
        
        // get input from input file
        InputFile inp (cmd.inputfile);
        
        // is there something to compute?
        if (inp.Etot.empty() or inp.instates.empty() or inp.outstates.empty())
        {
            std::cout << "Nothing to compute." << std::endl;
            return EXIT_SUCCESS;
        }
    
    //
    // Setup B-spline basis
    //
    
        // B-spline bases for all panels
        std::vector<Bspline> bspline;
        
        // number of B-splines (for each basis) that do not use shared knots
        std::vector<int> Nspline_fixed;
        
        // create all panels' bases
        rArray rknots = inp.rknots;
        Nspline_fixed.push_back(rknots.size() - inp.order - 2);
        rknots.append(rknots.back() + inp.overlap_knots.slice(1, inp.overlap_knots.size() - 1));
        for (int ipanel = 0; ipanel < cmd.panels; ipanel++)
        {
            // update B-spline knots for this panel
            if (ipanel > 0)
            {
                rknots.append(rknots.back() + inp.rknots_next.slice(1, inp.rknots_next.size() - 1));
                Nspline_fixed.push_back(rknots.size() - inp.order - 2);
                rknots.append(rknots.back() + inp.overlap_knots.slice(1, inp.overlap_knots.size() - 1));
            }
            
            // create the new B-spline object
            bspline.emplace_back
            (
                Bspline
                (
                    inp.order,
                    rknots,
                    inp.ecstheta,
                    rknots.back() + inp.cknots
                )
            );
            
            // print some information
            std::cout << "B-spline count for panel " << ipanel << ": " << bspline.back().Nspline() << std::endl;
        }
        std::cout << std::endl;
    
    //
    // Setup angular data
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
            return EXIT_SUCCESS;
        
        std::cout << std::endl;
    
    //
    // Zip solution file into VTK geometry if told so
    //
    /*
        if (cmd.zipfile.size() != 0 and par.IamMaster())
        {
            zip_solution (cmd, bspline, coupled_states);
            std::cout << std::endl << "Done." << std::endl << std::endl;
            return EXIT_SUCCESS;
        }
    */
    //
    // Write grid into VTK file if told so.
    //
    /*
        if (cmd.writegrid and par.IamMaster())
        {
            write_grid(bspline);
            std::cout << std::endl << "Done." << std::endl << std::endl;
            return EXIT_SUCCESS;
        }
    */
    
    //
    // Solve the equations
    //
    
        // create the solver instance
        Solver solver (cmd, inp, par, coupled_states, bspline, Nspline_fixed);
        
        // for all solver & propagator panels
        for (int ipanel = 0; ipanel < cmd.panels; ipanel++)
        {
            // pick preconditioner according to the user preferences
            solver.choose_preconditioner(ipanel);
            
            // initialize preconditioner
            solver.setup_preconditioner();
            
            // find the solution of the scattering equations
            solver.solve();
        }

    //
    // Extract amplitudes
    //

        if (cmd.itinerary & CommandLine::StgExtract)
        {
            // extract amplitudes
            Amplitudes ampl (bspline.front(), bspline.back(), inp, par, cmd, coupled_states);
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

    //
    // End
    //

        std::cout << std::endl << "Done." << std::endl << std::endl;
        
        return EXIT_SUCCESS;
}
