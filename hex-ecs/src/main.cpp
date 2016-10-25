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

#include "hex-arrays.h"
#include "hex-misc.h"
#include "hex-version.h"

#include "amplitudes.h"
#include "ang.h"
#include "bspline.h"
#include "io.h"
#include "parallel.h"
#include "radial.h"
#include "solver.h"

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
    
        // some factorizers need MPI environment
        if (not cmd.parallel)
        {
            if
            (
                cmd.factorizer == "superlu_dist" or
                cmd.factorizer == "mumps"
            )
            {
                HexException("You need to run the program using MPI launcher and with --mpi option to use the distributed SuperLU.");
            }
        }
        
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
        # pragma omp parallel
        # pragma omp master
        {
            int nthreads = omp_get_num_threads();
            bool nested = omp_get_nested();
            
            std::cout << "OpenMP environment" << std::endl;
            std::cout << "\tthreads: " << nthreads << std::endl;
            std::cout << "\tnesting: " << (nested ? "on" : "off") << std::endl;
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
    
        // create new B-spline basis for inner problem
        Bspline bspline_inner
        (
            inp.order,
            inp.rknots,
            inp.ecstheta,
            inp.rknots_ext.empty() ? inp.rknots.back() + inp.cknots : rArray{ inp.rknots.back() }
        );
        
        // create new B-spline basis for inner and outer problem (will be the same as 'bspline_atom' if 'rknots_ext' is empty)
        Bspline bspline_full
        (
            inp.order,
            inp.rknots_ext.empty() ? inp.rknots : concatenate(inp.rknots, inp.rknots.back() + inp.rknots_ext.slice(1, inp.rknots_ext.size())),
            inp.ecstheta,
            inp.rknots_ext.empty() ? inp.rknots.back() + inp.cknots : inp.rknots.back() + inp.rknots_ext.back() + inp.cknots
        );
    
    //
    // Setup angular data
    //
    
        AngularBasis ang (inp);
        
        std::cout << "\t-> The matrix of the set contains " << ang.states().size()
                << " diagonal blocks." << std::endl;
        
        // skip if there is nothing to compute
        if (ang.states().empty())
            return EXIT_SUCCESS;
        
        std::cout << std::endl;
    
    //
    // TODO Map solution file between two B-spline bases.
    //
    
        /*if (not cmd.map_solution.empty())
        {
            if (cmd.map_solution_target.empty())
                HexException("The option --map-solution-target is required.");
            
            // read the target basis description
            std::ifstream ifs (cmd.map_solution_target.c_str());
            if (not ifs.good())
                HexException("Failed to open target setup file \"%s\".", cmd.map_solution_target.c_str());
            InputFile target (ifs);
            ifs.close();
            
            // compose knot sequence
            if (not inp.overlap_knots.empty())
            {
                rArray overlap_knots = target.rknots.back() + target.overlap_knots.slice(1, target.overlap_knots.size());
                target.rknots.append(overlap_knots.begin(), overlap_knots.end());
            }
            
            // prepare target basis
            Bspline bspline_target
            (
                target.order,
                target.rknots,
                target.ecstheta,
                target.rknots.back() + target.cknots
            );
            
            // prepare one-electron radial integrals
            RadialIntegrals rad (bspline_panel.front(), bspline_target, bspline_target, 0);
            rad.verbose(false);
            rad.setupOneElectronIntegrals(par, cmd);
            
            // factorize target overlap matrix
            CsrMatrix<LU_int_t,Complex> tgtS = rad.S_proj().tocoo<LU_int_t>().tocsr();
            std::shared_ptr<LUft<LU_int_t,Complex>> tgtSlu = tgtS.factorize();
            
            // multiply inter-basis overlaps by the inverse target overlap matrix
            ColMatrix<Complex> matrix (bspline_target.Nspline(), bspline_panel.front().Nspline());
            ColMatrix<Complex> S21 = rad.S12().T();
            tgtSlu->solve
            (
                S21.data(),
                matrix.data(),
                bspline_panel.front().Nspline()
            );
            RowMatrix<Complex> rmatrix (matrix);
            
            // map all solutions
            for (unsigned i = 0; i < cmd.map_solution.size(); i++)
            {
                // load the solution
                cArray solution;
                solution.hdfload(cmd.map_solution[i]);
                
                // map the solution
                kron_dot(rmatrix, rmatrix, solution).hdfsave("mapped-" + cmd.map_solution[i]);
                std::cout << "Solution \"" << cmd.map_solution[i] << "\" mapped to a new basis." << std::endl;
            }
            
            return EXIT_SUCCESS;
        }*/
    
    //
    // Zip solution file into VTK geometry if told so
    //
    
        if (cmd.zipdata.file.size() != 0 and par.IamMaster())
        {
            zip_solution(cmd, inp, par, bspline_inner, bspline_full, ang.states());
            std::cout << std::endl << "Done." << std::endl << std::endl;
            return EXIT_SUCCESS;
        }
    
    //
    // Write grid into VTK file if told so.
    //
    
        if (cmd.writegrid and par.IamMaster())
        {
            write_grid(bspline_inner, "grid-inner");
            write_grid(bspline_full, "grid-full");
            std::cout << std::endl << "Done." << std::endl << std::endl;
            return EXIT_SUCCESS;
        }
    
    //
    // Solve the equations
    //
    
        // create the solver instance
        Solver solver (cmd, inp, par, ang, bspline_inner, bspline_full);
        
        std::cout << std::endl;
        std::cout << "Bspline basis summary" << std::endl;
        if (bspline_inner.Nspline() != bspline_full.Nspline())
        {
            std::cout << "\t- inner basis" << std::endl;
            std::cout << "\t\t- number of splines: " << bspline_inner.Nspline() << std::endl;
            std::cout << "\t\t- real knots : " << bspline_inner.rknots().front() << " to " << bspline_inner.rknots().back() << std::endl;
            std::cout << "\t\t- complex knots : " << bspline_inner.cknots().front() << " to " << bspline_inner.cknots().back() << std::endl;
        }
        std::cout << "\t- full basis" << std::endl;
        std::cout << "\t\t- number of splines: " << bspline_full.Nspline() << std::endl;
        std::cout << "\t\t- real knots : " << bspline_full.rknots().front() << " to " << bspline_full.rknots().back() << std::endl;
        std::cout << "\t\t- complex knots : " << bspline_full.cknots().front() << " to " << bspline_full.cknots().back() << std::endl;
        std::cout << std::endl;
        
        // pick preconditioner according to the user preferences
        solver.choose_preconditioner();
        
        // initialize preconditioner
        solver.setup_preconditioner();
        
        // find the solution of the scattering equations
        solver.solve();
        
        // release resources
        solver.finish();

    //
    // Extract amplitudes
    //

        if (cmd.itinerary & CommandLine::StgExtract)
        {
            // extract amplitudes
            Amplitudes ampl (bspline_inner, bspline_full, inp, par, cmd, ang.states());
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
