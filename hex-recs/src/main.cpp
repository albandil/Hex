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
#include "angular.h"
#include "arrays.h"
#include "bspline.h"
#include "complex.h"
#include "hydrogen.h"
#include "input.h"
#include "itersolve.h"
#include "misc.h"
#include "parallel.h"
#include "preconditioners.h"
#include "radial.h"
#include "version.h"

void zip_solution (CommandLine & cmd, Bspline const & bspline, AngularBasis const & ll)
{
    // use whole grid if not restricted
    if (cmd.zipmax < 0)
        cmd.zipmax = bspline.Rmax();
    
    cArray sol;     // stored solution expansion
    cArray ev;      // evaluated solution
    rArray grid;    // real evaluation grid
    
    // size of solution (l,l)-segment
    size_t N = bspline.Nspline() * bspline.Nspline();
    
    std::cout << "Zipping B-spline expansion of the solution: \"" << cmd.zipfile << "\"" << std::endl;
    
    // load the requested file
    if (not sol.hdfload(cmd.zipfile.c_str()))
        throw exception("Cannot load file %s.", cmd.zipfile.c_str());
    
    // evaluation grid
    grid = linspace (0., cmd.zipmax, cmd.zipcount);
    
    // for all coupled angular momentum pairs
    for (unsigned ill = 0; ill < ll.size(); ill++)
    {
        std::cout << "\t- partial wave l1 = " << ll[ill].l1 << ", l2 = " << ll[ill].l2 << "\n";
        
        // write to file
        std::ofstream out (format("%s_(%d,%d,%d,%d).vtk", cmd.zipfile.c_str(), ll[ill].L, ll[ill].S, ll[ill].l1, ll[ill].l2));
        writeVTK_points
        (
            out,
            bspline.zip
            (
                sol.slice(ill * N, (ill + 1) * N),
                grid, grid
            ),
            grid, grid, rArray({0.})
        );
    }
}

int main (int argc, char* argv[])
{
    // Preparations ------------------------------------------------------- //
    //
    
    // display logo
    std::cout << logo_raw() << std::endl;
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
        exit(0);
    }
    // --------------------------------------------------------------------- //
    
    
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
    // --------------------------------------------------------------------- //
    
    
    // Setup angular data -------------------------------------------------- //
    //
    std::cout << "Setting up the angular expansion basis..." << std::endl;
    
    // coupled angular momentum pairs
    AngularBasis coupled_states;
    
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
                coupled_states.push_back(state);
            }
            std::cout << std::endl;
        }
    }
    
    std::cout << "    The matrix of the set contains " << coupled_states.size()
              << " diagonal blocks." << std::endl;
    
    // skip if there is nothing to compute
    if (coupled_states.empty())
        exit(0);
    
    std::cout << "\n";
    // --------------------------------------------------------------------- //
    
    
    // zip file if told so
    if (cmd.zipfile.size() != 0 and par.IamMaster())
    {
        zip_solution (cmd, bspline, coupled_states);
        goto End;
    }
    
// StgSolve:
{
    // skip stage 2 if told so
    if (not (cmd.itinerary & CommandLine::StgSolve))
    {
        std::cout << "Skipped solution of the equation." << std::endl;
        goto StgExtract;
    }
    
    // create and initialize the preconditioner
    PreconditionerBase * prec = Preconditioners::choose(par, inp, coupled_states, bspline, cmd);
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
            int ni = std::get<0>(inp.instates[instate]);
            int li = std::get<1>(inp.instates[instate]);
            int two_ji = std::get<2>(inp.instates[instate]); double ji = 0.5 * two_ji;
            int two_mi = std::get<3>(inp.instates[instate]); double mi = 0.5 * two_mi;
            
            // get total energy of the system
            prevE = E;
            E = 0.5 * (inp.Ei[ie] - 1./(ni * ni));
            
            // add spin-orbital interaction energy contribution to non-S-states
            if (li > 0)
                E += special::constant::alpha_sqr * (ji*(ji+1) - li*(li+1) - 0.75) * Hydrogen::moment(-3, ni, li); // spin-orbital energy
            
            // update the preconditioner, if this is the first energy to compute or it changed from previous iteration
            if (not (E == prevE)) // does not use 'if (E != prevE)' to work with initial Nan values
                prec->update(E);
            
            // we may have already computed solution for this state and energy... is it so?
            SolutionIO reader (inp.J, ni, li, two_ji, two_mi, inp.Ei[ie]);
            if (reader.load(current_solution))
                continue;
            
            // create right hand side
            std::cout << "\tCreate RHS for ni = " << ni << ", li = " << li << ", ji = " << ji << ", mi = " << mi << ", σi = " << 0.5 * inp.two_M - mi << std::endl;
            cArray chi (coupled_states.size() * Nspline * Nspline);
            prec->rhs(chi, ie, instate);
            if (chi.norm() == 0.)
            {
                std::cout << "\t! Right-hand-side is zero (probably due to incompatible angular settings)." << std::endl;
                computations_done++;
                continue;
            }
            
            // we may have already computed the previous solution - it will serve as initial guess
            current_solution = cArray(chi.size());
            if (ie > 0)
            {
                SolutionIO prev_reader (inp.J, ni, li, two_ji, two_mi, inp.Ei[ie-1]);
                prev_reader.load(current_solution);
            }
            
            // custom conjugate gradients callback-based solver
            std::cout << "\tStart CG callback with tolerance " << cmd.itertol << std::endl;
            unsigned iterations = cg_callbacks<cArray,cArrayView>
            (
                chi,                      // right-hand side
                current_solution,         // on input, the initial guess, on return, the solution
                cmd.itertol,              // requested precision, |A·x - b|² < ε·|b|²
                0,                        // minimal iteration count
                (inp.maxell+1) * Nspline, // maximal iteration count
                apply_preconditioner,     // preconditioner callback
                matrix_multiply           // matrix multiplication callback
            );
            std::cout << "\tEnd CG callback\n";
            
            // update progress
            iterations_done += iterations;
            computations_done++;
            
            // save solution to disk
            reader.save(current_solution);
            
        } // end of For Spin, instate
        
    } // end of For ie = 0, ..., inp.Ei.size() - 1
    
    std::cout << "\rSolving the systems... ok                                                            \n";
    std::cout << "\t(typically " << iterations_done/inp.Ei.size() << " CG iterations per energy)\n";
}

StgExtract:
{
    // skip stage 3 if told so
    if (not (cmd.itinerary & CommandLine::StgExtract))
    {
        std::cout << "Skipped extraction of amplitudes." << std::endl;
        goto End;
    }
    
    // extract amplitudes
    Amplitudes ampl (bspline, inp, par, coupled_states);
    ampl.extract();
    
    // write T-matrices to a text file as SQL statements 
    ampl.writeSQL_files();
    
    // write integral cross sections to a text file
    ampl.writeICS_files();
}

End:
{
    std::cout << std::endl << "Done." << std::endl << std::endl;
}    
    return 0;
}
