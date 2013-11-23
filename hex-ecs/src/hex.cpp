/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2013                              *
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
#include <H5Cpp.h>

#include "amplitudes.h"
#include "arrays.h"
#include "bspline.h"
#include "complex.h"
#include "input.h"
#include "itersolve.h"
#include "misc.h"
#include "parallel.h"
#include "preconditioners.h"
#include "radial.h"
#include "version.h"

void zip_solution (CommandLine & cmd, Bspline const & bspline, std::vector<std::pair<int,int>> const & ll)
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
    grid = linspace(0., cmd.zipmax, cmd.zipcount);
    
    // for all coupled angular momentum pairs
    for (unsigned ill = 0; ill < ll.size(); ill++)
    {
        // angular momenta
        int l1 = ll[ill].first;
        int l2 = ll[ill].second;
        
        std::cout << "\t- partial wave l1 = " << l1 << ", l2 = " << l2 << "\n";
        
        // write to file
        std::ofstream out (format("%s_(%d,%d).vtk", cmd.zipfile.c_str(), l1, l2));
        bspline.writeVTK (
            out,
            sol.slice (ill * N, (ill + 1) * N),
            grid,
            grid
        );
    }
}

int main (int argc, char* argv[])
{
    // Preparations ------------------------------------------------------- //
    //
    
    // display logo
    std::cout << logo_raw();
    
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
        cmd.inputfile.open("hex.inp");
        if (not cmd.inputfile.good())
            throw exception("Input error: Cannot open the file \"hex.inp\".\n");
    }
    
    // get input from input file
    InputFile inp (cmd.inputfile);
    
    // is there something to compute?
    if (inp.Ei.empty() or inp.instates.empty() or inp.outstates.empty())
    {
        std::cout << "Nothing to compute.\n";
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
    std::cout << "Setting up the coupled angular states...\n";
    
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
    
    std::cout << "\t-> The matrix of the set contains " << coupled_states.size() << " diagonal blocks.\n";
    
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
        std::cout << "Skipped solution of the equation.\n";
        goto StgExtract;
    }
    
    // create and initialize the preconditioner
    PreconditionerBase * prec = Preconditioners::choose (par, inp, coupled_states, bspline, cmd);
    if (prec == nullptr)
        throw exception ("Preconditioner %d not implemented.", cmd.preconditioner);
    prec->setup();
    
    //
    // Solve the equations
    //
    
    std::cout << "Hamiltonian size: " << Nspline * Nspline * coupled_states.size() << "\n";
    int iterations_done = 0, computations_done = 0;
    for (unsigned ie = 0; ie < inp.Ei.size(); ie++)
    {
        // print progress information
        std::cout << "\nSolving the system for Ei[" << ie << "] = " << inp.Ei[ie] << " ("
                  << int(trunc(ie * 100. / inp.Ei.size() + 0.5)) << " % finished, typically "
                  << (computations_done == 0 ? 0 : iterations_done / computations_done)
                  << " CG iterations per energy)\n";
        
        cArray current_solution, previous_solution;
        
        // we may have already computed all solutions for this energy... is it so?
        bool all_done = true;
        for (auto instate : inp.instates)
        {
            int li = std::get<1>(instate);
            int mi = std::get<2>(instate);
            
            // skip angular forbidden states
            bool allowed = false;
            for (int l = abs(li - inp.L); l <= li + inp.L; l++)
                allowed = allowed or ClebschGordan(li,mi,l,0,inp.L,mi);
            if (not allowed)
            {
                std::cout << "\tThe initial state li=" << li << ", mi=" << mi << " is not allowed within the total angular variables.\n";
                continue;
            }
            
            // compose filename of the output file for this solution
            std::ostringstream oss;
            oss << "psi-" << inp.L << "-" << inp.Spin << "-" << inp.Pi << "-" << inp.ni << "-" << li << "-" << mi << "-" << inp.Ei[ie] << ".hdf";
            
            // check if there is some precomputed solution on the disk
            if ( not current_solution.hdfload(oss.str().c_str()) )
                all_done = false;
        }
        if (all_done)
        {
            std::cout << "\tAll solutions for Ei[" << ie << "] = " << inp.Ei[ie] << " loaded.\n";
            continue;
        }
        
        // get total energy of the system
        double E = 0.5 * (inp.Ei[ie] - 1./(inp.ni * inp.ni));
        
        // update the preconditioner
        prec->update(E);
        
        // for all initial states
        for (unsigned instate = 0; instate < inp.instates.size(); instate++)
        {
            int li = std::get<1>(inp.instates[instate]);
            int mi = std::get<2>(inp.instates[instate]);
            
            // skip angular forbidden states
            bool allowed = false;
            for (int l = abs(li - inp.L); l <= li + inp.L; l++)
                allowed = allowed or ClebschGordan(li,mi,l,0,inp.L,mi);
            if (not allowed)
                continue;
            
            // we may have already computed solution for this state and energy... is it so?
            std::ostringstream cur_oss;
            cur_oss << "psi-" << inp.L << "-" << inp.Spin << "-" << inp.Pi << "-" << inp.ni << "-" << li << "-" << mi << "-" << inp.Ei[ie] << ".hdf";
            if ( current_solution.hdfload(cur_oss.str().c_str()) )
                continue;
            
            // create right hand side
            std::cout << "\tCreate RHS for li = " << li << ", mi = " << mi << "\n";
            cArray chi (coupled_states.size() * Nspline * Nspline);
            prec->rhs(chi, ie, instate);
            if (chi.norm() == 0.)
            {
                std::cout << "\t! Right-hand-side is zero (probably due to incompatible angular settings).\n";
                computations_done++;
                continue;
            }
            
            // we may have already computed the previous solution - it will serve as initial guess
            current_solution = cArray(chi.size());
            if (ie > 0)
            {
                std::ostringstream prev_oss;
                prev_oss << "psi-" << inp.L << "-" << inp.Spin << "-" << inp.Pi << "-" << inp.ni << "-" << li << "-" << mi << "-" << inp.Ei[ie-1] << ".hdf";
                if (previous_solution.hdfload(prev_oss.str().c_str()))
                    current_solution = previous_solution;
            }
            
            // CG preconditioner callback
            auto apply_preconditioner = [ & ](cArray const & r, cArray & z) -> void { prec->precondition(r, z); };
            
            // CG matrix multiplication callback
            auto matrix_multiply = [ & ](cArray const & p, cArray & q) -> void { prec->multiply(p, q); };
            
            // custom conjugate gradients callback-based solver
            double tolerance = 1e-10;
            std::cout << "\tStart CG callback with tolerance " << tolerance << "\n";
            unsigned iterations = cg_callbacks<cArray,cArrayView>
            (
                chi,                      // right-hand side
                current_solution,         // on input, the initial guess, on return, the solution
                tolerance,                // requested precision, |A·x - b|² < ε·|b|²
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
            current_solution.hdfsave(cur_oss.str().c_str(), true /* = with compression */);
            
        } // end of For li
        
    } // end of For ie = 0, ..., inp.Ei.size() - 1
    
    std::cout << "\rSolving the systems... ok                                                            \n";
    std::cout << "\t(typically " << iterations_done/inp.Ei.size() << " CG iterations per energy)\n";
}

StgExtract:
{
    // skip stage 3 if told so
    if (not (cmd.itinerary & CommandLine::StgExtract))
    {
        std::cout << "Skipped extraction of amplitudes.\n";
        goto End;
    }
    
    // radial integrals
    RadialIntegrals rad (bspline);
    
    // compose output filename
    std::ostringstream ossfile;
    if (par.active())
        ossfile << inp.ni << "-" << inp.L << "-" << inp.Spin << "-" << inp.Pi << "-(" << par.iproc() << ").sql";
    else
        ossfile << inp.ni << "-" << inp.L << "-" << inp.Spin << "-" << inp.Pi << ".sql";
    
    // Create SQL batch file
    std::ofstream fsql(ossfile.str().c_str());
    
    // set exponential format for floating point output
    fsql.setf(std::ios_base::scientific);
    
    // write header
    fsql << "BEGIN TRANSACTION;\n";
    
    //
    // Extract the cross sections
    //
    
    std::cout << "\nExtracting T-matrices...\n";
    
    std::vector<std::tuple<int,int,int,int,int>> transitions;
    for (auto instate  : inp.instates)
    for (auto outstate : inp.outstates)
        transitions.push_back (
            std::make_tuple (
                inp.Spin,
                /*li*/ std::get<1>(instate),
                /*mi*/ std::get<2>(instate),
                /*nf*/ std::get<0>(outstate),
                /*lf*/ std::get<1>(outstate)
            )
        );
    
    int finished = 0;
    
    for (unsigned i = 0; i < transitions.size(); i++)
    {
        // skip other process's work
        if (not par.isMyWork(i))
            continue;
        
        // get quantum numbers
        int Spin = std::get<0>(transitions[i]);
        int li   = std::get<1>(transitions[i]);
        int mi   = std::get<2>(transitions[i]);
        int nf   = std::get<3>(transitions[i]);
        int lf   = std::get<4>(transitions[i]);
        
        // skip angular forbidden states
        bool allowed = false;
        for (int l = abs(li - inp.L); l <= li + inp.L; l++)
            allowed = allowed or ClebschGordan(li,mi,l,0,inp.L,mi);
        if (not allowed)
            continue;
        
        if (nf > 0)
        {
            //
            // Discrete transition
            //
            
            // precompute hydrogen function overlaps
            cArray Pf_overlaps = rad.overlapP(nf,lf,weightEndDamp(bspline));
            
            // compute radial integrals
            cArrays Lambda(2 * lf + 1);
            for (int mf = -lf; mf <= lf; mf++)
            {
                // final projectile momenta
                rArray kf = sqrt(inp.Ei - 1./(inp.ni*inp.ni) + 1./(nf*nf) + (mf-mi) * inp.B);
                
                // compute Λ for transitions to (nf,lf,mf); it will depend on [ie,ℓ]
                Lambda[mf+lf] = std::move (
                    computeLambda (bspline, kf, inp.ki, inp.maxell, inp.L, Spin, inp.Pi, inp.ni, li, mi, inp.Ei, lf, Pf_overlaps, coupled_states)
                );
            }
            
            // save the data
            for (int mf = -lf; mf <= lf; mf++)
            {
                // final projectile momenta
                rArray kf = sqrt(inp.Ei - 1./(inp.ni*inp.ni) + 1./(nf*nf) + (mf-mi) * inp.B);
                
                // compute Tℓ
                cArray T_ell(Lambda[mf+lf].size());
                for (unsigned i = 0; i < T_ell.size(); i++)
                {
                    int ie  = i / (inp.maxell + 1);
                    int ell = i % (inp.maxell + 1);
                    
                    T_ell[i] = Lambda[mf+lf][i] * 4. * M_PI / kf[ie] * pow(Complex(0.,1.), -ell)
                                    * ClebschGordan(lf, mf, ell, mi - mf, inp.L, mi) / sqrt(2.);
                }
                
                //
                // print out SQL
                //
                
                for (unsigned i = 0; i < T_ell.size(); i++)
                {
                    int ie  = i / (inp.maxell + 1);
                    int ell = i % (inp.maxell + 1);
                    
                    if (finite(T_ell[i].real()) and finite(T_ell[i].imag()))
                    if (T_ell[i].real() != 0. or T_ell[i].imag() != 0.)
                    {
                        fsql << "INSERT OR REPLACE INTO \"tmat\" VALUES ("
                             << inp.ni << "," << li << "," << mi << ","
                             << nf << "," << lf << "," << mf << ","
                             << inp.L  << "," << Spin << ","
                             << inp.Ei[ie] << "," << ell << "," 
                             << T_ell[i].real() << "," << T_ell[i].imag()
                             << ");\n";
                    }
                }
                
                //
                // print out the total cross section for quick overview
                //
                
                std::ostringstream sigmaname;
                sigmaname << "sigma-" 
                          << inp.ni << "-" << li << "-" << mi << "-"
                          << nf << "-" << lf << "-" << mf << "-"
                          << inp.L << "-" << Spin << "-" << inp.Pi << ".dat";

                std::ofstream ftxt(sigmaname.str().c_str());
                
                ftxt << "# Ei [Ry] sigma [a0^2]\n";
                for (unsigned ie = 0; ie < inp.Ei.size(); ie++)
                {
                    double sigma = 0.;
                    for (int ell = 0; ell <= inp.maxell; ell++)
                    {
                        double Re_f_ell = -T_ell[ie * (inp.maxell + 1) + ell].real() / (2 * M_PI);
                        double Im_f_ell = -T_ell[ie * (inp.maxell + 1) + ell].imag() / (2 * M_PI);
                        sigma += 0.25 * (2*Spin + 1) * kf[ie] / inp.ki[ie] * (Re_f_ell * Re_f_ell + Im_f_ell * Im_f_ell);
                    }
                    if (finite(sigma))
                    {
                        ftxt << inp.Ei[ie] << "\t" << sigma << "\n";
                    }
                }
                ftxt.close();
            }
        }
        else
        {
            //
            // Ionization
            //
            
            rArray ics;
            cArrays data = std::move (
                computeXi(bspline, inp.maxell, inp.L, Spin, inp.Pi, inp.ni, li, mi, inp.Ei, ics, coupled_states)
            );
            
            for (size_t ie = 0; ie < inp.Ei.size(); ie++)
            for (unsigned ill = 0; ill < coupled_states.size(); ill++) //??? or triangular
            {
                // save data as BLOBs
                fsql << "INSERT OR REPLACE INTO \"ionf\" VALUES ("
                     << inp.ni << "," << li << "," << mi << ","
                     << inp.L  << "," << Spin << ","
                     << inp.Ei[ie] << "," << coupled_states[ill].first << ","
                     << coupled_states[ill].second << ","
                     << data[ie * coupled_states.size() + ill].toBlob() << ");\n";
            }
            
            // print ionization cross section
            std::ostringstream fname;
            fname << "isigma-" << inp.L << "-" << Spin << "-" << inp.Pi << "-" << inp.ni << "-" << li << "-" << mi << ".dat";
            write_array(inp.Ei, ics, fname.str().c_str());
        }
        
        finished++;
        
        std::cout << "\rExtracting T-matrices... " 
                  << std::setw(3) << int(trunc(finished * 100. / transitions.size() + 0.5))
                  << " %        ";
    }
    
    fsql << "COMMIT;\n";
    fsql.close();
}

End:
{
    std::cout << "\nDone.\n\n";
}    
    return 0;
}
