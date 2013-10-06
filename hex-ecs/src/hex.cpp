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

#include <algorithm>
#include <cstdio>
#include <cmath>
#include <chrono>
#include <cstring>
#include <vector>
#include <iostream>
#include <string>
#include <tuple>

#include <gsl/gsl_errno.h>
#include <omp.h>
#include <H5Cpp.h>

#include "amplitudes.h"
#include "arrays.h"
#include "bspline.h"
#include "complex.h"
#include "input.h"
#include "itersolve.h"
#include "misc.h"
#include "moments.h"
#include "parallel.h"
#include "spmatrix.h"
#include "version.h"

#ifndef SuiteSparse_time
double SuiteSparse_time()
{
    return 0.;
}
#endif

int main(int argc, char* argv[])
{
    std::cout <<
        "\n"
        "       / /   / /    __    \\ \\  / /\n"
        "      / /__ / /   / _ \\    \\ \\/ /\n"
        "     /  ___  /   | |/_/    / /\\ \\\n"
        "    / /   / /    \\_\\      / /  \\ \\\n"
        "\n"
        "             UK MFF (c) 2013\n"
        "\n"
        "          version hash: " << commit_hash << "\n"
        "\n";
    
    // Preparations ------------------------------------------------------- //
    //
    
    gsl_set_error_handler_off();
    H5::Exception::dontPrint();
    
    // disable buffering of the standard output 
    // (so that all text messages are immediatelly visible)
    setvbuf(stdout, nullptr, _IONBF, 0);
    
    // variables that can be set by the user from the command line
    std::ifstream inputfile;    // input file
    std::string zipfile;        // HDF solution expansion file to zip
    int  zipcount = 0;          // zip sample count
    double zipmax = -1;         // zip bounding box
    bool parallel = false;      // whether to use OpenMPI
    double droptol = 1e-15;        // LU decomposition drop tolerance
    int preconditioner = ilu_prec; // preconditioner
    
    // which stages to run (default: all)
    int itinerary = StgNone;
    
    // get input from command line
    parse_command_line (
        argc, argv,
        inputfile,
        zipfile, zipcount, zipmax,
        parallel, preconditioner, droptol,
        itinerary
    );
    
    // run the whole sequence if nothing specified
    if (itinerary == StgNone)
        itinerary = StgRadial | StgSolve | StgExtract;
    
    // setup MPI
    Parallel par(parallel);
    
    // check input file
    if (not inputfile.is_open())
    {
        inputfile.open("hex.inp");
        if (not inputfile.good())
            throw exception("Input error: Cannot open the file \"hex.inp\".\n");
    }
    
    // create main variables
    int order;              // B-spline order
    double ecstheta;        // ECS rotation
    rArray real_knots;      // real knot sequence
    rArray complex_knots;   // complex knot sequence
    int ni;                 // atomic state parameters
    int L, Spin, Pi;        // global conserved variables
    int levels;             // angular momentum levels
    rArray Ei;              // energies in Rydbergs
    double B = 0;           // magnetic field in a.u.
    
    // initial and final atomic states
    std::vector<std::tuple<int,int,int>> instates, outstates;
    
    // get input from input file
    parse_input_file (
        inputfile, order, ecstheta,
        real_knots, complex_knots, ni, 
        instates, outstates,
        L, Spin, Pi, levels, Ei, B
    );
    
    // is there something to compute?
    if (Ei.empty() or instates.empty() or outstates.empty())
    {
        std::cout << "Nothing to compute.\n";
        exit(0);
    }
    
    // shorthand for total number of energies
    unsigned Nenergy = Ei.size();
    
    // projectile momenta, initial and final
    std::vector<double> ki(Nenergy), kf(Nenergy);        // in a.u.
    std::transform (
        Ei.begin(), Ei.end(),
        ki.begin(),
        [ = ](double E) -> double { return sqrt(E); }
    );
    // --------------------------------------------------------------------- //
    
    
    // Setup B-spline environment ------------------------------------------ //
    //
    double R0 = real_knots.back();          // end of real grid
    double Rmax = complex_knots.back();     // end of complex grid
    
    if (R0 >= Rmax)
        throw exception("ERROR: Rmax = %g (end of grid) must be greater than R0 = %g (end of real grid)!", Rmax, R0);

    Bspline bspline(order, real_knots, ecstheta, complex_knots);
    
    // shortcuts
    int Nspline = bspline.Nspline();
    
    // info
    std::cout << "B-spline total count: " << Nspline << "\n\n";
    //
    // --------------------------------------------------------------------- //
    
    
    // Setup angular data -------------------------------------------------- //
    //
    std::vector<std::pair<int,int>> coupled_states;
    std::vector<int> workers;
    
    std::cout << "Setting up the coupled angular states...\n";
    
    // for given L, Π and levels list all available (ℓ₁ℓ₂) pairs
    for (int ell = 0; ell <= levels; ell++)
    {
        std::cout << "\t-> [" << ell << "] ";
        
        // get sum of the angular momenta for this angular level
        int sum = 2 * ell + L + Pi;
        
        // for all angular momentum pairs that do compose L
        for (int l1 = ell; l1 <= sum - ell; l1++)
        {
            std::cout << "(" << l1 << "," << sum - l1 << ") ";
            coupled_states.push_back(std::make_pair(l1, sum - l1));
        }
        std::cout << "\n";
    }
    
    std::cout << "\t-> The matrix of the set contains " << coupled_states.size() << " diagonal blocks.\n";
    
    // skip if nothing to compute
    if (coupled_states.empty())
        exit(0);
    
    // store maximal angular momentum
    int maxell = levels + L;
    
    std::cout << "\n";
    // --------------------------------------------------------------------- //
    
    
    // --------------------------------------------------------------------- //
    // initialize radial integrals
    
    RadialIntegrals rad(bspline);
    
    // zip file if told so
    if (zipfile.size() != 0 and par.IamMaster())
    {
        if (zipmax < 0)
            zipmax = Rmax;
        
        cArray sol;     // stored solution expansion
        cArray ev;      // evaluated solution
        rArray grid;    // real evaluation grid
        
        std::cout << "Zipping B-spline expansion of the solution: \"" << zipfile << "\"" << std::endl;
        
        if (not sol.hdfload(zipfile.c_str()))
            throw exception("Cannot load file %s.", zipfile.c_str());
        
        grid = linspace(0., zipmax, zipcount);
        
        for (unsigned ill = 0; ill < coupled_states.size(); ill++)
        {
            int l1 = coupled_states[ill].first;
            int l2 = coupled_states[ill].second;
            
            std::cout << "\t- partial wave l1 = " << l1 << ", l2 = " << l2 << "\n";
            
            // zip this partial wave
            ev = bspline.zip (
                cArrayView (
                    sol,
                    ill * Nspline * Nspline,
                    Nspline * Nspline
                ),
                grid, grid
            );
            
            // write VTK header
            std::ostringstream ofname;
            ofname << zipfile << "_(" << l1 << "," << l2 << ").vtk";
            std::ofstream out(ofname.str().c_str());
            out << "# vtk DataFile Version 3.0\n";
            out << "Hex-ecs wave function partial waves\n";
            out << "ASCII\n";
            out << "DATASET RECTILINEAR_GRID\n";
            out << "DIMENSIONS " << zipcount << " " << zipcount << " 1\n";
            out << "X_COORDINATES " << zipcount << " float\n";
            out << grid.string() << "\n";
            out << "Y_COORDINATES " << zipcount << " float\n";
            out << grid.string() << "\n";
            out << "Z_COORDINATES 1 float\n";
            out << "0\n";
            out << "POINT_DATA " << zipcount * zipcount << "\n";
            out << "FIELD wavefunction 2\n";
            
            // save real part
            out << "pw_" << l1 << "_" << l2 << "_re 1 " << zipcount * zipcount << " float\n";
            for (int i = 0; i < zipcount; i++)
            {
                for (int j = 0; j < zipcount; j++)
                    out << ev[i * zipcount + j].real() << " ";
                out << "\n";
            }
            
            // save imaginary part
            out << "pw_" << l1 << "_" << l2 << "_im 1 " << zipcount * zipcount << " float\n";
            for (int i = 0; i < zipcount; i++)
            {
                for (int j = 0; j < zipcount; j++)
                    out << ev[i * zipcount + j].imag() << " ";
                out << "\n";
            }
        }
        
        goto End;
    }
    
// Stg1:
{
    // precompute one-electron integrals
    rad.setupOneElectronIntegrals();
    
    // precompute two-electron integrals
    if (not (itinerary & StgRadial))
    {
        std::cout << "Skipped computation of two-electron integrals.\n";
        goto Stg2;
    }
    rad.setupTwoElectronIntegrals(par, L + 2 * levels);
}
Stg2:
{
    // skip stage 2 if told so
    if (not (itinerary & StgSolve))
    {
        std::cout << "Skipped solution of the equation.\n";
        goto Stg3;
    }
    
    // Prepare B-spline overlaps and expansions of Ric-Bess functions ------ //
    //
    std::cout << "Computing B-spline expansions... ";
    
    //  j-overlaps of shape [Nenergy × Nangmom × Nspline]
    cArray ji_overlaps = rad.overlapj(maxell,ki,weightEdgeDamp(bspline));
    ji_overlaps.hdfsave("ji_overlaps_damp.hdf"); // just for debugging
    
    //  compute expansions; solve the system
    //      S * B_spline_expansion = B_spline_overlap
    unsigned ji_expansion_count = ji_overlaps.size()/Nspline;
    cArray ji_expansion = rad.S().tocoo().tocsr().solve(ji_overlaps, ji_expansion_count);
    ji_expansion.hdfsave("ji_expansion.hdf"); // just for debugging
    
    std::cout << "ok\n\n";
    // --------------------------------------------------------------------- //
    
    
    // Precompute some accelerators ---------------------------------------- //
    //
    // Kronecker producs
    std::cout << "Creating Kronecker products... ";
    SymDiaMatrix S_kron_S, S_kron_Mm1_tr, S_kron_Mm2, Mm1_tr_kron_S,
                 Mm2_kron_S, half_D_minus_Mm1_tr, half_D_minus_Mm1_tr_kron_S,
                 S_kron_half_D_minus_Mm1_tr;
    # pragma omp parallel sections
    {
        # pragma omp section
        S_kron_S   = rad.S().kron(rad.S());
        # pragma omp section
        S_kron_Mm1_tr = rad.S().kron(rad.Mm1_tr());
        # pragma omp section
        S_kron_Mm2 = rad.S().kron(rad.Mm2());
        # pragma omp section
        Mm1_tr_kron_S = rad.Mm1_tr().kron(rad.S());
        # pragma omp section
        Mm2_kron_S = rad.Mm2().kron(rad.S());
        # pragma omp section
        half_D_minus_Mm1_tr = 0.5 * rad.D() - rad.Mm1_tr();
    }
    # pragma omp parallel sections
    {
        # pragma omp section
        half_D_minus_Mm1_tr_kron_S = half_D_minus_Mm1_tr.kron(rad.S());
        # pragma omp section
        S_kron_half_D_minus_Mm1_tr = rad.S().kron(half_D_minus_Mm1_tr);
    }
    std::cout << "ok\n\n";
    // --------------------------------------------------------------------- //
    
    
    // Distribute LU factorizations among processes ------------------------ //
    //
    std::map<int,int> LUs;
    std::vector<int> info(par.Nproc());
    int worker = 0;
    std::cout << "Balancing " << coupled_states.size()
              << " diagonal blocks among " << par.Nproc() 
              << " worker processes...\n";
    for (unsigned ill = 0; ill < coupled_states.size(); ill++)
    {
        info[worker]++;                       // add work to the process 'worker'
        
        LUs[ill] = worker;                  // assign this block to worker
        worker = (worker + 1) % par.Nproc();        // move to next worker
    }
    // print statistics
    std::cout << "\t-> average " << coupled_states.size()/double(par.Nproc()) << " blocks/process\n";
    std::cout << "\t-> min " << *std::min_element(info.begin(), info.end()) << " blocks/process\n";
    std::cout << "\t-> max " << *std::max_element(info.begin(), info.end()) << " blocks/process\n";
    // --------------------------------------------------------------------- //
    
    
    // For all right hand sides -------------------------------------------- //
    //
    std::cout << "Hamiltonian properties:\n";
    std::cout << "\t-> hamiltonian size: " << Nspline * Nspline * coupled_states.size() << "\n\n";
    int iterations_done = 0, computations_done = 0;
    for (unsigned ie = 0; ie < Nenergy; ie++)
    {
        // print progress information
        std::cout << "\nSolving the system for Ei[" << ie << "] = " << Ei[ie] << " ("
                  << int(trunc(ie * 100. / Nenergy + 0.5)) << " % finished, typically "
                  << (computations_done == 0 ? 0 : iterations_done / computations_done)
                  << " CG iterations per energy)\n";
        
        cArray current_solution, previous_solution;
        
        // we may have already computed all solutions for this energy... is it so?
        bool all_done = true;
        for (auto instate : instates)
        {
            int li = std::get<1>(instate);
            int mi = std::get<2>(instate);
            
            // skip angular forbidden states
            bool allowed = false;
            for (int l = abs(li - L); l <= li + L; l++)
                allowed = allowed or ClebschGordan(li,mi,l,0,L,mi);
            if (not allowed)
            {
                std::cout << "\tThe initial state li=" << li << ", mi=" << mi << " is not allowed within the total angular variables.\n";
                continue;
            }
            
            // compose filename of the output file for this solution
            std::ostringstream oss;
            oss << "psi-" << L << "-" << Spin << "-" << Pi << "-" << ni << "-" << li << "-" << mi << "-" << Ei[ie] << ".hdf";
            
            // check if there is some precomputed solution on the disk
            if ( not current_solution.hdfload(oss.str().c_str()) )
                all_done = false;
        }
        if (all_done)
        {
            std::cout << "\tAll solutions for Ei[" << ie << "] = " << Ei[ie] << " loaded.\n";
            continue;
        }
        
        // get total energy of the system
        double E = 0.5 * (Ei[ie] - 1./(ni*ni));
        
        // diagonal blocks in DIA format
        // - these will be used in matrix multiplication
        std::vector<SymDiaMatrix> dia_blocks(coupled_states.size());
        
        // incomplete Cholesky factorization of the diagonal blocks
        //   L + Lt - I
        std::vector<SymDiaMatrix> icholL(coupled_states.size());
        //   D⁻¹
        std::vector<cArray> icholD(coupled_states.size());
        
        // diagonal incomplete Cholesky (DIC) preconditioner
        std::vector<SymDiaMatrix> DIC(coupled_states.size());
        
        // SSOR preconditioner
        std::vector<SymDiaMatrix> SSOR(coupled_states.size());
        
        // incomplete LU factorizations of the diagonal blocks
        std::vector<CsrMatrix> csr_blocks(coupled_states.size());        
        std::vector<CsrMatrix::LUft> iLU(coupled_states.size());
        std::vector<std::vector<CsrMatrix>> scsr_blocks(coupled_states.size());
        std::vector<std::vector<CsrMatrix::LUft>> siLU(coupled_states.size());
        
        // block incomplete D-ILU
        std::vector<std::vector<SymDiaMatrix>> bcks(coupled_states.size());
        std::vector<std::vector<CsrMatrix>> bcks_csr(coupled_states.size());
        std::vector<std::vector<CsrMatrix::LUft>> bcks_lufts(coupled_states.size());
        
        // sparse approximate inverse
        std::vector<SymDiaMatrix> spai(coupled_states.size());
        
        // setup the preconditioner - the diagonal block iChol-factorizations
        std::cout << "\tSetup preconditioner blocks... " << std::flush;
        for (unsigned ill = 0; ill < coupled_states.size(); ill++)
        {
            // skip computation of unwanted blocks for this process
            if (LUs.find(ill) == LUs.end() or LUs[ill] != par.iproc())
                continue;
            
            int l1 = coupled_states[ill].first;
            int l2 = coupled_states[ill].second;
            
            // one-electron parts
            SymDiaMatrix Hdiag =
                half_D_minus_Mm1_tr_kron_S
                + (0.5*l1*(l1+1)) * Mm2_kron_S
                + S_kron_half_D_minus_Mm1_tr
                + (0.5*l2*(l2+1)) * S_kron_Mm2;
            
            // two-electron part
            for (unsigned lambda = 0; lambda <= rad.maxlambda(); lambda++)
            {
                Complex f = computef(lambda,l1,l2,l1,l2,L);
                if (f != 0.)
                    Hdiag += f * rad.R_tr_dia(lambda);
            }
            
            // finalize the matrix
            dia_blocks[ill] = E*S_kron_S - Hdiag;
            csr_blocks[ill] = dia_blocks[ill].tocoo().tocsr();
            
            dia_blocks[ill].hdfsave(format("dia-%d-%d.hdf",l1,l2));
//             csr_blocks[ill].plot(format("csr-%d-%d.png",l1,l2));
        }
        std::cout << "ok\n";
        
        // setup the preconditioner
        std::cout << "\tCompose preconditioner matrices..." << std::flush;
        for (unsigned ill = 0; ill < coupled_states.size(); ill++)
        {
            // skip computation of unwanted blocks for this process
            if (LUs.find(ill) == LUs.end() or LUs[ill] != par.iproc())
                continue;
            
            // timer info
            std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
            
            int l1 = coupled_states[ill].first;
            int l2 = coupled_states[ill].second;
            
            // drop-tolerance-incomplete LU factorization
            if (preconditioner == ilu_prec)
            {
                // log output
                std::cout << "\n\t\t-> [" << par.iproc() << "] iLU factorization " 
                          << ill << " of (" << l1 << "," << l2 << ") block started\n";
                
                // drop-tolerance incomplete LU factorization
                iLU[ill] = csr_blocks[ill].factorize(droptol);
                
                // log output
                std::chrono::duration<int> sec = std::chrono::duration_cast<std::chrono::duration<int>>(std::chrono::steady_clock::now()-start);
                std::cout << "\t\t   [" << par.iproc() << "] "
                          << "droptol " << droptol << ", "
                          << "time " << sec.count() / 60 << ":" << std::setfill('0') << std::setw(2) << sec.count() % 60 << ", "
                          << "mem " << iLU[ill].size() / 1048576 << " MiB";
            }
            
            // block-diagonal iLU
            if (preconditioner == bilu_prec)
            {
                // log output
                std::cout << "\n\t\t-> [" << par.iproc() << "] block D-ILU factorization "
                          << ill << " of (" << l1 << "," << l2 << ") block started\n";
                
                // allocate space
                std::vector<SymDiaMatrix> ibcks(Nspline);
                bcks[ill].resize(Nspline);
                bcks_csr[ill].resize(Nspline);
                bcks_lufts[ill].resize(Nspline);
                
                std::cout << "\t\t   [" << par.iproc() << "] ";
                std::cout << "computing blocks and SPAIs\n";
                
                // for all diagonal blocks
                for (int iblock = 0; iblock < Nspline; iblock++)
                {
                    // - compute the block
                    bcks[ill][iblock] = E * rad.S().main_diagonal()[iblock] * rad.S()
                            - 0.5 * rad.D().main_diagonal()[iblock] * rad.S()
                            - 0.5 * rad.S().main_diagonal()[iblock] * rad.D()
                            - 0.5 * l1 * (l1 + 1.) * rad.Mm2().main_diagonal()[iblock] * rad.S()
                            - 0.5 * l2 * (l2 + 1.) * rad.S().main_diagonal()[iblock] * rad.Mm2()
                            + rad.Mm1_tr().main_diagonal()[iblock] * rad.S()
                            + rad.S().main_diagonal()[iblock] * rad.Mm1_tr();
                    
                    // - invert using the SPAI
                    ibcks[iblock] = SymDiaMatrix (
                        Nspline,
                        iArray(1), // = (int[]){ 0 }
                        cArray(Nspline,1.)/bcks[ill][iblock].main_diagonal()
                    );
                }
                
                std::cout << "\t\t   [" << par.iproc() << "] ";
                std::cout << "preconditioning the diagonal\n";
                
                // for all diagonals of overlap matrix
                for (int idiag = 1; idiag <= order; idiag++)
                {
                    // for all elements of this diagonal
                    for (int irow = 0; irow < Nspline - idiag; irow++)
                    {
                        // - construct corresponding block of Kronecker product
                        SymDiaMatrix bck = E * rad.S().dptr(idiag)[irow] * rad.S()
                                - 0.5 * rad.D().dptr(idiag)[irow] * rad.S()
                                - 0.5 * rad.S().dptr(idiag)[irow] * rad.D()
                                - 0.5 * l1 * (l1 + 1.) * rad.Mm2().dptr(idiag)[irow] * rad.S()
                                - 0.5 * l2 * (l2 + 1.) * rad.S().dptr(idiag)[irow] * rad.Mm2()
                                + rad.Mm1_tr().dptr(idiag)[irow] * rad.S()
                                + rad.S().dptr(idiag)[irow] * rad.Mm1_tr();
                        
                        /*if (irow + idiag == 1)
                        {
                            std::cout << "idiag = " << idiag << ", irow = " << irow << "\n";
                            std::cout << "\tPřed: bcks[ill][irow+idiag].diag().size() = " << bcks[ill][irow+idiag].diag().size() << "\n";
                            std::cout << "\t      bck.diag().size() = " << bck.diag().size() << "\n";
                            std::cout << "bcks[ill][irow+idiag]\n" << bcks[ill][irow+idiag] << "\n";
                            std::cout << "bck * (ibcks[irow] * bck)\n" << bck * (ibcks[irow] * bck) << "\n";
                        }*/
                        // - update pivot
                        bcks[ill][irow+idiag] -= bck * (ibcks[irow] * bck);
                        /*if (irow + idiag == 1)
                        {
                            std::cout << "\tPo: bcks[ill][irow+idiag].diag().size() = " << bcks[ill][irow+idiag].diag().size() << "\n";
                            std::cout << bcks[ill][irow+idiag] << "\n";
                        }*/
                    }
                }
                
                std::cout << "\t\t   [" << par.iproc() << "] ";
                std::cout << "factorization of the pivots\n";
                
                // convert and factorize the pivots
                size_t size = 0;
                for (int iblock = 0; iblock < Nspline; iblock++)
                {
                    /*std::cout << bcks[ill][iblock] << "\n";*/
                    
                    bcks_csr[ill][iblock] = bcks[ill][iblock].tocoo().tocsr();
                    bcks_lufts[ill][iblock] = bcks_csr[ill][iblock].factorize(droptol);
                    size += bcks_lufts[ill][iblock].size();
                }
                
                // log output
                std::chrono::duration<int> sec = std::chrono::duration_cast<std::chrono::duration<int>>(std::chrono::steady_clock::now()-start);
                std::cout << "\t\t   [" << par.iproc() << "] ";
                std::cout << "droptol " << droptol << ", ";
                std::cout << "time " << sec.count() / 60 << ":" << std::setfill('0') << std::setw(2) << sec.count() % 60 << ", ";
                std::cout << "average " << (sec.count() / Nspline) / 60 << ":" << std::setfill('0') << std::setw(2) << (sec.count() / Nspline) % 60 << ", ";
                std::cout << "mem " << size/1048576 << " MiB";
            }
            
            // drop-tolerance-incomplete single-electron LU factorization
            if (preconditioner == silu_prec)
            {
                // log output
                std::cout << "\n\t\t-> [" << par.iproc() << "] s-iLU factorization "
                          << ill << " of (" << l1 << "," << l2 << ") block started\n";
                
                // allocate space
                scsr_blocks[ill].resize(Nspline);
                siLU[ill].resize(Nspline);
                
                unsigned size = 0;
                
                // for all single-electron blocks
                for (int iblock = 0; iblock < Nspline; iblock++)
                {
                    // setup the single-electron block
                    scsr_blocks[ill][iblock] = (
                        E * rad.S().main_diagonal()[iblock] * rad.S()
                      - 0.5 * rad.D().main_diagonal()[iblock] * rad.S()
                      - 0.5 * rad.S().main_diagonal()[iblock] * rad.D()
                      - 0.5 * l1 * (l1 + 1.) * rad.Mm2().main_diagonal()[iblock] * rad.S()
                      - 0.5 * l2 * (l2 + 1.) * rad.S().main_diagonal()[iblock] * rad.Mm2()
                      + rad.Mm1_tr().main_diagonal()[iblock] * rad.S()
                      + rad.S().main_diagonal()[iblock] * rad.Mm1_tr()
                    ).tocoo().tocsr();
                    
//                     scsr_blocks[ill][iblock].plot(format("scsr-%d-%.03d.png", ill, iblock), 1.);
//                     scsr_blocks[ill][iblock].hdfsave(format("scsr-%d-%.03d.hdf", ill, iblock));
                    
                    // factorize the block
                    siLU[ill][iblock] = scsr_blocks[ill][iblock].factorize(droptol);
                    size += siLU[ill][iblock].size();
                }
                
                // log output
                std::chrono::duration<int> sec = std::chrono::duration_cast<std::chrono::duration<int>>(std::chrono::steady_clock::now()-start);
                std::cout << "\t\t   [" << par.iproc() << "] ";
                std::cout << "droptol " << droptol << ", ";
                std::cout << "time " << sec.count() / 60 << ":" << std::setfill('0') << std::setw(2) << sec.count() % 60 << ", ";
                std::cout << "average " << (sec.count() / Nspline) / 60 << ":" << std::setfill('0') << std::setw(2) << (sec.count() / Nspline) % 60 << ", ";
                std::cout << "mem " << size/1048576 << " MiB";
            }
            
            // sparse approximate inverse preconditioner
            if (preconditioner == spai_prec)
            {
                // use diagonal preconditioner
                spai[ill] = SymDiaMatrix (
                    Nspline * Nspline,
                    iArray(1),
                    cArray(Nspline * Nspline, 1.) / dia_blocks[ill].main_diagonal()
                );
            }
            
            // multi-resolution preconditioner
            if (preconditioner == res_prec)
            {
                
            }
            
            // diagonal incomplete Cholesky factorization
            if (preconditioner == dic_prec)
            {
                DIC[ill] = DIC_preconditioner(dia_blocks[ill]);
                cArray(DIC[ill].main_diagonal()).hdfsave(format("DIC-%d.hdf",ill));
            }
            
            // Jacobi or SSOR preconditioner (Jacobi will use only the diagonal)
            if (preconditioner == jacobi_prec or preconditioner == ssor_prec)
            {
                SSOR[ill] = SSOR_preconditioner(dia_blocks[ill]);
            }
        }
        std::cout << ((preconditioner != ilu_prec and preconditioner != silu_prec) ? "ok\n" : "\n");
        
        // For all initial states ------------------------------------------- //
        //
        
        for (auto instate : instates)
        {
            int li = std::get<1>(instate);
            int mi = std::get<2>(instate);
            
            // skip angular forbidden states
            bool allowed = false;
            for (int l = abs(li - L); l <= li + L; l++)
                allowed = allowed or ClebschGordan(li,mi,l,0,L,mi);
            if (not allowed)
                continue;
            
            // compute P-overlaps and P-expansion
            cArray Pi_overlaps, Pi_expansion;
            Pi_overlaps = rad.overlapP(ni,li,weightEndDamp(bspline));
            Pi_expansion = rad.S().tocoo().tocsr().solve(Pi_overlaps);
            
            // we may have already computed solution for this state and energy... is it so?
            std::ostringstream cur_oss;
            cur_oss << "psi-" << L << "-" << Spin << "-" << Pi << "-" << ni << "-" << li << "-" << mi << "-" << Ei[ie] << ".hdf";
            if ( current_solution.hdfload(cur_oss.str().c_str()) )
                continue;
            
            // create right hand side
            cArray chi ( coupled_states.size()*Nspline*Nspline );
            
            std::cout << "\tCreate RHS for li = " << li << ", mi = " << mi << "\n";
            
            // for all segments constituting the RHS
            # pragma omp parallel for 
            for (unsigned ill = 0; ill < coupled_states.size(); ill++)
            {
                int l1 = coupled_states[ill].first;
                int l2 = coupled_states[ill].second;
                
                // setup storage
                cArrayView chi_block(chi, ill * Nspline * Nspline, Nspline * Nspline);
                chi_block.fill(0);
                
                // for all allowed angular momenta (by momentum composition) of the projectile
                for (int l = abs(li - L); l <= li + L; l++)
                {
                    // skip wrong parity
                    if ((L + li + l) % 2 != Pi)
                        continue;
                    
                    // (anti)symmetrization
                    int Sign = ((Spin + Pi) % 2 == 0) ? 1. : -1.;
                    
                    // compute energy- and angular momentum-dependent prefactor
                    Complex prefactor = pow(Complex(0.,1.),l) * sqrt(2*M_PI*(2*l+1)) / Complex(ki[ie]); 
                    prefactor *= ClebschGordan(li,mi,l,0,L,mi);
                    if (prefactor == 0.)
                        continue;
                    
                    // pick the correct Bessel function expansion
                    cArrayView Ji_expansion (
                        ji_expansion,
                        Nspline * (ie * (maxell + 1) + l),
                        Nspline
                    );
                    
                    // compute outer products of B-spline expansions
                    cArray Pj1 = outer_product(Pi_expansion, Ji_expansion);
                    cArray Pj2 = outer_product(Ji_expansion, Pi_expansion);
                    
                    // skip angular forbidden right hand sides
                    for (unsigned lambda = 0; lambda <= rad.maxlambda(); lambda++)
                    {
                        Complex f1 = computef(lambda, l1, l2, li, l, L);
                        Complex f2 = computef(lambda, l1, l2, l, li, L);
                        
                        if (f1 != 0.)
                        {
                            chi_block += (prefactor * f1) * rad.R_tr_dia(lambda).dot(Pj1);
                        }
                        
                        if (f2 != 0.)
                        {
                            if (Sign > 0)
                                chi_block += (prefactor * f2) * rad.R_tr_dia(lambda).dot(Pj2);
                            else
                                chi_block -= (prefactor * f2) * rad.R_tr_dia(lambda).dot(Pj2);
                        }
                    }
                    
                    if (li == l1 and l == l2)
                    {
                        // direct contribution
                        chi_block -= prefactor * S_kron_Mm1_tr.dot(Pj1);
                    }
                    
                    if (li == l2 and l == l1)
                    {
                        // exchange contribution with the correct sign
                        if (Sign > 0)
                            chi_block -= prefactor * Mm1_tr_kron_S.dot(Pj2);
                        else
                            chi_block += prefactor * Mm1_tr_kron_S.dot(Pj2);
                    }
                }
            }
            
            // Solve the Equations ------------------------------------------------- //
            //
            
            std::cout << "\tSolve the equations.\n";
            
            // we may have already computed the previous solution - it will serve as initial guess
            current_solution = cArray(chi.size());
            if (ie > 0)
            {
                std::ostringstream prev_oss;
                prev_oss << "psi-" << L << "-" << Spin << "-" << Pi << "-" << ni << "-" << li << "-" << mi << "-" << Ei[ie-1] << ".hdf";
                if (previous_solution.hdfload(prev_oss.str().c_str()))
                    current_solution = previous_solution;
            }
            
            // CG preconditioner callback
            auto apply_preconditioner = [ & ](cArray const & r, cArray & z) -> void
            {
                // apply a block inversion preconditioner (parallel for ILU)
                # pragma omp parallel for if (preconditioner == ilu_prec)
                for (unsigned ill = 0; ill < coupled_states.size(); ill++)
                {
                    // skip computation of unwanted blocks for this process
                    if (LUs.find(ill) == LUs.end() or LUs[ill] != par.iproc())
                        continue;
                    
                    // create copy-to view of "z"
                    cArrayView zview(z, ill * Nspline * Nspline, Nspline * Nspline);
                    
                    // create copy-from view of "r"
                    cArrayView rview(r, ill * Nspline * Nspline, Nspline * Nspline);
                    
                    // preconditioner of the nested CG
                    auto apply_inner_preconditioner = [ & ](cArray const & r, cArray & z) -> void
                    {
                        // Incomplete LU factorization
                        if (preconditioner == ilu_prec)
                            z = iLU[ill].solve(r);
                        
                        // single-electron Incomplete LU factorization
                        // TODO Needs some modification to work!
                        if (preconditioner == silu_prec)
                        {
                            // for all single-electron blocks
                            # pragma omp parallel for
                            for (int iblock = 0; iblock < Nspline; iblock++)
                            {
                                // precondition by inverting a single diagonal block
                                cArrayView rview (r, iblock * Nspline, Nspline);
                                cArrayView zview (z, iblock * Nspline, Nspline);
                                zview = siLU[ill][iblock].solve(rview);
                            }
                        }
                        
                        // block D-ILU
                        if (preconditioner == bilu_prec)
                        {
                            // for all single-electron blocks
                            # pragma omp parallel for
                            for (int iblock = 0; iblock < Nspline; iblock++)
                            {
                                // precondition by inverting a single diagonal block
                                cArrayView rview (r, iblock * Nspline, Nspline);
                                cArrayView zview (z, iblock * Nspline, Nspline);
                                zview = bcks_lufts[ill][iblock].solve(rview);
                            }
                        }
                        
                        // Diagonal Incomplete Cholesky factorization
                        // TODO Needs to implement pivoting to work!
                        if (preconditioner == dic_prec)
                        {
                            z = DIC[ill].upperSolve( DIC[ill].dot( DIC[ill].lowerSolve(r), diagonal ) );
                            std::cout << z.norm()/r.norm() << "\n";
                            std::cout << cArray(DIC[ill].main_diagonal()).norm() << "\n";
                        }
                        
                        // Sparse Approximate Inverse preconditioner
                        if (preconditioner == spai_prec)
                            z = spai[ill].dot(r);
                        
                        // Symmetric Successive Over-Relaxation
                        // NOTE seems slower than Jacobi
                        if (preconditioner == ssor_prec)
                            z = SSOR[ill].upperSolve( SSOR[ill].dot( SSOR[ill].lowerSolve(r), diagonal ) );
                        
                        // Jacobi preconditioning
                        // NOTE seems the fastest
                        if (preconditioner == jacobi_prec)
                            z = SSOR[ill].dot(r,diagonal);
                        
                        // no preconditioning
                        // NOTE seems slower than Jacobi
                        if (preconditioner == no_prec)
                            z = r;
                    };
                    
                    // multiply by matrix block
                    auto inner_matrix_multiply = [ & ](cArray const & p, cArray & q) -> void
                    {
                        q = dia_blocks[ill].dot(p);
                    };
                    
                    // solve using the CG solver
                    cg_callbacks
                    (
                        rview,      // rhs
                        zview,      // solution
                        1e-11,      // tolerance
                        0,          // min. iterations
                        Nspline*Nspline,    // max. iteration
                        apply_inner_preconditioner,
                        inner_matrix_multiply,
                        true       // verbose
                    );
                }
                
#ifndef NO_MPI
                if (parallel)
                {
                    // synchronize across processes
                    for (unsigned ill = 0; ill < coupled_states.size(); ill++)
                    {
                        // relevant process will broadcast this segment's data
                        MPI_Bcast (
                            &z[0] + ill * Nspline * Nspline,
                            Nspline * Nspline,
                            MPI_DOUBLE_COMPLEX,
                            LUs[ill],
                            MPI_COMM_WORLD
                        );
                    }
                }
#endif
            };
            
            // CG matrix multiplication callback
            auto matrix_multiply = [ & ](cArray const & p, cArray & q) -> void
            {
                // clear all output segments that are going to be referenced by this process
                for (unsigned ill = 0; ill < coupled_states.size(); ill++)
                    if (LUs.find(ill) != LUs.end() and LUs[ill] == par.iproc())
                        cArrayView(q, ill * Nspline * Nspline, Nspline * Nspline).fill(0);
                
                // multiply "q" by the matrix of the system
                # pragma omp parallel for schedule (dynamic,1) collapse(2)
                for (unsigned ill = 0; ill < coupled_states.size(); ill++)
                for (unsigned illp = 0; illp < coupled_states.size(); illp++)
                if (LUs.find(ill) != LUs.end() and LUs[ill] == par.iproc())
                {
                    // row multi-index
                    int l1 = coupled_states[ill].first;
                    int l2 = coupled_states[ill].second;
                    
                    // column multi-index
                    int l1p = coupled_states[illp].first;
                    int l2p = coupled_states[illp].second;
                    
                    // product segment contribution
                    cArray q_contrib (Nspline * Nspline);
                    
                    // copy-from segment of "p"
                    cArrayView p_block(p, illp * Nspline * Nspline, Nspline * Nspline);
                    
                    // multiply by hamiltonian terms
                    if (ill == illp)
                    {
                        // reuse the diagonal block
                        q_contrib += dia_blocks[ill].dot(p_block);
                    }
                    else
                    {
                        // compute the offdiagonal block
                        for (unsigned lambda = 0; lambda <= rad.maxlambda(); lambda++)
                        {
                            Complex f = computef(lambda, l1, l2, l1p, l2p, L);
                            if (f != 0.)
                                q_contrib -= f * rad.R_tr_dia(lambda).dot(p_block);
                        }
                    }
                    
                    // safely update shared output array "q"
                    # pragma omp critical
                    cArrayView (q, ill * Nspline * Nspline, Nspline * Nspline) += q_contrib;
                }
                
#ifndef NO_MPI
                if (parallel)
                {
                    // synchronize across processes
                    for (unsigned ill = 0; ill < coupled_states.size(); ill++)
                    {
                        // relevant process will broadcast this segment's data
                        MPI_Bcast (
                            &q[0] + ill * Nspline * Nspline,
                            Nspline * Nspline,
                            MPI_DOUBLE_COMPLEX,
                            LUs[ill],
                            MPI_COMM_WORLD
                        );
                    }
                }
#endif
            };
                
            if (chi.norm() == 0.)
            {
                std::cout << "\t! Right-hand-side is zero (probably due to incompatible angular settings).\n";
                current_solution.clear();
            }
            else
            {
                // custom conjugate gradients callback-based solver
                double tolerance = 1e-10;
                std::cout << "\tStart CG callback with tolerance " << tolerance << "\n";
                unsigned iterations = cg_callbacks (
                    chi,                    // right-hand side
                    current_solution,        // on input, the initial guess, on return, the solution
                    tolerance,                // requested precision, |A·x - b|² < ε·|b|²
                    0,                        // minimal iteration count
                    (maxell+1) * Nspline,    // maximal iteration count
                    apply_preconditioner,    // preconditioner callback
                    matrix_multiply            // matrix multiplication callback
                );
                std::cout << "\tEnd CG callback\n";
                
                // update progress
                iterations_done += iterations;
            }
            
            computations_done++;
            
            // save solution to disk
            current_solution.hdfsave(cur_oss.str().c_str(), true /* = with compression */);
            
        } // end of For li
        
    } // end of For ie = 0, ..., Nenergy - 1
    
    std::cout << "\rSolving the systems... ok                                                            \n";
    std::cout << "\t(typically " << iterations_done/Nenergy << " CG iterations per energy)\n";
    
    // --------------------------------------------------------------------- //
    
}
Stg3:
{
    // skip stage 3 if told so
    if (not (itinerary & StgExtract))
    {
        std::cout << "Skipped extraction of amplitudes.\n";
        goto End;
    }
    
    // compose output filename
    std::ostringstream ossfile;
    if (parallel)
        ossfile << ni << "-" << L << "-" << Spin << "-" << Pi << "-(" << par.iproc() << ").sql";
    else
        ossfile << ni << "-" << L << "-" << Spin << "-" << Pi << ".sql";
    
    // Create SQL batch file
    std::ofstream fsql(ossfile.str().c_str());
    
    // set exponential format for floating point output
    fsql.setf(std::ios_base::scientific);
    
    // write header
    fsql << "BEGIN TRANSACTION;\n";
    
    
    // Extract the cross sections ------------------------------------------ //
    //
    
    std::cout << "\nExtracting T-matrices...\n";
    
    std::vector<std::tuple<int,int,int,int,int>> transitions;
    for (auto instate  : instates)
    for (auto outstate : outstates)
        transitions.push_back (
            std::make_tuple (
                Spin,
                /*li*/ std::get<1>(instate),
                /*mi*/ std::get<2>(instate),
                /*nf*/ std::get<0>(outstate),
                /*lf*/ std::get<1>(outstate)
            )
        );
    
    int finished = 0;
    
    for (int i = 0; i < (int)transitions.size(); i++)
    {
#ifndef NO_MPI
        // if MPI is active, compute only a subset of the transitions, corresponding to iproc
        if (parallel and i % par.Nproc() != par.iproc())
            continue;
#endif
        
        int Spin = std::get<0>(transitions[i]);
        int li   = std::get<1>(transitions[i]);
        int mi   = std::get<2>(transitions[i]);
        int nf   = std::get<3>(transitions[i]);
        int lf   = std::get<4>(transitions[i]);
        
        // skip angular forbidden states
        bool allowed = false;
        for (int l = abs(li - L); l <= li + L; l++)
            allowed = allowed or ClebschGordan(li,mi,l,0,L,mi);
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
                // compute expansions of final wave function for this final state (nf,lf)
                std::transform (
                    Ei.begin(), Ei.end(),
                    kf.begin(),
                    [ = ](double E) -> double { return sqrt(E - 1./(ni*ni) + 1./(nf*nf) + (mf-mi)*B); }
                );
                
                // compute Λ for transitions to (nf,lf,mf); it will depend on [ie,ℓ]
                Lambda[mf+lf] = std::move (
                    computeLambda (bspline, kf, ki, maxell, L, Spin, Pi, ni, li, mi, Ei, lf, Pf_overlaps, coupled_states)
                );
            }
            
            // save the data
            for (int mf = -lf; mf <= lf; mf++)
            {
                // compute Tℓ
                cArray T_ell(Lambda[mf+lf].size());
                for (unsigned i = 0; i < T_ell.size(); i++)
                {
                    int ie  = i / (maxell + 1);
                    int ell = i % (maxell + 1);
                    
                    T_ell[i] = Lambda[mf+lf][i] * 4. * M_PI / kf[ie] * pow(Complex(0.,1.), -ell)
                                    * ClebschGordan(lf, mf, ell, mi - mf, L, mi) / sqrt(2.);
                }
                
                //
                // print out SQL
                //
                
                for (unsigned i = 0; i < T_ell.size(); i++)
                {
                    int ie  = i / (maxell + 1);
                    int ell = i % (maxell + 1);
                    
                    if (finite(T_ell[i].real()) and finite(T_ell[i].imag()))
                    if (T_ell[i].real() != 0. or T_ell[i].imag() != 0.)
                    {
                        fsql << "INSERT OR REPLACE INTO \"tmat\" VALUES ("
                             << ni << "," << li << "," << mi << ","
                             << nf << "," << lf << "," << mf << ","
                             << L  << "," << Spin << ","
                             << Ei[ie] << "," << ell << "," 
                             << T_ell[i].real() << "," << T_ell[i].imag()
                             << ");\n";
                    }
                }
                
                //
                // print out the total cross section for quick overview
                //
                
                std::ostringstream sigmaname;
                sigmaname << "sigma-" 
                          << ni << "-" << li << "-" << mi << "-"
                          << nf << "-" << lf << "-" << mf << "-"
                          << L << "-" << Spin << "-" << Pi << ".dat";

                std::ofstream ftxt(sigmaname.str().c_str());
                
                ftxt << "# Ei [Ry] sigma [a0^2]\n";
                for (unsigned ie = 0; ie < Nenergy; ie++)
                {
                    double sigma = 0.;
                    for (int ell = 0; ell <= maxell; ell++)
                    {
                        double Re_f_ell = -T_ell[ie * (maxell + 1) + ell].real() / (2 * M_PI);
                        double Im_f_ell = -T_ell[ie * (maxell + 1) + ell].imag() / (2 * M_PI);
                        sigma += 0.25 * (2*Spin + 1) * kf[ie] / ki[ie] * (Re_f_ell * Re_f_ell + Im_f_ell * Im_f_ell);
                    }
                    if (finite(sigma))
                    {
                        ftxt << Ei[ie] << "\t" << sigma << "\n";
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
                computeXi(bspline, maxell, L, Spin, Pi, ni, li, mi, Ei, ics, coupled_states)
            );
            
            for (size_t ie = 0; ie < Ei.size(); ie++)
            for (unsigned ill = 0; ill < coupled_states.size(); ill++) //??? or triangular
            {
                // save data as BLOBs
                fsql << "INSERT OR REPLACE INTO \"ionf\" VALUES ("
                     << ni << "," << li << "," << mi << ","
                     << L  << "," << Spin << ","
                     << Ei[ie] << "," << coupled_states[ill].first << ","
                     << coupled_states[ill].second << ","
                     << data[ie * coupled_states.size() + ill].toBlob() << ");\n";
            }
            
            // print ionization cross section
            std::ostringstream fname;
            fname << "isigma-" << L << "-" << Spin << "-" << Pi << "-" << ni << "-" << li << "-" << mi << ".dat";
            write_array(Ei, ics, fname.str().c_str());
        }
        
        finished++;
        
        std::cout << "\rExtracting T-matrices... " 
                  << std::setw(3) << int(trunc(finished * 100. / transitions.size() + 0.5))
                  << " %        ";
    }
    // --------------------------------------------------------------------- //
    
    fsql << "COMMIT;\n";
    fsql.close();
}
End:
{    
#ifndef NO_MPI
    if (parallel)
        MPI_Finalize();
#endif
    
    std::cout << "\nDone.\n\n";
}    
    return 0;
}
