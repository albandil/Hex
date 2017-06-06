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
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <string>
#include <tuple>

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"
#include "hex-cmdline.h"
#include "hex-csrmatrix.h"
#include "hex-matrix.h"
#include "hex-vtkfile.h"

// --------------------------------------------------------------------------------- //

#include "io.h"
#include "preconditioners.h"

// --------------------------------------------------------------------------------- //

#ifdef WITH_OPENCL
    #include <CL/cl.h>
#endif

// --------------------------------------------------------------------------------- //

const std::string sample_input =
    "# This is a sample input file for the program hex-ecs.\n"
    "# Lines introduced by the hash symbol are comments, which document the options,\n"
    "# and can be omitted. Otherwise there are just a few numbers and characters.\n"
    "\n"
    "# ------------ B-spline parameters ----------------\n"
    "\n"
    "# B-spline order.\n"
    "# Recommended value is 4, higher values lead to larger matrices.\n"
    "# Smaller values can lead to insufficent description in coordinate origin.\n"
    "   4\n"
    "\n"
    "# ECS rotation angle in radians.\n"
    "# Typically around 45 degrees, the calculation is totally immune to small changes\n"
    "# in this parameter.\n"
    "   0.63\n"
    "\n"
    "# B-spline knots.\n"
    "# For each section use any of the following options and terminate by -1.\n"
    "#     L[inear] <Start> <End> <Samples>\n"
    "#     G[eometric] <Start> <End> <FirstInterval> <Quotient>\n"
    "# There are one or two sections; the projectile asymptotic extension is optional.\n"
    "#\n"
    "# a) Real knots of the basis that is common to atomic and projectile electron.\n"
    "  L  0.0  0.0   4\n"
    "  G  0.1 10.0  0.1  1.01\n"
    "  L   11  100  90\n"
    " -1\n"
    "# b) Real knots that are exclusive to the projectile, if any. (Start from zero.)\n"
    " -1\n"
    "# c) Complex region knots. (Start from zero.)\n"
    "  G    0   50   1  1.02\n"
    " -1\n"
    "\n"
    "# --------------- Atomic states -------------------\n"
    "\n"
    "# Initial atomic states (ni, li, mi).\n"
    "# Specified as vertical triplets terminated by -1 on the first line.\n"
    "# Computation of all angular quantum numbers can be requested by asterisk.\n"
    "  1 -1\n"
    "  *\n"
    "  *\n"
    "\n"
    "# Final atomic states (nf, lf).\n"
    "# Specified by vertical doublets terminated by -1 on the first line.\n"
    "  1  -1\n"
    "  *\n"
    "\n"
    "# Maximal energy (Ry) of states included in the asymptotic (outer) region.\n"
    "  -1\n"
    "\n"
    "# --------------- Other conditions ----------------\n"
    "\n"
    "# Angular momenta.\n"
    "#   L  ... total orbital momentum\n"
    "#   S  ... total spin\n"
    "#   Pi ... total parity\n"
    "#   nL, limit, exchange ... parameters controlling the number of coupled angular states\n"
    "# The angular basis is composed of coupled angular states (l1,l2) and looks like this:\n"
    "#     (Pi, L),     (Pi+1, L-1), ..., (L, Pi) \n"
    "#     (Pi+1, L+1), (Pi+2, L),   ..., (L+1, Pi+1) \n"
    "#     ...\n"
    "# The number of columns is equal to L + Pi - 1. The number of rows is equal to nL + 1.\n"
    "# When 'limit > 0' then all pairs with both l1 and l2 > limit are discarded.\n"
    "# When 'exchange = 0' then all pairs with l1 > l2 are discarded.\n"
    "# The options 'limit' and 'exchange' are useful for large angular momenta, where the projectile\n"
    "# is distinguishable from the atomic electron.\n"
    "# L   S   Pi  nL  limit exchange\n"
    "  0   *   0   4   -1    1\n"
    "\n"
    "# Projectile charge (+/-1)\n"
    "  -1\n"
    "\n"
    "# Atom + projectile total energies in Rydbergs.\n"
    "# Use any of the following options:\n"
    "#     L[inear] <Start> <End> <Samples>\n"
    "#     G[eometric] <Start> <End> <FirstInterval> <Quotient>\n"
    "#     E[xplicit] <Sample-1> <Sample-2> ... <Sample-last> -1\n"
    "  L  -0.35  -0.05  3\n"
    " -1\n"
    "\n"
    "# Weak magnetic field in atomic units.\n"
    " 0\n";


void CommandLine::parse (int argc, char* argv[])
{
    ParseCommandLine
    (
        argc, argv,
        
        "example", "e", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                std::cout << "Writing sample input file to \"example.inp\".\n\n";
                
                // produce sample input file
                std::ofstream out("example.inp");
                if (out.bad())
                    HexException("Error: Cannot write to \"example.inp\"\n");
                
                out << sample_input;
                    
                out.close();
                std::exit(EXIT_SUCCESS);
            },
        "input", "i", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // set custom input file
                inputfile.open(optargs[0]);
                if (not inputfile.good())
                    HexException("Error: Input file \"%s\" not found.\n", optargs[0].c_str());
                return true;
            },
        "help", "h", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // get all factorizers
                std::ostringstream f;
                for (LUft const * lu : *LUft::RTS_Table)
                    f << " '" << lu->name() << "'";
                
                // print usage information
                std::cout << "\n"
                    "Available switches (short forms in parentheses):                                                                                                          \n"
                    "                                                                                                                                                          \n"
                    "Basic usage                                                                                                                                               \n"
                    "\t--example                  (-e)  Create sample input file.                                                                                              \n"
                    "\t--help                     (-h)  Display this help.                                                                                                     \n"
                    "\t--input <filename>         (-i)  Use custom input file (other than \"ecs.inp\").                                                                        \n"
                    "\t--zip <parameters>         (-z)  Solution file to zip (i.e. evaluate in B-spline basis and produce VTK datafile).                                       \n"
                    "\t                                 The '<parameters>' stands for '<filename> <Xmin> <Ymin> <Xmax> <Ymax> <Xn> <Yn>'.                                      \n"
                    "\t--write-grid               (-g)  Write grid layout to a VTK file.                                                                                       \n"
                    "\t--write-intermediate-solutions   Write all intermediate solution (after every iteration of the PCOCG solver).                                           \n"
                    "\t--carry-initial-guess            Whether to use previous-energy solution as an initial guess for the new energy.                                        \n"
                    "\t--refine-solution                Load existing solutions and check that they are within tolerance, update if needed.                                    \n"
                    "                                                                                                                                                          \n"
#ifdef WITH_MPI
                    "MPI setup                                                                                                                                                 \n"
                    "\t--mpi                      (-m)  Use MPI (assuming that the program has been launched by mpiexec).                                                      \n"
                    "\t--groupsize <number>       (-G)  How many processes factorize single LU (used by 'superlu_dist' and 'mumps' preconditioners).                           \n"
                    "                                                                                                                                                          \n"
#endif
                    "LU decomposition                                                                                                                                          \n"
                    "\t--lu <name>                (-F)  Factorization library (one of" + f.str() + "). Default is 'umfpack'.\n"
#ifdef WITH_MUMPS
                    "\t--mumps-out-of-core              Use out-of-core capability of MUMPS (this is independent on --out-of-core option).                                     \n"
                    "\t--mumps-verbose <number>         Verbosity level of the MUMPS library. Zero ('0') means no output, higher numbers increase the verbosity.               \n"
#endif
                    "                                                                                                                                                          \n"
                    "Stage selection                                                                                                                                           \n"
                    "\t--stg-integ                (-a)  Only calculate needed radial integrals.                                                                                \n"
                    "\t--stg-integ-solve          (-b)  Only calculate integrals and the solution.                                                                             \n"
                    "\t--stg-extract              (-c)  Only extract amplitudes (assumes that the solution files exist).                                                       \n"
                    "                                                                                                                                                          \n"
                    "Extended grid parameters                                                                                                                                  \n"
                    "\t--channel-max-E <number>         Maximal energy (Ry) of states considered in the outer region.                                                          \n"
                    "                                                                                                                                                          \n"
                    "Right-hand side                                                                                                                                           \n"
                    "\t--analytic-eigenstates           Use analytic formulae for initial/final states instead of diagonalization.                                             \n"
                    "\t--fast-bessel                    Use faster Bessel function evaluation routine (not the Steed/Barnett) when calculating RHS.                            \n"
                    "                                                                                                                                                          \n"
                    "Disk access                                                                                                                                               \n"
                    "\t--own-radial-cache         (-w)  Keep two-electron radial integrals not referenced by preconditioner only on disk (slows down only the initialization). \n"
                    "\t--no-radial-cache          (-r)  Keep all two-electron radial integrals only on disk (slows down also the solution process).                            \n"
                    "\t--out-of-core              (-o)  Use hard disk drive to store most of intermediate data and thus to save RAM (considerably slower).                     \n"
                    "\t--out-of-core-continue     (-O)  Start solution from the existing OOC files.                                                                            \n"
                    "\t--whole-matrix             (-W)  In the above three cases: Load whole matrix from scratch file when calculating dot product (speeds them up a little).  \n"
                    "\t--shared-scratch           (-s)  Let every MPI process calculate only a subset of shared radial integrals (assume shared output directory).             \n"
//                     "\t--lightweight-radial-cache (-l)  Do not precalculate two-electron integrals and only apply them on the fly (slower, but saves RAM).                     \n"
                    "\t--lightweight-full         (-L)  Avoid precalculating all large matrices and only apply them on the fly (only available for KPA preconditioner).        \n"
                    "                                                                                                                                                          \n"
                    "Preconditioners (general)                                                                                                                                 \n"
                    "\t--preconditioner <name>    (-p)  Preconditioner to use (default: ILU).                                                                                  \n"
                    "\t--list-preconditioners     (-P)  List available preconditioners with short description of each.                                                         \n"
                    "\t--ssor <number>                  Apply SSOR coupling.                                                                                                   \n"
                    "\t--tolerance <number>       (-T)  Set tolerance for the conjugate gradients solver (default: 1e-8).                                                      \n"
                    "\t--prec-tolerance <number>  (-t)  Set tolerance for the conjugate gradients preconditioner (default: 1e-8).                                              \n"
                    "\t--drop-tolerance <number>  (-d)  Set drop tolerance for the ILU preconditioner (default: 1e-15).                                                        \n"
                    "\t--sub-prec-iter <number>         Maximal number of inner preconditioner iterations before giving up. Default is number of unknowns.                     \n"
                    "\t--sub-prec-nofail                Makes the solver continue with the calculation even when some of the inner preconditioners fails to converge.          \n"
                    "\t--sub-prec-verbose               Display convergence information for inner preconditioners.                                                             \n"
#ifndef DISABLE_PARALLEL_PRECONDITION
                    "\t--parallel-precondition          Apply multiple block preconditioners in parallel.                                                                      \n"
#endif
                    "                                                                                                                                                          \n"
                    "ILU preconditioner                                                                                                                                        \n"
                    "\t--parallel-factorization         Factorize multiple blocks simultaneously.                                                                              \n"
                    "\t--no-lu-update                   Do not recalculate LU factorization for different energies, use the first factorization for all of them.               \n"
                    "\t--ilu-max-iter <number>          Maximal number of iterations of the nested ILU preconditioner. When the number is exceeded, exception is thrown.       \n"
                    "\t--scratch <path>                 Scratch directory for out-of-core factorizers (currently only MUMPS). Also read from the $SCRATCHDIR env variable.     \n"
                    "                                                                                                                                                          \n"
                    "KPA preconditioner                                                                                                                                        \n"
                    "\t--kpa-simple-rad           (-R)  Use simplified radial integral matrix for nested KPA iterations (experimental).                                        \n"
                    "                                                                                                                                                          \n"
                    "HYB preconditioner                                                                                                                                        \n"
                    "\t--hyb-additional-levels <number> When using the HYB preconditioner: precondition more blocks with ILU, useful close below an excitation threshold.      \n"
                    "                                                                                                                                                          \n"
                    "Coupled preconditioner                                                                                                                                    \n"
                    "\t--coupling-limit                 Maximal multipole to be considered by the coupled preconditioner.                                                      \n"
                    "                                                                                                                                                          \n"
#ifdef WITH_OPENCL
                    "GPU preconditioner                                                                                                                                        \n"
                    "\t--cl-list                        List all OpenCL platforms and devices available.                                                                       \n"
                    "\t--cl-platform <index>            Use given OpenCL platform for GPU preconditioner (default is 0, i.e. the first platform found).                        \n"
                    "\t--cl-device <index>              Use given OpenCL device for GPU preconditioner (default is 0, i.e. the first device found).                            \n"
                    "\t--cl-use-host-memory             Keep large data in RAM instead of copying everything to the compute device. This will slow down the solution.          \n"
                    "\t--cl-multiply                    Do the sparse matrix multiplication on the OpenCL device (memory intensive!).                                          \n"
                    "\t--cl-host-multiply               Keep vectors in host memory when doing matrix multiplication on OpenCL device.                                         \n"
#endif
                    "                                                                                                                                                          \n"
                    "Multigrid preconditioner                                                                                                                                  \n"
                    "\t--multigrid-depth <number>       Depth of the geometric multigrid (= maximal refinement level of the original B-spline knot sequence).                  \n"
                    "\t--multigrid-coarse-prec <name>   What preconditioner to use for preconditioning of the solution of the coarse problem (specified in input file).        \n"
                    "                                                                                                                                                          \n"
                    "Domain decomposition preconditioner                                                                                                                       \n"
                    "\t--dom-xpanels <number>           Number of domain decomposition panels along x axis (default: 1).                                                       \n"
                    "\t--dom-ypanels <number>           Number of domain decomposition panels along y axis (default: 1).                                                       \n"
                    "\t--dom-preconditioner <name>      Panel preconditioner for domain decomposition (default: ILU).                                                          \n"
                    "                                                                                                                                                          \n"
                    "Post-processing                                                                                                                                           \n"
                    "\t--no-parallel-extraction         Disallow parallel extraction of T-matrices (e.g. when the whole solution does not fit into the memory).                \n"
                    "\t--extract-rho-begin              Where to start averaging / extrapolating the T-matrix.                                                                 \n"
                    "\t--extract-rho[-end]              Radius for T-matrix extraction.                                                                                        \n"
                    "\t--extract-samples                Number of evaluations of the T-matrix between --extract-rho-begin and --extract-rho.                                   \n"
                    "\t--extract-extrapolate            Radially extrapolate the extracted T-matrices instead of simple averaging.                                             \n"
                    "\t--runtime-postprocess            Evaluate T-matrices after every iteration.                                                                             \n"
                    "\n"
                ;
                std::exit(EXIT_SUCCESS);
            },
        "write-grid", "g", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // write grid to VTK
                writegrid = true;
                return true;
            },
        "zip", "z", 7, [&](std::vector<std::string> const & optargs) -> bool
            {
                // zip B-spline expansion file
                zipdata.file = optargs[0];
                zipdata.Xmin = std::stod(optargs[1]);
                zipdata.Ymin = std::stod(optargs[2]);
                zipdata.Xmax = std::stod(optargs[3]);
                zipdata.Ymax = std::stod(optargs[4]);
                zipdata.nX = std::stoi(optargs[5]);
                zipdata.nY = std::stoi(optargs[6]);
                return true;
            },
        "lu", "F", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // choose factorizer
                factorizer = optargs[0];
                return true;
            },
        "groupsize", "G", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // set the groupsize
                groupsize = std::atoi(optargs[0].c_str());
                return true;
            },
#ifdef WITH_MPI
        "mpi", "m", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // use MPI
                parallel = true;
                return true;
            },
#endif
        "stg-integ", "a", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // run only the first part (computation of radial integrals)
                itinerary = StgRadial;
                return true;
            },
        "stg-integ-solve", "b", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // run only the first part (computation of radial integrals)
                itinerary = StgRadial | StgSolve;
                return true;
            },
        "stg-extract", "c", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // run only the third part (extraction of amplitudes)
                itinerary = StgExtract;
                return true;
            },
        "own-radial-cache", "w", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // do not cache un-'owned' two-electron radial integrals in memory
                cache_own_radint = true;
                cache_all_radint = false;
                return true;
            },
        "no-radial-cache", "r", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // do not cache any two-electron radial integrals in memory at all
                cache_own_radint = false;
                cache_all_radint = false;
                return true;
            },
        "scratch", "", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // scratch directory for the out-of-core mode factorizers (currently MUMPS only)
                scratch = optargs[0];
                return true;
            },
        "out-of-core", "o", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // use full out-of-core functionality: store also diagonal blocks (and factorizations) on disk
                cache_all_radint = false;
                cache_own_radint = false;
                outofcore = true;
                return true;
            },
        "out-of-core-continue", "O", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // continue OOC calculation
                cache_all_radint = false;
                cache_own_radint = false;
                reuse_dia_blocks = true;
                outofcore = true;
                cont = true;
                return true;
            },
        "whole-matrix", "W", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // whether to load whole matrix from scratch disk when calculating dot product (etc.)
                wholematrix = true;
                return true;
            },
        "drop-tolerance", "d", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // drop tolerance for iLU-factorization
                droptol = std::atof(optargs[0].c_str());
                return true;
            },
        "tolerance", "T", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // iteration tolerance for terminating iteration solution
                itertol = std::atof(optargs[0].c_str());
                return true;
            },
        "prec-tolerance", "t", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // iteration tolerance for terminating iteration solution
                prec_itertol = std::atof(optargs[0].c_str());
                return true;
            },
        "preconditioner", "p", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // look-up the preconditioner
                std::vector<PreconditionerBase*>::const_iterator ip = std::find_if
                (
                    PreconditionerBase::RTS_Table->begin(),
                    PreconditionerBase::RTS_Table->end(),
                    [&](PreconditionerBase* ptr)
                    {
                        return ptr->name() == optargs[0];
                    }
                );
                
                if (ip == PreconditionerBase::RTS_Table->end())
                {
                    HexException("Unknown preconditioner \"%s\".", optargs[0].c_str());
                }
                
                preconditioner = optargs[0];
                return true;
            },
        "list-preconditioners", "P", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // look-up the preconditioners description
                std::cout << "\nPreconditioners description (\"ILU\" is the default):\n\n";
                for (PreconditionerBase const * ip : *PreconditionerBase::RTS_Table)
                {
                    std::cout << ip->name() << "\n";
                    std::cout << "\t" << ip->description() << "\n";
                }
                std::cout << std::endl;
                std::exit(EXIT_SUCCESS);
            },
#ifndef DISABLE_PARALLEL_PRECONDITION
        "parallel-precondition", "", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // parallelize preconditioning
                parallel_precondition = true;
                return true;
            },
#endif
        "lightweight-radial-cache", "l", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // do not precompute two-electron radial integral matrices but only apply them on the fly
                lightweight_radial_cache = true;
                return true;
            },
        "lightweight-full", "L", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // do not precompute large matrices but only apply them on the fly
                lightweight_full = lightweight_radial_cache = true;
                return true;
            },
        "write-intermediate-solutions", "", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // write intermediate solutions
                write_intermediate_solutions = true;
                return true;
            },
        "kpa-simple-rad", "R", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // use simplified radial matrix for KPA nested iterations
                kpa_simple_rad = true;
                return true;
            },
        "hyb-additional-levels", "", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // when using the HYB preconditioner: precondition more blocks with ILU, useful close below an excitation threshold
                hyb_additional_levels = std::stoi(optargs[0]);
                return true;
            },
        "no-lu-update", "", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // do not recalculate LU
                noluupdate = true;
                return true;
            },
        "ilu-max-iter", "", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // maximal number of ILU preconditioner iterations
                ilu_max_iter = std::atoi(optargs[0].c_str());
                return true;
            },
        "sub-prec-iter", "", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // maximal number of preconditioner iterations
                max_sub_iter = std::atoi(optargs[0].c_str());
                return true;
            },
        "sub-prec-nofail", "", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // continue calculation even when some of the preconditioners fails to converge
                fail_on_sub_iter = false;
                return true;
            },
        "sub-prec-verbose", "", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // verbose inner preconditioners
                sub_prec_verbose = true;
                return true;
            },
#ifdef WITH_MUMPS
        "coupling-limit", "", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // maximal multipole to be considered by the coupled preconditioner
                coupling_limit = std::atoi(optargs[0].c_str());
                return true;
            },
        "mumps-out-of-core", "", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // MUMPS out of core
                mumps_outofcore = true;
                return true;
            },
        "mumps-verbose", "", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // maximal multipole to be considered by the coupled preconditioner
                mumps_verbose = std::atoi(optargs[0].c_str());
                return true;
            },
#endif
        "shared-scratch", "s", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // precompute only the owned subset of radial integrals
                shared_scratch = true;
                return true;
            },
#ifdef WITH_OPENCL
        "cl-list", "", 0,  [&](std::vector<std::string> const & optargs) -> bool
            {
                // list all available OpenCL platforms and devices
                std::cout << "Available OpenCL devices" << std::endl;
                char platform_name[1024], platform_vendor[1024], platform_version[1024];
                char device_name[1024], device_vendor[1024];
                cl_platform_id platforms[10];
                cl_device_id devices[10];
                cl_uint nplatforms, ndevices;
                clGetPlatformIDs(10, platforms, &nplatforms);
                for (cl_uint i = 0; i < nplatforms; i++)
                {
                    clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, sizeof(platform_name), platform_name, nullptr);
                    clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR, sizeof(platform_vendor), platform_vendor, nullptr);
                    clGetPlatformInfo(platforms[i], CL_PLATFORM_VERSION, sizeof(platform_version), platform_version, nullptr);
                    std::cout << "\t- Platform " << i << ": " << platform_name << " (" << platform_vendor << ", " << platform_version << ")" << std::endl;
                    clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 10, devices, &ndevices);
                    for (cl_uint j = 0; j < ndevices; j++)
                    {
                        clGetDeviceInfo(devices[j], CL_DEVICE_NAME, sizeof(device_name), device_name, nullptr);
                        clGetDeviceInfo(devices[j], CL_DEVICE_VENDOR, sizeof(device_vendor), device_vendor, nullptr);
                        std::cout << "\t\t- Device " << j << ": " << device_name << " (" << device_vendor << ")" << std::endl;
                    }
                }
                std::cout << std::endl;
                std::exit(0);
            },
        "cl-platform", "", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // use given OpenCL platform
                ocl_platform = std::atol(optargs[0].c_str());
                return true;
            },
        "cl-device", "", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // use given OpenCL device
                ocl_device = std::atol(optargs[0].c_str());
                return true;
            },
        "cl-use-host-memory", "", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // keep large data in RAM
                gpu_large_data = true;
                gpu_host_multiply = true;
                return true;
            },
        "cl-multiply", "", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // do the lightweight multiplication on OpenCL device
                lightweight_radial_cache = true;
                gpu_multiply = true;
                return true;
            },
        "cl-host-multiply", "", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // keep vectors in host memory when doing matrix multiplication on GPU
                gpu_host_multiply = true;
                return true;
            },
#endif
        "parallel-factorization", "", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // allow multiple factorizations at a time
                parallel_factorization = true;
                return true;
            },
        "no-parallel-extraction", "", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // disallow parallel extraction
                parallel_extraction = false;
                return true;
            },
        "carry-initial-guess", "", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // use previous solution as an initial guess
                carry_initial_guess = true;
                return true;
            },
        "refine-solution", "", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // refine solution
                refine_solution = true;
                return true;
            },
        "extract-extrapolate", "", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // radially extrapolate extracted T-matrices
                extract_extrapolate = true;
                return true;
            },
        "extract-rho-begin", "", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // beginning of averaging/extrapolation window
                extract_rho_begin = std::atof(optargs[0].c_str());
                return true;
            },
        "extract-rho-end", "", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // end of averaging/extrapolation window
                extract_rho = std::atof(optargs[0].c_str());
                return true;
            },
        "extract-rho", "", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // end of averaging/extrapolation window
                extract_rho = std::atof(optargs[0].c_str());
                return true;
            },
        "extract-samples", "", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // number of extrapolation/averaging samples
                extract_samples = std::atoi(optargs[0].c_str());
                return true;
            },
        "ssor", "", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // apply SSOR coupling
                ssor = std::stod(optargs[0]);
                return true;
            },
        "fast-bessel", "", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // use faster Bessel function evaluation routine (not the Steed/Barnett) when calculating RHS
                fast_bessel = true;
                return true;
            },
        "multigrid-depth", "", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // multigrid levels
                multigrid_depth = std::stoi(optargs[0]);
                return true;
            },
        "multigrid-coarse-prec", "", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // look-up the preconditioner
                std::vector<PreconditionerBase*>::const_iterator ip = std::find_if
                (
                    PreconditionerBase::RTS_Table->begin(),
                    PreconditionerBase::RTS_Table->end(),
                    [&](PreconditionerBase* ptr)
                    {
                        return ptr->name() == optargs[0];
                    }
                );
                
                // check that it exists
                if (ip == PreconditionerBase::RTS_Table->end())
                    HexException("Unknown coarse preconditioner");
                
                return true;
            },
        "dom-xpanels", "", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // domain decomposition panels (x axis)
                dom_x_panels = std::stoi(optargs[0]);
                
                // check that the number make sense
                if (dom_x_panels < 1)
                    HexException("There must be at least one DOM panel in direction of X axis.");
                
                return true;
            },
        "dom-ypanels", "", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // domain decomposition panels (y axis)
                dom_y_panels = std::stoi(optargs[0]);
                
                // check that the number make sense
                if (dom_y_panels < 1)
                    HexException("There must be at least one DOM panel in direction of Y axis.");
                
                return true;
            },
        "dom-preconditioner", "", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // domain decomposition panels (y axis)
                dom_preconditioner = optargs[0];
                
                // look-up the preconditioner
                std::vector<PreconditionerBase*>::const_iterator ip = std::find_if
                (
                    PreconditionerBase::RTS_Table->begin(),
                    PreconditionerBase::RTS_Table->end(),
                    [&](PreconditionerBase* ptr)
                    {
                        return ptr->name() == dom_preconditioner;
                    }
                );
                
                // check that it exists
                if (ip == PreconditionerBase::RTS_Table->end())
                    HexException("Unknown domain decomposition panel preconditioner");
                
                return true;
            },
        "analytic-eigenstates", "", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // use analytic eigenstates instead of those obtained from the diagonalization
                analytic_eigenstates = true;
                return true;
            },
        "runtime-postprocess", "", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // use analytic eigenstates instead of those obtained from the diagonalization
                runtime_postprocess = true;
                return true;
            },
        
        [&] (std::string optname, std::vector<std::string> const & optargs) -> bool
        {
            HexException("Unknown switch \"%s\".", optname.c_str());
        }
    );
}

void ReadArrays (std::ifstream & inf, rArray & arr)
{
    arr.resize(0);
    std::string type;
    while ((type = ReadNext<std::string>(inf).val) != std::string("-1"))
    {
        if (type[0] == 'L')
        {
            Real begin = ReadNext<Real>(inf).val;
            Real end = ReadNext<Real>(inf).val;
            int samples = ReadNext<int>(inf).val;
            
            if (begin > end)
                HexException("Start of linear sequence is larger than its end (%g > %g).", begin, end);
            if (samples < 0)
                HexException("Invalid number of samples for linear sequence: %d.", samples);
            
            arr.append(linspace(begin, end, samples));
        }
        else if (type[0] == 'G')
        {
            Real begin = ReadNext<Real>(inf).val;
            Real end = ReadNext<Real>(inf).val;
            Real d = ReadNext<Real>(inf).val;
            Real quotient = ReadNext<Real>(inf).val;
            int samples = std::ceil(1 + std::log(1 + (end - begin) * (quotient - 1) / d) / std::log(quotient));
            
            if (begin > end)
                HexException("Start of geometric sequence is larger than its end (%g > %g).", begin, end);
            if (d <= 0)
                HexException("Initial interval size must be larger than zero (given: %g).", d);
            if (quotient <= 0)
                HexException("Quotient must be positive (given: %g).", quotient);
            
            arr.append(geomspace(begin, end, samples, quotient));
        }
        else if (type[0] == 'E')
        {
            // explicit list samples
            Real X;
            while ((X = ReadNext<Real>(inf).val) != -1)
                arr.push_back(X);
        }
        else
        {
            HexException("Unknown sequence type \"%s\".", type.c_str());
        }
    }
}

void InputFile::read (std::ifstream & inf)
{
    ReadItem<int> idata;
    
    // load B-spline parameters
    order = ReadNext<int>(inf).val;
    ecstheta = ReadNext<Real>(inf).val;
    
    std::cout << std::endl;
    std::cout << "B-spline basis" << std::endl;
    std::cout << "\torder = " << order << std::endl;
    std::cout << "\tecs angle = " << ecstheta << std::endl;
    
    //
    // load real atomic knot data
    //
    
    ReadArrays(inf, rknots);
    
    // print info
    std::cout << std::endl;
    if (rknots.size() < 2000)
    {
        std::cout << "Real knots (" << rknots.size() << ")" << std::endl;
        for (std::string line : rknots.lines(100))
            std::cout << '\t' << line << std::endl;
    }
    else
    {
        std::cout << "Real knots (first 1000 of " << rknots.size() << ")" << std::endl;
        for (std::string line : rknots.slice(0, 1000).lines(100))
            std::cout << '\t' << line << std::endl;
        std::cout << std::endl;
        std::cout << "Real knots (last 1000 of " << rknots.size() << ")" << std::endl;
        for (std::string line : rknots.slice(rknots.size() - 1000, rknots.size()).lines(100))
            std::cout << '\t' << line << std::endl;
    }
    
    // check order of knots
    for (unsigned i = 1; i < rknots.size(); i++)
        if (rknots[i] < rknots[i-1])
            HexException("The real knot sequence is not monotonous.");
    
    //
    // load real projectile extension knot data
    //
    
    ReadArrays(inf, rknots_ext);
    
    // print info
    std::cout << std::endl;
    if (rknots_ext.size() < 2000)
    {
        std::cout << "Extension knots (" << rknots_ext.size() << ")" << std::endl;
        for (std::string line : rknots_ext.lines(100))
            std::cout << '\t' << line << std::endl;
    }
    else
    {
        std::cout << "Extension knots (first 1000 of " << rknots_ext.size() << " total)" << std::endl;
        for (std::string line : rknots_ext.slice(0, 1000).lines(100))
            std::cout << '\t' << line << std::endl;
        std::cout << std::endl;
        std::cout << "Extension knots (last 1000 of " << rknots_ext.size() << " total)" << std::endl;
        for (std::string line : rknots_ext.slice(rknots_ext.size() - 1000, rknots_ext.size()).lines(100))
            std::cout << '\t' << line << std::endl;
    }
    
    // check that the first knot is zero
    if (rknots_ext.size() > 0 and rknots_ext[0] != 0.)
        HexException("The first knot in overlap region must be zero.");
    
    // check order of knots
    for (unsigned i = 1; i < rknots_ext.size(); i++)
        if (rknots_ext[i] < rknots_ext[i-1])
            HexException("The overlap knot sequence is not monotonous.");
    
    // determine whether only the inner problem is to be solved (i.e. no projectile extension grid)
    inner_only = rknots_ext.empty();
    
    //
    // load complex knot data
    //
    
    ReadArrays(inf, cknots);
    
    // print info
    std::cout << std::endl;
    if (cknots.size() < 2000)
    {
        std::cout << "Complex knots (before scaling; " << cknots.size() << ")" << std::endl;
        for (std::string line : cknots.lines(100))
            std::cout << '\t' << line << std::endl;
    }
    else
    {
        std::cout << "Complex knots (before scaling; first 1000 of " << cknots.size() << ")" << std::endl;
        for (std::string line : cknots.slice(0, 1000).lines(100))
            std::cout << '\t' << line << std::endl;
        std::cout << std::endl;
        std::cout << "Complex knots (before scaling; last 1000 of " << cknots.size() << ")" << std::endl;
        for (std::string line : cknots.slice(rknots.size() - 1000, rknots.size()).lines(100))
            std::cout << '\t' << line << std::endl;
    }
    
    // check that the first knot is zero
    if (cknots.size() > 0 and cknots[0] != 0.)
        HexException("The first knot in complex region must be zero.");
    
    // check order of knots
    for (unsigned i = 1; i < cknots.size(); i++)
        if (cknots[i] < cknots[i-1])
            HexException("The complex knot sequence is not monotonous.");
    
    //
    // load initial atomic quantum numbers
    //
    
    std::cout << std::endl;
    std::cout << "Initial atomic states" << std::endl;
    
    std::vector<ReadItem<int>> nis, lis, mis;
    
    // load initial principal quantum numbers
    while ((idata = ReadNext<int>(inf)).val != -1)
        nis.push_back(idata);
    
    // - orbital angular momentum
    for (std::size_t i = 0; i < nis.size(); i++)
        lis.push_back(ReadNext<int>(inf, ReadItem<int>::asterisk));
    
    // - magnetic quantum number
    for (std::size_t i = 0; i < nis.size(); i++)
        mis.push_back(ReadNext<int>(inf, ReadItem<int>::asterisk));
    
    for (unsigned i = 0; i < nis.size(); i++)
    for (int li = 0; li < nis[i].val; li++)
    for (int mi = -li; mi <= li; mi++)
    {
        // skip unused orbital angular momenta
        if (lis[i].val != li and not (lis[i].flags & ReadItem<int>::asterisk))
            continue;
        
        // skip unused angular momentum projections
        if (mis[i].val != mi and not (mis[i].flags & ReadItem<int>::asterisk))
            continue;
        
        // skip negative projections if asterisk was used (symmetry)
        if (mi < 0 and (mis[i].flags & ReadItem<int>::asterisk))
            continue;
        
        // add this initial state
        instates.push_back(std::make_tuple(nis[i].val,li,mi));
    }
    
    // print info
    std::cout << "\t[n l m]: ";
    for (auto state : instates)
    {
        std::cout << "["
                  << std::get<0>(state) << " "
                  << std::get<1>(state) << " "
                  << std::get<2>(state)
                  << "] ";
    }
    std::cout << std::endl;
    
    //
    // load final atomic quantum numbers
    //
    
    std::cout << std::endl;
    std::cout << "Final atomic states" << std::endl;
    std::vector<ReadItem<int>> nfs, lfs;
    
    // - principal quantum number
    while ((idata = ReadNext<int>(inf)).val != -1)
        nfs.push_back(idata);
    
    // - orbital angular momentum
    for (std::size_t i = 0; i < nfs.size(); i++)
        lfs.push_back(ReadNext<int>(inf, ReadItem<int>::asterisk));
    
    // - construct final states
    for (unsigned f = 0; f < nfs.size(); f++)
    for (int lf = 0; lf <= nfs[f].val; lf++)
    {
        // l=n only in ionization specification
        if (lf == nfs[f].val and nfs[f].val != 0)
            continue;
        
        // skip unused orbital angular momenta
        if (lfs[f].val != lf and not (lfs[f].flags & ReadItem<int>::asterisk))
            continue;
        
        // add this initial state
        outstates.push_back(std::make_tuple(nfs[f].val,lf,0));
    }
    
    // - asymptotic channel max energy
    channel_max_E = ReadNext<Real>(inf).val;
    
    // print info
    std::cout << "\t[n l m]: ";
    for (auto state : outstates)
    {
        std::cout << "["
                  << std::get<0>(state) << " "
                  << std::get<1>(state)
                  << " *] ";
    }
    std::cout << "\n" << std::endl;
    
    std::cout << "Asymptotic channels with energy up to " << channel_max_E << " Ry." << std::endl;
    
    //
    // load total quantum numbers etc.
    //
    
    std::cout << std::endl;
    std::cout << "Angular momentum limits" << std::endl;
    
    // total angular momentum
    L = ReadNext<int>(inf).val;
    
    // spin
    ReadItem<int> Sp = ReadNext<int>(inf, ReadItem<int>::asterisk);
    if ((Sp.flags & ReadItem<int>::asterisk) or Sp.val == 0) Spin.push_back(0);
    if ((Sp.flags & ReadItem<int>::asterisk) or Sp.val == 1) Spin.push_back(1);
    
    // parity
    Pi = ReadNext<int>(inf).val % 2;
    
    // number of angular momentum pairs per total angular momentum
    levels = ReadNext<int>(inf).val;
    
    // single-electron angular momentum limit
    limit = ReadNext<int>(inf).val;
    
    // whether to include also l1 > l2, or only l1 <= l2
    exchange = ReadNext<int>(inf).val;
    
    std::cout << "\tL = " << L << std::endl;
    std::cout << "\tS = " << Spin << std::endl;
    std::cout << "\tPi = " << Pi << std::endl;
    std::cout << "\tnL = " << levels << std::endl;
    std::cout << "\tlimit = " << limit << std::endl;
    
    Za = 1;
    Zp = ReadNext<int>(inf).val;
    
    if (Zp == 0)
        HexException("Invalid projectile charge 0. Use +1 or -1 for positron and electron respectively.");
    
    Zp = Zp / std::abs(Zp);
    
    std::cout << "\tZp = " << Zp << std::endl;
    
    //
    // load initial energies
    //
    
    ReadArrays(inf, Etot);
    
    // print info
    std::cout << std::endl << "Total projectile + atom energies [Ry]" << std::endl;
    std::cout << "\tcount: "     << Etot.size()  << std::endl;
    std::cout << "\tfull list: " << Etot         << std::endl;
    
    // check that all energies are allowed by the asymptotic expansion
    max_Etot = -1;
    if (not inner_only and not Etot.empty())
    {
        max_Etot = *std::max_element(Etot.begin(), Etot.end());
        
        if (max_Etot > channel_max_E)
        {
            std::cout << std::endl;
            std::cout << "Warning: The maximal asymptotic channel energy Easy = " << channel_max_E << " Ry " << std::endl;
            std::cout << "         is too low to cover all possible open channels for the given maximum total" << std::endl;
            std::cout << "         energy Etot = " << max_Etot << " Ry. Program will use all energetically allowed" << std::endl;
            std::cout << "         asymptotic channels in cases with Etot > Easy." << std::endl;

        }
    }
    
    //
    // load some other optional data
    //
    
    std::cout << std::endl;
    std::cout << "Other parameters" << std::endl;
    B = ReadNext<Real>(inf).val;
    
    // print info
    std::cout << "\tmagnetic field: " << B << " a.u." << std::endl;
    std::cout << std::endl;
}

void zip_solution
(
    CommandLine const & cmd,
    InputFile const & inp,
    Parallel const & par, 
    Bspline const & bspline_inner,
    Bspline const & bspline_full,
    std::vector<std::pair<int,int>> const & ll
)
{
    cArray sol;     // stored solution expansion
    cArray ev;      // evaluated solution
    rArray grid_x;  // real evaluation grid (atomic electron)
    rArray grid_y;  // real evaluation grid (projectile electron)
    
    // shorthands
    std::size_t Nspline_inner = bspline_inner.Nspline();
    std::size_t Nspline_full  = bspline_full .Nspline();
    std::size_t Nspline_outer = Nspline_full - Nspline_inner;
    
    std::cout << "Zipping B-spline expansion of the solution: \"" << cmd.zipdata.file << "\"" << std::endl;
    
    // evaluation grid
    Real Xmin = cmd.zipdata.Xmin < 0 ? 0 : cmd.zipdata.Xmin;
    Real Ymin = cmd.zipdata.Ymin < 0 ? 0 : cmd.zipdata.Ymin;
    Real Xmax = cmd.zipdata.Xmax < 0 ? bspline_full.Rmax() : cmd.zipdata.Xmax;
    Real Ymax = cmd.zipdata.Ymax < 0 ? bspline_full.Rmax() : cmd.zipdata.Ymax;
    grid_x = linspace(Xmin, Xmax, cmd.zipdata.nX);
    grid_y = linspace(Ymin, Ymax, cmd.zipdata.nY);
    
    // load the requested file
    HDFFile hdf (cmd.zipdata.file, HDFFile::readonly);
    if (not hdf.valid())
        HexException("Cannot load file %s.", cmd.zipdata.file.c_str());
    
    // get solution information
    int l1, l2, Nchan1 = 0, Nchan2 = 0; Real E;
    if (not hdf.read("l1", &l1, 1) or
        not hdf.read("l2", &l2, 1) or
        not hdf.read("E", &E, 1) or
        not sol.hdfload(cmd.zipdata.file))
        HexException("This is not a valid solution file.");
    
    // get number of asymptotic channels
    hdf.read("Nchan1", &Nchan1, 1);
    hdf.read("Nchan2", &Nchan2, 1);
    
    // prepare radial integrals structure
    RadialIntegrals r (bspline_inner, bspline_full, 0);
    r.verbose(false);
    r.setupOneElectronIntegrals(par, cmd);
    
    // prepare quadrature structure
    GaussLegendre g_inner;
    
    // factorize the overlap matrix
    CsrMatrix<LU_int_t,Complex> S_csr = r.S_x().tocoo<LU_int_t>().tocsr();
    std::shared_ptr<LUft> S_lu;
    S_lu.reset(LUft::Choose("lapack"));
    S_lu->factorize(S_csr);
    
    // compute all needed bound states
    cArrays Xp1 (Nchan2), Sp1 (Nchan2), Xp2 (Nchan1), Sp2 (Nchan1);
    for (int n1 = l1 + 1; n1 <= l1 + Nchan2; n1++)
    {
        Sp1[n1 - l1 - 1] = r.overlapP(bspline_inner, g_inner, inp.Za, n1, l1);
        Xp1[n1 - l1 - 1] = S_lu->solve(Sp1[n1 - l1 - 1]);
    }
    for (int n2 = l2 + 1; n2 <= l2 + Nchan1; n2++)
    {
        Sp2[n2 - l2 - 1] = r.overlapP(bspline_inner, g_inner, inp.Za, n2, l2);
        Xp2[n2 - l2 - 1] = S_lu->solve(Sp2[n2 - l2 - 1]);
    }
    
    // expand the solution
    cArray full_solution (Nspline_full * Nspline_full, 0.0_z);
    for (std::size_t i = 0; i < Nspline_full; i++)
    for (std::size_t j = 0; j < Nspline_full; j++)
    {
        if (i < Nspline_inner and j < Nspline_inner)
        {
            full_solution[i * Nspline_full + j] = sol[i * Nspline_inner + j];
        }
        else if (i >= Nspline_inner and j >= Nspline_inner)
        {
            full_solution[i * Nspline_full + j] = 0;
        }
        else if (i >= Nspline_inner)
        {
            // all channels r1 -> inf
            for (int n = 0; n < Nchan1; n++)
                full_solution[i * Nspline_full + j] += sol[Nspline_inner * Nspline_inner + n * Nspline_outer + i - Nspline_inner] * Xp2[n][j];
        }
        else /* if (j >= Nspline_inner) */
        {
            // all channels r2 -> inf
            for (int n = 0; n < Nchan2; n++)
                full_solution[i * Nspline_full + j] += Xp1[n][i] * sol[Nspline_inner * Nspline_inner + (Nchan1 + n) * Nspline_outer + j - Nspline_inner];
        }
    }
    
    // evaluate the solution
    cArray evPsi = Bspline::zip
    (
        bspline_full,
        bspline_full,
        full_solution,
        grid_x,
        grid_y
    );
    
    // evaluate the solution differentiated with respect to the first coordinate
    cArray evDxPsi = Bspline::zip
    (
        bspline_full,
        bspline_full,
        full_solution,
        grid_x,
        grid_y,
        &Bspline::dspline,
        &Bspline::bspline
    );
    
    // evaluate the solution differentiated with respect to the second coordinate
    cArray evDyPsi = Bspline::zip
    (
        bspline_full,
        bspline_full,
        full_solution,
        grid_x,
        grid_y,
        &Bspline::bspline,
        &Bspline::dspline
    );
    
    // write full solution to file
    VTKRectGridFile vtk;
    vtk.setGridX(grid_x);
    vtk.setGridY(grid_y);
    vtk.setGridZ(rArray{0.});
    vtk.appendScalarAttribute("RePsi", realpart(evPsi));
    vtk.appendScalarAttribute("ImPsi", imagpart(evPsi));
    vtk.appendScalarAttribute("probDensity", sqrabs(evPsi));
    vtk.appendVector3DAttribute("probFlux", imagpart(evPsi.conj() * evDxPsi), imagpart(evPsi.conj() * evDyPsi), rArray(evPsi.size(), 0.0_r));
    vtk.write(cmd.zipdata.file + "-full.vtk");
    
    // expand the channel functions
    for (int n = 0; n < Nchan1; n++)
    {
        cArray full_channel_function (Nspline_full, 0.0_z);
        
        for (std::size_t i = 0; i < Nspline_full; i++)
        {
            if (i < Nspline_inner)
            {
                for (std::size_t j = 0; j < Nspline_inner; j++)
                    full_channel_function[i] += sol[i * Nspline_inner + j] * Sp2[n][j];
            }
            else
            {
                full_channel_function[i] = sol[Nspline_inner * Nspline_inner + n * Nspline_outer + i - Nspline_inner];
            }
        }
        
        // write to file
        std::ofstream out (format("%s-channel-X-%d.vtk", cmd.zipdata.file.c_str(), n).c_str());
        writeVTK_points
        (
            out,
            bspline_full.zip(full_channel_function, grid_x),
            grid_x, rArray({0.}), rArray({0.})
        );
    }
    for (int n = 0; n < Nchan2; n++)
    {
        cArray full_channel_function (Nspline_full);
        
        for (std::size_t i = 0; i < Nspline_full; i++)
        {
            if (i < Nspline_inner)
            {
                for (std::size_t j = 0; j < Nspline_inner; j++)
                    full_channel_function[i] += Sp1[n][j] * sol[j * Nspline_inner + i];
            }
            else
            {
                full_channel_function[i] = sol[Nspline_inner * Nspline_inner + (Nchan1 + n) * Nspline_outer + i - Nspline_inner];
            }
        }
        
        // write to file
        std::ofstream out (format("%s-channel-%d-X.vtk", cmd.zipdata.file.c_str(), n).c_str());
        writeVTK_points
        (
            out,
            bspline_full.zip(full_channel_function, grid_x),
            grid_x, rArray({0.}), rArray({0.})
        );
    }
}

void write_grid (Bspline const & bspline, std::string const & basename)
{
    // get atomic grid
    rArray knots = bspline.rknots();
    knots.pop_back();
    knots.append(bspline.cknots2());
    
    // output file
    std::ofstream out (basename + ".vtk");
    
    // write knots (write header only for the first time)
    writeVTK_points(out, cArray(), knots, knots, rArray({0.}));
}
