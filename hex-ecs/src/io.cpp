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
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <string>
#include <tuple>

#include "hex-arrays.h"
#include "hex-cmdline.h"
#include "hex-matrix.h"

#include "io.h"
#include "preconditioners.h"

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
    "# The sections of the example are graphically demonstrated below.\n"
    "# Letters used here are the same as in numbering of the user input given\n"
    "# after this comment.\n"
    "#\n"
    "#         0      60   100   150 a.u.\n"
    "#  solver |aaaaaa|bbbb|ccccc|\n"
    "#                            160  200   250 a.u.\n"
    "#  1. propagator |bbbb|dddddd|bbbb|ccccc|\n"
    "#                                        260  300   350 a.u.\n"
    "#  2. propagator             |bbbb|dddddd|bbbb|ccccc|\n"
    "#                                                    360  400   450 a.u.\n"
    "#  3. propagator                         |bbbb|dddddd|bbbb|ccccc|\n"
    "#\n"
    "#  etc.\n"
    "#\n"
    "#   1) Solver grid (0 a.u. - 100 a.u. / 150 a.u.)\n"
    "#      - multiple knots at origin (4×)\n"
    "#      - geometric grid to 10 a.u.\n"
    "#      - uniform grid to distance (60 a.u. + 40 a.u. =) 100 a.u.\n"
    "#      - complex absorbtion layer to distance 150 a.u.\n"
    "#   2) First propagation grid (60 a.u. - 200 a.u. / 250 a.u.)\n"
    "#      - uniform grid of length (40 a.u. + 60 a.u. + 40 a.u. =) 140 a.u. to distance 200 a.u.\n"
    "#      - complex absorbtion layer of length 50 a.u. to distance 250 a.u.\n"
    "#   3+) Further propagation grids, always extending the distance by 100 a.u.\n"
    "# It is a good idea to visualise the grid before the calculation. You can use the command\n"
    "# \"hex-ecs --write-grid\" and view the result in ParaView, or any other program supporting VTK.\n"
    "#\n"
    "# a) Real knots of the basis that is common to atomic and projectile electron.\n"
    "  L  0.0  0.0   4\n"
    "  G  0.1 10.0  0.1  1.1\n"
    "  L   11   60  50\n"
    " -1\n"
    "# b) Real knots of the panel overlap, if any.\n"
    "  L    0   40  41\n"
    " -1\n"
    "# c) Complex region knots.\n"
    "  G    0   50   1  1.02\n"
    " -1\n"
    "# d) Knots of other panels (propagator projectile basis).\n"
    "  L    0   60  61\n"
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
    "# --------------- Other conditions ----------------\n"
    "\n"
    "# Angular momenta.\n"
    "# L  S  Pi limit\n"
    "  0  *  0  4\n"
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
                // print usage information
                std::cout << "\n"
                    "Available switches (short forms in parentheses):                                                                                                          \n"
                    "                                                                                                                                                          \n"
                    "\t--example                  (-e)  Create sample input file.                                                                                              \n"
                    "\t--help                     (-h)  Display this help.                                                                                                     \n"
                    "\t--input <filename>         (-i)  Use custom input file (other than \"ecs.inp\").                                                                        \n"
                    "\t--zip <parameters>         (-z)  Solution file to zip (i.e. evaluate in B-spline basis and produce VTK datafile).                                       \n"
                    "\t                                 The '<parameters>' stands for '<filename> <Xmin> <Ymin> <Xmax> <Ymax> <Xn> <Yn>'.                                      \n"
                    "\t--write-grid               (-g)  Write grid layout to a VTK file.                                                                                       \n"
                    "\t--map-solution <filename>        Convert solution between B-spline bases, requires target basis (--map-solution-target).                                \n"
                    "\t--map-solution-target <filename> Target B-spline basis in another input file (like ecs.inp).                                                            \n"
#ifdef WITH_MPI
                    "\t--mpi                      (-m)  Use MPI (assuming that the program has been launched by mpiexec).                                                      \n"
#endif
                    "\t--stg-integ                (-a)  Only calculate needed radial integrals.                                                                                \n"
                    "\t--stg-integ-solve          (-b)  Only calculate integrals and the solution.                                                                             \n"
                    "\t--stg-extract              (-c)  Only extract amplitudes (assumes that the solution files exist).                                                       \n"
                    "\t--preconditioner <name>    (-p)  Preconditioner to use (default: ILU).                                                                                  \n"
                    "\t--list-preconditioners     (-P)  List available preconditioners with short description of each.                                                         \n"
                    "\t--ssor <number>                  Apply SSOR coupling.                                                                                                   \n"
                    "\t--tolerance <number>       (-T)  Set tolerance for the conjugate gradients solver (default: 1e-8).                                                      \n"
                    "\t--prec-tolerance <number>  (-t)  Set tolerance for the conjugate gradients preconditioner (default: 1e-8).                                              \n"
                    "\t--drop-tolerance <number>  (-d)  Set drop tolerance for the ILU preconditioner (default: 1e-15).                                                        \n"
                    "\t--lu <name>                (-F)  Factorization library (one of 'umfpack', 'superlu' and 'superlu_dist'). Default is 'umfpack' (if available).           \n"
                    "\t--no-lu-update                   Do not recalculate LU factorization for different energies, use the first factorization for all of them.               \n"
                    "\t--parallel-factorization         Factorize multiple blocks simultaneously.                                                                              \n"
                    "\t--no-parallel-extraction         Disallow parallel extraction of T-matrices (e.g. when the whole solution does not fit into the memory).                \n"
                    "\t--extract-rho-begin              Where to start averaging / extrapolating the T-matrix.                                                                 \n"
                    "\t--extract-rho[-end]              Radius for T-matrix extraction.                                                                                        \n"
                    "\t--extract-samples                Number of evaluations of the T-matrix between --extract-rho-begin and --extract-rho.                                   \n"
                    "\t--extract-extrapolate            Radially extrapolate the extracted T-matrices instead of simple averaging.                                             \n"
                    "\t--groupsize <number>       (-G)  How many processes factorize single LU (only used for 'superlu_dist').                                                 \n"
                    "\t--own-radial-cache         (-w)  Keep two-electron radial integrals not referenced by preconditioner only on disk (slows down only the initialization). \n"
                    "\t--no-radial-cache          (-r)  Keep all two-electron radial integrals only on disk (slows down also the solution process).                            \n"
                    "\t--out-of-core              (-o)  Use hard disk drive to store most of intermediate data and thus to save RAM (considerably slower).                     \n"
                    "\t--out-of-core-continue     (-O)  Start solution from the existing OOC files.                                                                            \n"
                    "\t--whole-matrix             (-W)  In the above three cases: Load whole matrix from scratch file when calculating dot product (speeds them up a little).  \n"
                    "\t--shared-scratch           (-s)  Let every MPI process calculate only a subset of shared radial integrals (assume shared output directory).             \n"
                    "\t--lightweight-radial-cache (-l)  Do not precalculate two-electron integrals and only apply them on the fly (slower, but saves RAM).                     \n"
                    "\t--lightweight-full         (-L)  Avoid precalculating all large matrices and only apply them on the fly (only available for KPA preconditioner).        \n"
                    "\t--kpa-simple-rad           (-R)  Use simplified radial integral matrix for nested KPA iterations (experimental).                                        \n"
                    "\t--kpa-max-iter                   Maximal KPA iterations for hybrid preconditioner.                                                                      \n"
                    "\t--ilu-max-blocks                 Maximal number of ILU preconditioned blocks (per MPI node) for hybrid preconditioner.                                  \n"
#ifdef WITH_MUMPS
                    "\t--coupling-limit                 Maximal multipole to be considered by the coupled preconditioner.                                                      \n"
                    "\t--mumps-in-core                  Try to keep all factorization data in memory when using MUMPS.                                                         \n"
                    "\t--mumps-verbose                  Verbosity level of the MUMPS library.                                                                                  \n"
#endif
#ifndef DISABLE_PARALLEL_PRECONDITION
                    "\t--parallel-precondition          Apply multiple block preconditioners in parallel.                                                                      \n"
#endif
                    "\t--panels <number>                Propagate solution through given number of panels.                                                                     \n"
                    "\t--carry-initial-guess            Whether to use previous-energy solution as an initial guess for the new energy.                                        \n"
                    "\t--refine-solution                Load existing solutions and check that they are within tolerance, update if needed.                                    \n"
#ifdef WITH_OPENCL
                    "\t--cl-list                        List all OpenCL platforms and devices available.                                                                       \n"
                    "\t--cl-platform <index>            Use given OpenCL platform for GPU preconditioner (default is 0, i.e. the first platform found).                        \n"
                    "\t--cl-device <index>              Use given OpenCL device for GPU preconditioner (default is 0, i.e. the first device found).                            \n"
                    "\t--cl-use-host-memory             Keep large data in RAM instead of copying everything to the compute device. This will slow down the solution.          \n"
                    "\t--cl-multiply                    Do the sparse matrix multiplication on the OpenCL device (memory intensive!).                                          \n"
                    "\t--cl-host-multiply               Keep vectors in host memory when doing matrix multiplication on OpenCL device.                                         \n"
#endif
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
        "map-solution", "", -1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // solution file to map
                map_solution = optargs;
                return true;
            },
        "map-solution-target", "", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // solution file to map
                map_solution_target = optargs[0];
                return true;
            },
        "lu", "F", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // choose factorizer
                if (optargs[0] == "umfpack")
                    factorizer = LUFT_UMFPACK;
                else if (optargs[0] == "superlu")
                    factorizer = LUFT_SUPERLU;
                else if (optargs[0] == "superlu_dist")
                    factorizer = LUFT_SUPERLU_DIST;
                else
                    HexException("Unknown LU-factorizer '%s'.", optargs[0].c_str());
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
                // preconditioner
                if ((preconditioner = Preconditioners::findByName(optargs[0])) == -1)
                    HexException("Unknown preconditioner \"%s\".", optargs[0].c_str());
                return true;
            },
        "list-preconditioners", "P", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // preconditioners description
                std::cout << "\nPreconditioners description (first one is default):\n\n";
                for (unsigned i = 0; i < Preconditioners::size(); i++)
                {
                    std::cout << Preconditioners::name(i) << "\n";
                    std::cout << "\t" << Preconditioners::description(i) << "\n";
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
        "kpa-simple-rad", "R", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // use simplified radial matrix for KPA nested iterations
                kpa_simple_rad = true;
                return true;
            },
        "kpa-max-iter", "", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // number of KPA preconditioner iterations to trigger ILU preconditioning
                kpa_max_iter = std::atoi(optargs[0].c_str());
                return true;
            },
        "ilu-max-blocks", "", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // maximal number of ILU preconditioned blocks for hybrid preconditioner
                ilu_max_blocks = std::atoi(optargs[0].c_str());
                return true;
            },
        "no-lu-update", "", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // do not recalculate LU
                noluupdate = true;
                return true;
            },
#ifdef WITH_MUMPS
        "coupling-limit", "", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // maximal multipole to be considered by the coupled preconditioner
                coupling_limit = std::atoi(optargs[0].c_str());
                return true;
            },
        "mumps-in-core", "", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // MUMPS out of core
                mumps_outofcore = false;
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
        "reuse-dia-blocks", "", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // load dia blocks from scratch files
                reuse_dia_blocks = true;
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
        "panels", "", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // propagate solution along the projectile axis
                panels = std::atoi(optargs[0].c_str());
                return true;
            },
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
            double begin = ReadNext<double>(inf).val;
            double end = ReadNext<double>(inf).val;
            int samples = ReadNext<int>(inf).val;
            arr.append(linspace(begin, end, samples));
        }
        else if (type[0] == 'G')
        {
            double begin = ReadNext<double>(inf).val;
            double end = ReadNext<double>(inf).val;
            double d = ReadNext<double>(inf).val;
            double quotient = ReadNext<double>(inf).val;
            int samples = std::ceil(1 + std::log(1 + (end - begin) * (quotient - 1) / d) / std::log(quotient));
            arr.append(geomspace(begin, end, samples, quotient));
        }
        else if (type[0] == 'E')
        {
            // explicit list samples
            double X;
            while ((X = ReadNext<double>(inf).val) != -1)
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
    ecstheta = ReadNext<double>(inf).val;
    
    std::cout << std::endl;
    std::cout << "B-spline basis" << std::endl;
    std::cout << "\torder = " << order << std::endl;
    std::cout << "\tecs angle = " << ecstheta << std::endl;
    
    //
    // load real solver knot data
    //
    
    ReadArrays(inf, rknots);
    
    // print info
    std::cout << std::endl;
    std::cout << "Real knots (" << rknots.size() << ")" << std::endl;
    for (std::string line : rknots.lines(100))
        std::cout << '\t' << line << std::endl;
    
    // check order of knots
    for (unsigned i = 1; i < rknots.size(); i++)
        if (rknots[i] < rknots[i-1])
            HexException("The real knot sequence is not monotonous.");
    
    //
    // load basis overlap knot data
    //
    
    ReadArrays(inf, overlap_knots);
    
    // print info
    std::cout << std::endl;
    std::cout << "Overlap knots (" << overlap_knots.size() << ")" << std::endl;
    for (std::string line : overlap_knots.lines(100))
        std::cout << '\t' << line << std::endl;
    
    // check that the first knot is zero
    if (overlap_knots.size() > 0 and overlap_knots[0] != 0.)
        HexException("The first knot in overlap region must be zero.");
    
    // check order of knots
    for (unsigned i = 1; i < overlap_knots.size(); i++)
        if (overlap_knots[i] < overlap_knots[i-1])
            HexException("The overlap knot sequence is not monotonous.");
    
    //
    // load complex knot data
    //
    
    ReadArrays(inf, cknots);
    
    // print info
    std::cout << std::endl;
    std::cout << "Complex knots (before scaling; " << cknots.size() << ")" << std::endl;
    for (std::string line : cknots.lines(100))
        std::cout << '\t' << line << std::endl;
    
    // check that the first knot is zero
    if (cknots.size() > 0 and cknots[0] != 0.)
        HexException("The first knot in complex region must be zero.");
    
    // check order of knots
    for (unsigned i = 1; i < cknots.size(); i++)
        if (cknots[i] < cknots[i-1])
            HexException("The complex knot sequence is not monotonous.");
    
    //
    // load real propagator knot data
    //
    
    ReadArrays(inf, rknots_next);
    
    // print info
    std::cout << std::endl;
    std::cout << "Propagator real knots (" << rknots_next.size() << ")" << std::endl;
    for (std::string line : rknots_next.lines(100))
        std::cout << '\t' << line << std::endl;
    
    // check order of knots
    for (unsigned i = 1; i < rknots_next.size(); i++)
        if (rknots_next[i] < rknots_next[i-1])
            HexException("The propagator real knot sequence is not monotonous.");
    
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
    
    // print info
    std::cout << "\t[n l m]: ";
    for (auto state : outstates)
    {
        std::cout << "["
                  << std::get<0>(state) << " "
                  << std::get<1>(state)
                  << " *] ";
    }
    std::cout << std::endl;
    
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
    
    std::cout << "\tL = " << L << std::endl;
    std::cout << "\tS = " << Spin << std::endl;
    std::cout << "\tPi = " << Pi << std::endl;
    std::cout << "\tnL = " << levels << std::endl;
    
    //
    // load initial energies
    //
    
    ReadArrays(inf, Etot);
    
    // print info
    std::cout << std::endl << "Total projectile + atom energies [Ry]" << std::endl;
    std::cout << "\tcount: "     << Etot.size()  << std::endl;
    std::cout << "\tfull list: " << Etot         << std::endl;
    
    //
    // load some other optional data
    //
    
    std::cout << std::endl;
    std::cout << "Other parameters" << std::endl;
    B = ReadNext<double>(inf).val;
    
    // print info
    std::cout << "\tmagnetic field: " << B << " a.u." << std::endl;
    std::cout << std::endl;
}

void zip_solution (CommandLine & cmd, std::vector<Bspline> const & bspline, std::vector<std::pair<int,int>> const & ll)
{
    cArray sol;     // stored solution expansion
    cArray ev;      // evaluated solution
    rArray grid_x;  // real evaluation grid (atomic electron)
    rArray grid_y;  // real evaluation grid (projectile electron)
    
    std::cout << "Zipping B-spline expansion of the solution: \"" << cmd.zipdata.file << "\"" << std::endl;
    
    // load the requested file
    if (not sol.hdfload(cmd.zipdata.file.c_str()))
        HexException("Cannot load file %s.", cmd.zipdata.file.c_str());
    
    // determine which B-spline basis to use
    unsigned i;
    for (i = 0; i < bspline.size(); i++)
    {
        if (sol.size() == (std::size_t)bspline[0].Nspline() * (std::size_t)bspline[i].Nspline())
            break;
    }
    
    // was some basis appropriate?
    if (i == bspline.size())
        HexException("The solution file of size %ld is not compatible with defined B-spline basis. Did you specify the same number of panels?", sol.size());
    
    // evaluation grid
    double Xmin = cmd.zipdata.Xmin < 0 ? 0 : cmd.zipdata.Xmin;
    double Ymin = cmd.zipdata.Ymin < 0 ? 0 : cmd.zipdata.Ymin;
    double Xmax = cmd.zipdata.Xmax < 0 ? bspline[0].Rmax() : cmd.zipdata.Xmax;
    double Ymax = cmd.zipdata.Ymax < 0 ? bspline[i].Rmax() : cmd.zipdata.Ymax;
    grid_x = linspace(Xmin, Xmax, cmd.zipdata.nX);
    grid_y = linspace(Ymin, Ymax, cmd.zipdata.nY);
    
    // write to file
    std::ofstream out ((cmd.zipdata.file + ".vtk").c_str());
    writeVTK_points
    (
        out,                                // output file stream
        Bspline::zip
        (
            bspline[0], bspline[i],         // B-spline bases (for x and y)
            sol,                            // function expansion in those two bases
            grid_x, grid_y                  // evaluation grids
        ),
        grid_x, grid_y, rArray({0.})        // x,y,z
    );
}

void write_grid (std::vector<Bspline> const & bspline, std::string const & basename)
{
    // get atomic grid
    rArray knots0 = bspline[0].rknots();
    knots0.pop_back();
    knots0.append(bspline[0].cknots());
    
    // for all grids
    for (unsigned i = 0; i < bspline.size(); i++)
    {
        // output file
        std::ofstream out (format("%s-%d.vtk", basename.c_str(), i).c_str());
        
        // get knots
        rArray knots = bspline[i].rknots();
        knots.pop_back();
        knots.append(bspline[i].cknots());
        
        // write knots (write header only for the first time)
        writeVTK_points(out, cArray(), knots0, knots, rArray({0.}));
    }
}
