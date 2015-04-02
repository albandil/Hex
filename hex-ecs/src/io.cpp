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

#include "arrays.h"
#include "cmdline.h"
#include "io.h"
#include "matrix.h"
#include "preconditioners.h"

const std::string sample_input =
    "# B-spline parameters \n"
    "# order      θ\n"
    "      4   0.63\n"
    "\n"
    "# real knot sequences\n"
    "# L[inear] <Start> <End> <Sample>\n"
    "# G[eometric] <Start> <End> <FirstInterval> <Quotient>\n"
    "  L  0.0  0.0   4\n"
    "  L  0.1  2.0  20\n"
    "  L    3   60  58\n"
    " -1\n"
    "\n"
    "# complex knot sequences\n"
    "  G   60  100   1  1.02\n"
    " -1\n"
    "\n"
    "# initial atomic states (ni, li, mi)\n"
    "  1 -1\n"
    "  *\n"
    "  *\n"
    "\n"
    "# final atomic states (nf, lf)\n"
    "  1  -1\n"
    "  *\n"
    "\n"
    "# angular momenta\n"
    "# L  Pi limit\n"
    "  0  0  4\n"
    "\n"
    "# atom + projectile total energies in Rydbergs\n"
    "  L  -0.35  -0.05  3\n"
    " -1\n"
    "\n"
    "# magnetic field\n"
    " 0\n";


void CommandLine::parse (int argc, char* argv[])
{
    ParseCommandLine
    (
        argc, argv,
        
        "example", "e", 0, [&](std::string optarg) -> bool
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
        "input", "i", 1, [&](std::string optarg) -> bool
            {
                // set custom input file
                inputfile.open(optarg);
                if (not inputfile.good())
                    HexException("Error: Input file \"%s\" not found.\n", optarg.c_str());
                return true;
            },
        "help", "h", 0, [&](std::string optarg) -> bool
            {
                // print usage information
                std::cout << "\n"
                    "Available switches (short forms in parentheses):                                                                                                          \n"
                    "                                                                                                                                                          \n"
                    "\t--example                  (-e)  Create sample input file.                                                                                              \n"
                    "\t--help                     (-h)  Display this help.                                                                                                     \n"
                    "\t--input <filename>         (-i)  Use custom input file (other than \"ecs.inp\").                                                                        \n"
                    "\t--zipfile <filename>       (-z)  Solution file to zip (i.e. evaluate in B-spline basis and produce VTK datafile).                                       \n"
                    "\t--zipcount <number>        (-n)  Zip sample count (how many points along r1 and r2).                                                                    \n"
                    "\t--zipmax <number>          (-R)  Maximal radius to use for solution zipping.                                                                            \n"
#ifndef NO_MPI
                    "\t--mpi                      (-m)  Use MPI (assuming that the program has been launched by mpiexec).                                                      \n"
#endif
                    "\t--stg-integ                (-a)  Only calculate needed radial integrals.                                                                                \n"
                    "\t--stg-integ-solve          (-b)  Only calculate integrals and the solution.                                                                             \n"
                    "\t--stg-extract              (-c)  Only extract amplitudes (assumes that the solution files exist).                                                       \n"
                    "\t--preconditioner <name>    (-p)  Preconditioner to use (default: ILU).                                                                                  \n"
                    "\t--list-preconditioners     (-P)  List available preconditioners with short description of each.                                                         \n"
                    "\t--tolerance <number>       (-T)  Set tolerance for the conjugate gradients solver (default: 1e-8).                                                      \n"
                    "\t--prec-tolerance <number>  (-t)  Set tolerance for the conjugate gradients preconditioner (default: 1e-8).                                              \n"
                    "\t--drop-tolerance <number>  (-d)  Set drop tolerance for the ILU preconditioner (default: 1e-15).                                                        \n"
                    "\t--own-radial-cache         (-w)  Keep two-electron radial integrals not referenced by preconditioner only on disk (slows down only the initialization). \n"
                    "\t--no-radial-cache          (-r)  Keep all two-electron radial integrals only on disk (slows down also the solution process).                            \n"
                    "\t--out-of-core              (-o)  Use hard disk drive to store most of intermediate data and thus to save RAM (considerably slower).                     \n"
                    "\t--out-of-core-continue     (-O)  Start solution from the existing OOC files.                                                                            \n"
                    "\t--whole-matrix             (-W)  In the above three cases: Load whole matrix from scratch file when calculating dot product (speeds them up a little).  \n"
                    "\t--shared-scratch           (-s)  Let every MPI process calculate only a subset of shared radial integrals (assume shared output directory).             \n"
                    "\t--lightweight-radial-cache (-l)  Do not precalculate two-electron integrals and only apply them on the fly (slower, but saves RAM).                     \n"
#ifndef NO_LAPACK
                    "\t--lightweight-full         (-L)  Avoid precalculating all large matrices and only apply them on the fly (only available for KPA preconditioner).        \n"
                    "\t--kpa-simple-rad           (-R)  Use simplified radial integral matrix for nested KPA iterations (experimental).                                        \n"
#endif
                    "\t--parallel-dot                   OpenMP-parallelize SpMV operations.                                                                                    \n"
                    "\t--parallel-block                 Enable concurrent handling of matrix blocks by OpenMP (e.g. in preconditioning and multiplication).                    \n"
#if (!defined(NO_OPENCL) && !defined(NO_LAPACK))
                    "\t--cl-list                        List all OpenCL platforms and devices available.                                                                       \n"
                    "\t--cl-platform <index>            Use given OpenCL platform for GPU preconditioner (default is 0, i.e. the first platform found).                        \n"
                    "\t--cl-device <index>              Use given OpenCL device for GPU preconditioner (default is 0, i.e. the first device found).                            \n"
#endif
                    "                                                                                                                                                          \n"
                    "There are also some environment variables that control the execution.                                                                                     \n"
                    "                                                                                                                                                          \n"
                    "\tOMP_NUM_THREADS     Number of OpenMP threads to use.                                                                                                    \n"
                    "\t                    If empty, the physical number of hardware threads set by system is used.                                                            \n"
                    "                                                                                                                                                          \n"
                    "\tHEX_RHO             Scattering amplitude extraction distance (must be on the real part of the grid).                                                    \n"
                    "\t                    If empty, the end of the real grid is used.                                                                                         \n"
                    "                                                                                                                                                          \n"
                    "\tHEX_SAMPLES         How many times to evaluate the amplitude in the vicinity of the evaluation point. The evaluations will be uniformly                 \n"
                    "\t                    spread over one wave-length of the scattered electron and averaged to suppress numerical errors.                                    \n"
                    "\t                    If empty, 10 evaluations are used.                                                                                                  \n"
                    "                                                                                                                                                          \n"
                ;
                std::exit(EXIT_SUCCESS);
            },
        "zipfile", "z", 1, [&](std::string optarg) -> bool
            {
                // zip B-spline expansion file
                zipfile = optarg;
                return true;
            },
        "zipcount", "n", 1, [&](std::string optarg) -> bool
            {
                // zip samples
                zipcount = std::atol(optarg.c_str());
                return true;
            },
        "zipmax", "R", 1, [&](std::string optarg) -> bool
            {
                // zip bounding box
                zipmax = std::atof(optarg.c_str());
                return true;
            },
#ifndef NO_MPI
        "mpi", "m", 0, [&](std::string optarg) -> bool
            {
                // use MPI
                parallel = true;
                return true;
            },
#endif
        "stg-integ", "a", 0, [&](std::string optarg) -> bool
            {
                // run only the first part (computation of radial integrals)
                itinerary = StgRadial;
                return true;
            },
        "stg-integ-solve", "b", 0, [&](std::string optarg) -> bool
            {
                // run only the first part (computation of radial integrals)
                itinerary = StgRadial | StgSolve;
                return true;
            },
        "stg-extract", "c", 0, [&](std::string optarg) -> bool
            {
                // run only the third part (extraction of amplitudes)
                itinerary = StgExtract;
                return true;
            },
        "own-radial-cache", "w", 0, [&](std::string optarg) -> bool
            {
                // do not cache un-'owned' two-electron radial integrals in memory
                cache_own_radint = true;
                cache_all_radint = false;
                return true;
            },
        "no-radial-cache", "r", 0, [&](std::string optarg) -> bool
            {
                // do not cache any two-electron radial integrals in memory at all
                cache_own_radint = false;
                cache_all_radint = false;
                return true;
            },
        "out-of-core", "o", 0, [&](std::string optarg) -> bool
            {
                // use full out-of-core functionality: store also diagonal blocks (and factorizations) on disk
                cache_all_radint = false;
                cache_own_radint = false;
                outofcore = true;
                return true;
            },
        "out-of-core-continue", "O", 0, [&](std::string optarg) -> bool
            {
                // continue OOC calculation
                cache_all_radint = false;
                cache_own_radint = false;
                reuse_dia_blocks = true;
                outofcore = true;
                cont = true;
                return true;
            },
        "whole-matrix", "W", 0, [&](std::string optarg) -> bool
            {
                // whether to load whole matrix from scratch disk when calculating dot product (etc.)
                wholematrix = true;
                return true;
            },
        "drop-tolerance", "d", 1, [&](std::string optarg) -> bool
            {
                // drop tolerance for iLU-factorization
                droptol = std::atof(optarg.c_str());
                return true;
            },
        "tolerance", "T", 1, [&](std::string optarg) -> bool
            {
                // iteration tolerance for terminating iteration solution
                itertol = std::atof(optarg.c_str());
                return true;
            },
        "prec-tolerance", "t", 1, [&](std::string optarg) -> bool
            {
                // iteration tolerance for terminating iteration solution
                prec_itertol = std::atof(optarg.c_str());
                return true;
            },
        "preconditioner", "p", 1, [&](std::string optarg) -> bool
            {
                // preconditioner
                if ((preconditioner = Preconditioners::findByName(optarg)) == -1)
                    HexException("Unknown preconditioner \"%s\".", optarg.c_str());
                return true;
            },
        "list-preconditioners", "P", 0, [&](std::string optarg) -> bool
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
        "parallel-dot", "", 0, [&](std::string optarg) -> bool
            {
                // parallelize SpMV
                parallel_dot = true;
                return true;
            },
        "parallel-block", "", 0, [&](std::string optarg) -> bool
            {
                // parallelize preconditioning
                parallel_block = true;
                return true;
            },
        "gpu-slater", "", 0, [&](std::string optarg) -> bool
            {
                // compute diagonal two-electron integrals using OpenCL
                gpu_slater = true;
                return true;
            },
        "lightweight-radial-cache", "l", 0, [&](std::string optarg) -> bool
            {
                // do not precompute two-electron radial integral matrices but only apply them on the fly
                lightweight_radial_cache = true;
                return true;
            },
#ifndef NO_LAPACK
        "lightweight-full", "L", 0, [&](std::string optarg) -> bool
            {
                // do not precompute large matrices but only apply them on the fly
                lightweight_full = lightweight_radial_cache = true;
                return true;
            },
        "kpa-simple-rad", "R", 0, [&](std::string optarg) -> bool
            {
                // use simplified radial matrix for KPA nested iterations
                kpa_simple_rad = true;
                return true;
            },
#endif
        "shared-scratch", "s", 0, [&](std::string optarg) -> bool
            {
                // precompute only the owned subset of radial integrals
                shared_scratch = true;
                return true;
            },
        "reuse-dia-blocks", "", 0, [&](std::string optarg) -> bool
            {
                // load dia blocks from scratch files
                reuse_dia_blocks = true;
                return true;
            },
#if (!defined(NO_OPENCL) && !defined(NO_LAPACK))
        "cl-list", "", 0,  [&](std::string optarg) -> bool
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
                    for (cl_uint j = 0; j < nplatforms; j++)
                    {
                        clGetDeviceInfo(devices[j], CL_DEVICE_NAME, sizeof(device_name), device_name, nullptr);
                        clGetDeviceInfo(devices[j], CL_DEVICE_VENDOR, sizeof(device_vendor), device_vendor, nullptr);
                        std::cout << "\t\t- Device " << j << ": " << device_name << " (" << device_vendor << ")" << std::endl;
                    }
                }
                std::cout << std::endl;
                std::exit(0);
            },
        "cl-platform", "", 1, [&](std::string optarg) -> bool
            {
                // use given OpenCL platform
                ocl_platform = std::atol(optarg.c_str());
                return true;
            },
        "cl-device", "", 1, [&](std::string optarg) -> bool
            {
                // use given OpenCL device
                ocl_device = std::atol(optarg.c_str());
                return true;
            },
#endif
        
        [&] (std::string optname, std::string optarg) -> bool
        {
            HexException("Unknown switch \"%s\".", optname.c_str());
        }
    );
}

void InputFile::read (std::ifstream & inf)
{
    ReadItem<int> idata;
    
    // load B-spline parameters
    order = ReadNext<int>(inf).val;
    ecstheta = ReadNext<double>(inf).val;
    
    std::cout << std::endl;
    std::cout << "B-spline environment" << std::endl;
    std::cout << "\torder = " << order << std::endl;
    std::cout << "\tecs angle = " << ecstheta << std::endl;
    
    //
    // load real knot data
    //
    
    std::string type;
    while ((type = ReadNext<std::string>(inf).val) != std::string("-1"))
    {
        if (type[0] == 'L')
        {
            double begin = ReadNext<double>(inf).val;
            double end = ReadNext<double>(inf).val;
            int samples = ReadNext<int>(inf).val;
            rknots.append(linspace(begin, end, samples));
        }
        else if (type[0] == 'G')
        {
            double begin = ReadNext<double>(inf).val;
            double end = ReadNext<double>(inf).val;
            double d = ReadNext<double>(inf).val;
            double quotient = ReadNext<double>(inf).val;
            int samples = std::ceil(1 + std::log(1 + (end - begin) * (quotient - 1) / d) / std::log(quotient));
            rknots.append(geomspace(begin, end, samples, quotient));
        }
        else
        {
            HexException("Unknown sequence type \"%s\".", type.c_str());
        }
    }
    
    // print info
    std::cout << std::endl;
    std::cout << "Real knots" << std::endl;
    for (std::string line : rknots.lines(100))
        std::cout << '\t' << line << std::endl;
    
    // check order of knots
    for (unsigned i = 1; i < rknots.size(); i++)
        if (rknots[i] < rknots[i-1])
            HexException("The real knot sequence is not monotonous.");
    
    //
    // load complex knot data
    //
    
    while ((type = ReadNext<std::string>(inf).val) != std::string("-1"))
    {
        if (type[0] == 'L')
        {
            double begin = ReadNext<double>(inf).val;
            double end = ReadNext<double>(inf).val;
            int samples = ReadNext<int>(inf).val;
            cknots.append(linspace(begin, end, samples));
        }
        else if (type[0] == 'G')
        {
            double begin = ReadNext<double>(inf).val;
            double end = ReadNext<double>(inf).val;
            double d = ReadNext<double>(inf).val;
            double quotient = ReadNext<double>(inf).val;
            int samples = std::ceil(1 + std::log(1 + (end - begin) * (quotient - 1) / d) / std::log(quotient));
            cknots.append(geomspace(begin, end, samples, quotient));
        }
        else
        {
            HexException("Complex knot sequence: Unknown sequence type \"%s\".", type.c_str());
        }
    }
    
    // print info
    std::cout << std::endl;
    std::cout << "Complex knots (before scaling)" << std::endl;
    for (std::string line : cknots.lines(100))
        std::cout << '\t' << line << std::endl;
    
    // check continuity
    if (cknots.front() < rknots.back())
        HexException("Complex part of the grid overlaps with the real part (%g < %g).", cknots.front(), rknots.back());
    
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
    
    L = ReadNext<int>(inf).val;
    Pi = ReadNext<int>(inf).val % 2;
    levels = ReadNext<int>(inf).val;
    
    std::cout << "\tL = " << L << std::endl;
    std::cout << "\tPi = " << Pi << std::endl;
    std::cout << "\tnL = " << levels << std::endl;
    
    //
    // load initial energies
    //
    
    std::cout << std::endl << "Total projectile + atom energies [Ry]" << std::endl;
    while ((type = ReadNext<std::string>(inf).val) != std::string("-1"))
    {
        if (type[0] == 'L')
        {
            double begin = ReadNext<double>(inf).val;
            double end = ReadNext<double>(inf).val;
            int samples = ReadNext<int>(inf).val;
            Etot.append(linspace(begin, end, samples));
        }
        else if (type[0] == 'G')
        {
            double begin = ReadNext<double>(inf).val;
            double end = ReadNext<double>(inf).val;
            double d = ReadNext<double>(inf).val;
            double quotient = ReadNext<double>(inf).val;
            int samples = std::ceil(1 + std::log(1 + (end - begin) * (quotient - 1) / d) / std::log(quotient));
            Etot.append(geomspace(begin, end, samples, quotient));
        }
        else
        {
            HexException("Energy list: Unknown sequence type \"%s\".", type.c_str());
        }
    }
    
    // print info
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

void zip_solution (CommandLine & cmd, Bspline const & bspline, std::vector<std::pair<int,int>> const & ll)
{
    // use whole grid if not restricted
    if (cmd.zipmax < 0)
        cmd.zipmax = bspline.Rmax();
    
    cArray sol;     // stored solution expansion
    cArray ev;      // evaluated solution
    rArray grid;    // real evaluation grid
    
    // size of solution (l,l)-segment
    std::size_t N = bspline.Nspline() * bspline.Nspline();
    
    std::cout << "Zipping B-spline expansion of the solution: \"" << cmd.zipfile << "\"" << std::endl;
    
    // load the requested file
    if (not sol.hdfload(cmd.zipfile.c_str()))
        HexException("Cannot load file %s.", cmd.zipfile.c_str());
    
    // evaluation grid
    grid = linspace (0., cmd.zipmax, cmd.zipcount);
    
    // for all coupled angular momentum pairs
    for (unsigned ill = 0; ill < ll.size(); ill++)
    {
        // angular momenta
        int l1 = ll[ill].first;
        int l2 = ll[ill].second;
        
        std::cout << "\t- partial wave l1 = " << l1 << ", l2 = " << l2 << std::endl;
        
        // write to file
        std::ofstream out (format("%s_(%d,%d).vtk", cmd.zipfile.c_str(), l1, l2));
        writeVTK_points
        (
            out,
            bspline.zip
            (
                sol.slice (ill * N, (ill + 1) * N),
                grid, grid
            ),
            grid, grid, rArray({0.})
        );
    }
}
