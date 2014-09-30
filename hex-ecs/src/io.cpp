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
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <string>
#include <tuple>

#include <omp.h>

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
    " 0.0  0.1   3   -1\n"
    " 0.0  2.0  60\n"
    "   4   20  58\n"
    "\n"
    "# complex knot sequences\n"
    "  60    -1\n"
    " 100\n"
    "  41\n"
    "\n"
    "# initial atomic states\n"
    "# ni\n"
    "  1\n"
    "# angular states (li, mi)\n"
    "  0  -1\n"
    "  0\n"
    "\n"
    "# final atomic states (nf, lf)\n"
    "  1  -1\n"
    "  0\n"
    "\n"
    "# angular momenta\n"
    "# L  Pi limit\n"
    "  0  0  4\n"
    "\n"
    "# initial energies in Rydbergs\n"
    " 0.65   -1\n"
    " 0.95\n"
    "    3\n"
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
                    throw exception ("Error: Cannot write to \"example.inp\"\n");
                
                out << sample_input;
                    
                out.close();
                std::exit(0);
            },
        "input", "i", 1, [&](std::string optarg) -> bool
            {
                // set custom input file
                inputfile.open(optarg);
                if (not inputfile.good())
                    throw exception ("Error: Input file \"%s\" not found.\n", optarg.c_str());
                return true;
            },
        "help", "h", 0, [&](std::string optarg) -> bool
            {
                // print usage information
                std::cout << "\n"
                    "Available switches (short forms in parentheses):                                                                  \n"
                    "                                                                                                                  \n"
                    "\t--example                 (-e)  create sample input file                                                        \n"
                    "\t--help                    (-h)  display this help                                                               \n"
                    "\t--input <filename>        (-i)  use custom input file                                                           \n"
                    "\t--zipfile <filename>      (-z)  solution file to zip                                                            \n"
                    "\t--zipcount <number>       (-n)  zip samples                                                                     \n"
                    "\t--zipmax <number>         (-R)  maximal radius to use for solution zipping                                      \n"
#ifndef NO_MPI
                    "\t--mpi                     (-m)  use MPI                                                                         \n"
#endif
                    "\t--stg-integ               (-a)  only do radial integrals                                                        \n"
                    "\t--stg-integ-solve         (-b)  only do integrals & solve                                                       \n"
                    "\t--stg-extract             (-c)  only extract amplitudes                                                         \n"
                    "\t--preconditioner <name>   (-p)  preconditioner to use (default: ILU)                                            \n"
                    "\t--list-preconditioners    (-P)  list of available preconditioners with short description of each                \n"
                    "\t--tolerance <number>      (-T)  tolerance for the conjugate gradients solver                                    \n"
                    "\t--prec-tolerance <number> (-t)  tolerance for the conjugate gradients preconditioner                            \n"
                    "\t--drop-tolerance <number> (-d)  drop tolerance for the ILU preconditioner (default: 1e-15)                      \n"
                    "\t--out-of-core             (-O)  use hard disk drive to store intermediate results and thus to save RAM (slower) \n"
                    "\t--parallel-dot                  OpenMP-parallelize SpMV operations                                              \n"
                    "\t--no-parallel-block             disable simultaneous preconditioning of multiple blocks by OpenMP               \n"
                    "\t--gpu-slater                    compute diagonal two-electron interals using OpenCL                             \n"
                    "                                                                                                                  \n"
                ;
                exit(0);
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
        "out-of-core", "O", 0, [&](std::string optarg) -> bool
            {
                // use out-of-core functionality: store diagonal blocks on disk
                outofcore = true;
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
                    throw exception("Unknown preconditioner \"%s\".", optarg.c_str());
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
                std::exit(0);
            },
        "parallel-dot", "", 0, [&](std::string optarg) -> bool
            {
                // parallelize SpMV
                parallel_dot = true;
                return true;
            },
        "no-parallel-block", "", 0, [&](std::string optarg) -> bool
            {
                // un-parallelize preconditioning
                parallel_block = false;
                return true;
            },
        "gpu-slater", "", 0, [&](std::string optarg) -> bool
            {
                // compute diagonal two-electron integrals using OpenCL
                gpu_slater = true;
                return true;
            },
        
        [&] (std::string optname, std::string optarg) -> bool
        {
            throw exception ("Unknown switch \"%s\".", optname.c_str());
        }
    );
}

void InputFile::read (std::ifstream & inf)
{
    double x; int y;
    ReadItem<int> idata;
    
    // load B-spline parameters
    order = ReadNext<int>(inf).val;
    ecstheta = ReadNext<double>(inf).val;
    
    std::cout << "\n-----   B-spline environment  -------\n";
    std::cout << "order = " << order << "\n";
    std::cout << "ecsθ = " << ecstheta << "\n";
    
    //
    // load real knot data
    //
    
    std::vector<double> rknots_begin, rknots_end, rknots_samples;
    while ((x = ReadNext<double>(inf).val) != -1.)
        rknots_begin.push_back(x);
    for (std::size_t i = 0; i < rknots_begin.size(); i++)
        rknots_end.push_back(ReadNext<double>(inf).val);
    for (std::size_t i = 0; i < rknots_begin.size(); i++)
        rknots_samples.push_back(ReadNext<double>(inf).val);
    
    // construct real knot sequence
    for (unsigned i = 0; i < rknots_begin.size(); i++)
    {
        if (rknots_begin[i] > rknots_end[i])
        {
            std::cout << "\t" << rknots_begin[i] << " > " << rknots_end[i] << "\n";
            throw exception ("Inconsistent knot specification!");
        }
        
        rArray new_knots = linspace(rknots_begin[i], rknots_end[i], rknots_samples[i]);
        rknots = concatenate(rknots, new_knots);
    }
    
    // print info
    std::cout << "\n----------   Real knots  ------------\n";
    std::cout << rknots << std::endl;
    
    //
    // load complex knot data
    //
    
    std::vector<double> cknots_begin, cknots_end, cknots_samples;
    while ((x = ReadNext<double>(inf).val) != -1.)
        cknots_begin.push_back(x);
    for (std::size_t i = 0; i < cknots_begin.size(); i++)
        cknots_end.push_back(ReadNext<double>(inf).val);
    for (std::size_t i = 0; i < cknots_begin.size(); i++)
        cknots_samples.push_back(ReadNext<double>(inf).val);
    
    // construct complex(-to-become) knot sequence
    for (unsigned i = 0; i < cknots_begin.size(); i++)
    {
        cknots = concatenate
        (
            cknots,
            linspace
            (
                cknots_begin[i],
                cknots_end[i],
                cknots_samples[i]
            )
        );
    }
    
    // print info
    std::cout << "\n---------  Complex knots  ------------\n";
    std::cout << cknots << std::endl;
    
    //
    // load initial atomic quantum numbers
    //
    
    std::cout << "\n----------  Initial atomic states  -------------\n";
    
    std::vector<ReadItem<int>> lis, mis;
    
    // load initial principal quantum number
    ni = ReadNext<int>(inf).val;
    
    // - orbital angular momentum
    while ((idata = ReadNext<int>(inf)).val != -1)
        lis.push_back(idata);
    
    // - magnetic quantum number
    for (std::size_t i = 0; i < lis.size(); i++)
        mis.push_back(ReadNext<int>(inf, ReadItem<int>::asterisk));
    
    for (unsigned i = 0; i < lis.size(); i++)
    for (int li = 0; li < ni; li++)
    for (int mi = -li; mi <= li; mi++)
    {
        // skip unused orbital angular momenta
        if (lis[i].val != li and not (lis[i].flags & ReadItem<int>::asterisk))
            continue;
        
        // skip unused angular momentum projections
        if (mis[i].val != mi and not (mis[i].flags & ReadItem<int>::asterisk))
            continue;
        
        // add this initial state
        instates.push_back(std::make_tuple(ni,li,mi));
    }
    
    // print info
    std::cout << "[n l m]: ";
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
    
    std::cout << "\n----------  Final atomic states  -------------\n";
    std::vector<int> nfs;
    std::vector<ReadItem<int>> lfs;
    
    // - principal quantum number
    while ((y = ReadNext<int>(inf).val) != -1)
        nfs.push_back(y);
    
    // - orbital angular momentum
    for (std::size_t i = 0; i < nfs.size(); i++)
        lfs.push_back(ReadNext<int>(inf, ReadItem<int>::asterisk));
    
    // - construct final states
    for (unsigned f = 0; f < nfs.size(); f++)
    for (int lf = 0; lf <= nfs[f]; lf++)
    {
        // l=n only in ionization specification
        if (lf == nfs[f] and nfs[f] != 0)
            continue;
        
        // skip unused orbital angular momenta
        if (lfs[f].val != lf and not (lfs[f].flags & ReadItem<int>::asterisk))
            continue;
        
        // add this initial state
        outstates.push_back(std::make_tuple(nfs[f],lf,0));
    }
    
    // print info
    std::cout << "[n l m]: ";
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
    
    std::cout << "\n----------  Angular momentum limits  -------------\n";
    
    L = ReadNext<int>(inf).val;
    Pi = ReadNext<int>(inf).val % 2;
    levels = ReadNext<int>(inf).val;
    
    std::cout << "L = " << L << "\n";
    std::cout << "Π = " << Pi << "\n";
    std::cout << "ℓ = " << levels << "\n";
    
    //
    // load initial energies
    //
    
    std::cout << "\n---  Initial projectile energies  ----\n";
    std::vector<double> Ei_begin, Ei_end, Ei_samples;
    
    while ((x = ReadNext<double>(inf).val) != -1.)
        Ei_begin.push_back(x);
    for (std::size_t i = 0; i < Ei_begin.size(); i++)
        Ei_end.push_back(ReadNext<double>(inf).val);
    for (std::size_t i = 0; i < Ei_begin.size(); i++)
        Ei_samples.push_back(ReadNext<double>(inf).val);
    
    // construct energy sequence
    for (unsigned i = 0; i < Ei_begin.size(); i++)
    {
        Ei = concatenate
        (
            Ei,
            linspace
            (
                Ei_begin[i],
                Ei_end[i],
                Ei_samples[i]
            )
        );
    }
    
    // print info
    std::cout << "lowest energy: "    << Ei.front() << std::endl;
    std::cout << "highest energy: "   << Ei.back()  << std::endl;
    std::cout << "total enegies: "    << Ei.size()  << std::endl;
    std::cout << "full energy list: " << Ei         << std::endl;
    
    //
    // load some other optional data
    //
    
    std::cout << "\n---------- Other parameters -----------\n";
    B = ReadNext<double>(inf).val;
    
    // print info
    std::cout << "magnetic field: " << B << " a.u." << std::endl;
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
        throw exception("Cannot load file %s.", cmd.zipfile.c_str());
    
    // evaluation grid
    grid = linspace (0., cmd.zipmax, cmd.zipcount);
    
    // for all coupled angular momentum pairs
    for (unsigned ill = 0; ill < ll.size(); ill++)
    {
        // angular momenta
        int l1 = ll[ill].first;
        int l2 = ll[ill].second;
        
        std::cout << "\t- partial wave l1 = " << l1 << ", l2 = " << l2 << "\n";
        
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
