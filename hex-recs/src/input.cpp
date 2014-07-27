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
#include <iostream>
#include <string>
#include <tuple>

#include "arrays.h"
#include "cmdline.h"
#include "input.h"
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
    "# angular momenta\n"
    "# J  2M limit\n"
    "  0  0  4\n"
    "\n"
    "# initial atomic states (ni, li, 2ji, 2mi)\n"
    "  1  -1\n"
    "  *\n"
    "  *\n"
    "  *\n"
    "\n"
    "# final atomic states (nf, lf, 2jf)\n"
    "  1  -1\n"
    "  *\n"
    "  *\n"
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
                exit(0);
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
                    "\t--tolerance <number>      (-t)  tolerance for the conjugate gradients solver                                    \n"
                    "\t--drop-tolerance <number> (-d)  drop tolerance for the ILU preconditioner (default: 1e-15)                      \n"
                    "\t--out-of-core             (-O)  use hard disk drive to store intermediate results and thus to save RAM (slower) \n"
                    "\t--parallel-dot                  OpenMP-parallelize SpMV operations                                              \n"
                    "\t--no-parallel-block             disable simultaneous preconditioning of multiple blocks by OpenMP               \n"
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
                zipcount = atol(optarg.c_str());
                return true;
            },
        "zipmax", "R", 1, [&](std::string optarg) -> bool
            {
                // zip bounding box
                zipmax = atof(optarg.c_str());
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
                droptol = atof(optarg.c_str());
                return true;
            },
        "tolerance", "t", 1, [&](std::string optarg) -> bool
            {
                // iteration tolerance for terminating iteration solution
                itertol = atof(optarg.c_str());
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
                std::cout << "\n";
                exit (0);
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
    
    //
    // load B-spline parameters
    //
    
    order = ReadNext<int>(inf).val;
    ecstheta = ReadNext<double>(inf).val;
    
    // print info
    std::cout << "\n-----   B-spline environment  -------\n";
    std::cout << "order = " << order << std::endl;
    std::cout << "ecsθ = " << ecstheta << std::endl;
    
    //
    // load real knot data
    //
    
    std::vector<double> rknots_begin, rknots_end, rknots_samples;
    while ((x = ReadNext<double>(inf).val) != -1.)
        rknots_begin.push_back(x);
    for (size_t i = 0; i < rknots_begin.size(); i++)
        rknots_end.push_back(ReadNext<double>(inf).val);
    for (size_t i = 0; i < rknots_begin.size(); i++)
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
    for (size_t i = 0; i < cknots_begin.size(); i++)
        cknots_end.push_back(ReadNext<double>(inf).val);
    for (size_t i = 0; i < cknots_begin.size(); i++)
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
    // load total quantum numbers atc.
    //
    
    std::cout << "\n----------  Angular momentum limits  -------------\n";
    
    J = ReadNext<int>(inf).val;
    two_M = ReadNext<int>(inf).val;
    levels = ReadNext<int>(inf).val;
    
    std::cout << "J = " << J << std::endl;
    std::cout << "M = " << 0.5 * two_M << std::endl;
    std::cout << "ℓ = " << levels << std::endl;
    
    //
    // load initial atomic quantum numbers
    //
    
    std::cout << "\n----------  Initial atomic states  -------------\n";
    std::vector<ReadItem<int>> nis, lis, two_jis, two_mis, sis;
    
    // - principal quantum number
    while ((idata = ReadNext<int>(inf)).val != -1)
        nis.push_back(idata);
    
    // - orbital angular momentum
    for (size_t i = 0; i < nis.size(); i++)
        lis.push_back(ReadNext<int>(inf, ReadItem<int>::asterisk));
    
    // - spin-orbital angular momentum
    for (size_t i = 0; i < nis.size(); i++)
        two_jis.push_back(ReadNext<int>(inf, ReadItem<int>::asterisk));
    
    // - magnetic quantum number
    for (size_t i = 0; i < nis.size(); i++)
        two_mis.push_back(ReadNext<int>(inf, ReadItem<int>::asterisk));
    
    // - construct initial states
    for (unsigned i = 0; i < nis.size(); i++)
    for (int li = 0; li < nis[i].val; li++)
    for (int two_ji = std::abs(2*li-1); two_ji <= 2*li+1; two_ji += 2)
    for (int two_mi = -two_ji; two_mi <= two_ji; two_mi += 2)
    {
        // skip unused orbital angular momenta
        if (lis[i].val != li and not (lis[i].flags & ReadItem<int>::asterisk))
            continue;
        
        // skip unused spin-orbital angular momenta
        if (two_jis[i].val != two_ji and not (two_jis[i].flags & ReadItem<int>::asterisk))
            continue;
        
        // skip unused angular momentum projections
        if (two_mis[i].val != two_mi and not (two_mis[i].flags & ReadItem<int>::asterisk))
            continue;
        
        // add this initial state
        instates.push_back(std::make_tuple(nis[i].val,li,two_ji,two_mi));
    }
    
    // print info
    std::cout << "[n l 2j 2m]: ";
    for (auto state : instates)
    {
        std::cout << "["
                  << std::get<0>(state) << " "
                  << std::get<1>(state) << " "
                  << std::get<2>(state) << " "
                  << std::get<3>(state)
                  << "] ";
    }
    std::cout << std::endl;
    
    //
    // load final atomic quantum numbers
    //
    
    std::cout << "\n----------  Final atomic states  -------------\n";
    std::vector<int> nfs;
    std::vector<ReadItem<int>> lfs, two_jfs;
    
    // - principal quantum number
    while ((y = ReadNext<int>(inf).val) != -1)
        nfs.push_back(y);
    
    // - orbital angular momentum
    for (size_t i = 0; i < nfs.size(); i++)
        lfs.push_back(ReadNext<int>(inf, ReadItem<int>::asterisk));
    
    // - spin-orbital angular momentum
    for (size_t i = 0; i < nfs.size(); i++)
        two_jfs.push_back(ReadNext<int>(inf, ReadItem<int>::asterisk));
    
    // - construct final states
    for (unsigned f = 0; f < nfs.size(); f++)
    for (int lf = 0; lf < nfs[f]; lf++)
    for (int two_jf = std::abs(2*lf-1); two_jf <= 2*lf+1; two_jf += 2)
    for (int two_mf = -two_jf; two_mf <= two_jf; two_mf += 2)
    {
        // skip unused orbital angular momenta
        if (lfs[f].val != lf and not (lfs[f].flags & ReadItem<int>::asterisk))
            continue;
        
        // skip unused spin-orbital angular momenta
        if (two_jfs[f].val != two_jf and not (two_jfs[f].flags & ReadItem<int>::asterisk))
            continue;
        
        // add this initial state
        outstates.push_back(std::make_tuple(nfs[f],lf,two_jf,two_mf));
    }
    
    // print info
    std::cout << "[n l 2j 2m]: ";
    for (auto state : outstates)
    {
        std::cout << "["
                  << std::get<0>(state) << " "
                  << std::get<1>(state) << " "
                  << std::get<2>(state) << " "
                  << std::get<3>(state)
                  << "] ";
    }
    std::cout << std::endl;
    
    //
    // load initial energies
    //
    
    std::cout << "\n---  Initial projectile energies  ----\n";
    std::vector<double> Ei_begin, Ei_end, Ei_samples;
    
    while ((x = ReadNext<double>(inf).val) != -1.)
        Ei_begin.push_back(x);
    for (size_t i = 0; i < Ei_begin.size(); i++)
        Ei_end.push_back(ReadNext<double>(inf).val);
    for (size_t i = 0; i < Ei_begin.size(); i++)
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
