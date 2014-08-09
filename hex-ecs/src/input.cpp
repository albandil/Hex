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
                    "\t--concurrent-factorizations <number>   how many LU preconditioner factorizations to run simultaneously          \n"
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
        "concurrent-factorizations", "", 0, [&](std::string optarg) -> bool
            {
                // parallelize LU factorizations
                concurrent_factorizations = std::max(1,std::atoi(optarg.c_str()));
                return true;
            },
        
        [&] (std::string optname, std::string optarg) -> bool
        {
            throw exception ("Unknown switch \"%s\".", optname.c_str());
        }
    );
}

long read_int (std::ifstream& f)
{
    // text buffer
    std::string s;
    
    while (not f.eof())
    {
        // read string
        f >> s;
        
        // check length
        if (s.size() == 0)
            continue;
        
        // check if it is a beginning of a comment
        if (s[0] == '#')
        {
            // get the rest of the line
            std::getline(f, s);
            continue;
        }
        
        break;
    }
    
    // convert to long
    char* tail;
    long val = strtol(s.c_str(), &tail, 10);
    if (*tail != 0)
    {
        if (s == "*")
            throw false;
        else
            throw exception ("Can't read int.\n");
    }
    else
    {
        return val;
    }
}

double read_dbl (std::ifstream& f)
{
    // text buffer
    std::string s;
    
    while (not f.eof())
    {
        // read string
        f >> s;
        
        // check length
        if (s.size() == 0)
            continue;
        
        // check if it is a beginning of a comment
        if (s[0] == '#')
        {
            // get the rest of the line
            std::getline(f, s);
            continue;
        }
        
        break;
    }
    
    // convert to double
    char* tail;
    double val = strtod(s.c_str(), &tail);
    if (*tail != 0)
    {
        if (s == "*")
            throw false;
        else
            throw exception ("Can't read double.\n");
    }
    else
    {
        return val;
    }
}

void InputFile::read (std::ifstream & inputfile)
{
    double x;
    
    // load B-spline parameters
    try {
        order = read_int(inputfile);
        ecstheta = read_dbl(inputfile);
    } catch (std::exception e) {
        std::cerr << e.what() << std::endl;
        throw exception("Input error: Check B-spline parameters.\n");
    } catch (bool b) {
        throw exception("Wildcard not accepted here.\n");
    }
    
    std::cout << "\n-----   B-spline environment  -------\n";
    std::cout << "order = " << order << "\n";
    std::cout << "ecsθ = " << ecstheta << "\n";
    
    // load real knot data
    std::vector<double> rknots_begin, rknots_end, rknots_samples;
    try {
        while ((x = read_dbl(inputfile)) != -1.)
            rknots_begin.push_back(x);
        for (size_t i = 0; i < rknots_begin.size(); i++)
            rknots_end.push_back(read_dbl(inputfile));
        for (size_t i = 0; i < rknots_begin.size(); i++)
            rknots_samples.push_back(read_dbl(inputfile));
    } catch (std::exception e) {
        std::cerr << e.what() << std::endl;
        throw exception("Input error: Check real knot data.\n");
    } catch (bool b) {
        throw exception("Wildcard not accepted here.\n");
    }
    
    // construct real knot sequence
    for (unsigned i = 0; i < rknots_begin.size(); i++)
    {
        if (rknots_begin[i] > rknots_end[i])
        {
            std::cout << "\t" << rknots_begin[i] << " > " << rknots_end[i] << "\n";
            throw exception("Inconsistent knot specification!");
        }
        
        auto new_knots = linspace(rknots_begin[i], rknots_end[i], rknots_samples[i]);
        rknots = concatenate(rknots, new_knots);
    }
    
    std::cout << "\n----------   Real knots  ------------\n";
    for (auto knot = rknots.begin(); knot != rknots.end(); knot++)
        std::cout << *knot << " ";
    std::cout << std::endl;
    
    // load complex knot data
    std::vector<double> cknots_begin, cknots_end, cknots_samples;
    try {
        while ((x = read_dbl(inputfile)) != -1.)
            cknots_begin.push_back(x);
        for (size_t i = 0; i < cknots_begin.size(); i++)
            cknots_end.push_back(read_dbl(inputfile));
        for (size_t i = 0; i < cknots_begin.size(); i++)
            cknots_samples.push_back(read_int(inputfile));
    } catch (std::exception e) {
        std::cerr << e.what() << std::endl;
        throw exception("Input error: Check complex knot data.\n");
    } catch (bool b) {
        throw exception("Wildcard not accepted here.\n");
    }
    
    // construct complex(-to-be) knot sequence
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
    
    std::cout << "\n---------  Complex knots  ------------\n";
    for (auto knot = cknots.begin(); knot != cknots.end(); knot++)
        std::cout << *knot << " ";
    std::cout << std::endl;
    
    // load initial principal quantum number
    try {
        ni = read_int(inputfile);
    } catch (std::exception e) {
        std::cerr << e.what() << std::endl;
        throw exception("Input error: Check \"ni\".\n");
    } catch (bool b) {
        throw exception("Wildcard not accepted here.\n");
    }
    
    // load initial atomic angular states
    std::vector<int> lis, mis;
    int maxli = 0;
    try {
        while ((x = read_int(inputfile)) != -1.)
        {
            lis.push_back(x);
            if (lis.back() >= ni)
                throw exception("Input error: Angular momentum greater than \"ni\".\n");
            if (lis.back() > maxli)
                maxli = lis.back();
        }
    } catch (std::exception e) {
        std::cerr << e.what() << std::endl;
        throw exception("Input error: Check initial atomic state data.\n");
    } catch (bool b) {
        throw exception("Wildcard not accepted here.\n");
    }
    
    for (size_t i = 0; i < lis.size(); i++)
    {
        try {
            
            mis.push_back(read_int(inputfile));
            if (std::abs(mis[i]) > lis[i])
                throw exception("Input error: Magnetic quantum number greater than \"li\".\n");
            
            instates.push_back(std::make_tuple(ni,lis[i],mis[i]));
            
        } catch (std::exception e) {
            
            std::cerr << e.what() << std::endl;
            throw exception("Input error: Check initial atomic state data.\n");
            
        } catch (bool b) {
            
            // wildcard "*" found
            for (int j = -lis[i]; j <= lis[i]; j++)
            {
                mis.push_back(j);
                instates.push_back(std::make_tuple(ni,lis[i],j));
            }
            
        }
    }
    
    // load final atomic quantum numbers
    std::vector<int> nfs, lfs;
    int maxlf = 0;
    try {
        while ((x = read_int(inputfile)) != -1.)
            nfs.push_back(x);
    } catch (std::exception e) {
        std::cerr << e.what() << std::endl;
        throw exception("Input error: Check final atomic state data.\n");
    } catch (bool b) {
        throw exception("Wildcard not accepted here.\n");
    }
    
    for (size_t i = 0; i < nfs.size(); i++)
    {
        try {
            
            lfs.push_back(read_int(inputfile));
            
            if (nfs[i] == 0)
            {
                outstates.push_back(std::make_tuple(0,0,0));
                continue;
            }
            
            if (lfs[i] > nfs[i])
                throw exception("Input error: Angular momentum greater than \"nf\".\n");
            if (lfs[i] > maxlf)
                maxlf = lfs[i];
            
            outstates.push_back(std::make_tuple(nfs[i],lfs[i],0));
            
        } catch (std::exception e) {
            
            std::cerr << e.what() << std::endl;
            throw exception("Input error: Check final atomic state data.\n");
            
        } catch (bool b) {
            
            if (nfs[i] == 0)
            {
                outstates.push_back(std::make_tuple(0,0,0));
                continue;
            }
            
            // wildcard "*" found, add all allowed angular momenta
            for (int j = 0; j < nfs[i]; j++)
            {
                lfs.push_back(j);
                outstates.push_back(std::make_tuple(nfs[i],j,0));
            }
            
        }
    }
    
    
    // load total quantum numbers
    try {
        
        L = read_int(inputfile);
        Pi = read_int(inputfile) % 2;
        levels = read_int(inputfile);
        
        if (L + L%2 + levels < maxli)
            throw exception("Input error: ℓ is smaller than some initial angular momenta.\n");
        
        if (L + L%2 + levels < maxlf)
            throw exception("Input error: ℓ is smaller than some final angular momenta.\n");
        
    } catch (std::exception e) {
        
        std::cerr << e.what() << std::endl;
        throw exception("Input error: Check angular momentum data.\n");
        
    } catch (bool b) {
        throw exception("Wildcard not accepted here.\n");
    }
    
    std::cout << "\n----------  Angular momentum limits  -------------\n";
    std::cout << "L = " << L << "\n";
    std::cout << "Π = " << Pi << "\n";
    std::cout << "ℓ = " << levels << "\n";
    
    std::cout << "\n----------  Initial atomic states  -------------\n";
    for (auto state : instates)
        std::cout << "[" << std::get<0>(state) << " " << std::get<1>(state) << " " << std::get<2>(state) << "] ";
    std::cout << "\n";
    
    std::cout << "\n----------  Final atomic states  -------------\n";
    for (auto state : outstates)
        std::cout << "[" << std::get<0>(state) << " " << std::get<1>(state) << " " << std::get<2>(state) << "] ";
    std::cout << "\n";
    
    // load initial energies
    std::vector<double> Ei_begin, Ei_end, Ei_samples;
    try {
        while ((x = read_dbl(inputfile)) != -1.)
            Ei_begin.push_back(x);
        for (size_t i = 0; i < Ei_begin.size(); i++)
            Ei_end.push_back(read_dbl(inputfile));
        for (size_t i = 0; i < Ei_begin.size(); i++)
            Ei_samples.push_back(read_int(inputfile));
    } catch (std::exception e) {
        std::cerr << e.what() << std::endl;
        throw exception("Input error: Check energy data.\n");
    } catch (bool b) {
        throw exception("Wildcard not accepted here.\n");
    }
    
    // construct energy sequence
    for (unsigned i = 0; i < Ei_begin.size(); i++)
        Ei = concatenate(Ei, linspace(Ei_begin[i], Ei_end[i], Ei_samples[i]));
    
    std::cout << "\n---  Initial projectile energies  ----\n";
    std::cout << "lowest energy: " << Ei.front() << "\n";
    std::cout << "highest energy: " << Ei.back() << "\n";
    std::cout << "total enegies: " << Ei.size() << "\n";
    std::cout << "full energy list: " << Ei << "\n";
    
    try {
        B = read_dbl(inputfile);
    } catch (std::exception e) {
        std::cerr << e.what() << std::endl;
        throw exception("Input error: Check magnetic field data.\n");
    } catch (bool b) {
        throw exception("Wildcard not accepted here.\n");
    }
    
    std::cout << "\n---------- Other parameters -----------\n";
    std::cout << "magnetic field: " << B << " a.u.\n\n";
}
