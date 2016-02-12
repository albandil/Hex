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

#include "hex-cmdline.h"
#include "hex-special.h"

#include "io.h"

const std::string sample_input =
    "# This is a sample input file for the program hex-pecs.\n"
    "# Lines introduced by the hash symbol are comments, which document the options,\n"
    "# and can be omitted. Otherwise there are just a few numbers and characters.\n"
    "\n"
    "# ECS rotation angle in radians\n"
    "  0.78\n"
    "\n"
    "# Radial grid (real part).\n"
    "  L  0.00   0.20   21\n"
    "  L  0.26   2.00   30\n"
    "  L  2.20   20.0   90\n"
    "  L  20.4    200  450\n"
    " -1\n"
    "\n"
    "# Radial grid (complex part).\n"
    "  L   0.0    2.0   21\n"
    "  L   2.5   20.0   36\n"
    " -1\n"
    "\n"
    "# --------------- Atomic states -------------------\n"
    "\n"
    "# Initial atomic states (ni, li, mi).\n"
    "# Specified as vertical triplets terminated by -1 on the first line.\n"
    "# Computation of all angular quantum numbers can be requested by asterisk.\n"
    "  1  -1\n"
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
    "  0  0  0  4\n"
    "\n"
    "# Atom + projectile total energy in Rydbergs.\n"
    "  1\n"
    "\n"
    "# Nuclear charge.\n"
    "  1\n"
    ;

CommandLine::CommandLine (int argc, char* argv[])
{
    // default values
    fully_coupled = false;
    group_coupled = false;
    prepare_grid = true;
    propagate_states = true;
    max_iter = 100;
    inputfile = "pecs.inp";
    itertol = 1e-8;
    Epert.push_back(0.);
    
    // custom values
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
        "help", "h", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                // print usage information
                std::cout << "\n"
                    "Available switches (short forms in parentheses):                                                                        \n"
                    "                                                                                                                        \n"
                    "\t--example                  (-e)  Create sample input file.                                                            \n"
                    "\t--help                     (-h)  Display this help.                                                                   \n"
                    "\t--input <filename>         (-i)  Use custom input file (other than \"pecs.inp\").                                     \n"
                    "\t--propagate-only           (-p)  Skip preparation of propagation matrices and only propagate solutions.               \n"
                    "\t--extract-only             (-x)  Only extract T-matrices and cross sections.                                          \n"
                    "\t--fully-coupled            (-c)  Solve fully angularly coupled system.                                                \n"
                    "\t--group-coupled            (-g)  Solve partially coupled system (group by nL).                                        \n"
                    "\t--max-iter <number>              Maximal number of iterative coupling iterations.                                     \n"
                    "\t--tolerance                      Iterative coupling relative tolerance.                                               \n"
                    "\t--energy-pert <list>       (-E)  Solve several close energie in one run. Perturbations are in Ry, separated by spaces.\n"
                    "\t                                 If no perturbation is given, they will be read from the standard input.              \n"
                    "\n"
                ;
                std::exit(EXIT_SUCCESS);
            },
        "input", "i", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                // check custom input file
                if (not std::ifstream(optargs[0]).good())
                    HexException("Error: Input file \"%s\" not found.\n", optargs[0].c_str());
                inputfile = optargs[0];
                return true;
            },
        "propagate-only", "p", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                prepare_grid = false;
                return true;
            },
        "extract-only", "x", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                prepare_grid = false;
                propagate_states = false;
                return true;
            },
        "fully-coupled", "c", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                fully_coupled = true;
                max_iter = 1;
                return true;
            },
        "group-coupled", "g", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                group_coupled = true;
                return true;
            },
        "max-iter", "", 0, [&](std::vector<std::string> const & optargs) -> bool
            {
                max_iter = std::stoi(optargs[0]);
                return true;
            },
        "tolerance", "", 1, [&](std::vector<std::string> const & optargs) -> bool
            {
                itertol = std::stod(optargs[0]);
                return true;
            },
        "energy-pert", "E", -1, [&](std::vector<std::string> const & optargs) -> bool
            {
                if (optargs.size() == 0)
                {
                    // read from STDIN
                    double e;
                    while (std::cin >> e)
                        Epert.push_back(e);
                }
                else
                {
                    // read from given list
                    for (std::string const & e : optargs)
                        Epert.push_back(std::stod(e));
                }
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
            
            if (begin > end)
                HexException("Start of linear sequence is larger than its end (%g > %g).", begin, end);
            if (samples < 0)
                HexException("Invalid number of samples for linear sequence: %d.", samples);
            
            arr.append(linspace(begin, end, samples));
        }
        else if (type[0] == 'G')
        {
            double begin = ReadNext<double>(inf).val;
            double end = ReadNext<double>(inf).val;
            double d = ReadNext<double>(inf).val;
            double quotient = ReadNext<double>(inf).val;
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

InputFile::InputFile (const CommandLine& cmd)
{
    std::ifstream inf (cmd.inputfile);
    
    if (not inf.good())
        HexException("Failed to open file \"%s\".", cmd.inputfile.c_str());
    
    ReadItem<int> idata;
    
    ecstheta = ReadNext<double>(inf).val;
    
    ReadArrays(inf, rgrid);
    ReadArrays(inf, cgrid);
    
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
        istates.push_back(HState{nis[i].val,li,mi});
    }
    
    // print info
    std::cout << "\t[n l m]: ";
    for (HState const & state  : istates)
    {
        std::cout << "["
                  << state.n << " "
                  << state.l << " "
                  << state.m
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
        fstates.push_back(HState{nfs[f].val,lf,0});
    }
    
    // print info
    std::cout << "\t[n l m]: ";
    for (HState const & state : fstates)
    {
        std::cout << "["
                  << state.n << " "
                  << state.l
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
    S = ReadNext<int>(inf).val;
    
    // parity
    Pi = ReadNext<int>(inf).val % 2;
    
    // number of angular momentum pairs per total angular momentum
    nL = ReadNext<int>(inf).val;
    
    std::cout << "\tL = " << L << std::endl;
    std::cout << "\tS = " << S << std::endl;
    std::cout << "\tPi = " << Pi << std::endl;
    std::cout << "\tnL = " << nL << std::endl;
    
    //
    // load initial energies
    //
    
    Etot = ReadNext<double>(inf).val;
    
    // print info
    std::cout << std::endl << "Total projectile + atom energy [Ry]: " << Etot << std::endl;
    
    //
    // load the nuclear charge
    //
    
    Z = ReadNext<double>(inf).val;
    
    // print info
    std::cout << std::endl << "Nuclear charge: " << Z << std::endl << std::endl;
}
