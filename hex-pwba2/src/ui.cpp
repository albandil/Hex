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

#include <clocale>
#include <cmath>
#include <iostream>

#include <gsl/gsl_errno.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "hex-arrays.h"
#include "hex-cmdline.h"
#include "hex-special.h"
#include "hex-version.h"

#include "pwba2.h"

const std::string sample_input =
    "# ------------------------\n"
    "# Quantum state numbers\n"
    "# ------------------------\n"
    "\n"
    "# L  Pi\n"
    "  0  0\n"
    "\n"
    "# initial state\n"
    "# ni li\n"
    "   1  0\n"
    "\n"
    "# final state\n"
    "# nf lf\n"
    "   1  0\n"
    "\n"
    "# impact energy\n"
    "# Ei\n"
    "   4\n"
    "\n"
    "# ------------------------\n"
    "# Grid parameters\n"
    "# ------------------------\n"
    "\n"
    "# maximal radius\n"
    "# Rmax\n"
    "   100\n"
    "\n"
    "# linear samples\n"
    "# N\n"
    "  1000\n"
    "\n"
    "# ------------------------\n"
    "# Intermediate states\n"
    "# ------------------------\n"
    "\n"
    "# maximal quantum numbers\n"
    "# maxNn  nL     maxEn\n"
    "  8       3     20\n"
    "\n"
    "# continuum integration\n"
    "# allowed forbidden\n"
    "  1       1\n"
;

const std::string help_text = 
    "                                                                                                                  \n"
    "Available switches (short forms in parentheses):                                                                  \n"
    "                                                                                                                  \n"
    "\t--example                 (-e)  create sample input file                                                        \n"
    "\t--help                    (-h)  display this help                                                               \n"
    "\t--verbose                 (-v)  make the program richly comment all steps                                       \n"
    "\t--input <filename>        (-i)  use custom input file                                                           \n"
    "\t--direct-integrate        (-d)  compute exact second-order T-matrix for given angle                             \n"
    "\t--partial-wave            (-w)  compute only contribution of single partial wave                                \n"
    "                                                                                                                  \n"
;

void parse_input_file
(
    std::ifstream & inf,
    int & L, int & Pi,
    int & Ni, int & Li, int & Nf, int & Lf, double & Ei,
    double & Rmax, int & N,
    int & maxNn, int & maxLn, double & Enmax,
    bool & integrate_allowed, bool & integrate_forbidden
)
{
    // read all data
    L = ReadNext<int>(inf).val;  Pi = ReadNext<int>(inf).val;
    Ni = ReadNext<int>(inf).val; Li = ReadNext<int>(inf).val;
    Nf = ReadNext<int>(inf).val; Lf = ReadNext<int>(inf).val;
    Ei = ReadNext<double>(inf).val;
    Rmax = ReadNext<double>(inf).val;
    N = ReadNext<int>(inf).val;
    maxNn = ReadNext<int>(inf).val;
    maxLn = ReadNext<int>(inf).val;
    Enmax = ReadNext<double>(inf).val;
    integrate_allowed = ReadNext<int>(inf).val;
    integrate_forbidden = ReadNext<int>(inf).val;
}

typedef std::vector<std::string> const & Args;

int main (int argc, char* argv[])
{
    // set proper locale
    std::setlocale(LC_ALL, "en_GB.utf8");
    
    // print program logo
    std::cout << logo(" ") << std::endl;
    std::cout << "=== Plane wave second Born approximation ===" << std::endl << std::endl;
    
    // disable fatal GSL errors
    gsl_set_error_handler_off();
    
    // disable buffering of the standard output (-> immediate logging)
    setvbuf(stdout, nullptr, _IONBF, 0);
    
    // input file
    std::ifstream inputfile;
    
    // computation mode
    bool partial_wave = false, direct_integrate = false;
    
    // verbosity level
    bool verbose = false;
    
    // parse command line
    ParseCommandLine
    (
        argc, argv,
        
        "example", "e", 0, [&](Args optargs) -> bool
            {
                std::cout << "Writing sample input file to \"example.inp\"." << std::endl << std::endl;
                
                // produce sample input file
                std::ofstream out("example.inp");
                if (not out.good())
                    throw exception ("Error: Cannot write to \"example.inp\".");
                
                out << sample_input;
                    
                out.close();
                std::exit(EXIT_SUCCESS);
            },
        "help", "h", 0, [](Args optargs) -> bool
            {
                // print usage information
                std::cout << help_text << std::endl;
                std::exit(EXIT_SUCCESS);
            },
        "input", "i", 1, [&](Args optargs) -> bool
            {
                // set custom input file
                inputfile.open(optargs[0]);
                if (not inputfile.good())
                    throw exception ("Error: Input file \"%s\" not found.", optargs[0].c_str());
                std::cout << "Using input file \"" << optargs[0] << "\"." << std::endl << std::endl;
                return true;
            },
        "partial-wave", "w", 0, [&](Args optargs) -> bool
            {
                // compute only contribution of a single partial wave
                partial_wave = true;
                return true;
            },
        "direct-integrate", "d", 0, [&](Args optargs) -> bool
            {
                // compute multidimensional integral
                direct_integrate = true;
                return true;
            },
        "verbose", "v", 0, [&](Args optargs) -> bool
            {
                // verbose output
                verbose = true;
                return true;
            },
        
        [](std::string opt, Args optargs) -> bool
            {
                HexException("Unknown option \"%s\".", opt.c_str());
            }
    );
    
    // check mode
    if (direct_integrate and partial_wave)
    {
        std::cout << "Please select either --partial-wave, or --direct-integrate, not both." << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (not direct_integrate and not partial_wave)
    {
        // print usage information
        std::cout << help_text << std::endl;
        std::exit(EXIT_SUCCESS);
    }
    
    // check input file
    if (not inputfile.is_open())
    {
        const char * filename = "pwba2.inp";
        std::cout << "Using input file \"" << filename << "\"." << std::endl << std::endl;
        inputfile.open(filename);
        if (not inputfile.good())
        {
            std::cout << "Cannot open the file \"" << filename << "\"." << std::endl;
            std::cout << std::endl;
            std::cout << "Either (1) provide input settings in the file \"pwba2.inp\", " << std::endl;
            std::cout << "    or (2) give another name using the '--input' command line option." << std::endl;
            std::cout << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    
    // grid parameters
    int N;
    double Rmax;
    
    // quantum numbers
    int Pi, L, maxNn, nL;
    int Ni, Li, Nf, Lf;
    double Ei, Enmax;
    bool integrate_allowed;
    bool integrate_forbidden;
    
    // parse input file
    parse_input_file
    (
        inputfile,
        L, Pi,
        Ni, Li, Nf, Lf, Ei,
        Rmax, N,
        maxNn, nL, Enmax,
        integrate_allowed, integrate_forbidden
    );
    
    // compute other variables from input
    double ki = std::sqrt(Ei);
    double Etot = ki*ki - 1./(Ni*Ni);
    rArray grid = linspace(0., Rmax, N + 1);
    
    // echo input data
    if (partial_wave)
    {
        std::cout << "Quantum state parameters:" << std::endl;
        std::cout << "\t- total angular momentum: L = " << L << std::endl;
        std::cout << "\t- total parity: Π = " << Pi << std::endl;
        std::cout << "\t- initial atomic state: Ni = " << Ni << ", Li = " << Li << std::endl;
        std::cout << "\t- final atomic state: Nf = " << Nf << ", Lf = " << Lf << std::endl;
        std::cout << "\t- impact energy: Ei = " << ki * ki << std::endl;
        std::cout << "\t- total energy: Etot = " << Etot << std::endl;
        std::cout << std::endl;
        std::cout << "Grid parameters:" << std::endl;
        std::cout << "\t- grid length: Rmax = " << Rmax << std::endl;
        std::cout << "\t- grid total samples: N = " << N << std::endl;
        std::cout << "\t- grid spacing: h = " << Rmax / N << std::endl;
        std::cout << std::endl;
        std::cout << "Intermediate atomic states:" << std::endl;
        std::cout << "\t- maximal bound state principal quantum number: maxNn = " << maxNn << std::endl;
        std::cout << "\t- maximal intermediate angular momentum sum (- L): nL = " << nL << std::endl;
        std::cout << "\t- integrate allowed states: " << (integrate_allowed ? "yes" : "no") << std::endl;
        std::cout << "\t- integrate forbidden states: " << (integrate_forbidden ? "yes" : "no") << std::endl;
    }
    else
    {
        std::cout << "Quantum state parameters:" << std::endl;
        std::cout << "\t- initial atomic state: Ni = " << Ni << ", Li = " << Li << std::endl;
        std::cout << "\t- final atomic state: Nf = " << Nf << ", Lf = " << Lf << std::endl;
        std::cout << "\t- impact energy: Ei = " << ki * ki << std::endl;
        std::cout << "\t- total energy: Etot = " << Etot << std::endl;
        std::cout << std::endl;
        std::cout << "Intermediate atomic states:" << std::endl;
        std::cout << "\t- maximal bound state principal quantum number: maxNn = " << maxNn << std::endl;
        std::cout << "\t- maximal intermediate angular momentum sum (- L): nL = " << nL << std::endl;
        std::cout << "\t- integrate continuum states: " << (integrate_allowed ? "yes" : "no") << std::endl;
    }
    
    if (partial_wave)
    {
        // write all intermediate states' angular momenta pairs
        for (int ell = 0; ell <= nL; ell++)
        {
            std::cout << "\t- angular momenta [" << ell << "]: ";
            for (int Ln = ell; Ln <= ell + L + Pi; Ln++)
            {
                int ln = 2 * ell + L + Pi - Ln;
                std::cout << "(" << Ln << "," << ln << ") ";
            }
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;
    
    // final energy
    double ef = ki*ki - 1./(Ni*Ni) + 1./(Nf*Nf);
    
    // check energy
    if (ef < 0)
    {
        std::cout << "Excitation from Ni = " << Ni << " to Nf = " << Nf << " is not possible at given energy.";
        return EXIT_FAILURE;
    }
    
    // final momentum
    double kf = std::sqrt(ef);
    
    // check grid spacing
    if (partial_wave)
    {
        double min_wavelength = 2 * special::constant::pi / std::sqrt(Enmax);
        std::cout << "There are " << min_wavelength / (Rmax / N) << " grid samples per shortest wavelength." << std::endl;
        if (Rmax / N > min_wavelength / 10)
            std::cout << "Warning: Grid is not sufficiently fine!" << std::endl;
        std::cout << std::endl;
    }
    
    // write thread information
#ifdef _OPENMP
    # pragma omp parallel
    # pragma omp master
    std::cout << "Using " << omp_get_num_threads() << " OpenMP threads for the calculation." << std::endl << std::endl;
#endif
    
    Timer t;
    
    // outgoing electron partial T-matrices
    cArrays Tdir =
        partial_wave ?
        PWBA2::PartialWave_direct(grid, L, Pi, Ni, Li, ki, Nf, Lf, kf, nL, maxNn, Enmax, integrate_allowed, integrate_forbidden, verbose) :
        PWBA2::FullTMatrix_direct(grid, Ni, Li, ki, Nf, Lf, kf, maxNn, nL, Enmax, integrate_allowed, integrate_forbidden, verbose);
    
    std::cout << "Tdir = " << Tdir << std::endl;
    std::cout << "Tdir sums = " << sums(Tdir) << std::endl;
    std::cout << "Partial cross sections = " << sqrabs(sums(Tdir)) / (16 * special::constant::pi * special::constant::pi) << std::endl;
    std::cout << std::endl;
    
    std::cout << "Finished in " << t.nice_time() << std::endl;
    
    return EXIT_SUCCESS;
}
