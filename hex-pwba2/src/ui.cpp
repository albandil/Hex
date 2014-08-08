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

#include <clocale>
#include <cmath>
#include <iostream>

#include <gsl/gsl_errno.h>

#include "arrays.h"
#include "cmdline.h"
#include "pwba2.h"
#include "special.h"
#include "version.h"

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

int main (int argc, char* argv[])
{
    // set proper locale
    std::setlocale(LC_ALL, "en_GB.utf8");
    
    // print program logo
    std::cout << logo(" ") << std::endl;
    std::cout << "=== Plane wave second Born approximation ===" << std::endl << std::endl;
    
    // disable fatal GSL errors
    gsl_set_error_handler_off();
    
    // input file
    std::ifstream inputfile;
    
    // parse command line
    bool partial_wave = false;
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
        "help", "h", 0, [](std::string optarg) -> bool
            {
                // print usage information
                std::cout << "\n"
                    "Available switches (short forms in parentheses):                                                                  \n"
                    "                                                                                                                  \n"
                    "\t--example                 (-e)  create sample input file                                                        \n"
                    "\t--help                    (-h)  display this help                                                               \n"
                    "\t--input <filename>        (-i)  use custom input file                                                           \n"
                    "\t--partial-wave            (-w)  compute only contribution of single partial wave                                \n"
                    "                                                                                                                  \n"
                ;
                exit(0);
            },
        "input", "i", 1, [&](std::string optarg) -> bool
            {
                // set custom input file
                inputfile.open(optarg);
                if (not inputfile.good())
                    throw exception ("Error: Input file \"%s\" not found.\n", optarg.c_str());
                std::cout << "Using input file \"" << optarg << "\"." << std::endl << std::endl;
                return true;
            },
        "partial-wave", "w", 0, [&](std::string optarg) -> bool
            {
                // compute only contribution of a single partial wave
                partial_wave = true;
                return true;
            },
        
        [](std::string opt, std::string optarg) -> bool
            {
                throw exception
                (
                    "Unknown option \"%s\" with argument \"%s\".",
                    opt.c_str(), optarg.c_str()
                );
            }
    );
    
    // check input file
    if (not inputfile.is_open())
    {
        const char * filename = "pwba2.inp";
        std::cout << "Using input file \"" << filename << "\"." << std::endl << std::endl;
        inputfile.open(filename);
        if (not inputfile.good())
            throw exception ("Error: Input file \"%s\" not found.\n", filename);
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
    std::cout << "Quantum state parameters:" << std::endl;
    std::cout << "\t- total angular momentum: L = " << L << (partial_wave ? "" : " (not used)") << std::endl;
    std::cout << "\t- total parity: Π = " << Pi << (partial_wave ? "" : " (not used)") << std::endl;
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
    std::cout << "\t- maximal energy: Enmax = " << Enmax << std::endl;
    std::cout << "\t- integrate allowed states: " << (integrate_allowed ? "yes" : "no") << std::endl;
    std::cout << "\t- integrate forbidden states: " << (integrate_forbidden ? "yes" : "no") << std::endl;
    
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
        return 1;
    }
    
    // final momentum
    double kf = std::sqrt(ef);
    
    // check grid spacing
    double min_wavelength = 2 * special::constant::pi / std::sqrt(Enmax);
    std::cout << "There are " << min_wavelength / (Rmax / N) << " grid samples per shortest wavelength." << std::endl;
    if (Rmax / N > min_wavelength / 10)
        std::cout << "Warning: Grid is not sufficiently fine!" << std::endl;
    std::cout << std::endl;
    
    // outgoing electron partial T-matrices
    cArrays Tdir =
        partial_wave ?
        PWBA2::PartialWave_direct(grid, L, Pi, Ni, Li, ki, Nf, Lf, kf, nL, maxNn, Enmax, integrate_allowed, integrate_forbidden) :
        PWBA2::FullTMatrix_direct(grid, Ni, Li, ki, Nf, Lf, kf, maxNn, nL, Enmax, integrate_allowed, integrate_forbidden);
    
    std::cout << "Tdir = " << Tdir << std::endl;
    std::cout << "Tdir sums = " << sums(Tdir) << std::endl;
    std::cout << std::endl;
    
    return 0;
}
