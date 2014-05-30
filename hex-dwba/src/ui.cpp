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
#include <iostream>
#include <sstream>

#include <gsl/gsl_sf.h>

#include "arrays.h"
#include "dwba.h"
#include "pwba.h"
#include "special.h"
#include "version.h"

using special::constant::pi;

std::string help_text =
    "Usage:\n"
    "\thex-dwba [options] <ni> <li> <nf> <lf> <Ei> <L> [<rmax>]\n"
    "\n"
    "Available options:\n"
    "\t--help           display this help\n"
    "\t--nodistort      compute only plane wave Born approximation\n"
    "\t--nodirect       skip computation of direct T-matrix\n"
    "\t--noexchange     skip computation of exchange T-matrix\n";

int main (int argc, char *argv[])
{
    // draw package logo
    std::cout << logo_raw() << "\n";
    std::cout << "=== Distorted wave first Born approximation ===" << std::endl << std::endl;
    
    // echo command line
    std::cout << "Command line used:" << std::endl;
    std::cout << "\t";
    for (int iarg = 0; iarg < argc; iarg++)
        std::cout << argv[iarg] << " ";
    std::cout << std::endl << std::endl;
    
    // turn off GSL error jumps
    gsl_set_error_handler_off();
    
#ifndef NO_HDF
    H5::Exception::dontPrint();
#endif
    
    // disable STDOUT/STDERR buffering
    std::setvbuf(stdout, 0, _IONBF, 0);
    std::setvbuf(stderr, 0, _IONBF, 0);
    
    //
    // parse command line
    //
    
    bool distort = true, direct = true, exchange = true;
    std::vector<const char*> params;
    for (int iarg = 1; iarg < argc; iarg++)
    {
        // dashes introduce options
        if (argv[iarg][0] == '-')
        {
            // erase all leading dashes
            std::string param = argv[iarg];
            while (param.size() != 0 and param.front() == '-')
                param.erase(param.begin());
            
            // compare with known switches
            if (param == std::string("help"))
            {
                std::cout << std::endl << help_text << std::endl;
                exit (0);
            }
            else if (param == std::string("nodistort"))
            {
                std::cout << "Computing plane wave Born approximation." << std::endl;
                distort = false;
            }
            else if (param == std::string("nodirect"))
            {
                std::cout << "Not computing direct contribution." << std::endl;
                direct = false;
            }
            else if (param == std::string("noexchange"))
            {
                std::cout << "Not computing exchange contribution." << std::endl;
                exchange = false;
            }
            else
            {
                throw exception ("Unknown option \"%s\".", argv[iarg]);
            }
        }
        
        // otherwise it is a number paramater
        else
        {
            params.push_back(argv[iarg]);
        }
    }
    std::cout << std::endl;
    
    // add default rmax
    params.push_back("-1");
    
    if (params.size() < 7)
    {
        std::cout << std::endl << help_text << std::endl;
        exit(0);
    }
    
    // atomic quantum numbers
    int Ni = strtol(params[0], 0, 10);
    int Li = strtol(params[1], 0, 10);
    int Nf = strtol(params[2], 0, 10);
    int Lf = strtol(params[3], 0, 10);
    
    // energy of the projectile
    double Ei = strtod(params[4], 0);    // Ry
    double ki = std::sqrt(Ei);
    double Ef = Ei - 1./(Ni*Ni) + 1./(Nf*Nf);
    double kf = std::sqrt(Ef);
    
    // total angular momentum; can be set as range (e.g. "0-10")
    Range<int> Ls (params[5]);
    
    // grid size
    double rmax = strtod(params[6], 0);
    
    //
    // print input information
    //
    
    std::cout << "Total quantum numbers" << std::endl;
    if (Ls.first == Ls.last)
        std::cout << "\ttotal L = "   << Ls.first << std::endl << std::endl;
    else
        std::cout << "\ttotal L = "   << Ls.first << " ... " << Ls.last << std::endl << std::endl;
    std::cout << "Initial quantum numbers"  << std::endl;
    std::cout << "\thydrogen Ni = "   << Ni << std::endl;
    std::cout << "\thydrogen Li = "   << Li << std::endl;
    std::cout << "\tprojectile ki = " << ki << std::endl << std::endl;
    std::cout << "Final quantum numbers"    << std::endl;
    std::cout << "\thydrogen Nf = "   << Nf << std::endl;
    std::cout << "\thydrogen Lf = "   << Lf << std::endl;
    std::cout << "\tprojectile kf = " << kf << std::endl << std::endl;
    if (distort)
    {
        std::cout << "Numerical settings" << std::endl;
        if (rmax > 0)
            std::cout << "\tgrid rmax = " << rmax << std::endl;
        else
            std::cout << "\tgrid rmax = auto" << std::endl;
    }
    std::cout << std::endl << std::endl;
    
    //
    // main computational routine
    //
    
    std::cout << "Running the computation..." << std::endl << std::endl;
    
    // computed partial T-matrices
    std::vector<cArrays> Tdir, Texc;
    
    // for all total angular momenta
    for (int L = Ls.first; L <= Ls.last; L++)
    {
        // add space for new T-matrices
        Tdir.push_back({});
        Texc.push_back({});
        
        // compute the T-matrices
        if (distort)
            dwba (Ni, Li, ki, Nf, Lf, kf, L, Tdir.back(), Texc.back(), rmax, direct, exchange);
        else
            pwba (Ni, Li, ki, Nf, Lf, kf, L, Tdir.back(), Texc.back(), direct, exchange);
    }
    
    //
    // write the T-matrices, sum cross sections
    //
    
    // create the SQL batch file
    char sqlname[50];
    if (distort)
        std::sprintf(sqlname, "dwba-%d-%d-%d-%d-L%d-E%g.sql", Ni, Li, Nf, Lf, Ls.last, Ei);
    else
        std::sprintf(sqlname, "pwba-%d-%d-%d-%d-L%d-E%g.sql", Ni, Li, Nf, Lf, Ls.last, Ei);
    std::ofstream out(sqlname);
    out << "BEGIN TRANSACTION;\n";
    
    // write header to stdout
    std::cout << std::endl << "Cross sections for mi -> mf transitions (and sums)" << std::endl;
    std::cout << std::setw(5) << std::right << "E" << std::setw(5) << "L" << "    ";
    for (int Mi = -Li; Mi <= Li; Mi++)
    {
        for (int Mf = -Lf; Mf <= Lf; Mf++)
        {
            std::cout << std::left << std::setw(15) << format("%d -> %d", Mi, Mf);
        }
        std::cout << std::left << std::setw(15) << format("%d -> Sumf", Mi);
    }   
    std::cout << std::left << std::setw(15) << "Sumi -> Sumf" << std::endl;
    
    // for all L write T-matrices and cross sections
    for (int L = Ls.first; L <= Ls.last; L++)
    {
        std::cout << std::right << std::setw(5) << Ei << std::setw(5) << L << "    ";
        
        // summed cross section for [Ni,Li,*] -> [Nf,Lf,*]
        double sumsumsigma = 0.;
        
        // for all initial magnetic sublevels
        for (int Mi = -Li; Mi <= Li; Mi++)
        {
            // summed cross section for [Ni,Li,Mi] -> [Nf,Lf,*]
            double sumsigma = 0.;
            
            // for all final sublevels
            for (int Mf = -Lf; Mf <= Lf; Mf++)
            {
                // array index
                int idx = (Mi + Li)*(2*Lf + 1) + Mf + Lf;
                
                // cross section for [Ni,Li,Mi] -> [Nf,Lf,Mf]
                double sigma_singlet = 0., sigma_triplet = 0.;
                
                // for all outgoing partial waves
                for (int lf = std::abs(Lf - L); lf <= Lf + L; lf++)
                {
                    // compute parity of the partial wave
                    double parity = ((L + Lf + lf) % 2 == 0 ? 1 : -1);
                    
                    // assemble T-matrices for singlet and triplet
                    Complex T_singlet = Tdir[L-Ls.first][idx][lf-std::abs(Lf-L)] + parity * Texc[L-Ls.first][idx][lf-std::abs(Lf-L)];
                    Complex T_triplet = Tdir[L-Ls.first][idx][lf-std::abs(Lf-L)] - parity * Texc[L-Ls.first][idx][lf-std::abs(Lf-L)];
                    
                    // update cross sections
                    sigma_singlet += sqrabs(T_singlet);
                    sigma_triplet += sqrabs(T_triplet);
                    
                    // output SQL batch commands for singlet and triplet
                    if (T_singlet != 0.)
                    {
                        out << format
                        (
                            "INSERT OR REPLACE INTO \"tmat\" VALUES (%d,%d,%d, %d,%d,%d, %d,%d, %e, %d, %e, %e, 0, 0);\n",
                            Ni, Li, Mi, Nf, Lf, Mf, L, 0, Ei, lf, T_singlet.real(), T_singlet.imag()
                        );
                    }
                    if (T_triplet != 0.)
                    {
                        out << format
                        (
                            "INSERT OR REPLACE INTO \"tmat\" VALUES (%d,%d,%d, %d,%d,%d, %d,%d, %e, %d, %e, %e, 0, 0);\n",
                            Ni, Li, Mi, Nf, Lf, Mf, L, 1, Ei, lf, T_triplet.real(), T_triplet.imag()
                        );
                    }
                }
                
                // compute the cross section and write it to stdout
                double sigma = kf/ki * (sigma_singlet + 3*sigma_triplet) / (16*pi*pi);
                std::cout << std::setw(15) << std::left << sigma;
                
                // update the sum
                sumsigma += sigma;
            }
            
            // write the summed cross section to stdout
            std::cout << std::setw(15) << std::left << sumsigma;
            
            // update the sum of sums
            sumsumsigma += sumsigma;
        }
        
        // write the sum of sums
        std::cout << std::setw(15) << sumsumsigma << std::endl;
    }
    
    // finish writing the SQL batch file
    out << "COMMIT;" << std::endl;
    out.close();
    
    // finish the program
    std::cout << std::endl << "Done." << std::endl << std::endl;
    
    return 0;
}