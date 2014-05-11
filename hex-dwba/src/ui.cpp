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
#include "version.h"

std::string help_text =
    "Usage:\n"
    "\thex-dwba [options] <ni> <li> <nf> <lf> <L> <Ei>\n"
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
    
    // parse command line
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
    
    if (params.size() != 7)
    {
        std::cout << std::endl << help_text << std::endl;
        exit(0);
    }
    
    // atomic quantum numbers
    int Ni = strtol(params[0], 0, 10);
    int Li = strtol(params[1], 0, 10);
    int Nf = strtol(params[2], 0, 10);
    int Lf = strtol(params[3], 0, 10);
    int L  = strtol(params[4], 0, 10);
    
    // energy of the projectile
    double Ei = strtod(params[5], 0);    // Ry
    double ki = sqrt(Ei);
    double Ef = Ei - 1./(Ni*Ni) + 1./(Nf*Nf);
    double kf = sqrt(Ef);
    
    // grid size
    double rmax = strtod(params[6], 0);
    
    // computed partial T-matrices
    cArrays Tdir, Texc;
    
    // main computational routine
    if (distort)
        dwba (Ni, Li, ki, Nf, Lf, kf, L, Tdir, Texc, rmax, direct, exchange);
    else
        pwba (Ni, Li, ki, Nf, Lf, kf, L, Tdir, Texc, direct, exchange);
    
    std::cout << "Done." << std::endl << std::endl;
    
    //
    // write the T-matrices, sum cross sections
    //
    
    char sqlname[50];
    if (distort)
        sprintf(sqlname, "dwba-%d-%d-%d-%d-L%d-E%g.sql", Ni, Li, Nf, Lf, L, Ei);
    else
        sprintf(sqlname, "pwba-%d-%d-%d-%d-L%d-E%g.sql", Ni, Li, Nf, Lf, L, Ei);
    std::ofstream out(sqlname);
    out << "BEGIN TRANSACTION;\n";
    
    double sumsumsigma = 0.;
    std::cout << std::endl << "Cross sections for mi -> mf transitions (and sums)" << std::endl;
    std::cout << std::setw(5) << "E" << std::setw(5) << "L";
    for (int Mi = -Li; Mi <= Li; Mi++)
    {
        for (int Mf = -Lf; Mf <= Lf; Mf++)
        {
            std::cout << std::setw(15) << format("%d -> %d", Mi, Mf);
        }
        std::cout << std::setw(15) << format("%d -> Sumf", Mi);
    }   
    std::cout << std::setw(15) << "Sumi -> Sumf" << std::endl;
    std::cout << std::setw(5) << Ei << std::setw(5) << L;
    
    for (int Mi = -Li; Mi <= Li; Mi++)
    {
        double sumsigma = 0.;
        for (int Mf = -Lf; Mf <= Lf; Mf++)
        {
            int idx = (Mi + Li)*(2*Lf + 1) + Mf + Lf;
            double sigma_singlet = 0., sigma_triplet = 0.;
            
            for (int lf = std::abs(Lf - L); lf <= Lf + L; lf++)
            {
                // compute parity of the partial wave
                double parity = ((L + Lf + lf) % 2 == 0 ? 1 : -1);
                
                // assemble T-matrices for singlet and triplet
                Complex T_singlet = Tdir[idx][lf-std::abs(Lf - L)] + parity * Texc[idx][lf-std::abs(Lf - L)];
                Complex T_triplet = Tdir[idx][lf-std::abs(Lf - L)] - parity * Texc[idx][lf-std::abs(Lf - L)];
                
                // output SQL batch commands for singlet and triplet
                if (T_singlet != 0.)
                {
                    out << format
                    (
                        "INSERT OR REPLACE INTO \"tmat\" VALUES (%d,%d,%d, %d,%d,%d, %d,%d, %e, %d, %e, %e);\n",
                        Ni, Li, Mi, Nf, Lf, Mf, L, 0, Ei, lf, T_singlet.real(), T_singlet.imag()
                    );
                    
                    sigma_singlet += 0.25*sqrabs(T_singlet)/(4.*M_PI*M_PI);
                }
                if (T_triplet != 0.)
                {
                    out << format
                    (
                        "INSERT OR REPLACE INTO \"tmat\" VALUES (%d,%d,%d, %d,%d,%d, %d,%d, %e, %d, %e, %e);\n",
                        Ni, Li, Mi, Nf, Lf, Mf, L, 1, Ei, lf, T_triplet.real(), T_triplet.imag()
                    );
                    
                    sigma_triplet += 0.75*sqrabs(T_triplet)/(4.*M_PI*M_PI);
                }
            }
            
            double sigma = kf/ki * (sigma_singlet + sigma_triplet);
            std::cout << std::setw(15) << sigma;
            sumsigma += sigma;
        }
        
        std::cout << std::setw(15) << sumsigma;
        sumsumsigma += sumsigma;
    }
    std::cout << std::setw(15) << sumsumsigma << std::endl;
    
    out << "COMMIT;" << std::endl;
    out.close();
        
    return 0;
}
