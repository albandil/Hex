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

#include <complex>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include <gsl/gsl_sf.h>

#include "arrays.h"
#include "pwba.h"
#include "special.h"
#include "version.h"

std::string help_text =
    "Usage:\n"
    "\thex-pwba [options] <ni> <li> <nf> <lf> <L> <Ei>\n"
    "\n"
    "Available options:\n"
    "\t--help           display this help\n"
    "\t--nodirect       skip computation of direct T-matrix\n"
    "\t--noexchange     skip computation of exchange T-matrix\n";

int main (int argc, char* argv[])
{
    // write program logo
    std::cout << logo_raw() << std::endl;
    std::cout << "=== Plane wave first Born approximation ===" << std::endl << std::endl;
    
    // echo command line
    std::cout << "Command line used:" << std::endl;
    std::cout << "\t";
    for (int iarg = 0; iarg < argc; iarg++)
        std::cout << argv[iarg] << " ";
    std::cout << std::endl << std::endl;
    
    // disable GSL error handler
    gsl_set_error_handler_off();
    
    // disable buffering of the standard output (-> immediate logging)
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    
    // parse command line
    bool direct = true, exchange = true;
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
            else if (param == std::string("nodirect"))
            {
                direct = false;
            }
            else if (param == std::string("noexchange"))
            {
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
    
    if (params.size() != 6)
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
    
    // computed partial T-matrices
    cArrays Tdir, Texc;
    
    // main computational routine
    pwba (Ni, Li, ki, Nf, Lf, kf, L, Tdir, Texc, direct, exchange);
    std::cout << "Done." << std::endl << std::endl;
    
    //
    // write the T-matrices, sum cross sections
    //
    
    char sqlname[50];
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
                if (abs(T_singlet) != 0.)
                {
                    out << format
                    (
                        "INSERT INTO \"tmat\" VALUES (%d,%d,%d, %d,%d,%d, %d,%d, %e, %d, %e, %e);\n",
                        Ni, Li, Mi, Nf, Lf, Mf, L, 0, Ei, lf, T_singlet.real(), T_singlet.imag()
                    );
                    
                    sigma_singlet += 0.25*sqrabs(T_singlet)/(4.*M_PI*M_PI);
                }
                if (abs(T_triplet) != 0.)
                {
                    out << format
                    (
                        "INSERT INTO \"tmat\" VALUES (%d,%d,%d, %d,%d,%d, %d,%d, %e, %d, %e, %e);\n",
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
    
    return EXIT_SUCCESS;
}
