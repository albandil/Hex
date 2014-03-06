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
#include "specf.h"
#include "version.h"

int main (int argc, char* argv[])
{
    // disable GSL error handler
    gsl_set_error_handler_off();
    
    // disable buffering of the standard output (-> immediate logging)
    setvbuf(stdout, nullptr, _IONBF, 0);
    
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
            if (param == std::string("nodirect"))
                direct = false;
            else if (param == std::string("noexchange"))
                exchange = false;
            else
                throw exception ("Unknown option \"%s\".", argv[iarg]);
        }
        
        // otherwise it is a number paramater
        else
        {
            params.push_back(argv[iarg]);
        }
    }
    
    if (params.size() != 6)
    {
        printf("\nUsage:\n\thex-pwba [--nodirect] [--noexchange] <ni> <li> <nf> <lf> <L> <Ei>\n\n");
        exit(0);
    }
    
    // write program logo
    std::cout << logo_raw() << "\n";
    
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
    
    //
    // write the T-matrices, sum cross sections
    //
    
    char sqlname[50];
    sprintf(sqlname, "pwba-%d-%d-%d-%d-L%d-E%g.sql", Ni, Li, Nf, Lf, L, Ei);
    std::ofstream out(sqlname);
    out << "BEGIN TRANSACTION;\n";
    
    double sumsumsigma = 0.;
    std::cout << "# E, L, sigma ..." << std::endl;
    std::cout << Ei << " " << L << " ";
    
    for (int Mi = -Li; Mi <= Li; Mi++)
    {
        double sumsigma = 0.;
        for (int Mf = -Lf; Mf <= Lf; Mf++)
        {
            int idx = (Mi + Li)*(2*Lf + 1) + Mf + Lf;
            double sigma_singlet = 0., sigma_triplet = 0.;
            
            for (int lf = std::abs(Lf - L); lf <= Lf + L; lf++)
            {
                Complex T_singlet = Tdir[idx][lf-std::abs(Lf - L)] + Texc[idx][lf-std::abs(Lf - L)];
                Complex T_triplet = Tdir[idx][lf-std::abs(Lf - L)] - Texc[idx][lf-std::abs(Lf - L)];
                
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
            
            double sigma = sigma_singlet + sigma_triplet;
            std::cout << sigma << " ";
            sumsigma += sigma;
        }
        
        std::cout << sumsigma << " ";
        sumsumsigma += sumsigma;
    }
    std::cout << sumsumsigma << std::endl;
    
    out << "COMMIT;\n";
    out.close();
    
    return EXIT_SUCCESS;
}
