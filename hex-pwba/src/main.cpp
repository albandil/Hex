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
#include "complex.h"
#include "radial.h"
#include "misc.h"
#include "specf.h"
#include "version.h"

int main (int argc, char* argv[])
{
    // disable GSL error handler
    gsl_set_error_handler_off();
    
    if (argc != 7)
    {
        printf("\nUsage:\n\thex-pwba <ni> <li> <nf> <lf> <L> <Ei>\n\n");
        exit(0);
    }
    
    // write program logo
    std::cout << logo_raw() << "\n";
    
    // atomic quantum numbers
    int Ni = strtol(argv[1], 0, 10);
    int Li = strtol(argv[2], 0, 10);
    int Nf = strtol(argv[3], 0, 10);
    int Lf = strtol(argv[4], 0, 10);
    int L = strtol(argv[5], 0, 10);
    
    // energy of the projectile
    double Ei = strtod(argv[6], 0);    // Ry
    double ki = sqrt(Ei);
    double Ef = Ei - 1./(Ni*Ni) + 1./(Nf*Nf);
    double kf = sqrt(Ef);
    
    /* 
     * Computes all contributions to the T-matrix for the specified total angular
     * momentum L. The contributions depend on all such "li" and "lf" that
     * the following conditions must be satisfied:
     *    |li - Li| <= L <= li + Li
     *    |lf - Lf| <= L <= lf + Lf
     * This results in the following bounds on "li" and "lf":
     *    |Li - L| <= li <= Li + L
     *    |Lf - L| <= lf <= Li + L
     */
    
    // computed partial T-matrices
    cArrays Tdir((2*Li+1)*(2*Lf+1)), Texc((2*Li+1)*(2*Lf+1));
    
    // for all outgoing partial waves
    for (int lf = std::abs(Lf - L); lf <= Lf + L; lf++)
    {
        // add new T-matrix for this outgoing partial wave
        for (cArray & T : Tdir)
            T.push_back(0.);
        for (cArray & T : Texc)
            T.push_back(0.);
        
        // for all incoming partial waves
        for (int li = std::abs(Li - L); li <= Li + L; li++)
        {
            // get multipole bounds
            int minlambda = std::min
            (
                std::max(abs(Lf-Li), abs(lf-li)),   // direct T lowest contribution
                std::max(abs(lf-Li), abs(Lf-li))    // exchange T lowest contribution
            );
            int maxlambda = std::max
            (
                std::min(Li+Lf, lf+li),   // direct T highest contribution
                std::min(lf+Li, Lf+li)    // exchange T highest contribution
            );
            
            // for all multipoles
            for (int lam = minlambda; lam <= maxlambda; lam++)
            {
                // compute the needed radial integrals
                double Vdir = compute_Idir (li, lf, lam, Ni, Li, ki, Nf, Lf, kf);
                double Vexc = compute_Iexc (li, lf, lam, Ni, Li, ki, Nf, Lf, kf);
                
                // compute complex prefactor
                Complex prefactor = pow((2*M_PI),3) * 8./(ki*kf)*pow(Complex(0.,1.),li-lf)/(2.*lam+1)*sqrt((2*li+1)/(4*M_PI));
                
                // for all projections of the initial/final angular momentum
                for (int Mi = -Li; Mi <= Li; Mi++)
                for (int Mf = -Lf; Mf <= Lf; Mf++)
                {
                    // compute index in the array of T-matrices
                    int idx = (Mi + Li)*(2*Lf + 1) + Mf + Lf;
                    
                    // compute angular integrals (Gaunt coefficients)
                    double Gaunts_dir = Gaunt(Li,Mi,lam,Mf-Mi,Lf,Mf) * Gaunt(lam,Mf-Mi,lf,Mi-Mf,li,0);
                    double Gaunts_exc = Gaunt(li,0,lam,Mf,Lf,Mf) * Gaunt(lam,Mf,lf,Mi-Mf,Li,Mi);
                    
                    // add the T-matrix contributions
                    Tdir[idx][lf-std::abs(Lf - L)] += prefactor * Gaunts_dir * Vdir;
                    Texc[idx][lf-std::abs(Lf - L)] += prefactor * Gaunts_exc * Vexc;
                }
            }
        }
    }
    
    //
    // write the T-matrices
    //
    
    char sqlname[50];
    sprintf(sqlname, "pwba-%d-%d-%d-%d-E%g.sql", Ni, Li, Nf, Lf, Ei);
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
