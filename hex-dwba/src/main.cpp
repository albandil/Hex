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
#include "dwba1.h"
#include "hydrogen.h"
#include "potential.h"
#include "version.h"
#include "wave_distort.h"

std::string help_text =
    "Usage:\n"
    "\thex-dwba <ni> <li> <nf> <lf> <L> <Ei> [<rmax>]\n";

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
    
    // print usage info if called in a wrong way
    if (argc != 7 and argc != 8)
    {
        std::cout << std::endl << help_text << std::endl;
        exit(0);
    }
    
    // extract parameters
    int Ni = strtol(argv[1], 0, 10);
    int Li = strtol(argv[2], 0, 10);
    int Nf = strtol(argv[3], 0, 10);
    int Lf = strtol(argv[4], 0, 10);
    int L = strtol(argv[5], 0, 10);
    double Ei = strtod(argv[6], 0);
    double ki = sqrt(Ei);
    double kf = sqrt(Ei - 1./(Ni*Ni) + 1./(Nf*Nf));
    double rmax = (argc == 8) ? strtod(argv[7], 0) : -1.;
    
    int MM = (2*Li+1)*(2*Lf+1);		// m⟶m"transition count
    
    cArrays Tdir;
    cArrays Texc;
    
    // initial and final atomic state
    HydrogenFunction psii(Ni, Li);
    HydrogenFunction psif(Nf, Lf);
    
    // distorting potentials
    DistortingPotential Uf(Nf,rmax);
    DistortingPotential Ui(Ni,rmax);
    
    for (int lf = std::abs(Lf - L); lf <= Lf + L; lf++)
    {
        cArray Tdir_lf(MM), Texc_lf(MM);
        
        std::cout << "lf = " << lf << std::endl;
        
        DistortedWave chif(kf,lf,Ui);
        
        // direct 1e
        if (Ni == Nf and Li == Lf)
        {
            Complex tmat = DWBA1::computeDirect1e(Uf,lf,ki);
            
            std::cout << "\tdirect 1e = " << tmat << std::endl;
            
            for (int Mi = -Li; Mi <= Li; Mi++)
                for (int Mf = -Lf; Mf <= Lf; Mf++)
                    Tdir_lf[(Mi+Li)*(2*Lf+1)+Mf+Lf] += tmat;
        }
        
        for (int li = std::abs(Li - L); li <= Li + L; li++)
        {
            // conserve angular momentum
            if (li < lf - Li - Lf)
                continue;
            
            // conserve angular momentum
            if (li > lf + Li + Lf)
                break;
            
            // conserve parity
            if ((li + Li) % 2 != (lf + Lf) % 2)
                continue;
            
            std::cout << "\tli = " << li << std::endl;
            
            cArray DD_lf_li(MM), DE_lf_li(MM), ED_lf_li(MM), EE_lf_li(MM);
            DistortedWave chii(ki,li,Ui);
            
            // exchange 1e
            if (Li == lf and Lf == li)
            {
                Complex tmat = DWBA1::computeExchange1e(Uf, Ni, Li, ki, Nf, Lf, kf);
                
                std::cout << "\t\texchange 1e = " << tmat << std::endl;
                
                for (int Mi = -Li; Mi <= Li; Mi++)
                    for (int Mf = -Lf; Mf <= Lf; Mf++)
                        Texc_lf[(Mi+Li)*(2*Lf+1)+Mf+Lf] += (Mf == 0) ? tmat : 0.;
            }
            
            // direct 2e
            for (int lambda = std::max(abs(Li-Lf),abs(li-lf)); lambda <= std::min(Li+Lf,li+lf); lambda++)
            {
                Complex tmat = DWBA1::computeDirect2e(Uf, lambda, Nf, Lf, kf, lf, Ni, Li, ki, li);
                
                std::cout << "\t\tdirect 2e = " << tmat << std::endl;
                
                for (int Mi = -Li; Mi <= Li; Mi++)
                {
                    for (int Mf = -Lf; Mf <= Lf; Mf++)
                    {
                        double Gaunts_dir = Gaunt(lambda, Mi - Mf, li, 0, lf, Mi - Mf) * Gaunt(lambda, Mi - Mf, Lf, Mf, Li, Mi);
                        
                        if (not std::isfinite(Gaunts_dir))
                            throw exception ("Gaunt failure!\n");
                        
                        Tdir_lf[(Mi+Li)*(2*Lf+1)+Mf+Lf] += tmat * Gaunts_dir;
                    }
                }
            }
            
            // exchange 2e
            for (int lambda = std::max(abs(Li-lf),abs(li-Lf)); lambda <= std::min(Li+lf,li+Lf); lambda++)
            {
                Complex tmat = DWBA1::computeExchange2e(Uf, lambda, Nf, Lf, kf, lf, Ni, Li, ki, li);
                
                std::cout << "\t\texchange 2e = " << tmat << std::endl;
                
                for (int Mi = -Li; Mi <= Li; Mi++)
                {
                    for (int Mf = -Lf; Mf <= Lf; Mf++)
                    {
                        double Gaunts_exc = Gaunt(lambda, -Mf, Li, Mi, lf, Mi - Mf) * Gaunt(lambda, -Mf, Lf, Mf, li, 0);
                        
                        if (not std::isfinite(Gaunts_exc))
                            throw exception ("Gaunt failure!\n");
                        
                        Texc_lf[(Mi+Li)*(2*Lf+1)+Mf+Lf] += tmat * Gaunts_exc;
                    }
                }
            }
        } // end for li
        
        // update T-matrices
        Tdir.push_back(Tdir_lf);
        Texc.push_back(Texc_lf);
    }
    std::cout << "Done." << std::endl << std::endl;
    
    // write SQL
    std::ostringstream oss;
    oss << Ni << "_" << Li << "-" << Nf << "_" << Lf << "-" << Ei << "-" << L << ".sql";
    std::ofstream fsql(oss.str().c_str());
    fsql << "BEGIN TRANSACTION;\n";
    for (int Mi = -Li; Mi <= Li; Mi++)
    for (int Mf = -Lf; Mf <= Lf; Mf++)
    for (int ell = 0; ell < (int)Tdir.size(); ell++)
    {
        Complex tdir = Tdir[ell][(Mi+Li)*(2*Lf+1)+Mf+Lf];
        Complex texc = Texc[ell][(Mi+Li)*(2*Lf+1)+Mf+Lf];
        
        for (int S = 0; S <= 1; S++)
        {
            Complex T_ell = (S == 0) ? tdir + texc : tdir - texc;
            
            fsql << "INSERT OR REPLACE INTO \"tmat\" VALUES ("
            << Ni << "," << Li << "," << Mi << ","
            << Nf << "," << Lf << "," << Mf << ","
            << L << "," << S << "," << Ei << "," << ell << ","
            << T_ell.real() << "," << T_ell.imag()
            << ");\n";
        }
    }
    fsql << "COMMIT;\n";
    fsql.close();
    
    // extract integral cross section for all transitions
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
    
    double sumsumsigma = 0;
    for (int Mi = -Li; Mi <= Li; Mi++)
    {
        double sumsigma = 0;
        for (int Mf = -Lf; Mf <= Lf; Mf++)
        {
            double sigma = 0;
            
            for (int ell = 0; ell < (int)Tdir.size(); ell++)
            {
                Complex tdir = Tdir[ell][(Mi+Li)*(2*Lf+1)+Mf+Lf];
                Complex texc = Texc[ell][(Mi+Li)*(2*Lf+1)+Mf+Lf];
                double tsinglet = abs(tdir + texc);
                double ttriplet = abs(tdir - texc);
                sigma += kf * (0.25 * tsinglet * tsinglet + 0.75 * ttriplet * ttriplet) / ki;
            }
            
            sigma /= 4 * M_PI * M_PI;
            
            std::cout << std::setw(15) << sigma;
            sumsigma += sigma;
        }
        std::cout << std::setw(15) << sumsigma;
        sumsumsigma += sumsigma;
    }
    std::cout << std::setw(15) << sumsumsigma << std::endl;
        
    return 0;
}
