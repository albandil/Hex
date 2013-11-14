/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2013                              *
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
#include "wave_distort.h"

int main (int argc, char *argv[])
{
    // 	gsl_set_error_handler_off();
    H5::Exception::dontPrint();
    
    // disable STDOUT/STDERR buffering
    std::setvbuf(stdout, 0, _IONBF, 0);
    std::setvbuf(stderr, 0, _IONBF, 0);
    
    // print usage info if called in a wrong way
    if (argc != 8)
    {
        std::cout << "\nUsage:\n";
        std::cout << "\thex-dwba <ni> <li> <nf> <lf> <Ei> <sigmaeps> <rmax>\n\n";
        exit(0);
    }
    
    // extract parameters
    int Ni = strtol(argv[1], 0, 10);
    int Li = strtol(argv[2], 0, 10);
    int Nf = strtol(argv[3], 0, 10);
    int Lf = strtol(argv[4], 0, 10);
    double Ei = strtod(argv[5], 0);
    double ki = sqrt(Ei);
    double kf = sqrt(Ei - 1./(Ni*Ni) + 1./(Nf*Nf));
    double sigmaeps = strtod(argv[6], 0);
    double rmax = strtod(argv[7], 0);
    
    int MM = (2*Li+1)*(2*Lf+1);		// m⟶m"transition count
    
    cArrays Tdir;
    cArrays Texc;
    
    // accelerators: if a per-lf contribution decays under "accelerator_eps",
    // such term will not be computed anymore for higher lf-contributions
    double accelerator_eps = 1e-8;
    bool compute_Tdir = true;
    bool compute_Texc = true;
    
    // initial and final atomic state
    HydrogenFunction psii(Ni, Li);
    HydrogenFunction psif(Nf, Lf);
    
    // distorting potentials
    DistortingPotential Uf(Nf,rmax);
    DistortingPotential Ui(Ni,rmax);
    
    for (int lf = 0; ; lf++)
    {
        cArray Tdir_lf(MM), Texc_lf(MM);
        
        std::cout << "\n--------------------------------------------------------\n";
        std::cout << "lf = " << lf << "\n";
        
        DistortedWave chif(kf,lf,Ui);
        
        // direct 1e
        if (Ni == Nf and Li == Lf and compute_Tdir)
        {
            Complex tmat = DWBA1::computeDirect1e(Uf,lf,ki);
            
            std::cout << "\tdirect 1e = " << tmat << "\n";
            
            for (int Mi = -Li; Mi <= Li; Mi++)
                for (int Mf = -Lf; Mf <= Lf; Mf++)
                    Tdir_lf[(Mi+Li)*(2*Lf+1)+Mf+Lf] += tmat;
        }
        
        for (int li = 0; ; li++)
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
            
            std::cout << "\tli = " << li << "\n";
            
            cArray DD_lf_li(MM), DE_lf_li(MM), ED_lf_li(MM), EE_lf_li(MM);
            DistortedWave chii(ki,li,Ui);
            
            // exchange 1e
            if (Li == lf and Lf == li and compute_Texc)
            {
                Complex tmat = DWBA1::computeExchange1e(Uf, Ni, Li, ki, Nf, Lf, kf);
                
                std::cout << "\t\texchange 1e = " << tmat << "\n";
                
                for (int Mi = -Li; Mi <= Li; Mi++)
                    for (int Mf = -Lf; Mf <= Lf; Mf++)
                        Texc_lf[(Mi+Li)*(2*Lf+1)+Mf+Lf] += (Mf == 0) ? tmat : 0.;
            }
            
            // direct 2e
            if (compute_Tdir)
            for (int lambda = std::max(abs(Li-Lf),abs(li-lf)); lambda <= std::min(Li+Lf,li+lf); lambda++)
            {
                Complex tmat = DWBA1::computeDirect2e(Uf, lambda, Nf, Lf, kf, lf, Ni, Li, ki, li);
                
                std::cout << "\t\tdirect 2e = " << tmat << "\n";
                
                for (int Mi = -Li; Mi <= Li; Mi++)
                {
                    for (int Mf = -Lf; Mf <= Lf; Mf++)
                    {
                        double Gaunts_dir = Gaunt(lambda, Mi - Mf, li, 0, lf, Mi - Mf) * Gaunt(lambda, Mi - Mf, Lf, Mf, Li, Mi);
                        
                        if (not finite(Gaunts_dir))
                            throw exception ("Gaunt failure!\n");
                        
                        Tdir_lf[(Mi+Li)*(2*Lf+1)+Mf+Lf] += tmat * Gaunts_dir;
                    }
                }
            }
                
            // exchange 2e
            if (compute_Texc)
            for (int lambda = std::max(abs(Li-lf),abs(li-Lf)); lambda <= std::min(Li+lf,li+Lf); lambda++)
            {
                Complex tmat = DWBA1::computeExchange2e(Uf, lambda, Nf, Lf, kf, lf, Ni, Li, ki, li);
                
                std::cout << "\t\texchange 2e = " << tmat << "\n";
                
                for (int Mi = -Li; Mi <= Li; Mi++)
                {
                    for (int Mf = -Lf; Mf <= Lf; Mf++)
                    {
                        double Gaunts_exc = Gaunt(lambda, -Mf, Li, Mi, lf, Mi - Mf) * Gaunt(lambda, -Mf, Lf, Mf, li, 0);
                        
                        if (not finite(Gaunts_exc))
                            throw exception ("Gaunt failure!\n");
                        
                        Texc_lf[(Mi+Li)*(2*Lf+1)+Mf+Lf] += tmat * Gaunts_exc;
                    }
                }
            }
        } // end for li
        
        // update T-matrices
        Tdir.push_back(Tdir_lf);
        Texc.push_back(Texc_lf);
        
        // convergence check
        rArray sigma_singlet_lf = pow(abs(Tdir_lf + Texc_lf), 2);
        rArray sigma_triplet_lf = pow(abs(Tdir_lf - Texc_lf), 2);
        rArray sigma_singlet = sums(pow(abs(Tdir + Texc), 2));
        rArray sigma_triplet = sums(pow(abs(Tdir - Texc), 2));
        rArray relcng_singlet = sigma_singlet_lf/sigma_singlet;
        rArray relcng_triplet = sigma_triplet_lf/sigma_triplet;
        std::cout << "\tδ (singlet) = " << relcng_singlet << "\n";
        std::cout << "\tδ (triplet) = " << relcng_triplet << "\n";
        if (std::max(max(relcng_singlet), max(relcng_triplet)) < sigmaeps)
            break;
        
        // update regulators ("do not unnecessarily refine epsilons")
        double Td_contrib = std::max(max(abs(Tdir_lf)/sigma_singlet), max(abs(Tdir_lf)/sigma_triplet));
        double Te_contrib = std::max(max(abs(Texc_lf)/sigma_singlet), max(abs(Texc_lf)/sigma_triplet));
        if (compute_Tdir and Td_contrib < accelerator_eps)
        {
            compute_Tdir = false;
            std::cout << "\tAbandoning Tdir part of DWBA-1\n";
        }
        if (compute_Texc and Te_contrib < accelerator_eps)
        {
            compute_Texc = false;
            std::cout << "\tAbandoning Texc part of DWBA-1\n";
        }
    }
    
    // write SQL
    std::ostringstream oss;
    oss << Ni << "_" << Li << "-" << Nf << "_" << Lf << "-" << Ei << ".sql";
    std::ofstream fsql(oss.str().c_str());
    fsql << "BEGIN TRANSACTION;\n";
    for (int Mi = -Li; Mi <= Li; Mi++)
    for (int Mf = -Lf; Mf <= Lf; Mf++)
    for (int ell = 0; ell < (int)Tdir.size(); ell++)
    {
        Complex tdir = Tdir[ell][(Mi+Li)*(2*Lf+1)+Mf+Lf];
        Complex texc = Texc[ell][(Mi+Li)*(2*Lf+1)+Mf+Lf];
        
        int L = 0; // ??? FIXME consistent with hex-ecs ???
        
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
    double sumsumsigma = 0;
    std::cout << "\n" << Ei << "\t";
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
            
            std::cout << sigma << " ";
            sumsigma += sigma;
        }
        std::cout << sumsigma << " ";
        sumsumsigma += sumsigma;
    }
    std::cout << sumsumsigma << " " << Tdir.size() << "\n";
        
    return 0;
}
