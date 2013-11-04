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

#include <complex>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <vector>

#include <gsl/gsl_sf.h>

#include "complex.h"
#include "gausskronrod.h"
#include "misc.h"
#include "specf.h"
#include "symbolic.h"

double compute_Idir (int li, int lf, int lam, int Ni, int Li, double ki, int Nf, int Lf, double kf)
{
    if (lam == 0)
    {
        auto outer_integrand = [ = ](double x) -> double {
            
            auto inner_integrand = [ = ](double y) -> double {
                try {
                    return (1./y - 1./x) * hydro_P(Ni,Li,y) * hydro_P(Nf,Lf,y);
                } catch (...) {
                    return 0.;
                }
            };
            
            GaussKronrod<decltype(inner_integrand)> Q(inner_integrand);
            Q.integrate(x, Inf);
            
            try {
                return ric_j(li,ki*x) * ric_j(lf,kf*x) * Q.result();
            } catch (...) {
                return 0.;
            }
            
        };
        
        GaussKronrod<decltype(outer_integrand)> Q(outer_integrand);
        Q.integrate(0, Inf);
        return Q.result();
    }
    else
    {
        // construct inner integrands
        SymbolicPoly inner_integrand_1 = HydrogenP(Ni,Li) * HydrogenP(Nf,Lf);
        SymbolicPoly inner_integrand_2 = HydrogenP(Ni,Li) * HydrogenP(Nf,Lf);
        for (SymbolicTerm & term : inner_integrand_1)
            term.a += -lam - 1;
        for (SymbolicTerm & term : inner_integrand_2)
            term.a += lam;
        
        // integrate
        SymbolicPoly ii1 = integrate_inf(inner_integrand_1);
        SymbolicPoly ii2 = integrate_low(inner_integrand_2);
        
        auto outer_integrand = [ = ] (double x) -> double {
            
            // evaluate inner integrals at "x"
            double i1 = eval(ii1, x);
            double i2 = eval(ii2, x);
            
            // evaluate outer integrand at "x"
            try {
                return ric_j(li,ki*x) * ric_j(lf,kf*x) * ( pow(x,lam) * i1 + pow(x,-lam-1) * i2 );
            } catch (...) {
                return 0.;
            }
        };
        
        // outer integrate
        GaussKronrod<decltype(outer_integrand)> Q(outer_integrand);
        Q.integrate(0., Inf);
        return Q.result();
    }
}

double compute_Iexc (int li, int lf, int lam, int Ni, int Li, double ki, int Nf, int Lf, double kf)
{
    if (lam == 0)
    {
        auto outer_integrand = [ = ](double x) -> double {
            
            auto inner_integrand = [ = ](double y) -> double {
                try {
                    return (1./x - 1./y) * ric_j(li,ki*y) * hydro_P(Nf,Lf,y);
                } catch (...) {
                    return 0.;
                }
            };
            
            GaussKronrod<decltype(inner_integrand)> Q(inner_integrand);
            Q.integrate(0, x);
            
            try {
                return hydro_P(Ni,Li,x) * ric_j(lf,kf*x) * Q.result();
            } catch (...) {
                return 0.;
            }
            
        };
        
        GaussKronrod<decltype(outer_integrand)> Q(outer_integrand);
        Q.integrate(0, Inf);
        return Q.result();
    }
    else
    {
        // construct inner integrands
        SymbolicPoly inner_integrand_1 = RiccatiBessel(li,ki) * HydrogenP(Nf,Lf);
        SymbolicPoly inner_integrand_2 = HydrogenP(Ni,Li) * RiccatiBessel(lf,kf);
        
        for (SymbolicTerm& term : inner_integrand_1)
            term.a += lam;
        for (SymbolicTerm& term : inner_integrand_2)
            term.a += lam;
        
        // integrate
        SymbolicPoly ii1 = integrate_low(inner_integrand_1);
        SymbolicPoly ii2 = integrate_low(inner_integrand_2);
        
        auto outer_integrand = [ = ] (double x) -> double {
            
            // evaluate inner integrals at "x"
            double i1 = eval(ii1, x);
            double i2 = eval(ii2, x);
            
            // evaluate outer integrand at "x"
            try {
                return (ric_j(lf,kf*x) * hydro_P(Ni,Li,x) * i1 + ric_j(li,ki*x) * hydro_P(Nf,Lf,x) * i2) * pow(x,-lam-1);
            } catch (...) {
                return 0.;
            }
        };
        
        // outer integrate
        GaussKronrod<decltype(outer_integrand)> Q(outer_integrand);
        Q.integrate(0., Inf);
        return Q.result();
    }
}

std::pair<double,double> cross_section (
    double ki, int Li, int Mi,
    double kf, int Lf, int Mf, int maxL,
    std::map<unsigned long long,double> Idir,
    std::map<unsigned long long,double> Iexc
) {
    double sigma_singlet = 0., sigma_triplet = 0.;
    for (int li = 0; li <= maxL; li++)
    for (int lit = 0; lit <= maxL; lit++)
    for (int lam = abs(Li-Lf); lam <= Li+Lf; lam++)
    for (int lamt = abs(Li-Lf); lamt <= Li+Lf; lamt++)
    for (int lf = std::max(abs(li-lam),abs(lit-lamt)); lf <= std::min(li+lam,lit+lamt); lf++)
    {
        unsigned long long li_lf_lam = ((((unsigned long long)li << 16) + lf) << 16) + lam;
        unsigned long long lit_lf_lamt = ((((unsigned long long)lit << 16) + lf) << 16) + lamt;
        
        Complex pref = 4*pow(4*M_PI,3)/(ki*ki*ki*kf) * pow(std::complex<double>(0.,1.),li-lit) *
                    sqrt((2*li+1)*(2*lit+1)) / ((2.*lam+1)*(2.*lamt+1));
        
        double G_G_Idir   = Gaunt(li,0,lam,Mi-Mf,lf,Mi-Mf) * Gaunt(Lf,Mf,lam,Mi-Mf,Li,Mi) * Idir[li_lf_lam];
        double G_G_Idir_t = Gaunt(lit,0,lamt,Mi-Mf,lf,Mi-Mf) * Gaunt(Lf,Mf,lamt,Mi-Mf,Li,Mi) * Idir[lit_lf_lamt];
//         double G_G_Iexc   = Gaunt(li,0,lam,Mf,Lf,Mf) * Gaunt(lam,Mf,lf,Mi-Mf,Li,Mi) * Iexc[li_lf_lam];
//         double G_G_Iexc_t = Gaunt(lit,0,lamt,Mf,Lf,Mf) * Gaunt(lamt,Mf,lf,Mi-Mf,Li,Mi) * Iexc[lit_lf_lamt];
        
//          sigma_singlet += 0.25 * (pref * (G_G_Idir + G_G_Iexc) * (G_G_Idir_t + G_G_Iexc_t)).real();
//          sigma_triplet += 0.75 * (pref * (G_G_Idir - G_G_Iexc) * (G_G_Idir_t - G_G_Iexc_t)).real();
         sigma_singlet += ( pref * G_G_Idir * G_G_Idir_t ).real();
    }
    return std::make_pair (sigma_singlet, sigma_triplet);
}

int main(int argc, char* argv[])
{
    if (argc != 8)
    {
        printf("\nUsage:\n\thex-pwba <ni> <li> <nf> <lf> <maxL> <Ei> <sigmaeps>\n\n");
        exit(0);
    }
    
    // atomic quantum numbers
    int Ni = strtol(argv[1], 0, 10);
    int Li = strtol(argv[2], 0, 10);
    int Nf = strtol(argv[3], 0, 10);
    int Lf = strtol(argv[4], 0, 10);
    int maxl_limit = strtol(argv[5], 0, 10);
    
    // energy of the projectile
    double Ei = strtod(argv[6], 0);    // Ry
    double ki = sqrt(Ei);
    double Ef = Ei - 1./(Ni*Ni) + 1./(Nf*Nf);
    double kf = sqrt(Ef);
    
    // sigma increment threshold
    double sigma_threshold = strtod(argv[7], 0);
    
    // radial integrals
    std::map<unsigned long long,double> Idir, Iexc;

    // cross section
    double sigma_singlet = 0;
    
    // partial wave count limit
    int maxL;
    
    // for all partial wave count limits
    for (maxL = 0; maxL <= maxl_limit or maxl_limit < 0; maxL++)
    {
        // compute radial integrals
        for (int li = 0; li <= maxL; li++)
        for (int lam = abs(Li - Lf); lam <= Li + Lf; lam++)
        for (int lf = abs(li - lam); lf <= li + lam and lf <= maxL; lf++)
        {
            unsigned long long li_lf_lam = ((((unsigned long long)li << 16) + lf) << 16) + lam;
            
            if (Idir.find(li_lf_lam) == Idir.end())
                Idir[li_lf_lam] = compute_Idir(li, lf, lam, Ni, Li, ki, Nf, Lf, kf);
            
//             if (Iexc.find(li_lf_lam) == Iexc.end())
//                 Iexc[li_lf_lam] = compute_Iexc(li, lf, lam, Ni, Li, ki, Nf, Lf, kf);
        }
        
        // compute total cross section
        double new_sigma_singlet = 0;
        for (int Mi = -Li; Mi <= Li; Mi++)
        for (int Mf = -Lf; Mf <= Lf; Mf++)
            new_sigma_singlet += cross_section(ki,Li,Mi,kf,Lf,Mf,maxL,Idir,Iexc).first;
        
        // check convergence
        if (fabs(new_sigma_singlet - sigma_singlet) / new_sigma_singlet < sigma_threshold)
            break;
        else
            sigma_singlet = new_sigma_singlet;
    }
    
    // check if hard limit was enough
    if (maxL >= maxl_limit)
    {
        std::cerr << "Warning: Partial wave limit " << maxl_limit
                  << " reached for E = " << Ei << " Ry, the result may lack the desired presicion.\n";
    }
    
    //
    // compute the cross sections
    //
    
    std::cout << Ei << "\t";
    for (int Mi = -Li; Mi <= Li; Mi++)
    for (int Mf = -Lf; Mf <= Lf; Mf++)
    {
        auto tcs = cross_section(ki,Li,Mi,kf,Lf,Mf,maxL,Idir,Iexc);
        std::cout << tcs.first + tcs.second << "\t" << tcs.first << "\t" << tcs.second;
    }
    std::cout << maxL << "\n";
    
    //
    // write the T-matrices
    //
    
    char sqlname[50];
    sprintf(sqlname, "pwba-%d-%d-%d-%d-E%g.sql", Ni, Li, Nf, Lf, Ei);
    std::ofstream out(sqlname);
    out << "BEGIN TRANSACTION;\n";
    
    for (int Mi = -Li; Mi <= Li; Mi++)
    for (int Mf = -Lf; Mf <= Lf; Mf++)
    for (int lf = 0; lf < maxL; lf++)
    {
        std::complex<double> Tmat_dir = 0.;
        std::complex<double> Tmat_exc = 0.;
        
        for (int li = 0; li < maxL; li++)
        {
            int min_lambda = std::min(
                std::max(abs(Lf-Li), abs(lf-li)),    // direct T lowest contribution
                std::max(abs(lf-Li), abs(Lf-li))    // exchange T lowest contribution
            );
            int max_lambda = std::max(
                std::min(Li+Lf, lf+li),    // direct T highest contribution
                std::min(lf+Li, Lf+li)    // exchange T highest contribution
            );
            for (int lam = min_lambda; lam <= max_lambda; lam++)
            {
                std::complex<double> pref = pow((2*M_PI),3) * 8./(ki*kf)*pow(std::complex<double>(0.,1.),li-lf)/(2.*lam+1)*sqrt((2*li+1)/(4*M_PI));
                double Gaunts_dir = Gaunt(Li,Mi,lam,Mf-Mi,Lf,Mf) * Gaunt(lam,Mf-Mi,lf,Mi-Mf,li,0);
                double Gaunts_exc = Gaunt(li,0,lam,Mf,Lf,Mf) * Gaunt(lam,Mf,lf,Mi-Mf,Li,Mi);
                unsigned long long idx = ((((unsigned long long)li << 16) + lf) << 16) + lam;
                Tmat_dir += pref * Gaunts_dir * Idir[idx];
                Tmat_exc += pref * Gaunts_exc * Iexc[idx];
            }
        }
        
        std::complex<double> T_singlet = Tmat_dir + Tmat_exc;
        std::complex<double> T_triplet = Tmat_dir - Tmat_exc;
        
        if (abs(T_singlet) != 0.)
        {
            out << format("INSERT INTO \"tmat\" VALUES (%d,%d,%d, %d,%d,%d, %d,%d, %e, %d, %e, %e);\n",
                Ni, Li, Mi, Nf, Lf, Mf, 0, 0, Ei, lf, T_singlet.real(), T_singlet.imag());
        }
        if (abs(T_triplet) != 0.)
        {
            out << format("INSERT INTO \"tmat\" VALUES (%d,%d,%d, %d,%d,%d, %d,%d, %e, %d, %e, %e);\n",
                Ni, Li, Mi, Nf, Lf, Mf, 0, 1, Ei, lf, T_triplet.real(), T_triplet.imag());
        }
    }
    out << "COMMIT;\n";
    out.close();
    
    //
    // end
    //
    
    return EXIT_SUCCESS;
}
