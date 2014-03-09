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

#include <cmath>

#include <gsl/gsl_sf.h>

#include "arrays.h"
#include "complex.h"
#include "radial.h"
#include "specf.h"

void pwba
(
    int Ni, int Li, double ki,
    int Nf, int Lf, double kf,
    int L,
    cArrays & Tdir, cArrays & Texc,
    bool direct, bool exchange
)
{    
    Tdir.resize((2*Li+1)*(2*Lf+1));
    Texc.resize((2*Li+1)*(2*Lf+1));

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
                double Vdir = direct ? compute_Idir (li, lf, lam, Ni, Li, ki, Nf, Lf, kf) : 0.;
                double Vexc = exchange ? compute_Iexc (li, lf, lam, Ni, Li, ki, Nf, Lf, kf) : 0.;
                
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
                    
                    std::cout << "Gaunts_dir = " << Gaunts_dir << std::endl;
                    std::cout << "Gaunts_exc = " << Gaunts_exc << std::endl;
                    
                    // add the T-matrix contributions
                    Tdir[idx][lf-std::abs(Lf - L)] += prefactor * Gaunts_dir * Vdir;
                    Texc[idx][lf-std::abs(Lf - L)] += prefactor * Gaunts_exc * Vexc;
                }
            }
        }
    }
}
