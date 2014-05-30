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

#include "arrays.h"
#include "complex.h"
#include "radial.h"
#include "special.h"

using special::constant::pi;

void pwba
(
    int Ni, int Li, double ki,
    int Nf, int Lf, double kf,
    int L,
    cArrays & Tdir, cArrays & Texc,
    bool direct, bool exchange
)
{
    // allocate memory
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
            // conserve parity
            if ((lf + Lf) % 2 != (li + Li) % 2)
                continue;
            
            //
            // compute direct contribution
            //
            
            // for all multipoles
            for (int lam = std::max(std::abs(Lf-Li), std::abs(lf-li)); direct and lam <= std::min(Li+Lf, lf+li); lam++)
            {
                // compute the needed radial integrals
                double Vdir = compute_Idir (li, lf, lam, Ni, Li, ki, Nf, Lf, kf);
                
                // compute complex prefactor
                Complex prefactor = std::pow(4*pi,2)/(ki*kf) * std::pow(Complex(0.,1.),li-lf) * std::sqrt((2*li+1)/(4*pi));
                
                // for all projections of the initial/final angular momentum
                for (int Mi = -Li; Mi <= Li; Mi++)
                for (int Mf = -Lf; Mf <= Lf; Mf++)
                {
                    // compute index in the array of T-matrices
                    int idx = (Mi + Li)*(2*Lf + 1) + Mf + Lf;
                    
                    // compute angular integrals (Gaunt coefficients)
                    double ang = ClebschGordan(Lf,Mf,lf,Mi-Mf,L,Mi) * ClebschGordan(Li,Mi,li,0,L,Mi) * computef(lam,Lf,lf,Li,li,L);
                    
                    // add the T-matrix contributions
                    Tdir[idx][lf-std::abs(Lf - L)] += prefactor * ang * Vdir;
                }
            }
            
            //
            // compute exchange contribution
            //
            
            // for all multipoles
            for (int lam = std::max(std::abs(lf-Li), std::abs(Lf-li)); exchange and lam <= std::min(lf+Li, Lf+li); lam++)
            {
                // compute the needed radial integrals
                double Vexc = compute_Iexc (li, lf, lam, Ni, Li, ki, Nf, Lf, kf);
                
                // compute complex prefactor
                Complex prefactor = std::pow(4*pi,2)/(ki*kf)*std::pow(Complex(0.,1.),li-lf)*std::sqrt((2*li+1)/(4*pi));
                
                // for all projections of the initial/final angular momentum
                for (int Mi = -Li; Mi <= Li; Mi++)
                for (int Mf = -Lf; Mf <= Lf; Mf++)
                {
                    // compute index in the array of T-matrices
                    int idx = (Mi + Li)*(2*Lf + 1) + Mf + Lf;
                    
                    // compute angular integrals (Gaunt coefficients)
                    double ang = ClebschGordan(Lf,Mf,lf,Mi-Mf,L,Mi) * ClebschGordan(Li,Mi,li,0,L,Mi) * computef(lam,lf,Lf,Li,li,L);
                    
                    // add the T-matrix contributions
                    Texc[idx][lf-std::abs(Lf - L)] += prefactor * ang * Vexc;
                }
            }
        }
    }
}
