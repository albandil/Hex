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
#include "misc.h"
#include "pwba2.h"
#include "radial.h"

#include "clenshawcurtis.h"
#include "gausskronrod.h"
#include "nodeintegrate.h"

cArrays PWBA2::PartialWave_direct
(
    rArray grid,
    int L, int Pi,
    int Ni, int Li, double ki,
    int Nf, int Lf, double kf,
    int nL, int maxNn, double Enmax,
    int maxlevel_allowed, int maxlevel_forbidden
)
{
    cArrays Tdir;
    double Etot = ki*ki - 1./(Ni*Ni);
    
    // compute all angular contributions to the T-matrix
    for (int lf = std::abs(L - Lf); lf <= L + Lf; lf++)
    {
        // conserve parity
        if ((L + Lf + lf) % 2 != Pi)
        {
            std::cout << "Skipping lf = " << lf << " due to parity conservation." << std::endl;
            continue;
        }
        
        std::cout << "---------- lf = " << lf << " ----------" << std::endl << std::endl;
        
        cArray Tdir_lf ((2*Li+1)*(2*Lf+1));
        
        for (int li = std::abs(L - Li); li <= L + Li; li++)
        {
            Complex Tdir_lf_li = 0;
            
            // conserve parity
            if ((L + Li + li) % 2 != Pi)
            {
                std::cout << "Skipping li = " << li << " due to parity conservation." << std::endl;
                continue;
            }
            
            for (int ell = 0; ell <= nL; ell++)
            for (int Ln = ell; Ln <= ell + L + Pi; Ln++)
            {
                int ln = 2 * ell + L + Pi - Ln;
                
                std::cout << "\nli = " << li << ", Ln = " << Ln << ", ln = " << ln << std::endl << std::endl;
                
                // sum over bound states
                std::cout << "\tBound intermediate states" << std::endl;
                for (int Nn = Ln + 1; Nn <= maxNn; Nn++)
                {
                    // compute energy of the intermediate projectile state
                    double en = ki*ki - 1./(Ni*Ni) + 1./(Nn*Nn);
                    
                    // check energy
                    if (en > 0)
                    {
                        double kn = std::sqrt(en);
                        
                        // integrate
                        Tdir_lf_li += Idir_nBound_allowed
                        (
                            grid, L,
                            Nf, Lf, kf, lf,
                            Nn, Ln, kn, ln,
                            Ni, Li, ki, li
                        ) * (-2. / kn);
                    }
                    else
                    {
                        double kappan = std::sqrt(-en);
                        
                        // integrate
                        Tdir_lf_li += Idir_nBound_forbidden
                        (
                            grid, L,
                            Nf, Lf, kf, lf,
                            Nn, Ln, kappan, ln,
                            Ni, Li, ki, li
                        ) * (-2. / kappan);
                    }
                }
                
                //
                // integrate over allowed free states (Kn^2 < ki^2 - 1/Ni^2)
                //
                
                if (maxlevel_allowed != 0)
                {
                    std::cout << "\n\tAllowed intermediate states" << std::endl;
                    auto allowed_energy_contribution = [&](double En) -> Complex
                    {
                        if (En == 0 or En == Etot)
                            return 0.;
                        
                        // get momentum of the intermediate hydrogen continuum state
                        double Kn = std::sqrt(En);
                            
                        // get momentum of the projectile
                        double kn = std::sqrt(ki*ki - 1./(Ni*Ni) - En);
                        
                        // compute the radial integral
                        return En * Idir_nFree_allowed
                        (
                            grid, L,
                            Nf, Lf, kf, lf,
                            Kn, Ln, kn, ln,
                            Ni, Li, ki, li
                        ) * (-2. / kn);
                    };
                    
                    ClenshawCurtis<decltype(allowed_energy_contribution),Complex> CCa(allowed_energy_contribution);
                    CCa.setVerbose(true, "\t\tcc");
                    CCa.setEps(1e-5); // relative tolerance
//                     CCa.setSubdiv(6); // evaluation points
//                     CCa.setStack(5);  // subdivision limit
                    Tdir_lf_li += CCa.integrate(0., std::min(Enmax, Etot));
                }
                
                //
                // integrate over forbidden free states (Kn^2 > ki^2 - 1/Ni^2)
                //
                
                if (maxlevel_forbidden != 0 and Etot < Enmax)
                {
                    std::cout << "\n\tForbidden intermediate states" << std::endl;
                    auto forbidden_energy_contribution = [&](double En) -> Complex
                    {
                        if (En == Etot)
                            return 0.;
                        
                        // get momentum of the intermediate hydrogen continuum state
                        double Kn = std::sqrt(En);
                        
                        // get momentum of the projectile
                        double kappan = std::sqrt(En - ki*ki + 1./(Ni*Ni));
                        
                        // compute the radial integral
                        return En * Idir_nFree_forbidden
                        (
                            grid, L,
                            Nf, Lf, kf, lf,
                            Kn, Ln, kappan, ln,
                            Ni, Li, ki, li
                        ) * (-2. / kappan);
                            
                    };
                    ClenshawCurtis<decltype(forbidden_energy_contribution),Complex> CCf(forbidden_energy_contribution);
                    Tdir_lf_li += CCf.integrate(Etot, Enmax);
                }
            }
            
            Complex factor = std::pow(Complex(0.,1.),li-lf) * std::pow(4*special::constant::pi, 1.5) * std::sqrt(2*li + 1.);
            
            for (int Mi = -Li; Mi <= Li; Mi++)
            for (int Mf = -Lf; Mf <= Lf; Mf++)
            {
                double Cf = ClebschGordan(Lf, Mf, lf, Mi-Mf, L, Mi);
                double Ci = ClebschGordan(Li, Mi, li, 0, L, Mi);
                
                std::cout << "Contribution of li = " << li << ": " << factor * Cf * Ci * Tdir_lf_li << std::endl;
                
                Tdir_lf[(Mi + Li) * (2 * Lf + 1) + Mf + Lf] += factor * Cf * Ci * Tdir_lf_li;
            }
        }
        
        Tdir.push_back(Tdir_lf / (ki * kf));
    }
    
    return Tdir;
}

cArrays PWBA2::FullTMatrix_direct
(
    rArray grid,
    int Ni, int Li, double ki,
    int Nf, int Lf, double kf,
    int maxNn, int maxLn, double maxEn,
    int maxlevel_allowed, int maxlevel_forbidden
)
{
    cArray Tdir;
    
    // for all bound intermediate states
    std::cout << "\tBound intermediate states" << std::endl;
    for (int Nn = 1; Nn <= maxNn; Nn++)
    for (int Ln = 0; Ln <= std::min(maxLn, Nn-1); Ln++)
    {
        // compute energy of the intermediate projectile state
        Complex en = ki*ki - 1./(Ni*Ni) + 1./(Nn*Nn);
        
        // get momentum of the projectile in the intermediate state
        Complex kn = std::sqrt(en);
        
        // TODO
    }
    
    // for all free intermediate states
    // TODO
    
    return cArrays ({ Tdir });
}
