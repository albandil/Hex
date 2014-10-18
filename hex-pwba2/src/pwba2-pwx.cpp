//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2014, Jakub Benda, Charles University in Prague                    //
//                                                                                   //
// MIT License:                                                                      //
//                                                                                   //
//  Permission is hereby granted, free of charge, to any person obtaining a          //
// copy of this software and associated documentation files (the "Software"),        //
// to deal in the Software without restriction, including without limitation         //
// the rights to use, copy, modify, merge, publish, distribute, sublicense,          //
// and/or sell copies of the Software, and to permit persons to whom the             //
// Software is furnished to do so, subject to the following conditions:              //
//                                                                                   //
//  The above copyright notice and this permission notice shall be included          //
// in all copies or substantial portions of the Software.                            //
//                                                                                   //
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS          //
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       //
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE       //
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, //
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF         //
// OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.  //
//                                                                                   //
//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "arrays.h"
#include "clenshawcurtis.h"
#include "complex.h"
#include "misc.h"
#include "pwba2.h"
#include "radial.h"
#include "vec3d.h"

cArrays PWBA2::PartialWave_direct
(
    rArray grid,
    int L, int Pi,
    int Ni, int Li, double ki,
    int Nf, int Lf, double kf,
    int nL, int maxNn, double Enmax,
    bool integrate_allowed, bool integrate_forbidden,
    bool verbose
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
        
        std::cout << std::endl << "---------- lf = " << lf << " ----------" << std::endl << std::endl;
        
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
                
                // conserve parity
                if ((L + Ln + ln) % 2 != Pi)
                {
                    std::cout << "\tSkipping ln = " << ln << " due to parity conservation." << std::endl;
                    continue;
                }
                
                std::cout << "\nli = " << li << ", Ln = " << Ln << ", ln = " << ln << std::endl << std::endl;
                
                // sum over bound states
                std::cout << "\tBound intermediate states" << std::endl;
                Complex bound_contrib = 0;
                for (int Nn = Ln + 1; Nn <= maxNn; Nn++)
                {
                    // compute energy of the intermediate projectile state
                    double en = ki*ki - 1./(Ni*Ni) + 1./(Nn*Nn);
                    
                    // check energy
                    if (en > 0)
                    {
                        double kn = std::sqrt(en);
                        
                        // integrate
                        bound_contrib += 2. * Idir_nBound_allowed
                        (
                            grid, L,
                            Nf, Lf, kf, lf,
                            Nn, Ln, kn, ln,
                            Ni, Li, ki, li
                        );
                    }
                    else
                    {
                        double kappan = std::sqrt(-en);
                        
                        // integrate
                        bound_contrib += 2. * Idir_nBound_forbidden
                        (
                            grid, L,
                            Nf, Lf, kf, lf,
                            Nn, Ln, kappan, ln,
                            Ni, Li, ki, li
                        );
                    }
                }
                
                //
                // integrate over allowed free states (Kn^2 < ki^2 - 1/Ni^2)
                //
                
                Complex allowed_contrib = 0;
                std::cout << "\n\tAllowed intermediate states" << std::endl;
                if (integrate_allowed)
                {
                    auto allowed_energy_contribution = [&](double Kn) -> Complex
                    {
                        // get momentum of the projectile
                        double kn = std::sqrt(std::abs(Etot - Kn * Kn));
                        
                        // compute the radial integral
                        return Kn == 0. ? 0. : 2. * Kn * Kn * Idir_nFree_allowed
                        (
                            grid, L,
                            Nf, Lf, kf, lf,
                            Kn, Ln, kn, ln,
                            Ni, Li, ki, li
                        );
                    };
                    
                    ClenshawCurtis<decltype(allowed_energy_contribution),Complex> CCa(allowed_energy_contribution);
                    CCa.setVerbose(verbose, "\t\tcc");
                    CCa.setEps(1e-5); // relative tolerance
                    allowed_contrib += CCa.integrate(0., std::sqrt(std::min(Enmax, Etot)));
                }
                
                //
                // integrate over forbidden free states (Kn^2 > ki^2 - 1/Ni^2)
                //
                
                double forbidden_contrib = 0;
                std::cout << "\n\tForbidden intermediate states" << std::endl;
                if (integrate_forbidden and Etot < Enmax)
                {
                    auto forbidden_energy_contribution = [&](double Kn) -> double
                    {
                        // get momentum of the projectile
                        double kappan = std::sqrt(std::abs(Etot - Kn * Kn));
                        
                        // compute the radial integral
                        return Kn == 0. ? 0. : 2. * Kn * Kn * Idir_nFree_forbidden
                        (
                            grid, L,
                            Nf, Lf, kf, lf,
                            Kn, Ln, kappan, ln,
                            Ni, Li, ki, li
                        );
                    };
                    ClenshawCurtis<decltype(forbidden_energy_contribution),double> CCf(forbidden_energy_contribution);
                    CCf.setVerbose(verbose, "\t\tcc");
                    CCf.setEps(1e-5); // relative tolerance
                    forbidden_contrib += CCf.integrate(std::sqrt(Etot), std::sqrt(Enmax));
                }
                
                // update T-matrix
                Complex contrib = bound_contrib + allowed_contrib + forbidden_contrib;
                Tdir_lf_li += contrib;
                
                // convergence check for Ln-loop (and high angular momenta)
                if (Ln > 5 and Ln != ell + L + Pi and std::abs(contrib) < 1e-8 * std::abs(Tdir_lf_li))
                {
                    std::cout << "\t\tSkipping the following values of Ln due to negligible contribution: ";
                    for (int LLn = Ln + 1; LLn <= ell + L + Pi; LLn++)
                        std::cout << LLn << ' ';
                    std::cout << std::endl;
                    break;
                }
            }
            
            Complex factor = std::pow(Complex(0.,1.),li-lf) * std::pow(4*special::constant::pi, 1.5) * std::sqrt(2*li + 1.);
            
            for (int Mi = -Li; Mi <= Li; Mi++)
            for (int Mf = -Lf; Mf <= Lf; Mf++)
            {
                double Cf = special::ClebschGordan(Lf, Mf, lf, Mi-Mf, L, Mi);
                double Ci = special::ClebschGordan(Li, Mi, li, 0, L, Mi);
                
                std::cout << "Contribution of li = " << li << ": " << factor * Cf * Ci * Tdir_lf_li << std::endl;
                
                Tdir_lf[(Mi + Li) * (2 * Lf + 1) + Mf + Lf] += factor * Cf * Ci * Tdir_lf_li;
            }
        }
        
        Tdir.push_back(Tdir_lf / (ki * kf));
    }
    
    return Tdir;
}
