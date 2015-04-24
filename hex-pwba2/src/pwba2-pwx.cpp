//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2015, Jakub Benda, Charles University in Prague                    //
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
#include "gausskronrod.h"
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
            
            # pragma omp parallel for schedule (dynamic,1)
            for (int i = 0; i < (nL + 1) * (L + Pi + 1); i++)
            {
                // extract helper variables
                int ell = i / (L + Pi + 1); // i.e. 0 <= ell <= nL
                int iLn = i % (L + Pi + 1); // i.e. 0 <= iLn <= L + Pi;
                
                // compute intermediate angular momenta
                int Ln = ell + iLn; // i.e. ell <= L <= ell + L + Pi
                int ln = 2 * ell + L + Pi - Ln;
                
                // create unique output stream for this thread
                std::string filename = format("pwba2-pwx-L%d-Pi%d-(%d,%d).log", L, Pi, Ln, ln);
                std::ofstream log (0);
                if (verbose) log.open(filename);
                
                // conserve parity
                if ((L + Ln + ln) % 2 != Pi)
                {
                    # pragma omp critical
                    std::cout << "Skipping ln = " << ln << " due to parity conservation." << std::endl;
                    continue;
                }
                
                # pragma omp critical
                {
                    std::cout << "li = " << li << ", (Ln,ln) = (" << Ln << "," << ln << ")" << std::endl;
                    log << "li = " << li << ", Ln = " << Ln << ", ln = " << ln << std::endl << std::endl;
                    if (verbose)
                        std::cout << "\t-> writing debug output to file \"" << filename << "\"" << std::endl;
                    std::cout << std::endl;
                }
                
                // sum over bound states
                log << "\tBound intermediate states" << std::endl;
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
                            Ni, Li, ki, li,
                            log
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
                            Ni, Li, ki, li,
                            log
                        );
                    }
                }
                
                //
                // integrate over allowed free states (0 <= Kn^2 <= ki^2 - 1/Ni^2)
                //
                
                Complex allowed_contrib = 0;
                log << std::endl << "\tAllowed intermediate states" << std::endl;
                if (integrate_allowed)
                {
                    auto allowed_energy_contribution = [&](double Kn) -> Complex
                    {
                        // boundary values
                        if (Kn == 0)
                            return 0.;
                        
                        // get momentum of the projectile
                        double kn = std::sqrt(std::abs(Etot - Kn * Kn));
                        
                        // compute the radial integral
                        return 2. * Kn * Kn * Idir_nFree_allowed
                        (
                            grid, L,
                            Nf, Lf, kf, lf,
                            Kn, Ln, kn, ln,
                            Ni, Li, ki, li,
                            log
                        );
                    };
                    
                    // use complex Clenshaw-Curtis integrator
                    ClenshawCurtis<decltype(allowed_energy_contribution),Complex> CC(allowed_energy_contribution);
                    CC.setVerbose(verbose, "\t\tcc", log);
                    CC.setEps(1e-5); // relative tolerance
                    allowed_contrib += CC.integrate(0., std::sqrt(std::min(Enmax, Etot)));
                }
                
                //
                // integrate over forbidden free states (ki^2 - 1/Ni^2 <= Kn^2 <= +Inf)
                //
                
                double forbidden_contrib = 0;
                if (verbose) log << std::endl << "\tForbidden intermediate states" << std::endl;
                if (integrate_forbidden and Etot < Enmax)
                {
                    double far = special::constant::Inf;
                    
                    auto forbidden_energy_contribution = [&](double Kn) -> double
                    {
                        // boundary values
                        if (Kn == 0 or not std::isfinite(Kn) or Kn >= far)
                            return 0.;
                        
                        // get momentum of the projectile
                        double kappan = std::sqrt(std::abs(Etot - Kn * Kn));
                        
                        // compute the radial integral
                        double val = 2. * Kn * Kn * Idir_nFree_forbidden
                        (
                            grid, L,
                            Nf, Lf, kf, lf,
                            Kn, Ln, kappan, ln,
                            Ni, Li, ki, li,
                            log
                        );
                        
                        // shift numerical threshold
                        if (val == 0)
                            far = Kn;
                        
                        return val;
                    };
                    
                    // use real compactified Gauss-Kronrod integrator
                    GaussKronrod<decltype(forbidden_energy_contribution)> GK(forbidden_energy_contribution);
                    GK.integrate(std::sqrt(Etot), special::constant::Inf);
                    GK.setEpsRel(1e-5);
                    if (not GK.ok())
                    {
                        # pragma omp critical
                        std::cerr << "Integration of forbidden contribution failed (" << GK.status() << ").";
                        std::terminate();
                    }
                    forbidden_contrib += GK.result();
                }
                
                // update T-matrix
                Complex contrib = bound_contrib + allowed_contrib + forbidden_contrib;
                # pragma omp critical
                Tdir_lf_li += contrib;
                
                // convergence check for Ln-loop (and high angular momenta)
//                 if (Ln > 5 and Ln != ell + L + Pi and std::abs(contrib) < 1e-8 * std::abs(Tdir_lf_li))
//                 {
//                     std::cout << "\t\tSkipping the following values of Ln due to negligible contribution: ";
//                     for (int LLn = Ln + 1; LLn <= ell + L + Pi; LLn++)
//                         std::cout << LLn << ' ';
//                     std::cout << std::endl;
//                     break;
//                 }
            }
            
            // calculate angular factor
            Complex factor = std::pow(Complex(0.,1.),li-lf) * std::pow(4*special::constant::pi, 1.5) * std::sqrt(2*li + 1.);
            
            // evaluate T-matrix for all transitions between available magnetic sublevels
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
