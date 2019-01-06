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

#include "hex-arrays.h"
#include "hex-gausskronrod.h"
#include "hex-misc.h"
#include "hex-vec3d.h"

#include "pwba2.h"
#include "radial.h"

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

        // partial T-matrices
        cArray Tdir_lf ((2*Li+1)*(2*Lf+1));

        // for all initial partial waves
        for (int li = std::abs(L - Li); li <= L + Li; li++)
        {
            // contribution of this initial partial wave to the partial T-matrix
            Complex Tdir_lf_li = 0;

            // conserve parity
            if ((L + Li + li) % 2 != Pi)
            {
                std::cout << "Skipping li = " << li << " due to parity conservation." << std::endl;
                continue;
            }

            // contributions of individual intermediate angular states
            std::vector<std::tuple<Complex,Complex,double>> contributions ((nL + 1) * (L + Pi + 1));

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
                        bound_contrib += 2. * Idir_nBound_allowed(grid, L, Nf, Lf, kf, lf, Nn, Ln, kn, ln, Ni, Li, ki, li, log);
                    }
                    else
                    {
                        double kappan = std::sqrt(-en);
                        bound_contrib += 2. * Idir_nBound_forbidden(grid, L, Nf, Lf, kf, lf, Nn, Ln, kappan, ln, Ni, Li, ki, li, log);
                    }
                }

                //
                // integrate over allowed free states (0 <= Kn^2 <= ki^2 - 1/Ni^2)
                //

                Complex allowed_contrib = 0;
                log << std::endl << "\tAllowed intermediate states" << std::endl;
                if (integrate_allowed and Etot > 0)
                {
                    if (verbose) log << "\t\t(real)" << std::endl;

                    auto allowed_energy_contribution_re = [&](double Kn) -> double
                    {
                        double kn = std::sqrt(std::abs(Etot - Kn * Kn));
                        if (Kn == 0.) return 0.;
                        double inte = 2. * Kn * Kn * Idir_nFree_allowed_re(grid, L, Nf, Lf, kf, lf, Kn, Ln, kn, ln, Ni, Li, ki, li, log);
                        return std::isfinite(inte) ? inte : 0.;
                    };

                    GaussKronrod<decltype(allowed_energy_contribution_re)> GKre (allowed_energy_contribution_re);
                    GKre.setEpsAbs(1e-8);
                    GKre.setEpsRel(1e-5);
                    GKre.integrate(0, std::sqrt(std::min(Enmax, Etot)));

                    if (not GKre.ok())
                    {
                        # pragma omp critical
                        std::cerr << "Integration of real allowed contribution failed (" << GKre.status() << ").";
                        std::terminate();
                    }

                    if (verbose) log << "\t\t(imag)" << std::endl;

                    auto allowed_energy_contribution_im = [&](double Kn) -> double
                    {
                        double kn = std::sqrt(std::abs(Etot - Kn * Kn));
                        if (Kn == 0.) return 0.;
                        double inte = 2. * Kn * Kn * Idir_nFree_allowed_im(grid, L, Nf, Lf, kf, lf, Kn, Ln, kn, ln, Ni, Li, ki, li, log);
                        return std::isfinite(inte) ? inte : 0.;
                    };

                    GaussKronrod<decltype(allowed_energy_contribution_im)> GKim (allowed_energy_contribution_im);
                    GKim.setEpsAbs(1e-8);
                    GKim.setEpsRel(1e-5);
                    GKim.integrate(0, std::sqrt(std::min(Enmax, Etot)));

                    if (not GKim.ok())
                    {
                        # pragma omp critical
                        std::cerr << "Integration of imag allowed contribution failed (" << GKim.status() << ").";
                        std::terminate();
                    }

                    allowed_contrib += Complex(GKre.result(), GKim.result());
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
                        double val = 2. * Kn * Kn * Idir_nFree_forbidden(grid, L, Nf, Lf, kf, lf, Kn, Ln, kappan, ln, Ni, Li, ki, li, log);

                        // numerical threshold to avoid torturing the integrator
                        if (std::abs(val) < 1e-10)
                            far = Kn;

                        return std::isfinite(val) ? val : 0;
                    };

                    GaussKronrod<decltype(forbidden_energy_contribution)> GK (forbidden_energy_contribution);
                    GK.setEpsAbs(1e-8);
                    GK.setEpsRel(1e-5);
                    GK.integrate(std::sqrt(std::max(0., Etot)), std::sqrt(Enmax));
                    if (not GK.ok())
                    {
                        # pragma omp critical
                        std::cerr << "Integration of forbidden contribution failed (" << GK.status() << ").";
                        std::terminate();
                    }
                    forbidden_contrib += GK.result();
                }

                contributions[i] = std::make_tuple(bound_contrib, allowed_contrib, forbidden_contrib);
            }

            // print contributions of individual angular intermediate states
            OutputTable table (std::cout);
            table.setWidth(5, 5, 30, 30, 30);
            table.write("Ln", "ln", "bound", "free - allowed", "free - forbidden");
            table.write("--", "--", "-----", "--------------", "----------------");
            for (int i = 0; i < (nL + 1) * (L + Pi + 1); i++)
            {
                // get data
                Complex bound_contrib = std::get<0>(contributions[i]);
                Complex allowed_contrib = std::get<1>(contributions[i]);
                double forbidden_contrib = std::get<2>(contributions[i]);

                // extract helper variables
                int ell = i / (L + Pi + 1); // i.e. 0 <= ell <= nL
                int iLn = i % (L + Pi + 1); // i.e. 0 <= iLn <= L + Pi;

                // compute intermediate angular momenta
                int Ln = ell + iLn; // i.e. ell <= L <= ell + L + Pi
                int ln = 2 * ell + L + Pi - Ln;

                // print data
                table.write(Ln, ln, bound_contrib, allowed_contrib, forbidden_contrib);

                // update T-matrix
                Tdir_lf_li += bound_contrib + allowed_contrib + forbidden_contrib;
            }
            std::cout << std::endl;

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

        // add this partial T-matrix to storage container
        Tdir.push_back(Tdir_lf / (ki * kf));
    }

    return Tdir;
}
