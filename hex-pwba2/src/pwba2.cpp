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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "arrays.h"
#include "complex.h"
#include "misc.h"
#include "pwba2.h"
#include "radial.h"
#include "vec3d.h"

#include "clenshawcurtis.h"
#include "diophantine.h"
#include "gausskronrod.h"
#include "nodeintegrate.h"
#include "spgrid.h"

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
                        bound_contrib += Idir_nBound_allowed
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
                        bound_contrib += Idir_nBound_forbidden
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
                
                Complex allowed_contrib = 0;
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
                        return Idir_nFree_allowed
                        (
                            grid, L,
                            Nf, Lf, kf, lf,
                            Kn, Ln, kn, ln,
                            Ni, Li, ki, li
                        ) * (-Kn / kn);
                    };
                    
                    ClenshawCurtis<decltype(allowed_energy_contribution),Complex> CCa(allowed_energy_contribution);
                    CCa.setVerbose(true, "\t\tcc");
                    CCa.setEps(1e-5); // relative tolerance
                    allowed_contrib += CCa.integrate(0., std::min(Enmax, Etot));
                }
                
                //
                // integrate over forbidden free states (Kn^2 > ki^2 - 1/Ni^2)
                //
                
                Complex forbidden_contrib = 0;
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
                        return Idir_nFree_forbidden
                        (
                            grid, L,
                            Nf, Lf, kf, lf,
                            Kn, Ln, kappan, ln,
                            Ni, Li, ki, li
                        ) * (-Kn / kappan);
                            
                    };
                    ClenshawCurtis<decltype(forbidden_energy_contribution),Complex> CCf(forbidden_energy_contribution);
                    forbidden_contrib += CCf.integrate(Etot, Enmax);
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

inline Complex WA (double k, double nu, double q)
{
    return Complex (k*k + nu*nu - q*q, -2.*nu*q);
}

inline double WB (geom::vec3d k, double nu, geom::vec3d q)
{
    return nu*nu + (k.x+q.x)*(k.x+q.x) + (k.y+q.y)*(k.y+q.y) + (k.z+q.z)*(k.z+q.z);
}

Complex W_1s (geom::vec3d vk, geom::vec3d vq)
{
    double nu = 1.;
    double k = geom::vec3d::norm(vk);
    double q = geom::vec3d::norm(vq);
    
    Complex miq (0.,-1./q); // = -i/q
    Complex A = WA(k,nu,q);
    double B = WB(vk,nu,vq);
    
    Complex B_A = B / A;
    return 2.0 * (-miq) * std::pow(B_A,-miq) * (Complex(-1.,q) * B_A + Complex(1.,q)) / (B * B * k * k);
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
    // DEBUG output
    std::ostringstream table;
    OutputTable tab (table);
    tab.setWidth(15, 15, 30, 15, 30, 15, 30);
    tab.write("theta", "phi", "fU", "evaluations", "fW", "evaluations", "fU+fW");
    
    int Ntheta = 10, Nphi = 10;
//     for (int itheta = 0; itheta < Ntheta; itheta++)
    for (int iphi = 0; iphi < Nphi; iphi++)
    {
//         double theta = itheta * special::constant::pi / (Ntheta - 1);
        double theta = special::constant::pi_quart;
        double phi = iphi * special::constant::pi / (Nphi - 1);
        
        // total energy of the system
        double Etot = 0.5 * ki * ki - 0.5 / (Ni * Ni);
        
        // initial and final wave-vectors
        geom::vec3d vki = { 0, 0, ki };
        geom::vec3d vkf = {
            kf * std::sin(theta) * std::cos(phi),
            kf * std::sin(theta) * std::sin(phi),
            kf * std::cos(theta)
        };
        
        // maximal total (and on-shell) intermediate wave-number
        double Qmax = 10, Qon = std::sqrt(2 * Etot);
        
        // U-integrand (real part of propagator)
        auto integrand_U_wrap = [ki,kf,vki,vkf,Qmax,Qon,Etot](int n, double const * coords) -> Complex
        {
            // check dimensions
            assert(n == 6);
            
            // unpack polar coordinates
            double costheta1 = 2 * coords[0] - 1, sintheta1 = std::sqrt(1 - costheta1 * costheta1);
            double costheta2 = 2 * coords[1] - 1, sintheta2 = std::sqrt(1 - costheta2 * costheta2);
            double phi1  = special::constant::two_pi  * coords[2];
            double phi2  = special::constant::two_pi  * coords[3];
            double alpha = special::constant::pi_half * coords[4];
            double Q = Qmax * coords[5];
            
            // compute both off- and on-shell momentum magnitudes
            double q1 = Q * std::cos(alpha), q1on = Qon * std::cos(alpha);
            double q2 = Q * std::sin(alpha), q2on = Qon * std::sin(alpha);
            
            // compute both off- and on-shell momentum vectors
            geom::vec3d nq1 = { sintheta1 * std::cos(phi1), sintheta1 * std::sin(phi1), costheta1 };
            geom::vec3d nq2 = { sintheta2 * std::cos(phi2), sintheta2 * std::sin(phi2), costheta2 };
            geom::vec3d vq1 = q1 * nq1, vq1on = q1on * nq1;
            geom::vec3d vq2 = q2 * nq2, vq2on = q2on * nq2;
            
            // compute carthesian coordinates: "vql" = q<, "vqg" = q>
            double ql = std::min(q1,q2), qlon = std::min(q1on,q2on);
            geom::vec3d vql = (q1 < q2 ? vq1 : vq2), vqlon = (q1on < q2on ? vq1on : vq2on);
            geom::vec3d vqg = (q1 > q2 ? vq1 : vq2), vqgon = (q1on > q2on ? vq1on : vq2on);
            
            // the value of the off-shell integrand
            double norm = 4. / (special::constant::pi * ql * (1. - std::exp(-2.*special::constant::pi/ql)));
            Complex Wf = W_1s(vkf - vqg, vql), Wi = std::conj(W_1s(vki - vqg, vql));
            Complex integrand_U_off = Q * norm * Wf * Wi;
            
            // the value of the on-shell integrand
            double normon = 4. / (special::constant::pi * qlon * (1. - std::exp(-2.*special::constant::pi/qlon)));
            Complex Wfon = W_1s(vkf - vqgon, vqlon), Wion = std::conj(W_1s(vki - vqgon, vqlon));
            Complex integrand_U_on = Qon * normon * Wfon * Wion;
            
            // evaluate integrand
            return special::constant::two_inv_pi * Qmax * q1*q1 * q2*q2 * (integrand_U_off - integrand_U_on) / (Etot - 0.5*Q*Q);
        };
        
        // iW-integrand (imag part of propagator)
        auto integrand_W_wrap = [ki,kf,vki,vkf,Qmax,Etot](int n, double const * coords) -> Complex
        {
            // check dimensions
            assert(n == 5);
            
            // unpack polar coordinates
            double costheta1 = 2 * coords[1] - 1;
            double costheta2 = 2 * coords[0] - 1;
            double sintheta1 = std::sqrt(1 - costheta1*costheta1);
            double sintheta2 = std::sqrt(1 - costheta2*costheta2);
            double phi1 = special::constant::two_pi * coords[2];
            double phi2 = special::constant::two_pi * coords[3];
            double alpha = special::constant::pi_half * coords[4];
            double Q = std::sqrt(2 * Etot);
            double q1 = Q * std::cos(alpha);
            double q2 = Q * std::sin(alpha);
            
            // compute directions
            geom::vec3d vq1 = { q1 * sintheta1 * std::cos(phi1), q1 * sintheta1 * std::sin(phi1), q1 * costheta1 };
            geom::vec3d vq2 = { q2 * sintheta2 * std::cos(phi2), q2 * sintheta2 * std::sin(phi2), q2 * costheta2 };
            
            // compute carthesian coordinates: "vql" = q<, "vqg" = q>
            double ql = std::min(q1,q2);
            geom::vec3d vql = (q1 < q2 ? vq1 : vq2);
            geom::vec3d vqg = (q1 > q2 ? vq1 : vq2);
            
            // the value of the integrand
            double norm = 4. / (special::constant::pi * ql * (1. - std::exp(-2.*special::constant::pi/ql)));
            Complex Wf = W_1s(vkf - vqg, vql), Wi = std::conj(W_1s(vki - vqg, vql));
            Complex integrand_W = Complex(0.,-special::constant::pi) * norm * Wf * Wi;
            
            // evaluate integrand
            return special::constant::two_inv_pi * q1*q1 * q2*q2 * integrand_W;
        };
        
        // setup the sparse grid integrator
        spgrid::SparseGrid<Complex> G;
        G.setMaxLevel(6);
        G.setLocEpsAbs(1e-9);
        G.setLocEpsRel(1e-5);
        G.setGlobEpsAbs(1e-6);
        G.setGlobEpsRel(1e-6);
        G.setVerbose(true);
        
        // DEBUG output
        std::cout << underline(format("Sparse grid integration of U for θ = %g, φ = %g", theta, phi)) << std::endl;
        
        // integrate U on 6-dimensional sparse grid
        G.integrate_adapt_v2(integrand_U_wrap, Unit_6Cube, spgrid::d6l4n257, spgrid::d6l5n737);
        Complex fU = G.result();
        std::size_t nEvalU = G.evalcount();
        
        // DEBUG output
        std::cout << underline(format("Sparse grid integration of W for θ = %g, φ = %g", theta, phi)) << std::endl;
        
        // integrate W on 5-dimensional sparse grid
        G.integrate_adapt_v2(integrand_W_wrap, Unit_5Cube, spgrid::d5l4n151, spgrid::d5l5n391);
        Complex fW = G.result();
        std::size_t nEvalW = G.evalcount();
        
        // DEBUG output
        tab.write(theta, phi, fU, nEvalU, fW, nEvalW, fU+fW);
    }
    
    // write out the table with T-matrices
    std::cout << std::endl;
    std::cout << underline("Resulting T-matrices") << std::endl;
    std::cout << table.str() << std::endl;
    
    // for now, exit
    std::exit(0);
}
