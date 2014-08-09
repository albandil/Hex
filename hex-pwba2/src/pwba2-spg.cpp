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
#include <sstream>

#include "arrays.h"
#include "complex.h"
#include "misc.h"
#include "pwba2.h"
#include "vec3d.h"
#include "special.h"
#include "spgrid.h"

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

double Wb_1s (geom::vec3d vk, int Nn, int Ln)
{
    assert(Nn == 1);
    assert(Ln == 0);
    
    double k = geom::vec3d::norm(vk);
    
    return -(k*k + 8) / ((k*k + 4)*(k*k + 4));
}

cArrays PWBA2::FullTMatrix_direct
(
    rArray grid,
    int Ni, int Li, double ki,
    int Nf, int Lf, double kf,
    int maxNn, int maxLn, double maxEn,
    bool integrate_allowed, bool integrate_forbidden
)
{
    // output text table
    std::ostringstream table;
    OutputTable tab (table);
    tab.setWidth(15, 15, 30, 15, 30, 15, 30, 15, 30, 15, 30);
    tab.write("theta", "phi", "fUb", "evaluations", "fWb", "evaluations", "fU", "evaluations", "fW", "evaluations", "sum");
    
    int Ntheta = 1, Nphi = 10;
//     for (int itheta = 0; itheta < Ntheta; itheta++)
//     for (int iphi = 0; iphi < Nphi; iphi++)
    {
//         double theta = itheta * special::constant::pi / (Ntheta - 1);
        double theta = special::constant::pi_quart;
//         double theta = 0;
//         double phi = iphi * special::constant::pi / (Nphi - 1);
        double phi = 0;
        
        // total energy of the system
        double Etot = 0.5 * ki * ki - 0.5 / (Ni * Ni);
        
        // initial and final wave-vectors
        geom::vec3d vki = { 0, 0, ki };
        geom::vec3d vkf = {
            kf * std::sin(theta) * std::cos(phi),
            kf * std::sin(theta) * std::sin(phi),
            kf * std::cos(theta)
        };
        
        // on-shell intermediate wave-number
        double Qon = std::sqrt(2 * Etot);
        
        // the sparse grid integrator
        spgrid::SparseGrid<Complex> G;
        
        // criteria to stop integration cell currently being processed
        G.setLocEpsRel(1e-5);  // ... if the rules are equal up to 1e-5 (relative)
        G.setLocEpsAbs(1e-9);  // ... if the rules are equal up to 1e-9 (absolute)
        G.setGlobEpsRel(1e-6); // ... if the extrapolated relative contribution to the whole domain is within 1e-6
        G.setGlobEpsAbs(0);    // ... if the extrapolated contribution to the whole domain is within ... [not used]
        G.setVerbose(true);    // print detailed progress information
        G.setPrefix("  ");     // output formatting prefix
        
        // marching integration criteria
        double marching_epsrel = 1e-5;
        double marching_epsabs = 1e-9;
        
        // results (bound/continuum intermediate states, real/imag parts of propagator)
        Complex fUb = 0, fWb = 0, fU = 0, fW = 0;
        
        // evaluation counts for reference and debugging
        std::size_t nEvalUb = 0, nEvalWb = 0, nEvalU = 0, nEvalW = 0;
        
        // for all intermediate bound state contributions
        for (int Nn = 1; Nn <= maxNn; Nn++)
        for (int Ln = 0; Ln < Nn and Ln <= maxLn; Ln++)
        {
            // allow at most 70 integration cell bisections
            G.setMaxLevel(70);
            
            // marching integration bounds
            double Qmin = 0, Qmax = Qon;
            
            // marching index
            int step = 0;
            
            // contribution to fUb from the current state
            Complex fUb_contrib = 0;
            
            std::cout << underline(format("Sparse grid marching integration of bound U for θ = %g, φ = %g", theta, phi)) << std::endl;
            
            // start marching integration along Q
            do
            {
                // real bound integrand
                auto integrand_Ub_wrap = [ki,kf,vki,vkf,Nn,Ln,Qmin,Qmax,Etot](int n, double const * coords) -> Complex
                {
                    // check dimensions
                    assert(n == 3);
                    
                    // unpack polar coordinates
                    double costheta = 2 * coords[0] - 1;
                    double sintheta = std::sqrt(1 - costheta * costheta);
                    double phi = special::constant::two_pi  * coords[1];
                    double kn = Qmin + (Qmax - Qmin) * coords[2];
                    
                    // compute on-shell momentum magnitude
                    double Qn = std::sqrt(2*Etot + 1./(Nn*Nn));
                    
                    // compute both off- and on-shell momentum vectors
                    geom::vec3d vn = { sintheta * std::cos(phi), sintheta * std::sin(phi), costheta };
                    geom::vec3d vkn = kn * vn;
                    geom::vec3d vQn = Qn * vn;
                    
                    // the value of the off-shell integrand
                    Complex Wf = Wb_1s(vkf - vkn,Nn,Ln), Wi = std::conj(Wb_1s(vki - vkn,Nn,Ln));
                    Complex integrand_Ub_off = kn * kn * Wf * Wi;
                    
                    // the value of the on-shell integrand
                    Complex Wfon = Wb_1s(vkf - vQn,Nn,Ln), Wion = std::conj(Wb_1s(vki - vQn,Nn,Ln));
                    Complex integrand_Ub_on = Qn * Qn * Wfon * Wion;
                    
                    // Jacobian
                    double Jac = 2. // for cos theta
                            * special::constant::two_pi // for phi
                            * (Qmax - Qmin); // for Q
                    
                    // evaluate integrand
                    return special::constant::two_inv_pi * Jac * (integrand_Ub_off - integrand_Ub_on) / (0.5*Qn*Qn - 0.5*kn*kn);
                };
                
                std::cout << std::endl << "Marching: Step " << ++step << " from Q = " << Qmin << " to " << Qmax << std::endl;
                
                // integrate U on 6-dimensional sparse grid
//                 G.setWriteVTK(true, format("spgrid-U-%g-%g.vtk", Qmin, Qmax));
                G.integrate_adapt(integrand_Ub_wrap, Unit_3Cube, spgrid::d3l4n39, spgrid::d3l5n87);
                fUb_contrib += G.result();
                nEvalUb += G.evalcount();
                
                std::cout << std::endl << "  Relative change = " << std::abs(G.result()) / std::abs(fUb_contrib) << std::endl;
                
                Qmin = Qmax;
                Qmax *= 2;
            }
            while
            (
                std::abs(G.result()) > marching_epsrel * std::abs(fUb_contrib) and
                std::abs(G.result()) > marching_epsabs
            );
            
            fUb += fUb_contrib;
            std::cout << std::endl;
            
            // imag bound integrand
            auto integrand_Wb_wrap = [ki,kf,vki,vkf,Nn,Ln,Etot](int n, double const * coords) -> Complex
            {
                // check dimensions
                assert(n == 2);
                
                // unpack polar coordinates
                double costheta = 2 * coords[0] - 1;
                double sintheta = std::sqrt(1 - costheta * costheta);
                double phi = special::constant::two_pi  * coords[1];
                
                // compute on-shell momentum magnitude
                double Qn = std::sqrt(2*Etot + 1./(Nn*Nn));
                
                // compute both off- and on-shell momentum vectors
                geom::vec3d vn = { sintheta * std::cos(phi), sintheta * std::sin(phi), costheta };
                geom::vec3d vQn = Qn * vn;
                
                // the value of the on-shell integrand
                Complex Wf = Wb_1s(vkf - vQn,Nn,Ln), Wi = std::conj(Wb_1s(vki - vQn,Nn,Ln));
                Complex integrand_Wb = Complex(0.,-special::constant::pi) * Qn * Wf * Wi;
                
                // Jacobian
                double Jac = 2. // for cos theta
                           * special::constant::two_pi; // for phi
                
                // evaluate integrand
                return special::constant::two_inv_pi * Jac * integrand_Wb;
            };
            
            std::cout << underline(format("Sparse grid integration of bound W for θ = %g, φ = %g", theta, phi)) << std::endl;
            
            // integrate W on 5-dimensional sparse grid
//             G.setWriteVTK(true, "spgrid-W.vtk");
            G.integrate_adapt(integrand_Wb_wrap, Unit_2Cube, spgrid::d2l4n17, spgrid::d2l5n33);
//             G.setWriteVTK(false);
            fWb += G.result();
            nEvalWb += G.evalcount();
            
            std::cout << std::endl;
        }
        
        // shall we integrate continuum intermediate state contributions?
        if (integrate_allowed)
        {
            // allow at most 7 integration cell bisections
            G.setMaxLevel(7);
            G.setParallel(true);
            G.setGlobEpsAbs(0);
            
            std::cout << underline(format("Sparse grid integration of continuum U for θ = %g, φ = %g", theta, phi)) << std::endl;
            
            // marching integration bounds
            double Qmin = 0, Qmax = Qon;
            
            // marching index
            int step = 0;
            
            // start marching integration
            do
            {
                // U-integrand (real part of propagator)
                auto integrand_U_wrap = [ki,kf,vki,vkf,Qmin,Qmax,Qon,Etot](int n, double const * coords) -> Complex
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
                    
                    // Jacobian
                    double Jac = 2. // for costheta1
                            * 2. // for costheta2
                            * special::constant::two_pi // for phi1
                            * special::constant::two_pi // for phi2
                            * special::constant::pi_half // for alpha
                            * Qmax; // for Q
                    
                    // evaluate integrand
                    return special::constant::two_inv_pi * Jac * q1*q1 * q2*q2 * (integrand_U_off - integrand_U_on) / (Etot - 0.5*Q*Q);
                };
                
                std::cout << std::endl << "Marching: Step " << ++step << " from Q = " << Qmin << " to " << Qmax << std::endl;
                
                // integrate U on 6-dimensional sparse grid
                G.integrate_adapt(integrand_U_wrap, Unit_6Cube, spgrid::d6l4n257, spgrid::d6l5n737);
                fU += G.result();
                nEvalU += G.evalcount();
                
                std::cout << std::endl << "  Relative change = " << std::abs(G.result()) / std::abs(fU) << std::endl;
                
                Qmin = Qmax;
                Qmax *= 2;
                
                // do not integrate cells that are incomparable with the total marching estimate
                G.setGlobEpsAbs(1e-6 * std::abs(fU));
            }
            while
            (
                std::abs(G.result()) > marching_epsrel * std::abs(fU) and
                std::abs(G.result()) > marching_epsabs
            );
            
            std::cout << std::endl;
            G.setGlobEpsAbs(0);
            
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
                
                // Jacobian
                double Jac = 2. // for costheta1
                           * 2. // for costheta2
                           * special::constant::two_pi // for phi1
                           * special::constant::two_pi // for phi2
                           * special::constant::pi_half // for alpha
                           * Qmax; // for Q
                
                // evaluate integrand
                return special::constant::two_inv_pi * Jac * q1*q1 * q2*q2 * integrand_W;
            };
            
            std::cout << underline(format("Sparse grid integration of continuum W for θ = %g, φ = %g", theta, phi)) << std::endl;
            
            // integrate W on 5-dimensional sparse grid
            G.integrate_adapt(integrand_W_wrap, Unit_5Cube, spgrid::d5l4n151, spgrid::d5l5n391);
            fW = G.result();
            nEvalW = G.evalcount();
            
            std::cout << std::endl;
        }
        
        // write new row to the output table
        tab.write(theta, phi, fUb, nEvalUb, fWb, nEvalWb, fU, nEvalU, fW, nEvalW, fUb+fWb+fU+fW);
    }
    
    // write out the table with T-matrices
    std::cout << std::endl;
    std::cout << underline("Resulting T-matrices") << std::endl;
    std::cout << table.str() << std::endl;
    
    // for now, exit
    std::exit(0);
}
