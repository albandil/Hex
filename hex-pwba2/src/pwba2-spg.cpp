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
#include <map>
#include <sstream>

#include "arrays.h"
#include "born.h"
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
    
    Complex iq (0.,1./q); // = i/q
    Complex A = WA(k,nu,q);
    double B = WB(vk,nu,vq);
    
    Complex B_A = B / A;
    return special::constant::two_inv_sqrt_pi * iq * std::pow(B_A,iq) * (Complex(-1.,q) * B_A + Complex(1.,q)) / (B * B * k * k);
}

Complex W_cont (int n, int l, int m, geom::vec3d vk, geom::vec3d vq)
{
    if (n == 1 and l == 0 and m == 0)
        return W_1s(vk,vq);
    
    throw exception ("Matrix element ⟨χ⁻(q)|exp(-ik·r)|%d,%d,%d⟩ not implemented.", n, l, m);
}

cArrays PWBA2::FullTMatrix_direct
(
    rArray grid,
    int Ni, int Li, double ki,
    int Nf, int Lf, double kf,
    int maxNn, int maxLn, double maxEn,
    bool integrate_allowed, bool integrate_forbidden,
    bool verbose
)
{
    // TEMPORARY
    int Mi = 0;
    int Mf = 0;
    
    // output text table
    std::ostringstream table;
    OutputTable tab (table);
    tab.setWidth(15, 15, 30, 15, 30, 15, 30, 15, 30, 15, 30);
    tab.write("theta", "phi", "fUb", "evaluations", "fWb", "evaluations", "fU", "evaluations", "fW", "evaluations", "sum");
    
//     for (int itheta = 0; itheta <= 180; itheta += 45)
    int itheta = 0;
    {
        double theta = itheta * special::constant::pi / 180.;
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
        G.setMinLevel(1);       // ... do at least one subdivision
        G.setLocEpsRel(1e-4);   // ... if the rules are equal up to 1e-4 (relative)
        G.setLocEpsAbs(1e-8);   // ... if the rules are equal up to 1e-8 (absolute)
        G.setGlobEpsRel(1e-6);  // ... if the extrapolated relative contribution to the whole domain is within 1e-6
        G.setGlobEpsAbs(0);     // ... if the extrapolated contribution to the whole domain is within ... [not used]
        G.setVerbose(verbose);  // print detailed progress information
        G.setParallel(true);    // evaluate integration domains in parallel
        G.setPrefix("    ");    // output formatting prefix
        
        // marching integration tolerance
        double marchingEpsRel = 1e-4;   // ... stop if curr. section contributes less than 1e-4 of the curr. estimate
        double marchingEpsAbs = 0;      // ... allow any absolute contribution of a marching section
        
        // results (bound/continuum intermediate states, real/imag parts of propagator)
        Complex fUb = 0, fWb = 0, fU = 0, fW = 0;
        
        // evaluation counts for reference and debugging
        std::size_t nEvalUb = 0, nEvalWb = 0, nEvalU = 0, nEvalW = 0;
        
        // for all intermediate bound state contributions
        for (int Ln = 0; Ln <= maxLn; Ln++)
        for (int Mn = -Ln; Mn <= Ln; Mn++)
        for (int Nn = Ln + 1; Nn <= maxNn; Nn++)
        {
            Complex fUb_contrib = 0, fWb_contrib = 0;
            
            std::cout << underline(format("Intermediate state (%d,%d,%d)",Nn,Ln,Mn)) << std::endl << std::endl;
            
            // create polynomial expansion of the integrals
            auto Wbf = Wb_symb_in(Nf,Lf,Mf,Nn,Ln,Mn);
            auto Wbi = Wb_symb_in(Ni,Li,Mi,Nn,Ln,Mn);
            
            // marching integration bounds
            double Q1 = 0, Q2 = Qon;
            
            // real bound integrand
            auto integrand_Ub_wrap_lin = [ki,kf,vki,vkf,Ni,Nn,Nf,Ln,Qon,Etot,Wbf,Wbi,&Q1,&Q2](int npt, int dim, double const * origin, double range, double const * scale, Complex * eval) -> void
            {
                // check dimensions
                assert(dim == 3);
                
                // for all evaluation points
                for (int ipt = 0; ipt < npt; ipt++)
                {
                    double coords[3];
                    for (int icoo = 0; icoo < 3; icoo++)
                        coords[icoo] = origin[icoo] + range * scale[ipt * dim + icoo];
                    
                    // unpack polar coordinates
                    double costheta = 2 * coords[0] - 1;
                    double sintheta = std::sqrt(1 - costheta * costheta);
                    double phi = special::constant::two_pi  * coords[1];
                    double kn = Q1 + (Q2 - Q1) * coords[2];
                    
                    // compute on-shell momentum magnitude
                    double Qn = std::sqrt(2*Etot + 1./(Nn*Nn));
                    
                    // compute both off- and on-shell momentum vectors
                    geom::vec3d vn = { sintheta * std::cos(phi), sintheta * std::sin(phi), costheta };
                    geom::vec3d vkn = kn * vn;
                    geom::vec3d vQn = Qn * vn;
                    
                    // the value of the off-shell integrand
                    Complex Wf = eval_Wb(Wbf, 1./Nn + 1./Nf, vkf - vkn);
                    Complex Wi = eval_Wb(Wbi, 1./Nn + 1./Ni, vki - vkn);
                    Complex integrand_Ub_off = kn * kn * Wf * std::conj(Wi);
                    
                    // the value of the on-shell integrand
                    Complex Wfon = eval_Wb(Wbf, 1./Nn + 1./Nf, vkf - vQn);
                    Complex Wion = eval_Wb(Wbi, 1./Nn + 1./Ni, vki - vQn);
                    Complex integrand_Ub_on = Qn * Qn * Wfon * std::conj(Wion);
                    
                    // Jacobian
                    double Jac = 2. // for cos theta
                               * special::constant::two_pi // for phi
                               * (Q2 - Q1); // for Q
                    
                    // evaluate integrand
                    eval[ipt] = special::constant::two_inv_pi * Jac * (integrand_Ub_off - integrand_Ub_on) / (0.5*Qn*Qn - 0.5*kn*kn);
                }
            };
            
            std::cout << format("  Sparse grid marching integration of bound U for theta = %g, phi = %g", theta, phi) << std::endl;
            
            // allow at most 70 integration cell bisections
            G.setMaxLevel(70);
            
            // this state contribution
            Complex fUb_state_contrib = 0.;
            
            // marching integration
            do
            {
                // integrate U on 3-dimensional sparse grid
                std::cout << std::endl<< "  Linear integrand for Q = " << Q1 << " .. " << Q2 << std::endl;
                static int i = 0;
                G.setWriteVTK(true, format("%02d-%g-%g.vtk", i++, Q1, Q2));
                G.integrate_adapt<3>(integrand_Ub_wrap_lin, spgrid::d3l4n39, spgrid::d3l6n135);
                G.setWriteVTK(false);
                std::cout << std::endl;
                
                // update total variables
                fUb_state_contrib += G.result();
                nEvalUb += G.evalcount();
                
                // new Q1 = old Q2
                // new Q2 = old Q2 + 2 * (old Q2 - old Q1)
                double Qdiff = Q2 - Q1;
                Q1 = Q2;
                Q2 = Q2 + 2 * Qdiff;
            }
            while
            (
                std::abs(G.result()) > marchingEpsAbs or
                std::abs(G.result()) > marchingEpsRel * std::abs(fUb_state_contrib)
            );
            
            // update bound state contribution
            fUb_contrib = fUb_state_contrib;
            
            // imag bound integrand
            auto integrand_Wb_wrap = [ki,kf,vki,vkf,Ni,Nn,Nf,Ln,Qon,Etot,Wbf,Wbi](int npt, int dim, double const * origin, double range, double const * scale, Complex * eval) -> void
            {
                // check dimensions
                assert(dim == 2);
                
                // for all evaluation points
                for (int ipt = 0; ipt < npt; ipt++)
                {
                    double coords[2];
                    for (int icoo = 0; icoo < 2; icoo++)
                        coords[icoo] = origin[icoo] + range * scale[ipt * dim + icoo];
                    
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
                    Complex Wf = eval_Wb(Wbf, 1./Nn + 1./Nf, vkf - vQn);
                    Complex Wi = eval_Wb(Wbi, 1./Nn + 1./Ni, vki - vQn);
                    Complex integrand_Wb = Complex(0.,-special::constant::pi) * Qn * Wf * std::conj(Wi);
                    
                    // Jacobian
                    double Jac = 2. // for cos theta
                               * special::constant::two_pi; // for phi
                    
                    // evaluate integrand
                    eval[ipt] = special::constant::two_inv_pi * Jac * integrand_Wb;
                }
            };
            
            // integrate W on 2-dimensional sparse grid
            std::cout << format("  Sparse grid integration of bound W for theta = %g, phi = %g", theta, phi) << std::endl;
            G.integrate_adapt<2>(integrand_Wb_wrap, spgrid::d2l4n17, spgrid::d2l7n65);
            std::cout << std::endl;
            
            // update contributions
            fWb_contrib = G.result();
            nEvalWb += G.evalcount();
            
            // check convergence with respect to Nn
            fUb += fUb_contrib;
            fWb += fWb_contrib;
            if (std::abs(fUb_contrib) < 1e-6 * std::abs(fUb) and std::abs(fWb_contrib) < 1e-6 * std::abs(fWb))
                break;
        }
        
        // shall we integrate continuum intermediate state contributions?
        if (integrate_allowed)
        {
            std::cout << underline(format("Intermediate state CONTINUUM")) << std::endl << std::endl;
            
            auto Wcf = W_symb_in(Nf,Lf,Mf);
            auto Wci = W_symb_in(Ni,Li,Mi);
            
            // marching integration bounds
            double Q1 = 0., Q2 = Qon;
            
            // U-integrand (real part of propagator)
            auto integrand_U_wrap_lin = [Ni,Nf,ki,kf,vki,vkf,Qon,Etot,Wcf,Wci,&Q1,&Q2](int npt, int dim, double const * origin, double range, double const * scale, Complex * eval) -> void
            {
                // check dimensions
                assert(dim == 6);
                
                // for all evaluation points
                for (int ipt = 0; ipt < npt; ipt++)
                {
                    double coords[6];
                    for (int icoo = 0; icoo < 6; icoo++)
                        coords[icoo] = origin[icoo] + range * scale[ipt * dim + icoo];
                    
                    // unpack polar coordinates
                    double costheta1 = 2. * coords[0] - 1., sintheta1 = std::sqrt(1. - costheta1 * costheta1);
                    double costheta2 = 2. * coords[1] - 1., sintheta2 = std::sqrt(1. - costheta2 * costheta2);
                    double phi1  = special::constant::two_pi  * coords[2];
                    double phi2  = special::constant::two_pi  * coords[3];
                    double alpha = special::constant::pi_half * coords[4];
                    double Q = Q1 + (Q2 - Q1) * coords[5];
                    
                    // compute both off- and on-shell momentum magnitudes
                    double qn = Q * std::sin(alpha), qnon = Qon * std::sin(alpha);
                    double kn = Q * std::cos(alpha), knon = Qon * std::cos(alpha);
                    
                    // compute both off- and on-shell momentum vectors
                    geom::vec3d n1 = { sintheta1 * std::cos(phi1), sintheta1 * std::sin(phi1), costheta1 };
                    geom::vec3d n2 = { sintheta2 * std::cos(phi2), sintheta2 * std::sin(phi2), costheta2 };
                    geom::vec3d vqn = qn * n1, vqnon = qnon * n1;
                    geom::vec3d vkn = kn * n2, vknon = knon * n2;
                    
                    // the value of the off-shell integrand
                    double norm = 4. / (qn * (1. - std::exp(-special::constant::two_pi/qn)));
                    Complex Wf = eval_W(Wcf,1./Nf,vkf-vkn,vqn), Wi = std::conj(eval_W(Wci,1./Ni,vki-vkn,vqn));
                    Complex integrand_U_off = qn * qn * kn * kn * Q * norm * Wf * Wi;
                    
                    // the value of the on-shell integrand
                    double normon = 4. / (qnon * (1. - std::exp(-special::constant::two_pi/qnon)));
                    Complex Wfon = eval_W(Wcf,1./Nf,vkf-vknon,vqnon), Wion = std::conj(eval_W(Wci,1./Ni,vki-vknon,vqnon));
                    Complex integrand_U_on = qnon * qnon * knon * knon * Qon * normon * Wfon * Wion;
                    
                    // Jacobian
                    double Jac = 2. // for costheta1
                               * 2. // for costheta2
                               * special::constant::two_pi // for phi1
                               * special::constant::two_pi // for phi2
                               * special::constant::pi_half // for alpha
                               * (Q2 - Q1); // for Q
                    
                    // evaluate integrand
                    eval[ipt] = special::constant::two_inv_pi * Jac * (integrand_U_off - integrand_U_on) / (Etot - 0.5*Q*Q);
                }
            };
            
            // allow at most 10 integration cell bisections
            G.setMaxLevel(10);
            G.setGlobEpsAbs(0);
            
            // contribution from continuum U
            Complex fU_contrib = 0;
            
            // marching integration
            std::cout << format("  Sparse grid integration of continuum U for theta = %g, phi = %g", theta, phi) << std::endl;
            do
            {
                // integrate U on 6-dimensional sparse grid
                std::cout << std::endl << "  Linear integrand for Q = " << Q1 << " .. " << Q2 << std::endl;
                G.integrate_adapt<6>(integrand_U_wrap_lin, spgrid::d6l4n257, spgrid::d6l5n737);
                std::cout << std::endl;
                
                // update variables
                fU_contrib += G.result();
                nEvalU += G.evalcount();
                
                // new Q1 = old Q2
                // new Q2 = old Q2 + 2 * (old Q2 - old Q1)
                double Qdiff = Q2 - Q1;
                Q1 = Q2;
                Q2 = Q2 + 2 * Qdiff;
            }
            while
            (
                std::abs(G.result()) > marchingEpsAbs or
                std::abs(G.result()) > marchingEpsRel * std::abs(fU_contrib)
            );
            
            // update real contribution
            fU += fU_contrib;
            
            // iW-integrand (imag part of propagator)
            auto integrand_W_wrap = [Ni,Nf,ki,kf,vki,vkf,Etot,Qon,Wcf,Wci](int npt, int dim, double const * origin, double range, double const * scale, Complex * eval) -> void
            {
                // check dimensions
                assert(dim == 5);
                
                // for all evaluation points
                for (int ipt = 0; ipt < npt; ipt++)
                {
                    double coords[5];
                    for (int icoo = 0; icoo < 5; icoo++)
                        coords[icoo] = origin[icoo] + range * scale[ipt * dim + icoo];
                    
                    // unpack polar coordinates
                    double costheta1 = 2. * coords[0] - 1., sintheta1 = std::sqrt(1. - costheta1 * costheta1);
                    double costheta2 = 2. * coords[1] - 1., sintheta2 = std::sqrt(1. - costheta2 * costheta2);
                    double phi1  = special::constant::two_pi  * coords[2];
                    double phi2  = special::constant::two_pi  * coords[3];
                    double alpha = special::constant::pi_half * coords[4];
                    
                    // compute both off- and on-shell momentum magnitudes
                    double qnon = Qon * std::sin(alpha);
                    double knon = Qon * std::cos(alpha);
                    
                    // compute directions
                    geom::vec3d n1 = { sintheta1 * std::cos(phi1), sintheta1 * std::sin(phi1), costheta1 };
                    geom::vec3d n2 = { sintheta2 * std::cos(phi2), sintheta2 * std::sin(phi2), costheta2 };
                    geom::vec3d vqnon = qnon * n1;
                    geom::vec3d vknon = knon * n2;
                    
                    // the value of the integrand
                    double normon = 4. / (qnon * (1. - std::exp(-special::constant::two_pi/qnon)));
                    Complex Wfon = eval_W(Wcf,1./Nf,vkf-vknon,vqnon), Wion = std::conj(eval_W(Wci,1./Ni,vki-vknon,vqnon));
                    Complex integrand_W_on = qnon * qnon * knon * knon * normon * Wfon * Wion;
                    
                    // Jacobian
                    double Jac = 2. // for costheta1
                               * 2. // for costheta2
                               * special::constant::two_pi // for phi1
                               * special::constant::two_pi // for phi2
                               * special::constant::pi_half; // for alpha
                    
                    // evaluate integrand
                    eval[ipt] = Complex(0.,-2) * Jac * integrand_W_on;
                }
            };
            
            // integrate W on 5-dimensional sparse grid
            std::cout << format("  Sparse grid integration of continuum W for theta = %g, phi = %g", theta, phi) << std::endl;
            G.integrate_adapt<5>(integrand_W_wrap, spgrid::d5l4n151, spgrid::d5l5n391);
            std::cout << std::endl;
            
            // update contributions
            fW = G.result();
            nEvalW = G.evalcount();
        }
        
        // write new row to the output table
        tab.write(itheta, phi, fUb, nEvalUb, fWb, nEvalWb, fU, nEvalU, fW, nEvalW, fUb+fWb+fU+fW);
    }
    
    // write out the table with T-matrices
    std::cout << std::endl;
    std::cout << underline("Resulting T-matrices") << std::endl;
    std::cout << table.str() << std::endl;
    
    // for now, exit
    std::exit(EXIT_SUCCESS);
}
