//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2017, Jakub Benda, Charles University in Prague                    //
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
#include <iostream>
#include <cstdlib>
#include <cmath>

// --------------------------------------------------------------------------------- //

#ifdef __linux__
    #include <fenv.h>
#endif

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"
#include "hex-hydrogen.h"
#include "hex-symbandmatrix.h"
#include "hex-version.h"

// --------------------------------------------------------------------------------- //

#include "bspline.h"
#include "radial.h"

// --------------------------------------------------------------------------------- //

extern "C" void zgbtrf_ (int*, int*, int*, int*, Complex*, int*, int*, int*);
extern "C" void zgbtrs_ (char*, int*, int*, int*, int*, Complex*, int*, int*, Complex*, int*, int*);

// --------------------------------------------------------------------------------- //

int main (int argc, char * argv[])
{
    //
    // Program initialization
    //
    
        // display logo
        std::cout << logo(" ") << std::endl;
        std::cout << "=== Semi-classical trajectory Monte Carlo ===" << std::endl << std::endl;
        
        // echo command line
        std::cout << "Command line used" << std::endl;
        std::cout << "\t";
        for (int iarg = 0; iarg < argc; iarg++)
            std::cout << argv[iarg] << " ";
        std::cout << std::endl << std::endl;
        
        // turn off GSL exceptions
        gsl_set_error_handler_off();
        
        // disable buffering of the standard output (-> immediate logging)
        std::setvbuf(stdout, nullptr, _IONBF, 0);
        
        // get input from command line
        CommandLine cmd (argc, argv);
        
#ifdef __linux__
        if (cmd.fpe)
        {
            // abort on non-numerical values
            feenableexcept(FE_INVALID);
        }
#endif
    
    //
    // 1. Set up atomic electron radial basis.
    //
    
        //- B-spline knots
        int order = 4;
        rArray rknots = concatenate
        (
            linspace( 0.0,  0.0, order),
            linspace( 0.1, 10.0, 100),
            linspace(10.5, 100.0, 180)
        );
        rArray cknots = linspace(100., 150., 101);
        
        //- ECS basis
        Real theta = 0.63;
        Bspline bspline (order, theta, {}, rknots, cknots);
        int Nspline = bspline.Nspline();
        
        //- Radial integrals
        RadialIntegrals rad (bspline, bspline, 1);
        rad.setupOneElectronIntegrals();
        
        std::cout << "Calculating one-electron eigenstates" << std::endl << std::endl;
        
        //- Eigenstates and eigenenergies
        cArray omega (Nspline), D (Nspline);
        
        ColMatrix<Complex> CR (Nspline, Nspline), invCR (Nspline, Nspline), invsqrtS (Nspline, Nspline);
        
        std::cout << "S(0,0) = " << rad.S()(0,0) << std::endl;
        ColMatrix<Complex> S = rad.S().torow().T();
        ColMatrix<Complex> H = (0.5_z * rad.D() - rad.Mm1()).torow().T();
        
        std::cout << "\t- basis overlap matrix diagonalization" << std::endl;
        Timer timer;
        
        S.diagonalize(D, nullptr, &CR);
        CR.invert(invCR);
        
        // Now S = CR * (D * CR⁻¹)
        std::cout << "\t\t- time: " << timer.nice_time() << std::endl;
        for (std::size_t i = 0; i < std::size_t(Nspline) * std::size_t(Nspline); i++)
            invCR.data()[i] *= D[i % Nspline];
        
        // S = S - CR * invCR
        blas::gemm(-1., CR, invCR, 1., S);
        std::cout << "\t\t- residual: " << S.data().norm() << std::endl;
        
        // compute √S⁻¹
        for (std::size_t i = 0; i < std::size_t(Nspline) * std::size_t(Nspline); i++)
            invCR.data()[i] /= std::pow(D.data()[i % Nspline], 1.5);
        blas::gemm(1., CR, invCR, 0., invsqrtS);
        
        std::cout << "\t- one-electron Hamiltonian matrix diagonalization (Z = " << 1 << ", l = " << 0 << ")" << std::endl;
        timer.reset();
        
        blas::gemm(1., invsqrtS, H, 0., S);
        blas::gemm(1., S, invsqrtS, 0., H);
        
        H.diagonalize(omega, nullptr, &CR);
        CR.invert(invCR);
        
        ColMatrix<Complex> invsqrtS_Cl (Nspline, Nspline);
        blas::gemm(1., invsqrtS, CR, 0., invsqrtS_Cl);
        ColMatrix<Complex> const & P = invsqrtS_Cl;
        
        // Now Hl = ClR * D * ClR⁻¹
        std::cout << "\t\t- time: " << timer.nice_time() << std::endl;
        for (std::size_t i = 0; i < std::size_t(Nspline) * std::size_t(Nspline); i++)
            invCR.data()[i] *= omega[i % Nspline];
        
        blas::gemm(-1., CR, invCR, 1., H);
        std::cout << "\t\t- residual: " << H.data().norm() << std::endl;
        
        iArray indices (Nspline);
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), [&](int i, int j) -> bool { return omega[i].real() < omega[j].real(); });
        
        std::cout << "\t\t- eigen-energies: " << std::endl;
        for (int nr = 0; nr < Nspline; nr++)
        {
            Real E_exact = -1./(2*(nr+1)*(nr+1));
            Real rel_err = std::abs(E_exact - omega[indices[nr]].real()) / std::abs(E_exact);
            
            if (rel_err >= 0.01)
                break;
            
            std::cout << format("\t\t  E(%3s) = %7.6f Ry (err %4.3f %%)", Hydrogen::stateName(nr+1,0).c_str(), 2*omega[indices[nr]].real(), 100*rel_err) << std::endl;
        }
    
    //
    // 2. Initial atomic state.
    //
    
        int ni = 1;
        
        //- Dimension of atomic electron space to consider (1 <= N <= Nspline)
        int N = 1;
        
        //- Atomic state
        cArray psia (N);
        psia[ni - 1] = 1; // pure initial state
    
    //
    // 3. Initial projectile state.
    //
    
        Real Ei = 4; // impact energy [Ry]
        
        Real x2 = -bspline.R2();
        Real v2 = std::sqrt(Ei);
        Real r2 = std::abs(x2);
    
    //
    // 4. Time loop.
    //
    
        Real dt = 1e-6;
        
        //- Potential matrix
        SymBandMatrix<Complex> v (Nspline, order + 1);
        RowMatrix<Complex> V (N,N);
        ColMatrix<Complex> R (Nspline,N);
        
        //- Auxiliary workspaces
        cArray wrkb1 (Nspline), wrkb2 (Nspline), wrke (N);
        SymBandMatrix<Complex> M (Nspline, 2 * order + 1);
        cArray w (Nspline * (4 * order + 1));
        iArray piv (Nspline);
        
        //- Other helper variables
        int ld = 4 * order + 1, info = 0;
        char norm = 'N';
        int one = 1;
        
        std::cout << std::endl << "Starting time loop" << std::endl;
        
        std::cout << std::endl << "Initial populations: " << sqrabs(psia) << std::endl << std::endl;
        
        //- Partial integral moments (0 and -1 degree)
        Complex const * const Mi0  = rad.Mitr_L_x(0).data();
        Complex const * const Mim1 = rad.Mitr_mLm1_x(0).data();
        
        std::ofstream out ("out.txt");
        
        for (int it = 1; x2 < bspline.R2() and it < 10; it++)
        {
            Real t = dt * it;
            
            // 4a. Calculate potential matrix based on the projectile position.
            
                //- evaluate potential matrix in B-spline basis
                #pragma omp parallel for
                for (int ispline = 0; ispline < Nspline; ispline++)
                for (int jspline = ispline; jspline < Nspline and jspline <= ispline + order; jspline++)
                {
                    v(ispline,jspline) = 0;
                    
                    Complex const * const pMi0  = Mi0  + (ispline * (2 * order + 1) + jspline + order - ispline) * (order + 1);
                    Complex const * const pMim1 = Mim1 + (ispline * (2 * order + 1) + jspline + order - ispline) * (order + 1);
                    
                    for (int iknot = jspline; iknot <= ispline + order and iknot < bspline.Nreknot() - 1; iknot++)
                    {
                        if (r2 <= bspline.t(iknot).real())
                        {
                            // r1 > r2
                            v(ispline,jspline) += pMim1[iknot-ispline] / bspline.t(iknot + 1).real();
                        }
                        else if (bspline.t(iknot + 1).real() <= r2)
                        {
                            // r2 > r1
                            v(ispline,jspline) += pMi0[iknot-ispline] / r2;
                        }
                        else
                        {
                            // r2 ~ r1
                            // FIXME: This is just an approximation.
                            v(ispline,jspline) += pMim1[iknot-ispline] / bspline.t(iknot + 1).real();
                        }
                    }
                }
                std::cout << "|v| ~ " << v.data().norm() << std::endl;
                std::cout << "v00 = " << v(0,0) << std::endl;
                std::cout << " -> = " << rad.S()(0,0) / r2 << std::endl;
                
            // 4b. Update the atomic quantum state.
                
                //- free hamiltonian action, O(Nspline)
                #pragma omp parallel for
                for (int i = 0; i < N; i++)
                    psia[i] *= std::exp(-1.0_i * omega[indices[i]] * t);
                
                //- convert atomic state to the B-spline space
                #pragma omp parallel for
                for (int j = 0; j < Nspline; j++)
                {
                    wrkb1[j] = 0;
                    for (int i = 0; i < N; i++)
                        wrkb1[j] += P.col(indices[i])[j] * psia[i];
                }
                
                M.data().fill(0.);
/*
                //- Crank-Nicholson I. (forward Euler half-step)
                #pragma omp parallel for
                for (int i = 0; i < Nspline; i++)
                for (int j = i; j <= i + order and j < Nspline; j++)
                    M(i,j) = rad.S()(i,j) - 0.5_i * dt * v(i,j);
                M.dot(1., wrkb1, 0., wrkb2);
                
                //- Crank-Nicholson II. (backward Euler half-step)
                #pragma omp parallel for
                for (int i = 0; i < Nspline; i++)
                for (int j = i; j <= i + order and j < Nspline; j++)
                    M(i,j) = rad.S()(i,j) + 0.5_i * dt * v(i,j);
                M.toPaddedRows(w);
                zgbtrf_(&Nspline,&Nspline,&order,&order,w.data(),&ld,piv.data(),&info);
                zgbtrs_(&norm,&Nspline,&order,&order,&one,w.data(),&ld,piv.data(),wrkb2.data(),&Nspline,&info);
*/
/*
                // backward Euler
                #pragma omp parallel for
                for (int i = 0; i < Nspline; i++)
                for (int j = i; j <= i + order and j < Nspline; j++)
                    M(i,j) = rad.S()(i,j) - 1.0_i * dt * v(i,j);
                M.toPaddedRows(w);
                zgbtrf_(&Nspline,&Nspline,&order,&order,w.data(),&ld,piv.data(),&info);
                zgbtrs_(&norm,&Nspline,&order,&order,&one,w.data(),&ld,piv.data(),wrkb1.data(),&Nspline,&info);
*/
                //- forward Euler
                #pragma omp parallel for
                for (int i = 0; i < Nspline; i++)
                for (int j = i; j <= i + order and j < Nspline; j++)
                    M(i,j) = rad.S()(i,j) - 1.0_i * dt * v(i,j);
                M.dot(1., wrkb1, 0., wrkb2);
                
                //- apply overlap matrix
                rad.S().dot(1., wrkb2, 0., wrkb1);
                
                //- project back to eigenstate space
                #pragma omp parallel for
                for (int i = 0; i < N; i++)
                    psia[i] = (P.col(indices[i]) | wrkb1);
                
                //- free hamiltonian action, O(Nspline)
                #pragma omp parallel for
                for (int i = 0; i < N; i++)
                    psia[i] *= std::exp(+1.0_i * omega[indices[i]] * t);
                
                //- renormalize, O(Nspline)
                std::cout << std::endl;
                std::cout << "Time = " << t << std::endl;
                //std::cout << "Potential matrix: " << V1.data() << ", " << V2.data() << std::endl;
                std::cout << "Atomic state populations for first " << std::min(N,5) << " states: " << sqrabs(psia.slice(0,std::min(N,5))) << std::endl;
                std::cout << "Atomic state norm: " << psia.norm() << std::endl;
                std::cout << "Projectile kinematics: v2 = " << v2 << ", x2 = " << x2 << std::endl;
                
                out << t << '\t' << x2 << '\t' << psia.norm();
                for (int i = 0; i < std::min(N,5); i++)
                    out << '\t' << sqrabs(psia[i]);
                out << std::endl;
                
                psia /= psia.norm();
                
            // 4c. Calculate total mean charge below projectile radius.
                
                //- evaluate the mean atomic electron charge
//                 Real Qc = 0;
//                 #pragma omp parallel for reduction(+:Qc)
//                 for (int n  = 0; n  < N; n ++)
//                 for (int np = 0; np < N; np++)
//                     Qc += (V1(n,np) * psia[n] * psia[np] * std::exp(1.0_i * (omega[indices[n]] - omega[indices[np]]) * t)).real();
                
                //- get total real charge
//                 Real Qr = 1. - Qc;
//                 std::cout << "Effective charge felt by projectile: " << Qr << std::endl;
            
            // 4d. Update the projectile classical state.
            
                //- projectile acceletation
                Real a2 = 0; // -x2/(r2*r2*r2) * Qr;
                
                //- forward Euler scheme
                v2 += a2 * dt;
                x2 += v2 * dt;
                
                //- also calculate projectile distance from the origin
                r2 = std::abs(x2);
        }
    
    //
    // 5. Get resulting state populations.
    //
    
        std::cout << std::endl << "Final populations: " << sqrabs(psia) << std::endl;
    
    
    return EXIT_SUCCESS;
}
