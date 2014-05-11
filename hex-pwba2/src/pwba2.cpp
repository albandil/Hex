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

Complex W_1s (double kx, double ky, double kz, double px, double py, double pz)
{
    double alpha = 1;
    double k2 = kx*kx + ky*ky + kz*kz, k = std::sqrt(k2);
    double p2 = px*px + py*py + pz*pz, p = std::sqrt(p2);
    
    if (k == 0)
        return 0.;
    
    Complex i_over_k (0.,1/k);
    Complex A (k2 + alpha*alpha - p2, -2*alpha*p);
    double B = alpha*alpha + (kx-px)*(kx-px) + (ky-py)*(ky-py) + (kz-pz)*(kz-pz);
    Complex A_B = A/B;
    
    //           d   |
    // evaluate -- I |
    //          dα   |α=1
    
    Complex dA_da (2*alpha, -2*p);
    double dB_da = 2 * alpha;
    Complex dI_da = std::pow(A_B,i_over_k) * (i_over_k / (A*B) * dA_da - (i_over_k+1.) / (B*B) * dB_da);
    
    return -dI_da / special::constant::pi_sqrt;
}

Complex W_2s (double kx, double ky, double kz, double px, double py, double pz)
{
    double alpha = 0.5;
    double k2 = kx*kx + ky*ky + kz*kz, k = std::sqrt(k2);
    double p2 = px*px + py*py + pz*pz, p = std::sqrt(p2);
    
    if (k == 0)
        return 0.;
    
    Complex i_over_k (0.,1/k);
    Complex A = k2 + std::pow(Complex(alpha,-p), 2);
    Complex B = alpha*alpha + std::pow(kx-px,2) + std::pow(ky-py,2) + std::pow(kz-pz,2);
    
    //           d   |
    // evaluate -- I |
    //          dα   |α=1/2
    
    Complex dA_da (2*alpha, -2*p);
    Complex dB_da = 2*alpha;
    Complex dI_da = i_over_k * std::pow(A,i_over_k-1.) / std::pow(B,i_over_k+1.) * dA_da
                  -(i_over_k+1.) * std::pow(A,i_over_k) / std::pow(B,i_over_k+2.) * dB_da;
    
    //           d²   |
    // evaluate --  I |
    //          dα²   |α=1/2
    
    Complex d2A_da2 = 2;
    Complex d2B_da2 = 2;
    Complex d2I_da2 = i_over_k * (i_over_k-1.) * std::pow(A,i_over_k-2.) / std::pow(B,i_over_k+1.) * dA_da * dA_da
                    - 2. * i_over_k * (i_over_k+1.) * std::pow(A,i_over_k-1.) / std::pow(B,i_over_k) * dA_da * dB_da
                    + i_over_k * std::pow(A,i_over_k-1.) / std::pow(B,i_over_k-1.) * d2A_da2
                    + (i_over_k+1.) * (i_over_k+2.) * std::pow(A,i_over_k) / std::pow(B,i_over_k+3.) * dB_da * dB_da
                    - (i_over_k+1.) * std::pow(A,i_over_k) / std::pow(B,i_over_k+2.) * d2B_da2;
    
    return (-dI_da - 0.5*d2I_da2) / std::sqrt(8. * special::constant::pi);
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
    double Etot = 0.5 * ki * ki - 0.5 / (Ni * Ni);
    Complex Eplus (Etot, 1e-5);
    double En = 0.5 * Etot;
    double qn = std::sqrt(2 * En);
    
    double qx  = 0,  qy  = 0,  qz  = qn;
    double kix = ki, kiy = 0,  kiz = 0;
    double kfx = 0,  kfy = kf, kfz = 0;
    
    rArray kxgrid = linspace<double>(-2.5,2.5,101);
    rArray kygrid = linspace<double>(-2.5,2.5,101);
    rArray kzgrid = linspace<double>(-2.5,2.5,101);
    cArray Wi (kxgrid.size() * kygrid.size() * kzgrid.size());
    cArray Wf (kxgrid.size() * kygrid.size() * kzgrid.size());
    cArray Green (kxgrid.size() * kygrid.size() * kzgrid.size());
    
    for (unsigned ix = 0; ix < kxgrid.size()-1; ix++)
    for (unsigned iy = 0; iy < kygrid.size()-1; iy++)
    for (unsigned iz = 0; iz < kzgrid.size()-1; iz++)
    {
        double knx = 0.5*(kxgrid[ix]+kxgrid[ix+1]);
        double kny = 0.5*(kygrid[iy]+kygrid[iy+1]);
        double knz = 0.5*(kzgrid[iz]+kzgrid[iz+1]);
        
        Wi[ix + kxgrid.size() * (iy + kygrid.size() * iz)] = W_1s(kix-knx, kiy-kny, kiz-knz, qx, qy, qz)
            / ( std::pow(kix-knx,2) + std::pow(kiy-kny,2) + std::pow(kiz-knz,2) );
        
        Wf[ix + kxgrid.size() * (iy + kygrid.size() * iz)] = W_1s(kfx-knx, kfy-kny, kfz-knz, qx, qy, qz)
            / ( std::pow(kfx-knx,2) + std::pow(kfy-kny,2) + std::pow(kfz-knz,2) );
        
        Green[ix + kxgrid.size() * (iy + kygrid.size() * iz)] = 1.
            / (Eplus - 0.5*(qx*qx+qy*qy+qz*qz) - 0.5*(knx*knx+kny*kny+knz*knz));
    }
    
    cArray result = Green * Wi * Wf;
    for (Complex & x : result)
    {
        if (not std::isfinite(x.real()) or not std::isfinite(x.imag()))
            x = 0.;
    }
    std::ofstream vtk ("integrand.vtk");
    writeVTK_cells (vtk, result, kxgrid, kygrid, kzgrid);
    vtk.close();
    
    throw exception
    (
        "Computation of the full T-matrix element is not implemented yet."
    );
}
