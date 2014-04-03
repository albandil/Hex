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

#include "diophantine.h"
#include "complex.h"
#include "gausskronrod.h"
#include "misc.h"
#include "multipot.h"
#include "specf.h"

Complex Idir_nBound
(
    int lambdaf, int lambdai,
    int Nf, int Lf, int lf, double kf,
    int Nn, int Ln, int ln, double kn,
    int Ni, int Li, int li, double ki
)
{
    // construct potentials
    MultipolePotential Vfn (lambdaf, Nf, Lf, Nn, Ln);
    MultipolePotential Vni (lambdai, Nn, Ln, Ni, Li);
    
//     auto Vfn = [&](double r2) -> double
//     {
//         if (lambdaf == 0)
//         {
//             auto integrand = [&](double r1) -> double { return hydro_P(Nf,Lf,r1) * (1./r1 - 1./r2) * hydro_P(Nn,Ln,r1) ; };
//             GaussKronrod<decltype(integrand)> Q(integrand);
//             Q.integrate(r2, Inf);
//             
//             return Q.result();
//         }
//         else
//         {
//             auto integrand1 = [&](double r1) -> double { return hydro_P(Nf,Lf,r1) * hydro_P(Nn,Ln,r1) * std::pow(r1/r2,lambdaf); };
//             GaussKronrod<decltype(integrand1)> Q1(integrand1);
//             Q1.integrate(0., r2);
//             
//             auto integrand2 = [&](double r1) -> double { return hydro_P(Nf,Lf,r1) * hydro_P(Nn,Ln,r1) * std::pow(r2/r1,lambdaf+1); };
//             GaussKronrod<decltype(integrand2)> Q2(integrand2);
//             Q2.integrate(r2, Inf);
//             
//             return (Q1.result() + Q2.result()) / r2;
//         }
//     };
//     
//     auto Vni = [&](double r2) -> double
//     {
//         if (lambdaf == 0)
//         {
//             auto integrand = [&](double r1) -> double { return hydro_P(Nn,Ln,r1) * (1./r1 - 1./r2) * hydro_P(Ni,Li,r1) ; };
//             GaussKronrod<decltype(integrand)> Q(integrand);
//             Q.integrate(r2, Inf);
//             
//             return Q.result();
//         }
//         else
//         {
//             auto integrand1 = [&](double r1) -> double { return hydro_P(Nn,Ln,r1) * hydro_P(Ni,Li,r1) * std::pow(r1/r2,lambdaf); };
//             GaussKronrod<decltype(integrand1)> Q1(integrand1);
//             Q1.integrate(0., r2);
//             
//             auto integrand2 = [&](double r1) -> double { return hydro_P(Nn,Ln,r1) * hydro_P(Ni,Li,r1) * std::pow(r2/r1,lambdaf+1); };
//             GaussKronrod<decltype(integrand2)> Q2(integrand2);
//             Q2.integrate(r2, Inf);
//             
//             return (Q1.result() + Q2.result()) / r2;
//         }
//     };
    
    // setup the integral
    auto integrand_re = [&](double r2) -> double
    {
        // inner integrand 1
        // r2 < r2p
        auto iintegrand1 = [&](double r2p) -> double
        {
            return Vni(r2p) * ric_j(li,ki*r2p) * (-ric_n(ln,kn*r2p));
        };
        BesselNodeIntegrator<decltype(iintegrand1),GaussKronrod<decltype(iintegrand1)>> R1(iintegrand1, ki, li);
        R1.setEpsAbs(1e-15);
        R1.setEpsRel(1e-8);
        R1.integrate(r2, Inf);
        if (not R1.ok())
        {
            throw exception
            (
                "Can't integrate 1st inner integral.\nError: %s\n", R1.status().c_str()
            );
        }
        
        // inner integrand 2
        // r2 > r2p
        auto iintegrand2 = [&](double r2p) -> double
        {
            return Vni(r2p) * ric_j(li,ki*r2p) * ric_j(ln,kn*r2p);
        };
        BesselNodeIntegrator<decltype(iintegrand2),GaussKronrod<decltype(iintegrand2)>> R2(iintegrand2, ki, li);
        R2.setEpsAbs(1e-15);
        R2.setEpsRel(1e-8);
        R2.integrate(0, r2);
        if (not R2.ok())
        {
            throw exception
            (
                "Can't integrate 2nd inner integral.\nError: %s\n", R2.status().c_str()
            );
        }
        
        // evaluate outer integrand
        return Vfn(r2) * ric_j(lf,kf*r2) * (ric_j(ln,kn*r2) * R1.result() - ric_n(ln,kn*r2) * R2.result());
    };
    auto integrand_im = [&](double r2) -> double
    {
        // inner integrand 1
        // r2 < r2p
        auto iintegrand1 = [&](double r2p) -> double
        {
            return Vni(r2p) * ric_j(li,ki*r2p) * ric_j(ln,kn*r2p);
        };
        BesselNodeIntegrator<decltype(iintegrand1),GaussKronrod<decltype(iintegrand1)>> R1(iintegrand1, ki, li);
        R1.setEpsAbs(1e-15);
        R1.setEpsRel(1e-8);
        R1.integrate(r2, Inf);
        if (not R1.ok())
        {
            throw exception
            (
                "Can't integrate 1st inner integral.\nError: %s\n", R1.status().c_str()
            );
        }
        
        // inner integrand 2
        // r2 > r2p
        auto iintegrand2 = [&](double r2p) -> double
        {
            return Vni(r2p) * ric_j(li,ki*r2p) * ric_j(ln,kn*r2p);
        };
        BesselNodeIntegrator<decltype(iintegrand2),GaussKronrod<decltype(iintegrand2)>> R2(iintegrand2, ki, li);
        R2.setEpsAbs(1e-15);
        R2.setEpsRel(1e-8);
        R2.integrate(0, r2);
        if (not R2.ok())
        {
            throw exception
            (
                "Can't integrate 2nd inner integral.\nError: %s\n", R2.status().c_str()
            );
        }
        
        // evaluate outer integrand
        return Vfn(r2) * ric_j(lf,kf*r2) * (ric_j(ln,kn*r2) * R1.result() + ric_j(ln,kn*r2) * R2.result());
    };
    
    // integrate (outer)
    BesselNodeIntegrator<decltype(integrand_re),GaussKronrod<decltype(integrand_re)>> Re(integrand_re, kf, lf);
    Re.setEpsAbs(1e-14);
    Re.setEpsRel(1e-8);
    Re.integrate(0, Inf);
    if (not Re.ok())
    {
        throw exception
        (
            "Can't integrate real outer integral.\nError: %s\n", Re.status().c_str()
        );
    }
    BesselNodeIntegrator<decltype(integrand_im),GaussKronrod<decltype(integrand_im)>> Im(integrand_im, kf, lf);
    Im.setEpsAbs(1e-14);
    Im.setEpsRel(1e-8);
    Im.integrate(0, Inf);
    if (not Im.ok())
    {
        throw exception
        (
            "Can't integrate imag outer integral.\nError: %s\n", Re.status().c_str()
        );
    }
    
    return Complex (Re.result(),Im.result());
}

Complex Idir_nFree
(
    int lambdaf, int lambdai,
    int Nf, int Lf, int lf, double kf,
    double Kn, int Ln, int ln, double kn,
    int Ni, int Li, int li, double ki
)
{
    // construct potentials
//     MultipolePotential Vfn (lambdaf, Nf, Lf, Kn, Ln);
//     MultipolePotential Vni (lambdai, Kn, Ln, Ni, Li);
    
    auto Vfn = [&](double r2) -> double
    {
        if (lambdaf == 0)
        {
            auto integrand = [&](double r1) -> double { return hydro_P(Nf,Lf,r1) * (1./r1 - 1./r2) * hydro_F(Kn,Ln,r1) ; };
            GaussKronrod<decltype(integrand)> Q(integrand);
            Q.integrate(r2, Inf);
            
            return Q.result();
        }
        else
        {
            auto integrand1 = [&](double r1) -> double { return hydro_P(Nf,Lf,r1) * hydro_F(Kn,Ln,r1) * std::pow(r1/r2,lambdaf); };
            GaussKronrod<decltype(integrand1)> Q1(integrand1);
            Q1.integrate(0., r2);
            
            auto integrand2 = [&](double r1) -> double { return hydro_P(Nf,Lf,r1) * hydro_F(Kn,Ln,r1) * std::pow(r2/r1,lambdaf+1); };
            GaussKronrod<decltype(integrand2)> Q2(integrand2);
            Q2.integrate(r2, Inf);
            
            return (Q1.result() + Q2.result()) / r2;
        }
    };
    
    auto Vni = [&](double r2) -> double
    {
        if (lambdaf == 0)
        {
            auto integrand = [&](double r1) -> double { return hydro_F(Kn,Ln,r1) * (1./r1 - 1./r2) * hydro_P(Ni,Li,r1) ; };
            GaussKronrod<decltype(integrand)> Q(integrand);
            Q.integrate(r2, Inf);
            
            return Q.result();
        }
        else
        {
            auto integrand1 = [&](double r1) -> double { return hydro_F(Kn,Ln,r1) * hydro_P(Ni,Li,r1) * std::pow(r1/r2,lambdaf); };
            GaussKronrod<decltype(integrand1)> Q1(integrand1);
            Q1.integrate(0., r2);
            
            auto integrand2 = [&](double r1) -> double { return hydro_F(Kn,Ln,r1) * hydro_P(Ni,Li,r1) * std::pow(r2/r1,lambdaf+1); };
            GaussKronrod<decltype(integrand2)> Q2(integrand2);
            Q2.integrate(r2, Inf);
            
            return (Q1.result() + Q2.result()) / r2;
        }
    };
    
    // setup the integral
    auto integrand_re = [&](double r2) -> double
    {
        // inner integrand 1
        // r2 < r2p
        auto iintegrand1 = [&](double r2p) -> double
        {
            return Vni(r2p) * ric_j(li,ki*r2p) * (-ric_n(ln,kn*r2p));
        };
        BesselNodeIntegrator<decltype(iintegrand1),GaussKronrod<decltype(iintegrand1)>> R1(iintegrand1, ki, li);
        R1.setEpsAbs(1e-15);
        R1.setEpsRel(1e-8);
        R1.integrate(r2, Inf);
        if (not R1.ok())
        {
            throw exception
            (
                "Can't integrate 1st inner integral.\nError: %s\n", R1.status().c_str()
            );
        }
        
        // inner integrand 2
        // r2 > r2p
        auto iintegrand2 = [&](double r2p) -> double
        {
            return Vni(r2p) * ric_j(li,ki*r2p) * ric_j(ln,kn*r2p);
        };
        BesselNodeIntegrator<decltype(iintegrand2),GaussKronrod<decltype(iintegrand2)>> R2(iintegrand2, ki, li);
        R2.setEpsAbs(1e-15);
        R2.setEpsRel(1e-8);
        R2.integrate(0, r2);
        if (not R2.ok())
        {
            throw exception
            (
                "Can't integrate 2nd inner integral.\nError: %s\n", R2.status().c_str()
            );
        }
        
        // evaluate outer integrand
        return Vfn(r2) * ric_j(lf,kf*r2) * (ric_j(ln,kn*r2) * R1.result() - ric_n(ln,kn*r2) * R2.result());
    };
    auto integrand_im = [&](double r2) -> double
    {
        // inner integrand 1
        // r2 < r2p
        auto iintegrand1 = [&](double r2p) -> double
        {
            return Vni(r2p) * ric_j(li,ki*r2p) * ric_j(ln,kn*r2p);
        };
        BesselNodeIntegrator<decltype(iintegrand1),GaussKronrod<decltype(iintegrand1)>> R1(iintegrand1, ki, li);
        R1.setEpsAbs(1e-15);
        R1.setEpsRel(1e-8);
        R1.integrate(r2, Inf);
        if (not R1.ok())
        {
            throw exception
            (
                "Can't integrate 1st inner integral.\nError: %s\n", R1.status().c_str()
            );
        }
        
        // inner integrand 2
        // r2 > r2p
        auto iintegrand2 = [&](double r2p) -> double
        {
            return Vni(r2p) * ric_j(li,ki*r2p) * ric_j(ln,kn*r2p);
        };
        BesselNodeIntegrator<decltype(iintegrand2),GaussKronrod<decltype(iintegrand2)>> R2(iintegrand2, ki, li);
        R2.setEpsAbs(1e-15);
        R2.setEpsRel(1e-8);
        R2.integrate(0, r2);
        if (not R2.ok())
        {
            throw exception
            (
                "Can't integrate 2nd inner integral.\nError: %s\n", R2.status().c_str()
            );
        }
        
        // evaluate outer integrand
        return Vfn(r2) * ric_j(lf,kf*r2) * (ric_j(ln,kn*r2) * R1.result() + ric_j(ln,kn*r2) * R2.result());
    };
    
    // integrate (outer)
    BesselNodeIntegrator<decltype(integrand_re),GaussKronrod<decltype(integrand_re)>> Re(integrand_re, kf, lf);
    Re.setEpsAbs(1e-14);
    Re.setEpsRel(1e-8);
    Re.integrate(0, Inf);
    if (not Re.ok())
    {
        throw exception
        (
            "Can't integrate real outer integral.\nError: %s\n", Re.status().c_str()
        );
    }
    BesselNodeIntegrator<decltype(integrand_im),GaussKronrod<decltype(integrand_im)>> Im(integrand_im, kf, lf);
    Im.setEpsAbs(1e-14);
    Im.setEpsRel(1e-8);
    Im.integrate(0, Inf);
    if (not Im.ok())
    {
        throw exception
        (
            "Can't integrate imag outer integral.\nError: %s\n", Re.status().c_str()
        );
    }
    
    return Complex (Re.result(),Im.result());
}
