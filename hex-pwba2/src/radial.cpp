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
    
    // setup the integral
    auto integrand_re = [&](double r2) -> double
    {
        if (r2 == 0.)
            return 0;
        
        // inner integrand 1
        // r2 < r2p
        auto iintegrand1 = [&](double r2p) -> double
        {
            return Vni(r2p) * ric_j(li,ki*r2p) * (-ric_n(ln,kn*r2p));
        };
        BesselNodeIntegrator<decltype(iintegrand1),GaussKronrod<decltype(iintegrand1)>> R1(iintegrand1, ki, li);
        R1.setEpsAbs(1e-15);
        R1.setEpsRel(1e-6);
        R1.integrate(r2, special::constant::Inf);
        if (not R1.ok())
        {
            throw exception
            (
                "Can't integrate 1st real inner integral.\nError: %s\n", R1.status().c_str()
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
        R2.setEpsRel(1e-6);
        R2.integrate(0, r2);
        if (not R2.ok())
        {
            throw exception
            (
                "Can't integrate 2nd real inner integral.\nError: %s\n", R2.status().c_str()
            );
        }
        
        // evaluate outer integrand
        return Vfn(r2) * ric_j(lf,kf*r2) * (ric_j(ln,kn*r2) * R1.result() - ric_n(ln,kn*r2) * R2.result());
    };
    auto integrand_im = [&](double r2) -> double
    {
        if (r2 == 0.)
            return 0.;
        
        // inner integrand 1
        // r2 < r2p
        auto iintegrand1 = [&](double r2p) -> double
        {
            return Vni(r2p) * ric_j(li,ki*r2p) * ric_j(ln,kn*r2p);
        };
        BesselNodeIntegrator<decltype(iintegrand1),GaussKronrod<decltype(iintegrand1)>> R1(iintegrand1, ki, li);
        R1.setEpsAbs(1e-15);
        R1.setEpsRel(1e-6);
        R1.integrate(r2, special::constant::Inf);
        if (not R1.ok())
        {
            throw exception
            (
                "Can't integrate 1st imag inner integral.\nError: %s\n", R1.status().c_str()
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
        R2.setEpsRel(1e-6);
        R2.integrate(0, r2);
        if (not R2.ok())
        {
            throw exception
            (
                "Can't integrate 2nd imag inner integral.\nError: %s\n", R2.status().c_str()
            );
        }
        
        // evaluate outer integrand
        return Vfn(r2) * ric_j(lf,kf*r2) * (ric_j(ln,kn*r2) * R1.result() + ric_j(ln,kn*r2) * R2.result());
    };
    
    
    rArray rs = linspace(0., 10., 1001);
    
    write_1D_data(rs.size(), "Vni.dat", [&](int i) -> double { return Vni(rs[i]); });
    write_1D_data(rs.size(), "Vfn.dat", [&](int i) -> double { return Vfn(rs[i]); });
    
    // integrate (outer)
    BesselNodeIntegrator<decltype(integrand_re),GaussKronrod<decltype(integrand_re)>> Re(integrand_re, kf, lf);
    write_1D_data(rs.size(), "re.dat", [&](int i) -> double { return integrand_re(rs[i]); });
    Re.setEpsAbs(1e-14);
    Re.setEpsRel(1e-8);
//     Re.integrate(0, special::constant::Inf);
    if (not Re.ok())
    {
        throw exception
        (
            "Can't integrate real outer integral.\nError: %s\n", Re.status().c_str()
        );
    }
    BesselNodeIntegrator<decltype(integrand_im),GaussKronrod<decltype(integrand_im)>> Im(integrand_im, kf, lf);
    write_1D_data(rs.size(), "im.dat", [&](int i) -> double { return integrand_im(rs[i]); });
    Im.setEpsAbs(1e-14);
    Im.setEpsRel(1e-8);
//     Im.integrate(0, special::constant::Inf);
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
    MultipolePotential Vfn (lambdaf, Nf, Lf, Kn, Ln);
    MultipolePotential Vni (lambdai, Kn, Ln, Ni, Li);
    
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
        R1.integrate(r2, special::constant::Inf);
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
        R1.integrate(r2, special::constant::Inf);
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
    Re.integrate(0, special::constant::Inf);
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
    Im.integrate(0, special::constant::Inf);
    if (not Im.ok())
    {
        throw exception
        (
            "Can't integrate imag outer integral.\nError: %s\n", Re.status().c_str()
        );
    }
    
    return Complex (Re.result(),Im.result());
}
