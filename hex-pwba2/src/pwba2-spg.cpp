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

#include <ginac/ginac.h>

#include "arrays.h"
#include "complex.h"
#include "misc.h"
#include "pwba2.h"
#include "vec3d.h"
#include "special.h"
#include "spgrid.h"

#include "clenshawcurtis.h"
#include "gausskronrod.h"

GiNaC::ex rl_Y
(
    int l, int m,
    GiNaC::symbol const & rp, GiNaC::symbol const & rm,
    GiNaC::realsymbol const & rz
)
{
    // evaluate prefactor
    GiNaC::ex prefactor = GiNaC::sqrt
    (
        (2*l+1) / (4*GiNaC::Pi) * GiNaC::factorial(l+m) * GiNaC::factorial(l-m)
    );
    
    // assemble the series
    GiNaC::ex poly;
    for (int k = std::max(0,-m); k <= (l-m)/2; k++)
    {
        poly += (GiNaC::pow(-rp,k+m) * GiNaC::pow(rm,k) * GiNaC::pow(rz,l-m-2*k)) /
                (GiNaC::pow(2,2*k+m) * GiNaC::factorial(k+m) * GiNaC::factorial(k) * GiNaC::factorial(l-m-2*k));
    }
    
    // return final expression for the spherical harmonic
    return prefactor * poly;
}

GiNaC::ex psi_nlm_poly
(
    int n, int l, int m,
    GiNaC::possymbol const & r,
    GiNaC::symbol const & rp, GiNaC::symbol const & rm,
    GiNaC::realsymbol const & rz
)
{
    // construct spherical part r^l Y_lm
    GiNaC::ex Y = rl_Y(l, m, rp, rm, rz);
    
    // construct radial prefactor
    GiNaC::ex N_nl = GiNaC::sqrt(GiNaC::pow(GiNaC::numeric(2,n),3)*GiNaC::factorial(n-l-1)/(2*n*GiNaC::factorial(n+l)));
    
    // construct remaining radial part R_nl / (N_nl r^l)
    GiNaC::ex R_nl;
    for (int i = 0; i <= n-l-1; i++)
        R_nl += GiNaC::pow(-1,i) * GiNaC::binomial(n+l,n-l-1-i) * GiNaC::pow(2*r/n,i) / GiNaC::factorial(i);
    
    // return the function
    return N_nl * R_nl * GiNaC::pow(GiNaC::numeric(2,n),l) * Y;
}

GiNaC::ex Wb_symb
(
    int n1, int l1, int m1,
    int n2, int l2, int m2,
    GiNaC::realsymbol const & kx,
    GiNaC::realsymbol const & ky,
    GiNaC::realsymbol const & kz
)
{
    // This routine will compute integral
    //                   -ik·r
    //   W = Int ψf(r)* e     ψi(r) d³r
    //
    // There are four symbols that the basic integral
    //            -ik·r  -ν|r| d³r      4π
    //   J = Int e      e      ——— = —————————
    //                          r    υ² + |k|²
    // depends on:
    //     ν  ... the combined exponential factor (1/n₁ + 1/n₂)
    //     k+ ... component of momentum k₁ + ik₂
    //     k- ...                       k₁ - ik₂
    //     kz ... third component of momentum
    // Derivative with respect to any of these parameters can be used to construct
    // a different integrand.
    
    // the basic integral
    GiNaC::possymbol k("k"), r("r"), nu("nu");
    GiNaC::realsymbol rz("rz");
    GiNaC::symbol rp("rp"), rm("rm"), kp("kp"), km("km");
    GiNaC::ex J = 4*GiNaC::Pi / (nu*nu + kp*km + kz*kz);
    
    // calculate value of ν as a sum of two fractions
    GiNaC::ex nu_val = GiNaC::numeric(1,n1) + GiNaC::numeric(1,n2);
    
    // construct the initial and final state (drop exponential factor which is handled separately)
    GiNaC::ex psif = psi_nlm_poly(n1,l1,m1,r,rp,rm,rz);
    GiNaC::ex psii = psi_nlm_poly(n2,l2,m2,r,rp,rm,rz);
    
    // calculate the product of the wave functions
    GiNaC::ex poly = (psif.conjugate() * psii).expand();
    
    // this is the resulting integral
    GiNaC::ex integral;
    
    // for all powers of r appearing in the polynomial 'poly'
    for (int inu = poly.ldegree(r); inu <= poly.degree(r); inu++)
    {
        // extract factor in front of r^inu
        GiNaC::ex nu_poly = poly.collect(r).coeff(r,inu);
        
        // (inu + 1)-times differentiate the basic integral and substitute
        GiNaC::ex J_nu = J.diff(nu,inu+1).subs(nu == nu_val) * GiNaC::pow(-1,inu+1);
        
        // for all powers of r+ appearing in the polynomial 'nu_poly'
        for (int irp = nu_poly.ldegree(rp); irp <= nu_poly.degree(rp); irp++)
        {
            // extract factor in front of ν^inu rp^irp
            GiNaC::ex nu_rp_poly = nu_poly.coeff(rp,irp);
            
            // irp-times differentiate the basic integral
            GiNaC::ex J_nu_rp = J_nu.diff(km,irp) * GiNaC::pow(2*GiNaC::I,irp);
            
            // for all powers of r- appearing in the polynomial 'nu_rp_poly'
            for (int irm = nu_rp_poly.ldegree(rm); irm <= nu_rp_poly.degree(rm); irm++)
            {
                // extract factor in front of ν^inu rp^irp rm^irm
                GiNaC::ex nu_rp_rm_poly = nu_rp_poly.coeff(rm,irm);
                
                // irm-times differentiate the basic integral
                GiNaC::ex J_nu_rp_rm = J_nu_rp.diff(kp,irm) * GiNaC::pow(2*GiNaC::I,irm);
                
                // for all powers of r- appearing in the polynomial 'nu_rp_poly'
                for (int irz = nu_rp_rm_poly.ldegree(rz); irz <= nu_rp_rm_poly.degree(rz); irz++)
                {
                    // extract factor in front of ν^inu rp^irp rm^irm rz^irz (it is a number)
                    GiNaC::ex nu_rp_rm_rz_coef = nu_rp_rm_poly.coeff(rz,irz);
                    
                    // irz-times differentiate the basic integral
                    GiNaC::ex J_nu_rp_rm_rz = J_nu_rp_rm.diff(kz,irz) * GiNaC::pow(GiNaC::I,irz);
                    
                    // add new term to the result
                    integral += nu_rp_rm_rz_coef * J_nu_rp_rm_rz;
                }
            }
        }
    }
    
    // finally, correct the elastic case
    if (n1 == n2 and l1 == l2 and m1 == m2)
        integral -= 1;
    
//     Debug << format("⟨%d,%d,%d|W|%d,%d,%d⟩: ", n1, l1, m1, n2, l2, m2) << (integral / (k*k)).subs(kz == GiNaC::sqrt(k*k-kp*km)).normal() << std::endl;
    
    // return the integral in its normal form
    return (integral / (k*k))
           .subs(k == GiNaC::sqrt(kx * kx + ky * ky + kz * kz))
           .subs(kp == kx + GiNaC::I * ky)
           .subs(km == kx - GiNaC::I * ky)
           .normal();
}

std::map<std::tuple<int,int,int,int,int>,Complex> Wb_symb_in
(
    int n1, int l1, int m1,
    int n2, int l2, int m2
)
{
    // This routine will compute integral
    //                   -ik·r
    //   W = Int ψf(r)* e     ψi(r) d³r
    //
    // There are four symbols that the basic integral
    //            -ik·r  -ν|r| d³r      4π
    //   J = Int e      e      ——— = —————————
    //                          r    υ² + |k|²
    // depends on:
    //     ν  ... the combined exponential factor (1/n₁ + 1/n₂)
    //     k+ ... component of momentum k₁ + ik₂
    //     k- ...                       k₁ - ik₂
    //     kz ... third component of momentum
    // Derivative with respect to any of these parameters can be used to construct
    // a different integrand.
    
    std::map<std::tuple<int,int,int,int,int>,Complex> res;
    
    //
    // Construct the product of wave functions.
    //
    
    // necessary symbolic variables
    GiNaC::possymbol r("r");
    GiNaC::realsymbol rz("rz");
    GiNaC::symbol rp("rp"), rm("rm");
    
    // construct the initial and final state (drop exponential factor which is handled separately)
    GiNaC::ex psif = psi_nlm_poly(n1,l1,m1,r,rp,rm,rz);
    GiNaC::ex psii = psi_nlm_poly(n2,l2,m2,r,rp,rm,rz);
    
    // calculate the product of the wave functions
    GiNaC::ex poly = (psif.conjugate() * psii).expand();
    
    //
    // Translate the symbolical form to list of exponent tuples.
    //
    
    // list of exponent tuples (one tuple for every term of 'poly')
    std::map<std::tuple<int,int,int,int>,Complex> terms;
    
    // for all powers of r appearing in the polynomial 'poly'
    for (int inu = poly.ldegree(r); inu <= poly.degree(r); inu++)
    {
        // extract factor in front of r^inu
        GiNaC::ex nu_poly = poly.collect(r).coeff(r,inu);
        
        // for all powers of r₊ appearing in the polynomial 'nu_poly'
        for (int irp = nu_poly.ldegree(rp); irp <= nu_poly.degree(rp); irp++)
        {
            // extract factor in front of ν^inu r₊^irp
            GiNaC::ex nu_rp_poly = nu_poly.coeff(rp,irp);
            
            // for all powers of r- appearing in the polynomial 'nu_rp_poly'
            for (int irm = nu_rp_poly.ldegree(rm); irm <= nu_rp_poly.degree(rm); irm++)
            {
                // extract factor in front of ν^inu r₊^irp r₋^irm
                GiNaC::ex nu_rp_rm_poly = nu_rp_poly.coeff(rm,irm);
                
                // for all powers of r- appearing in the polynomial 'nu_rp_poly'
                for (int irz = nu_rp_rm_poly.ldegree(rz); irz <= nu_rp_rm_poly.degree(rz); irz++)
                {
                    // extract factor in front of ν^inu r₊^irp r₋^irm rz^irz (it is a number)
                    GiNaC::ex nu_rp_rm_rz_coef = 4 * GiNaC::Pi * nu_rp_rm_poly.coeff(rz,irz);
                    terms[std::make_tuple(inu,irp,irm,irz)] = Complex
                    (
                        GiNaC::ex_to<GiNaC::numeric>(GiNaC::real_part(nu_rp_rm_rz_coef.evalf())).to_double(),
                        GiNaC::ex_to<GiNaC::numeric>(GiNaC::imag_part(nu_rp_rm_rz_coef.evalf())).to_double()
                    );
                } // irz
            } // irm
        } // irp
    } // inu
    
    //
    // Compute derivative for every term of the product.
    //
    
    for (auto term : terms)
    {
        // decode exponents
        int inu, irp, irm, irz;
        std::tie(inu, irp, irm, irz) = term.first;
        
        // resulting derivative for this particular term
        std::map<std::tuple<int,int,int,int,int>,Complex> rest, tmp;
        rest[std::make_tuple(1,0,0,0,0)] = term.second;
        
        // (inu + 1)-times differentiate the basic integral (-∂/∂ν)
        for (int j = 0; j <= inu; j++)
        {
            // move current result to temporary storage
            std::swap(rest,tmp); rest.clear();
            
            // differentiate all terms
            for (auto entry : tmp)
            {
                // decode entry
                int n, a, b, c, d; std::tie(n,a,b,c,d) = entry.first;
                
                // differentiate with respect to ν
                if (a != 0)
                {
                    // create key for the new term
                    auto key = std::make_tuple(n,a-1,b,c,d);
                    
                    // look for the term
                    auto it = rest.find(key);
                    
                    // update the term
                    rest[key] = -double(a) * entry.second + (it == rest.end() ? 0. : it->second);
                }
                {
                    // create key for the new term
                    auto key = std::make_tuple(n+1,a+1,b,c,d);
                    
                    // look for the term
                    auto it = rest.find(key);
                    
                    // update the term
                    rest[key] = 2. * n * entry.second + (it == rest.end() ? 0. : it->second);
                }
            }
        }
        
        // irp times differentiate the integral (2i∂/∂k₋)
        for (int j = 0; j < irp; j++)
        {
            // move current result to temporary storage
            std::swap(rest,tmp); rest.clear();
            
            // differentiate all terms
            for (auto entry : tmp)
            {
                // decode entry
                int n, a, b, c, d; std::tie(n,a,b,c,d) = entry.first;
                
                // differentiate with respect to k₋
                if (c != 0)
                {
                    // create key for the new term
                    auto key = std::make_tuple(n,a,b,c-1,d);
                    
                    // look for the term
                    auto it = rest.find(key);
                    
                    // update the term
                    rest[key] = Complex(0.,2.*c) * entry.second + (it == rest.end() ? 0. : it->second);
                }
                {
                    // create key for the new term
                    auto key = std::make_tuple(n+1,a,b+1,c,d);
                    
                    // look for the term
                    auto it = rest.find(key);
                    
                    // update the term
                    rest[key] = Complex(0.,-2.*n) * entry.second + (it == rest.end() ? 0. : it->second);
                }
            }
        }
        
        // irm times differentiate the integral (2i∂/∂k₊)
        for (int j = 0; j < irm; j++)
        {
            // move current result to temporary storage
            std::swap(rest,tmp); rest.clear();
            
            // differentiate all terms
            for (auto entry : tmp)
            {
                // decode entry
                int n, a, b, c, d; std::tie(n,a,b,c,d) = entry.first;
                
                // differentiate with respect to k₊
                if (b != 0)
                {
                    // create key for the new term
                    auto key = std::make_tuple(n,a,b-1,c,d);
                    
                    // look for the term
                    auto it = rest.find(key);
                    
                    // update the term
                    rest[key] = Complex(0.,2.*b) * entry.second + (it == rest.end() ? 0. : it->second);
                }
                {
                    // create key for the new term
                    auto key = std::make_tuple(n+1,a,b,c+1,d);
                    
                    // look for the term
                    auto it = rest.find(key);
                    
                    // update the term
                    rest[key] = Complex(0.,-2.*n) * entry.second + (it == rest.end() ? 0. : it->second);
                }
            }
        }
        
        // irz times differentiate the integral (i∂/∂kz)
        for (int j = 0; j < irz; j++)
        {
            // move current result to temporary storage
            std::swap(rest,tmp); rest.clear();
            
            // differentiate all terms
            for (auto entry : tmp)
            {
                // decode entry
                int n, a, b, c, d; std::tie(n,a,b,c,d) = entry.first;
                
                // differentiate with respect to kz
                if (d != 0)
                {
                    // create key for the new term
                    auto key = std::make_tuple(n,a,b,c,d-1);
                    
                    // look for the term
                    auto it = rest.find(key);
                    
                    // update the term
                    rest[key] = Complex(0.,d) * entry.second + (it == rest.end() ? 0. : it->second);
                }
                {
                    // create key for the new term
                    auto key = std::make_tuple(n+1,a,b,c,d+1);
                    
                    // look for the term
                    auto it = rest.find(key);
                    
                    // update the term
                    rest[key] = Complex(0.,-2.*n) * entry.second + (it == rest.end() ? 0. : it->second);
                }
            }
        }
        
        // add all new terms to result
        for (auto t : rest)
        {
            // check if this term already exists
            auto pt = res.find(t.first);
            
            // update term coefficient
            res[t.first] = t.second + (pt != res.end() ? pt->second : 0.);
        }
    }
    
    // finally, correct the elastic case
    if (n1 == n2 and l1 == l2 and m1 == m2)
        res[std::make_tuple(0,0,0,0,0)] = -1;
    
    return res;
}

Complex eval_Wb (std::map<std::tuple<int,int,int,int,int>,Complex> const & poly, double nu, geom::vec3d const & k)
{
    double ksqr = k.x*k.x + k.y*k.y + k.z*k.z, Z = nu*nu + ksqr;
    Complex kp(k.x,k.y), km(k.x,-k.y), kz = k.z;
    
    Complex result = 0;
    for (auto term : poly)
    {
        int n,a,b,c,d;
        std::tie(n,a,b,c,d) = term.first;
        result += term.second * std::pow(Z,-n) * std::pow(nu,a) * std::pow(kp,b) * std::pow(km,c) * std::pow(kz,d);
    }
    
    return result / ksqr;
}

GiNaC::ex W_symb
(
    int n, int l, int m,
    GiNaC::realsymbol const & kx,
    GiNaC::realsymbol const & ky,
    GiNaC::realsymbol const & kz,
    GiNaC::realsymbol const & qx,
    GiNaC::realsymbol const & qy,
    GiNaC::realsymbol const & qz
)
{
    // This routine will compute integral
    //                  -ik·r
    //   W = Int χ(r)* e     ψ(r) d³r
    //
    // (without the normalization factor).
    // There are seven symbols that the basic integral
    //                   -ik·r  -ν|r| d³r   (k² + (ν-iq)²)^(-i/q)
    //   I = Int χ(q,r) e      e      ——— = —————————————————————
    //                                 r    (υ² + |k+q|²)^(1-i/q)
    // depends on:
    //     ν  ... the combined exponential factor (1/n₁ + 1/n₂)
    //     kx,ky,kz ... components of momentum of the plane wave
    //     qx,qy,qz ... components of momentum of the Coulomb wave
    // Derivative with respect to any of the parameters ν,kx,ky,kz can be used to construct
    // a different integrand.
    
    // the basic integral
    GiNaC::possymbol k("k"), r("r"), nu("nu"), q("q");
    GiNaC::realsymbol rz("rz");
    GiNaC::symbol rp("rp"), rm("rm"), kp("kp"), km("km"), qp("qp"), qm("qm");
    GiNaC::ex I = GiNaC::pow(k*k+GiNaC::pow(nu-GiNaC::I*q,2),-GiNaC::I/q) /
                  GiNaC::pow(nu*nu+(kp+qp)*(km+qm)+(kz+qz)*(kz+qz),1-GiNaC::I/q);
    
    // calculate value of ν as a sum of two fractions
    GiNaC::ex nu_val = GiNaC::numeric(1,n);
    
    // construct the initial and final state (drop exponential factor which is handled separately)
    GiNaC::ex psi = psi_nlm_poly(n,l,m,r, rp, rm, rz);
    
    // this is the resulting integral
    GiNaC::ex integral;
    
    // for all powers of r appearing in the polynomial 'poly'
    for (int inu = psi.ldegree(r); inu <= psi.degree(r); inu++)
    {
        // extract factor in front of r^inu
        GiNaC::ex nu_poly = psi.collect(r).coeff(r,inu);
        
        // (inu + 1)-times differentiate the basic integral and substitute
        GiNaC::ex I_nu = I.diff(nu,inu+1).subs(nu == nu_val) * GiNaC::pow(-1,inu+1);
        
        // for all powers of r+ appearing in the polynomial 'nu_poly'
        for (int irp = nu_poly.ldegree(rp); irp <= nu_poly.degree(rp); irp++)
        {
            // extract factor in front of ν^inu rp^irp
            GiNaC::ex nu_rp_poly = nu_poly.coeff(rp,irp);
            
            // irp-times differentiate the basic integral
            GiNaC::ex I_nu_rp = I_nu.diff(km,irp) * GiNaC::pow(2*GiNaC::I,irp);
            
            // for all powers of r- appearing in the polynomial 'nu_rp_poly'
            for (int irm = nu_rp_poly.ldegree(rm); irm <= nu_rp_poly.degree(rm); irm++)
            {
                // extract factor in front of ν^inu rp^irp rm^irm
                GiNaC::ex nu_rp_rm_poly = nu_rp_poly.coeff(rm,irm);
                
                // irm-times differentiate the basic integral
                GiNaC::ex I_nu_rp_rm = I_nu_rp.diff(kp,irm) * GiNaC::pow(2*GiNaC::I,irm);
                
                // for all powers of r- appearing in the polynomial 'nu_rp_poly'
                for (int irz = nu_rp_rm_poly.ldegree(rz); irz <= nu_rp_rm_poly.degree(rz); irz++)
                {
                    // extract factor in front of ν^inu rp^irp rm^irm rz^irz (it is a number)
                    GiNaC::ex nu_rp_rm_rz_coef = nu_rp_rm_poly.coeff(rz,irz);
                    
                    // irz-times differentiate the basic integral
                    GiNaC::ex I_nu_rp_rm_rz = I_nu_rp_rm.diff(kz,irz) * GiNaC::pow(GiNaC::I,irz);
                    
                    // add new term to the result
                    integral += nu_rp_rm_rz_coef * I_nu_rp_rm_rz;
                }
            }
        }
    }
    
//     Debug << format("⟨%d,%d,%d|W|χ(q)⟩: ", n, l, m) << (integral / (k*k)).subs(kz == GiNaC::sqrt(k*k-kp*km)).normal() << std::endl;
    
    // return the integral in its normal form
    return (integral / (k*k))
           .subs(k == GiNaC::sqrt(kx * kx + ky * ky + kz * kz))
           .subs(kp == kx + GiNaC::I * ky)
           .subs(km == kx - GiNaC::I * ky)
           .subs(q == GiNaC::sqrt(qx * qx + qy * qy + qz * qz))
           .subs(qp == qx + GiNaC::I * qy)
           .subs(qm == qx - GiNaC::I * qy)
           .normal();
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
    
    Complex iq (0.,1./q); // = i/q
    Complex A = WA(k,nu,q);
    double B = WB(vk,nu,vq);
    
    Complex B_A = B / A;
    return 2.0 * iq * std::pow(B_A,iq) * (Complex(-1.,q) * B_A + Complex(1.,q)) / (B * B * k * k);
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
    bool integrate_allowed, bool integrate_forbidden
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
        G.setLocEpsRel(1e-5);   // ... if the rules are equal up to 1e-5 (relative)
        G.setLocEpsAbs(1e-9);   // ... if the rules are equal up to 1e-9 (absolute)
        G.setGlobEpsRel(1e-6);  // ... if the extrapolated relative contribution to the whole domain is within 1e-6
        G.setGlobEpsAbs(0);     // ... if the extrapolated contribution to the whole domain is within ... [not used]
        G.setVerbose(true);     // print detailed progress information
        G.setParallel(true);    // evaluate integration domains in parallel
        G.setPrefix("    ");    // output formatting prefix
        
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
            
            // real bound integrand
            auto integrand_Ub_wrap_lin = [ki,kf,vki,vkf,Ni,Nn,Nf,Ln,Qon,Etot,Wbf,Wbi](int n, double const * coords) -> Complex
            {
                // check dimensions
                assert(n == 3);
                
                // unpack polar coordinates
                double costheta = 2 * coords[0] - 1;
                double sintheta = std::sqrt(1 - costheta * costheta);
                double phi = special::constant::two_pi  * coords[1];
                double kn = Qon * coords[2];
                
                // compute on-shell momentum magnitude
                double Qn = std::sqrt(2*Etot + 1./(Nn*Nn));
                
                // compute both off- and on-shell momentum vectors
                geom::vec3d vn = { sintheta * std::cos(phi), sintheta * std::sin(phi), costheta };
                geom::vec3d vkn = kn * vn;
                geom::vec3d vQn = Qn * vn;
                
                // the value of the off-shell integrand
                Complex Wf = eval_Wb(Wbf, 1./Nn + 1./Nf, vkf - vkn);
                Complex Wi = eval_Wb(Wbi, 1./Nn + 1./Ni, vki - vkn);
                Complex integrand_Ub_off = kn * kn * Wf * Wi;
                
                // the value of the on-shell integrand
                Complex Wfon = eval_Wb(Wbf, 1./Nn + 1./Nf, vkf - vQn);
                Complex Wion = eval_Wb(Wbi, 1./Nn + 1./Ni, vki - vQn);
                Complex integrand_Ub_on = Qn * Qn * Wfon * Wion;
                
                // Jacobian
                double Jac = 2. // for cos theta
                            * special::constant::two_pi // for phi
                            * Qon; // for Q
                
                // evaluate integrand
                return special::constant::two_inv_pi * Jac * (integrand_Ub_off - integrand_Ub_on) / (0.5*Qn*Qn - 0.5*kn*kn);
            };
            
            // real bound integrand
            auto integrand_Ub_wrap_log = [ki,kf,vki,vkf,Ni,Nn,Nf,Ln,Qon,Etot,Wbf,Wbi](int n, double const * coords) -> Complex
            {
                // check dimensions
                assert(n == 3);
                
                // compactification parameter
                const int m = 1;
                
                // unpack polar coordinates
                double costheta = 2 * coords[0] - 1;
                double sintheta = std::sqrt(1 - costheta * costheta);
                double phi = special::constant::two_pi  * coords[1];
                double kn = Qon / std::pow(1 - coords[2], m);
                
                // compute on-shell momentum magnitude
                double Qn = std::sqrt(2*Etot + 1./(Nn*Nn));
                
                // compute both off- and on-shell momentum vectors
                geom::vec3d vn = { sintheta * std::cos(phi), sintheta * std::sin(phi), costheta };
                geom::vec3d vkn = kn * vn;
                geom::vec3d vQn = Qn * vn;
                
                // the value of the off-shell integrand
                Complex Wf = eval_Wb(Wbf, 1./Nn + 1./Nf, vkf - vkn);
                Complex Wi = eval_Wb(Wbi, 1./Nn + 1./Ni, vki - vkn);
                Complex integrand_Ub_off = kn * kn * Wf * Wi;
                
                // the value of the on-shell integrand
                Complex Wfon = eval_Wb(Wbf, 1./Nn + 1./Nf, vkf - vQn);
                Complex Wion = eval_Wb(Wbi, 1./Nn + 1./Ni, vki - vQn);
                Complex integrand_Ub_on = Qn * Qn * Wfon * Wion;
                
                // Jacobian
                double Jac =  2. // for cos theta
                            * special::constant::two_pi // for phi
                            * m*Qon/std::pow(1 - coords[2],m+1); // for Q
                
                // evaluate integrand
                return special::constant::two_inv_pi * Jac * (integrand_Ub_off - integrand_Ub_on) / (0.5*Qn*Qn - 0.5*kn*kn);
            };
            
            std::cout << format("  Sparse grid marching integration of bound U for theta = %g, phi = %g", theta, phi) << std::endl;
            
            // allow at most 70 integration cell bisections
            G.setMaxLevel(70);
            
            // integrate U on 3-dimensional sparse grid
            std::cout << std::endl<< "  Linear integrand for Q = 0 .. " << Qon << std::endl;
            G.integrate_adapt<3>(integrand_Ub_wrap_lin, spgrid::d3l4n39, spgrid::d3l5n87);
            fUb_contrib = G.result();
            nEvalUb += G.evalcount();
            
            // integrate U on 3-dimensional sparse grid
            std::cout << std::endl<< "  Compactified integrand for Q = " << Qon << " .. ∞" << std::endl;
            G.integrate_adapt<3>(integrand_Ub_wrap_log, spgrid::d3l4n39, spgrid::d3l5n87);
            fUb_contrib += G.result();
            nEvalUb += G.evalcount();
            
            std::cout << std::endl;
            
            // imag bound integrand
            auto integrand_Wb_wrap = [ki,kf,vki,vkf,Ni,Nn,Nf,Ln,Qon,Etot,Wbf,Wbi](int n, double const * coords) -> Complex
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
                Complex Wf = eval_Wb(Wbf, 1./Nn + 1./Nf, vkf - vQn);
                Complex Wi = eval_Wb(Wbi, 1./Nn + 1./Ni, vki - vQn);
                Complex integrand_Wb = Complex(0.,-special::constant::pi) * Qn * Wf * Wi;
                
                // Jacobian
                double Jac = 2. // for cos theta
                        * special::constant::two_pi; // for phi
                
                // evaluate integrand
                return special::constant::two_inv_pi * Jac * integrand_Wb;
            };
            
            std::cout << format("  Sparse grid integration of bound W for theta = %g, phi = %g", theta, phi) << std::endl;
            
            // integrate W on 2-dimensional sparse grid
            G.integrate_adapt<2>(integrand_Wb_wrap, spgrid::d2l4n17, spgrid::d2l7n65);
            fWb_contrib = G.result();
            nEvalWb += G.evalcount();
            
            std::cout << std::endl;
            
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
            
            // U-integrand (real part of propagator)
            auto integrand_U_wrap_lin = [Ni,Li,Mi,ki,kf,vki,vkf,Qon,Etot](int n, double const * coords) -> Complex
            {
                // check dimensions
                assert(n == 6);
                
                // unpack polar coordinates
                double costheta1 = 2 * coords[0] - 1, sintheta1 = std::sqrt(1 - costheta1 * costheta1);
                double costheta2 = 2 * coords[1] - 1, sintheta2 = std::sqrt(1 - costheta2 * costheta2);
                double phi1  = special::constant::two_pi  * coords[2];
                double phi2  = special::constant::two_pi  * coords[3];
                double alpha = special::constant::pi_half * coords[4];
                double Q = Qon * coords[5];
                
                // compute both off- and on-shell momentum magnitudes
                double qn = Q * std::sin(alpha), qnon = Qon * std::sin(alpha);
                double kn = Q * std::cos(alpha), knon = Qon * std::cos(alpha);
                
                // compute both off- and on-shell momentum vectors
                geom::vec3d n1 = { sintheta1 * std::cos(phi1), sintheta1 * std::sin(phi1), costheta1 };
                geom::vec3d n2 = { sintheta2 * std::cos(phi2), sintheta2 * std::sin(phi2), costheta2 };
                geom::vec3d vqn = qn * n1, vqnon = qnon * n1;
                geom::vec3d vkn = kn * n2, vknon = knon * n2;
                
                // the value of the off-shell integrand
                double norm = 4. / (special::constant::pi * qn * (1. - std::exp(-special::constant::two_pi/qn)));
                Complex Wf = W_cont(Ni,Li,Mi,vkf-vkn,vqn), Wi = std::conj(W_cont(Ni,Li,Mi,vki-vkn,vqn));
                Complex integrand_U_off = qn * qn * kn * kn * Q * norm * Wf * Wi;
                
                // the value of the on-shell integrand
                double normon = 4. / (special::constant::pi * qnon * (1. - std::exp(-special::constant::two_pi/qnon)));
                Complex Wfon = W_cont(Ni,Li,Mi,vkf-vknon,vqnon), Wion = std::conj(W_cont(Ni,Li,Mi,vki-vknon,vqnon));
                Complex integrand_U_on = qnon * qnon * knon * knon * Qon * normon * Wfon * Wion;
                
                // Jacobian
                double Jac =  2. // for costheta1
                            * 2. // for costheta2
                            * special::constant::two_pi // for phi1
                            * special::constant::two_pi // for phi2
                            * special::constant::pi_half // for alpha
                            * Qon; // for Q
                
                // evaluate integrand
                return special::constant::two_inv_pi * Jac * (integrand_U_off - integrand_U_on) / (Etot - 0.5*Q*Q);
            };
            
            // U-integrand (real part of propagator)
            auto integrand_U_wrap_log = [Ni,Li,Mi,ki,kf,vki,vkf,Qon,Etot](int n, double const * coords) -> Complex
            {
                // check dimensions
                assert(n == 6);
                
                // compactification parameter
                const int m = 1;
                
                // unpack polar coordinates
                double costheta1 = 2 * coords[0] - 1, sintheta1 = std::sqrt(1 - costheta1 * costheta1);
                double costheta2 = 2 * coords[1] - 1, sintheta2 = std::sqrt(1 - costheta2 * costheta2);
                double phi1  = special::constant::two_pi  * coords[2];
                double phi2  = special::constant::two_pi  * coords[3];
                double alpha = special::constant::pi_half * coords[4];
                double Q = Qon / std::pow(1 - coords[5], m);
                
                // compute both off- and on-shell momentum magnitudes
                double qn = Q * std::sin(alpha), qnon = Qon * std::sin(alpha);
                double kn = Q * std::cos(alpha), knon = Qon * std::cos(alpha);
                
                // compute both off- and on-shell momentum vectors
                geom::vec3d n1 = { sintheta1 * std::cos(phi1), sintheta1 * std::sin(phi1), costheta1 };
                geom::vec3d n2 = { sintheta2 * std::cos(phi2), sintheta2 * std::sin(phi2), costheta2 };
                geom::vec3d vqn = qn * n1, vqnon = qnon * n1;
                geom::vec3d vkn = kn * n2, vknon = knon * n2;
                
                // the value of the off-shell integrand
                double norm = 4. / (special::constant::pi * qn * (1. - std::exp(-special::constant::two_pi/qn)));
                Complex Wf = W_cont(Ni,Li,Mi,vkf-vkn,vqn), Wi = std::conj(W_cont(Ni,Li,Mi,vki-vkn,vqn));
                Complex integrand_U_off = qn * qn * kn * kn * Q * norm * Wf * Wi;
                
                // the value of the on-shell integrand
                double normon = 4. / (special::constant::pi * qnon * (1. - std::exp(-special::constant::two_pi/qnon)));
                Complex Wfon = W_cont(Ni,Li,Mi,vkf-vknon,vqnon), Wion = std::conj(W_cont(Ni,Li,Mi,vki-vknon,vqnon));
                Complex integrand_U_on = qnon * qnon * knon * knon * Qon * normon * Wfon * Wion;
                
                // Jacobian
                double Jac =  2. // for costheta1
                            * 2. // for costheta2
                            * special::constant::two_pi // for phi1
                            * special::constant::two_pi // for phi2
                            * special::constant::pi_half // for alpha
                            * m * Qon / std::pow(1 - coords[4], m + 1); // for Q
                
                // evaluate integrand
                return special::constant::two_inv_pi * Jac * (integrand_U_off - integrand_U_on) / (Etot - 0.5*Q*Q);
            };
            
            // allow at most 10 integration cell bisections
            G.setMaxLevel(10);
            G.setGlobEpsAbs(0);
            
            std::cout << underline(format("Sparse grid integration of continuum U for theta = %g, phi = %g", theta, phi)) << std::endl;
            
            // integrate U on 6-dimensional sparse grid
            std::cout << std::endl << "Linear integrand for Q = 0 .. " << Qon << std::endl;
            G.integrate_adapt<6>(integrand_U_wrap_lin, spgrid::d6l4n257, spgrid::d6l5n737);
            fU += G.result();
            nEvalU += G.evalcount();
            
            // integrate U on 6-dimensional sparse grid
            std::cout << std::endl << "Compactified integrand for Q = " << Qon << " .. ∞" << std::endl;
            G.integrate_adapt<6>(integrand_U_wrap_log, spgrid::d6l4n257, spgrid::d6l5n737);
            fU += G.result();
            nEvalU += G.evalcount();
            std::cout << std::endl;
            G.setGlobEpsAbs(0);
            
            // iW-integrand (imag part of propagator)
            auto integrand_W_wrap = [Ni,Li,Mi,ki,kf,vki,vkf,Etot,Qon](int n, double const * coords) -> Complex
            {
                // check dimensions
                assert(n == 5);
                
                // unpack polar coordinates
                double costheta1 = 2 * coords[0] - 1, sintheta1 = std::sqrt(1 - costheta1 * costheta1);
                double costheta2 = 2 * coords[1] - 1, sintheta2 = std::sqrt(1 - costheta2 * costheta2);
                double phi1  = special::constant::two_pi  * coords[2];
                double phi2  = special::constant::two_pi  * coords[3];
                double alpha = special::constant::pi_half * coords[4];
                double qn = Qon * std::sin(alpha);
                double kn = Qon * std::cos(alpha);
                
                // compute directions
                geom::vec3d vqn = { qn * sintheta1 * std::cos(phi1), qn * sintheta1 * std::sin(phi1), qn * costheta1 };
                geom::vec3d vkn = { kn * sintheta2 * std::cos(phi2), kn * sintheta2 * std::sin(phi2), kn * costheta2 };
                
                // the value of the integrand
                double norm = 4. / (special::constant::pi * qn * (1. - std::exp(-special::constant::two_pi/qn)));
                Complex Wf = W_cont(Ni,Li,Mi,vkf-vkn,vqn), Wi = std::conj(W_cont(Ni,Li,Mi,vki-vkn,vqn));
                Complex integrand_W = Complex(0.,-special::constant::pi) * norm * Wf * Wi;
                
                // Jacobian
                double Jac = 2. // for costheta1
                           * 2. // for costheta2
                           * special::constant::two_pi // for phi1
                           * special::constant::two_pi // for phi2
                           * special::constant::pi_half; // for alpha
                
                // evaluate integrand
                return special::constant::two_inv_pi * Jac * qn*qn * kn*kn * integrand_W;
            };
            
            std::cout << format("  Sparse grid integration of continuum W for theta = %g, phi = %g", theta, phi) << std::endl;
            
            // integrate W on 5-dimensional sparse grid
            G.integrate_adapt<5>(integrand_W_wrap, spgrid::d5l4n151, spgrid::d5l5n391);
            fW = G.result();
            nEvalW = G.evalcount();
            std::cout << std::endl;
        }
        
        // write new row to the output table
        tab.write(itheta, phi, fUb, nEvalUb, fWb, nEvalWb, fU, nEvalU, fW, nEvalW, fUb+fWb+fU+fW);
    }
    
    // write out the table with T-matrices
    std::cout << std::endl;
    std::cout << underline("Resulting T-matrices") << std::endl;
    std::cout << table.str() << std::endl;
    
    // for now, exit
    std::exit(0);
}
