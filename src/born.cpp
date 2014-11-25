//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2014, Jakub Benda, Charles University in Prague                    //
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

#include <cmath>
#include <map>

#include <ginac/ginac.h>
#include <gsl/gsl_sf.h>

#include "complex.h"
#include "special.h"
#include "vec3d.h"

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
    ) * GiNaC::pow(-1, std::max(0,m));
    
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
    
    // resulting terms (ν²+k²)^(-n) ν^a k₊^b k₋^c kz^d
    std::map<std::tuple<int,int,int,int,int>,Complex> res;
    
    // for all terms
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
                    auto key = std::make_tuple(n,a-1,b,c,d); auto it = rest.find(key);
                    rest[key] = -double(a) * entry.second + (it == rest.end() ? 0. : it->second);
                }
                {
                    auto key = std::make_tuple(n+1,a+1,b,c,d); auto it = rest.find(key);
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
                    auto key = std::make_tuple(n,a,b,c-1,d); auto it = rest.find(key);
                    rest[key] = Complex(0.,2.*c) * entry.second + (it == rest.end() ? 0. : it->second);
                }
                {
                    auto key = std::make_tuple(n+1,a,b+1,c,d); auto it = rest.find(key);
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
                    auto key = std::make_tuple(n,a,b-1,c,d); auto it = rest.find(key);
                    rest[key] = Complex(0.,2.*b) * entry.second + (it == rest.end() ? 0. : it->second);
                }
                {
                    auto key = std::make_tuple(n+1,a,b,c+1,d); auto it = rest.find(key);
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
                    auto key = std::make_tuple(n,a,b,c,d-1); auto it = rest.find(key);
                    rest[key] = Complex(0.,d) * entry.second + (it == rest.end() ? 0. : it->second);
                }
                {
                    auto key = std::make_tuple(n+1,a,b,c,d+1); auto it = rest.find(key);
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
    // prepare some used variables
    double kxysqr = k.x*k.x + k.y*k.y, ksqr = kxysqr + k.z*k.z;
    double rho = std::sqrt(kxysqr), phi = std::atan2(k.y,k.x);
    double J = 1. / (nu*nu + ksqr);
    
    // evaluate all terms J^n ν^a k₊^b k₋^c kz^d
    Complex result = 0;
    for (auto term : poly)
    {
        // extract exponents
        int n,a,b,c,d; std::tie(n,a,b,c,d) = term.first;
        
        // contribution from this term
        Complex contrib = 4. * special::constant::pi * term.second * gsl_sf_pow_int(J,n);
        
        // add ν^a
        if (a != 0)
            contrib *= gsl_sf_pow_int(nu,a);
        
        // add k₊^b k₋^c
        if (b != 0 or c != 0)
            contrib *= gsl_sf_pow_int(rho,b+c) * Complex(std::cos((b-c)*phi),std::sin((b-c)*phi));
        
        // add kz^d
        if (d != 0)
            contrib *= gsl_sf_pow_int(k.z,d);
        
        // update sum
        result += contrib;
    }
    
    // return result
    return result / ksqr;
}

std::map<std::tuple<int,int,int,int,int,int,int,int,int,int>,Complex> W_symb_in
(
    int n, int l, int m
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
    GiNaC::possymbol r("r");
    GiNaC::realsymbol rz("rz");
    GiNaC::symbol rp("rp"), rm("rm");
    
    // construct the initial and final state (drop exponential factor which is handled separately)
    GiNaC::ex psi = psi_nlm_poly(n,l,m,r, rp, rm, rz);
    
    //
    // Translate the symbolical form to list of exponent tuples.
    //
    
    // list of exponent tuples (one tuple for every term of 'poly')
    std::map<std::tuple<int,int,int,int>,Complex> terms;
    
    // for all powers of r appearing in the polynomial 'poly'
    for (int inu = psi.ldegree(r); inu <= psi.degree(r); inu++)
    {
        // extract factor in front of r^inu
        GiNaC::ex nu_poly = psi.collect(r).coeff(r,inu);
        
        // for all powers of r+ appearing in the polynomial 'nu_poly'
        for (int irp = nu_poly.ldegree(rp); irp <= nu_poly.degree(rp); irp++)
        {
            // extract factor in front of ν^inu rp^irp
            GiNaC::ex nu_rp_poly = nu_poly.coeff(rp,irp);
            
            // for all powers of r- appearing in the polynomial 'nu_rp_poly'
            for (int irm = nu_rp_poly.ldegree(rm); irm <= nu_rp_poly.degree(rm); irm++)
            {
                // extract factor in front of ν^inu rp^irp rm^irm
                GiNaC::ex nu_rp_rm_poly = nu_rp_poly.coeff(rm,irm);
                
                // for all powers of r- appearing in the polynomial 'nu_rp_poly'
                for (int irz = nu_rp_rm_poly.ldegree(rz); irz <= nu_rp_rm_poly.degree(rz); irz++)
                {
                    // extract factor in front of ν^inu rp^irp rm^irm rz^irz (it is a number)
                    GiNaC::ex nu_rp_rm_rz_coef = nu_rp_rm_poly.coeff(rz,irz);
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
    
    // resulting terms A^(-i/q-m) B^(-i/q-n) ν^a k₊^b k₋^c kz^d q^t q₊^u q₋^v qz^w
    std::map<std::tuple<int,int,int,int,int,int,int,int,int,int>,Complex> res;
    
    // for all terms
    for (auto term : terms)
    {
        // decode exponents
        int inu, irp, irm, irz;
        std::tie(inu, irp, irm, irz) = term.first;
        
        // resulting derivative for this particular term
        std::map<std::tuple<int,int,int,int,int,int,int,int,int,int>,Complex> rest, tmp;
        rest[std::make_tuple(0,1, 0,0,0,0, 0,0,0,0)] = term.second;
        
        // (inu + 1)-times differentiate the basic integral (-∂/∂ν)
        for (int j = 0; j <= inu; j++)
        {
            // move current result to temporary storage
            std::swap(rest,tmp); rest.clear();
            
            // differentiate all terms
            for (auto entry : tmp)
            {
                // decode entry
                int m,n,a,b,c,d,t,u,v,w; std::tie(m,n,a,b,c,d,t,u,v,w) = entry.first;
                
                // differentiate with respect to ν
                if (a != 0)
                {
                    auto key = std::make_tuple(m,n, a-1,b,c,d, t,u,v,w); auto it = rest.find(key);
                    rest[key] = -double(a) * entry.second + (it == rest.end() ? 0. : it->second);
                }
                {
                    auto key = std::make_tuple(m+1,n, a+1,b,c,d, t-1,u,v,w); auto it = rest.find(key);
                    rest[key] = Complex(0.,2.) * entry.second + (it == rest.end() ? 0. : it->second);
                }
                if (m != 0)
                {
                    auto key = std::make_tuple(m+1,n, a+1,b,c,d, t,u,v,w); auto it = rest.find(key);
                    rest[key] = (2.*m) * entry.second + (it == rest.end() ? 0. : it->second);
                }
                {
                    auto key = std::make_tuple(m+1,n, a,b,c,d, t,u,v,w); auto it = rest.find(key);
                    rest[key] = 2. * entry.second + (it == rest.end() ? 0. : it->second);
                }
                if (m != 0)
                {
                    auto key = std::make_tuple(m+1,n, a,b,c,d, t+1,u,v,w); auto it = rest.find(key);
                    rest[key] = Complex(0.,-2.*m) * entry.second + (it == rest.end() ? 0. : it->second);
                }
                {
                    auto key = std::make_tuple(m,n+1, a+1,b,c,d, t-1,u,v,w); auto it = rest.find(key);
                    rest[key] = Complex(0.,-2.) * entry.second + (it == rest.end() ? 0. : it->second);
                }
                if (n != 0)
                {
                    auto key = std::make_tuple(m,n+1, a+1,b,c,d, t,u,v,w); auto it = rest.find(key);
                    rest[key] = (2.*n) * entry.second + (it == rest.end() ? 0. : it->second);
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
                int m,n,a,b,c,d,t,u,v,w; std::tie(m,n,a,b,c,d,t,u,v,w) = entry.first;
                
                // differentiate with respect to k₋
                if (c != 0)
                {
                    auto key = std::make_tuple(m,n, a,b,c-1,d, t,u,v,w); auto it = rest.find(key);
                    rest[key] = Complex(0.,2.*c) * entry.second + (it == rest.end() ? 0. : it->second);
                }
                {
                    auto key = std::make_tuple(m+1,n, a,b+1,c,d, t-1,u,v,w); auto it = rest.find(key);
                    rest[key] = 2. * entry.second + (it == rest.end() ? 0. : it->second);
                }
                if (m != 0)
                {
                    auto key = std::make_tuple(m+1,n, a,b+1,c,d, t,u,v,w); auto it = rest.find(key);
                    rest[key] = Complex(0.,-2.*m) * entry.second + (it == rest.end() ? 0. : it->second);
                }
                {
                    auto key = std::make_tuple(m,n+1, a,b+1,c,d, t-1,u,v,w); auto it = rest.find(key);
                    rest[key] = -2. * entry.second + (it == rest.end() ? 0. : it->second);
                }
                {
                    auto key = std::make_tuple(m,n+1, a,b,c,d, t-1,u+1,v,w); auto it = rest.find(key);
                    rest[key] = Complex(0.,2.) * entry.second + (it == rest.end() ? 0. : it->second);
                }
                if (n != 0)
                {
                    auto key = std::make_tuple(m,n+1, a,b+1,c,d, t,u,v,w); auto it = rest.find(key);
                    rest[key] = Complex(0.,-2.*n) * entry.second + (it == rest.end() ? 0. : it->second);
                }
                if (n != 0)
                {
                    auto key = std::make_tuple(m,n+1, a,b,c,d, t,u+1,v,w); auto it = rest.find(key);
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
                int m,n,a,b,c,d,t,u,v,w; std::tie(m,n,a,b,c,d,t,u,v,w) = entry.first;
                
                // differentiate with respect to k₊
                if (b != 0)
                {
                    auto key = std::make_tuple(m,n, a,b-1,c,d, t,u,v,w); auto it = rest.find(key);
                    rest[key] = Complex(0.,2.*b) * entry.second + (it == rest.end() ? 0. : it->second);
                }
                {
                    auto key = std::make_tuple(m+1,n, a,b,c+1,d, t-1,u,v,w); auto it = rest.find(key);
                    rest[key] = 2. * entry.second + (it == rest.end() ? 0. : it->second);
                }
                if (m != 0)
                {
                    auto key = std::make_tuple(m+1,n, a,b,c+1,d, t,u,v,w); auto it = rest.find(key);
                    rest[key] = Complex(0.,-2.*m) * entry.second + (it == rest.end() ? 0. : it->second);
                }
                {
                    auto key = std::make_tuple(m,n+1, a,b,c+1,d, t-1,u,v,w); auto it = rest.find(key);
                    rest[key] = -2. * entry.second + (it == rest.end() ? 0. : it->second);
                }
                {
                    auto key = std::make_tuple(m,n+1, a,b,c,d, t-1,u,v+1,w); auto it = rest.find(key);
                    rest[key] = Complex(0.,2.) * entry.second + (it == rest.end() ? 0. : it->second);
                }
                if (n != 0)
                {
                    auto key = std::make_tuple(m,n+1, a,b,c+1,d, t,u,v,w); auto it = rest.find(key);
                    rest[key] = Complex(0.,-2.*n) * entry.second + (it == rest.end() ? 0. : it->second);
                }
                if (n != 0)
                {
                    auto key = std::make_tuple(m,n+1, a,b,c,d, t,u,v+1,w); auto it = rest.find(key);
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
                int m,n,a,b,c,d,t,u,v,w; std::tie(m,n,a,b,c,d,t,u,v,w) = entry.first;
                
                // differentiate with respect to kz
                if (d != 0)
                {
                    auto key = std::make_tuple(m,n, a,b,c,d-1, t,u,v,w); auto it = rest.find(key);
                    rest[key] = Complex(0.,d) * entry.second + (it == rest.end() ? 0. : it->second);
                }
                {
                    auto key = std::make_tuple(m+1,n, a,b,c,d+1, t-1,u,v,w); auto it = rest.find(key);
                    rest[key] = 2. * entry.second + (it == rest.end() ? 0. : it->second);
                }
                if (m != 0)
                {
                    auto key = std::make_tuple(m+1,n, a,b,c,d+1, t,u,v,w); auto it = rest.find(key);
                    rest[key] = Complex(0.,-2.*m) * entry.second + (it == rest.end() ? 0. : it->second);
                }
                {
                    auto key = std::make_tuple(m,n+1, a,b,c,d+1, t-1,u,v,w); auto it = rest.find(key);
                    rest[key] = -2. * entry.second + (it == rest.end() ? 0. : it->second);
                }
                {
                    auto key = std::make_tuple(m,n+1, a,b,c,d, t-1,u,v,w+1); auto it = rest.find(key);
                    rest[key] = Complex(0.,2.) * entry.second + (it == rest.end() ? 0. : it->second);
                }
                if (n != 0)
                {
                    auto key = std::make_tuple(m,n+1, a,b,c,d+1, t,u,v,w); auto it = rest.find(key);
                    rest[key] = Complex(0.,-2.*n) * entry.second + (it == rest.end() ? 0. : it->second);
                }
                if (n != 0)
                {
                    auto key = std::make_tuple(m,n+1, a,b,c,d, t,u,v,w+1); auto it = rest.find(key);
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
    
    return res;
}

Complex eval_W
(
    std::map<std::tuple<int,int,int,int,int,int,int,int,int,int>,Complex> const & poly,
    double nu,
    geom::vec3d const & vk,
    geom::vec3d const & vq
)
{
    // prepare some used variables
    double kxysqr = vk.x*vk.x + vk.y*vk.y, ksqr = kxysqr + vk.z*vk.z,
           kxy = std::sqrt(kxysqr), /*k = std::sqrt(ksqr),*/ kphi = std::atan2(vk.y,vk.x);
    double qxysqr = vq.x*vq.x + vq.y*vq.y, qsqr = qxysqr + vq.z*vq.z,
           qxy = std::sqrt(qxysqr), q = std::sqrt(qsqr), qphi = std::atan2(vq.y,vq.x);
    
    Complex A (ksqr + nu*nu - qsqr, -2*nu*q);
    double absA = std::hypot(A.real(),A.imag()), argA = std::atan2(A.imag(),A.real());
    double B = nu*nu + (vk.x+vq.x)*(vk.x+vq.x) + (vk.y+vq.y)*(vk.y+vq.y) + (vk.z+vq.z)*(vk.z+vq.z);
    double log_abs_A_B_over_q = std::log(absA / B) / q;
    
    // evaluate all terms A^(-m) B^(-n) ν^a k₊^b k₋^c kz^d q^t q₊^u q₋^v qz^w
    Complex result = 0;
    for (auto term : poly)
    {
        // extract exponents
        int m,n,a,b,c,d,t,u,v,w; std::tie(m,n,a,b,c,d,t,u,v,w) = term.first;
        
        // contribution from this term
        Complex contrib = term.second;
        
        // add A^(-m)
        if (m != 0)
            contrib *= Complex(std::cos(m*argA),-std::sin(m*argA)) / gsl_sf_pow_int(absA,m);
        
        // add B^(-n)
        if (n != 0)
            contrib *= 1./gsl_sf_pow_int(B,n);
        
        // add ν^a
        if (a != 0)
            contrib *= gsl_sf_pow_int(nu,a);
        
        // add k₊^b k₋^c
        if (b != 0 or c != 0)
            contrib *= Complex(std::cos((b-c)*kphi),std::sin((b-c)*kphi)) * gsl_sf_pow_int(kxy,b+c);
        
        // add kz^d
        if (d != 0)
            contrib *= gsl_sf_pow_int(vk.z,d);
        
        // add q^t
        if (t != 0)
            contrib *= gsl_sf_pow_int(q,t);
        
        // add q₊^u q₋^v
        if (u != 0 or v != 0)
            contrib *= Complex(std::cos((u-v)*qphi),std::sin((u-v)*qphi)) * gsl_sf_pow_int(qxy,u+v);
        
        // add qz^w
        if (w != 0)
            contrib *= gsl_sf_pow_int(vq.z,w);
        
        // update sum
        result += contrib;
    }
    
    // add missing factors
//     result *= std::pow(A/B,Complex(0.,-1./q)) / ksqr;
    result *= Complex(std::cos(log_abs_A_B_over_q),-std::sin(log_abs_A_B_over_q)) * std::exp(argA/q) / ksqr;
    
    return result;
}
