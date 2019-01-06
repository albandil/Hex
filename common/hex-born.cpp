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

#include <cmath>
#include <map>

// --------------------------------------------------------------------------------- //

#ifdef WITH_GINAC
#include <ginac/ginac.h>
#endif

// --------------------------------------------------------------------------------- //

#include <gsl/gsl_sf.h>

// --------------------------------------------------------------------------------- //

#include "hex-gausskronrod.h"
#include "hex-hydrogen.h"
#include "hex-numbers.h"
#include "hex-special.h"
#include "hex-symbolic.h"
#include "hex-vec3d.h"

// --------------------------------------------------------------------------------- //

#ifdef WITH_GINAC
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
#endif

#ifdef WITH_GINAC
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
#endif

#ifdef WITH_GINAC
std::map<std::tuple<int,int,int,int,int>,Complex> Wb_symb_in
(
    int n1, int l1, int m1,
    int n2, int l2, int m2
)
{
    // This routine will compute integral
    //                    -ik·r
    //   W = Int ψf(r)* (e      - 1) ψi(r) d³r
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

    // known rule: "rp" and "rm" are complex conjugates of each other
    GiNaC::exmap rule = {{ rp.conjugate(), rm },{ rm.conjugate(), rp }};
    poly = poly.subs(rule);

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
#endif

#ifdef WITH_GINAC
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
#endif

#ifdef WITH_GINAC
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
#endif

#ifdef WITH_GINAC
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
#endif


template <class Functor, class Integrator>
class BesselNodeIntegrator1D
{
    private:

        Integrator Q_;

        double k_;
        int l_;

        double result_;
        bool ok_;
        std::string status_;

        double epsabs_;
        double epsrel_;
        int limit_;


    public:

        BesselNodeIntegrator1D (Functor f, double k, int l)
            : Q_(f), k_(k), l_(l), result_(0), ok_(true), status_(),
              epsabs_(1e-8), epsrel_(1e-5), limit_(100)
        {
        }

        double result () const { return result_; }
        bool ok () const { return ok_; }
        std::string const & status () const { return status_; }

        double epsabs () const { return epsabs_; }
        void setEpsAbs (double eps) { epsabs_ = eps; }

        double epsrel () const { return epsrel_; }
        void setEpsRel (double eps) { epsrel_ = eps; }

        int limit () const { return limit_; }
        void setLimit (double n) { limit_ = n; }

        bool integrate (double a, double b)
        {
            // set parcel integrator accuracy to one order more
            Q_.setEpsAbs(0.1 * epsabs_);
            Q_.setEpsRel(0.1 * epsrel_);

            // overall integral
            double integral = 0;

            // end of previous integration parcel
            double prevR = 0;

            // for all integration parcels (nodes of the Bessel function)
            for (int inode = 1; inode < limit_; inode++)
            {
                // get integration bounds
                gsl_sf_result res;
                int err = gsl_sf_bessel_zero_Jnu_e(l_+0.5,inode,&res);
                if (err != GSL_SUCCESS)
                {
                    throw exception
                    (
                        "Cannot find %d-th root of the Bessel function j[%d](%g r) -- %s.",
                        inode, l_, k_, gsl_strerror(err)
                    );
                }
                double rmin = prevR;
                double rmax = res.val / k_;

                // skip intervals that are below the lower limit
                if (rmax < a)
                    continue;

                // skip intervals that are above the upper limit
                if (rmin > b)
                    break;

                // shrink integration interval, if necessary
                rmin = std::max (a, rmin);
                rmax = std::min (rmax, b);

                // integrate and check success
                if (not Q_.integrate(rmin,rmax))
                {
                    ok_ = false;
                    status_ = Q_.status();
                    result_ = integral;
                    return ok_;
                }

                // update result
                integral += Q_.result();
                prevR = rmax;

                // check convergence
                if (std::abs(Q_.result()) < epsabs_ or std::abs(Q_.result()) < epsrel_ * std::abs(integral))
                    break;
            }

            ok_ = true;
            result_ = integral;
            status_ = "";
            return ok_;
        }
};

#if defined(WITH_CLN) || defined(WITH_GINAC)
double compute_Idir (int li, int lf, int lambda, int Ni, int Li, double ki, int Nf, int Lf, double kf)
{
//     std::cout << format
//     (
//         "Precompute Idir\n"
//         "\tlambda = %d\n"
//         "\tNi = %d, Li = %d, ki = %g, li = %d\n"
//         "\tNf = %d, Lf = %d, kf = %g, lf = %d\n",
//         lambda, Ni, Li, ki, li, Nf, Lf, kf, lf
//     );

    if (lambda == 0)
    {
        //
        // r1 > r2
        //

        auto integrand = [Ni,Li,Nf,Lf,li,ki,lf,kf](double r2) -> double
        {
            // inner integrand
            auto iintegrand = [Ni,Li,Nf,Lf,r2](double r1) -> double { return Hydrogen::P(Ni,Li,r1) * Hydrogen::P(Nf,Lf,r1) * (1./r1 - 1./r2); };

            // inner integrator
            GaussKronrod<decltype(iintegrand)> Qi(iintegrand);
            Qi.setEpsAbs(0);

            // integrate and check success
            if (not Qi.integrate(r2,special::constant::Inf))
            {
                throw exception
                (
                    "compute_Idir (inner) failed for λ=0, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d, r2=%g (\"%s\").\n\tresult = %g\n",
                    Ni, Li, ki, li, Nf, Lf, kf, lf, r2, Qi.status().c_str(), Qi.result()
                );
            }

            return Qi.result() * special::ric_j(li,ki*r2) * special::ric_j(lf,kf*r2);
        };

        // which Bessel function oscillates slowlier ?
        int    l = (ki < kf ? li : lf);
        double k = (ki < kf ? ki : kf);

        // outer integrator
        BesselNodeIntegrator1D<decltype(integrand),GaussKronrod<decltype(integrand)>> R(integrand, k, l);
        R.setEpsAbs(0);
        R.integrate(0,special::constant::Inf);

//         std::cout << "\tIdir = " << R.result() << std::endl;
        return R.result();
    }
    else
    {
        // compute normalization factors of the generalized Laguerre polynomials
        double Normi = std::sqrt(std::pow(2./Ni,3) * gsl_sf_fact(Ni-Li-1) / (2. * Ni * gsl_sf_fact(Ni+Li)));
        double Normf = std::sqrt(std::pow(2./Nf,3) * gsl_sf_fact(Nf-Lf-1) / (2. * Nf * gsl_sf_fact(Nf+Lf)));

        // compute the polynomials
        symbolic::poly Lagi = symbolic::GeneralizedLaguerre (Ni-Li-1, 2*Li+1);
        symbolic::poly Lagf = symbolic::GeneralizedLaguerre (Nf-Lf-1, 2*Lf+1);

        // multiply by the angular factor
        for (symbolic::term & pi : Lagi) pi.a += Li + 1;
        for (symbolic::term & pf : Lagf) pf.a += Lf + 1;

        // part of the angular factor goes to normalization
        Normi *= std::pow(2./Ni,Li);
        Normf *= std::pow(2./Nf,Lf);

        // compute the product of the polynomials
        symbolic::poly PP = Lagi * Lagf;

        // factor in the argument of the exponential
        double c = 1./Ni + 1./Nf;

        // outer integrand evaluated at "r"
        auto integrand = [PP,c,lambda,Normi,Normf,ki,kf,li,lf](double r) -> double
        {
            // inner integral
            double integral = 0;

            // integrate term by term
            for (symbolic::term const & p : PP)
            {
                // compute the high integral
                gsl_sf_result res;
                int err_high = gsl_sf_gamma_inc_e (p.a - lambda, c * r, &res);
                if (err_high != GSL_SUCCESS and err_high != GSL_EUNDRFLW)
                {
                    throw exception
                    (
                        "Unable to evaluate incomplete gamma-function Gamma(%d,%g) - %s.",
                        p.a - lambda, c * r, gsl_strerror(err_high)
                    );
                }
                double int_high = 0;
                if (err_high != GSL_EUNDRFLW)
                    int_high = gsl_sf_pow_int(c*r,lambda) * res.val;

                // compute the low integral
                int err_low = gsl_sf_gamma_inc_P_e (p.a + lambda + 1, c * r, &res);
                double scale = gsl_sf_gamma (p.a + lambda + 1);
                if (err_low != GSL_SUCCESS)
                {
                    throw exception
                    (
                        "Unable to evaluate scaled complementary incomplete gamma-function P(%d,%g) - %s.",
                        p.a + lambda + 1, c * r, gsl_strerror(err_low)
                    );
                }
                double int_low = gsl_sf_pow_int(c*r,-lambda-1) * res.val * scale;

                // sum both contributions
                integral += (int_low + int_high) * symbolic::double_approx(p.kr) / gsl_sf_pow_int(c,p.a);
            }

            return Normi * Normf * integral * special::ric_j(li,ki*r) * special::ric_j(lf,kf*r);
        };

        // which Bessel function oscillates slowlier ?
        int    l = std::max(li, lf);
        double k = (ki == kf ? ki : std::abs(ki - kf));

        // outer integrator
        BesselNodeIntegrator1D<decltype(integrand),GaussKronrod<decltype(integrand)>> R(integrand, k, l);
        R.integrate (0,special::constant::Inf);
//         std::cout << "\tIdir = " << R.result() << std::endl;

        return R.result();
    }
}

double compute_Iexc (int li, int lf, int lambda, int Ni, int Li, double ki, int Nf, int Lf, double kf)
{
//     std::cout << format
//     (
//         "Precompute Iexc\n"
//         "\tlambda = %d\n"
//         "\tNi = %d, Li = %d, ki = %g, li = %d\n"
//         "\tNf = %d, Lf = %d, kf = %g, lf = %d\n",
//         lambda, Ni, Li, ki, li, Nf, Lf, kf, lf
//     );

    if (lambda == 0)
    {
        //
        // r1 > r2
        //

        auto integrand = [Ni,Li,Nf,Lf,li,ki,lf,kf](double r2) -> double
        {
            // inner integrand
            auto iintegrand = [Ni,Li,kf,lf,r2](double r1) -> double
            {
                return Hydrogen::P(Ni,Li,r1) * special::ric_j(lf,kf*r1) * (1./r1 - 1./r2);
            };

            // inner integrator
            GaussKronrod<decltype(iintegrand)> Qi(iintegrand);

            // integrate and check success
            if (not Qi.integrate(r2,special::constant::Inf))
            {
                throw exception
                (
                    "compute_Iexc (inner) failed for λ=0, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d, r2=%g (\"%s\").\n\tresult = %g\n",
                    Ni, Li, ki, li, Nf, Lf, kf, lf, r2, Qi.status().c_str(), Qi.result()
                );
            }

            return Qi.result() * special::ric_j(li,ki*r2) * Hydrogen::P(Nf,Lf,r2);
        };

        // outer integrator
        GaussKronrod<decltype(integrand)> Q(integrand);

        // integrate and check success
        if (not Q.integrate(0.,special::constant::Inf))
        {
            throw exception
            (
                "compute_Iexc (outer) failed for λ=0, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d (\"%s\").\n\tresult = %g",
                Ni, Li, ki, li, Nf, Lf, kf, lf, Q.status().c_str(), Q.result()
            );
        }

//         std::cout << "\tIexc = " << Q.result() << std::endl;
        return Q.result();
    }
    else
    {
        // r1 < r2
        auto integrand1 = [Ni,Li,Nf,Lf,li,ki,lf,kf,lambda](double r2) -> double
        {
            // inner integrand
            auto iintegrand1 = [Ni,Li,kf,lf,r2,lambda](double r1) -> double
            {
                return Hydrogen::P(Ni,Li,r1) * special::ric_j(lf,kf*r1) * std::pow(r1/r2,lambda);
            };

            // inner integrator
            BesselNodeIntegrator1D<decltype(iintegrand1),GaussKronrod<decltype(iintegrand1)>> Qi1 (iintegrand1, kf, lf);

            // integrate and check success
            if (not Qi1.integrate(0.,r2))
            {
                throw exception
                (
                    "compute_Iexc (inner1) failed for λ=%d, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d, r2=%g (\"%s\").\n\tresult = %g\n",
                    lambda, Ni, Li, ki, li, Nf, Lf, kf, lf, r2, Qi1.status().c_str(), Qi1.result()
                );
            }

            return Qi1.result() / r2 * special::ric_j(li,ki*r2) * Hydrogen::P(Nf,Lf,r2);
        };

        BesselNodeIntegrator1D<decltype(integrand1),GaussKronrod<decltype(integrand1)>> Q1 (integrand1, ki, li);
        if (not Q1.integrate(0.,special::constant::Inf))
        {
            throw exception
            (
                "compute_Iexc (outer1) failed for λ=%d, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d (\"%s\").\n\tresult = %g",
                lambda, Ni, Li, ki, li, Nf, Lf, kf, lf, Q1.status().c_str(), Q1.result()
            );
        }

        // r2 < r1
        auto integrand2 = [Ni,Li,Nf,Lf,li,ki,lf,kf,lambda](double r1) -> double
        {
            // inner integrand
            auto iintegrand2 = [Nf,Lf,ki,li,r1,lambda](double r2) -> double
            {
                return Hydrogen::P(Nf,Lf,r2) * special::ric_j(li,ki*r2) * std::pow(r2/r1,lambda);
            };

            // inner integrator
            BesselNodeIntegrator1D<decltype(iintegrand2),GaussKronrod<decltype(iintegrand2)>> Qi2 (iintegrand2, ki, li);

            // integrate and check success
            if (not Qi2.integrate(0.,r1))
            {
                throw exception
                (
                    "compute_Iexc (inner2) failed for λ=%d, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d, r2=%g (\"%s\").\n\tresult = %g\n",
                    lambda, Ni, Li, ki, li, Nf, Lf, kf, lf, r1, Qi2.status().c_str(), Qi2.result()
                );
            }

            return Qi2.result() / r1 * special::ric_j(lf,kf*r1) * Hydrogen::P(Ni,Li,r1);
        };

        BesselNodeIntegrator1D<decltype(integrand2),GaussKronrod<decltype(integrand2)>> Q2 (integrand2, kf, lf);
        if (not Q2.integrate(0.,special::constant::Inf))
        {
            throw exception
            (
                "compute_Iexc (outer2) failed for λ=%d, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d (\"%s\").\n\tresult = %g",
                lambda, Ni, Li, ki, li, Nf, Lf, kf, lf, Q2.status().c_str(), Q2.result()
            );
        }

//         std::cout << "\tIexc = " << Q1.result() + Q2.result() << std::endl;
        return Q1.result() + Q2.result();
    }
}

void pwba
(
    int Ni, int Li, double ki,
    int Nf, int Lf, double kf,
    int L,
    cArrays & Tdir, cArrays & Texc,
    bool direct, bool exchange
)
{
    // re-allocate memory
    Tdir = cArrays((2*Li+1)*(2*Lf+1), cArray());
    Texc = cArrays((2*Li+1)*(2*Lf+1), cArray());

    // for all outgoing partial waves
    for (int lf = std::abs(Lf - L); lf <= Lf + L; lf++)
    {
        // add new T-matrix for this outgoing partial wave
        for (cArray & T : Tdir)
            T.push_back(0.);
        for (cArray & T : Texc)
            T.push_back(0.);

        // for all incoming partial waves
        for (int li = std::abs(Li - L); li <= Li + L; li++)
        {
            // conserve parity
            if ((lf + Lf) % 2 != (li + Li) % 2)
                continue;

            //
            // compute direct contribution
            //

            // for all multipoles
            for (int lam = std::max(std::abs(Lf-Li), std::abs(lf-li)); direct and lam <= std::min(Li+Lf, lf+li); lam++)
            {
                // compute the needed radial integrals
                double Vdir = compute_Idir(li, lf, lam, Ni, Li, ki, Nf, Lf, kf);

                // compute complex prefactor
                Complex prefactor = std::pow(4*special::constant::pi,2)/(ki*kf) * std::pow(Complex(0.,1.),li-lf) * std::sqrt((2*li+1)/(4*special::constant::pi));

                // for all projections of the initial/final angular momentum
                for (int Mi = -Li; Mi <= Li; Mi++)
                for (int Mf = -Lf; Mf <= Lf; Mf++)
                {
                    // compute index in the array of T-matrices
                    int idx = (Mi + Li)*(2*Lf + 1) + Mf + Lf;

                    // compute angular integrals (Gaunt coefficients)
                    double ang = special::ClebschGordan(Lf,Mf,lf,Mi-Mf,L,Mi)
                               * special::ClebschGordan(Li,Mi,li,0,L,Mi)
                               * special::computef(lam,Lf,lf,Li,li,L);

                    // add the T-matrix contributions
                    Tdir[idx][lf-std::abs(Lf - L)] += prefactor * ang * Vdir;
                }
            }

            //
            // compute exchange contribution
            //

            // for all multipoles
            for (int lam = std::max(std::abs(lf-Li), std::abs(Lf-li)); exchange and lam <= std::min(lf+Li, Lf+li); lam++)
            {
                // compute the needed radial integrals
                double Vexc = compute_Iexc (li, lf, lam, Ni, Li, ki, Nf, Lf, kf);

                // compute complex prefactor
                Complex prefactor = std::pow(4*special::constant::pi,2)/(ki*kf) * std::pow(Complex(0.,1.),li-lf) * std::sqrt((2*li+1)/(4*special::constant::pi));

                // for all projections of the initial/final angular momentum
                for (int Mi = -Li; Mi <= Li; Mi++)
                for (int Mf = -Lf; Mf <= Lf; Mf++)
                {
                    // compute index in the array of T-matrices
                    int idx = (Mi + Li)*(2*Lf + 1) + Mf + Lf;

                    // compute angular integrals (Gaunt coefficients)
                    double ang = special::ClebschGordan(Lf,Mf,lf,Mi-Mf,L,Mi)
                               * special::ClebschGordan(Li,Mi,li,0,L,Mi)
                               * special::computef(lam,lf,Lf,Li,li,L);

                    // add the T-matrix contributions
                    Texc[idx][lf-std::abs(Lf - L)] += prefactor * ang * Vexc;
                }
            }
        }
    }
}
#endif
