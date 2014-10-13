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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>
#include <vector>

#include <cln/cln.h>
#include <gsl/gsl_sf.h>

#include "complex.h"
#include "misc.h"
#include "special.h"
#include "symbolic.h"

symbolic::rational symbolic::onehalf = cln::cl_RA(1)/cln::cl_I(2);

symbolic::rational symbolic::combination (symbolic::rational const & alpha, int k)
{
    symbolic::rational res = 1;
    for (int w = 0; w < k; w++)
        res *= (alpha - w) / (k - w);
    return res;
}

const char* gfname (unsigned i)
{
    static const char * const names[] = {"none", "sin", "cos"};
    return names[i];
}

symbolic::term symbolic::operator + (symbolic::term const & A, symbolic::term const & B)
{
    if ((A.gf != B.gf) or (A.a != B.a) or (A.b != B.b) or (A.c != B.c))
        throw "Can't add non-simillar symbolic terms.";
    
    symbolic::term S;
    S.a = A.a;
    S.b = A.b;
    S.c = A.c;
    S.gf = A.gf;
    
    if (A.ki == B.ki)
    {
        S.ki = A.ki;
        S.kr = A.kr + B.kr;
    }
    else
    {
        S.ki = A.ki * cln::double_approx(A.kr) + B.ki * cln::double_approx(B.kr);
        S.kr = 1;
        
        // leave ki positive
        if (S.ki < 0)
        {
            S.ki = -S.ki;
            S.kr = -S.kr;
        }
    }
    
    return S;
}

symbolic::term symbolic::operator / (symbolic::term const & A, symbolic::term const & B) throw (exception)
{
    symbolic::term S = A;
    
    if (A.gf == B.gf and A.b == B.b)
    {
        S.gf = GF_NONE;
    }
    else if (B.gf != GF_NONE)
    {
        throw exception
        (
            "Can't divide non-similar symbolic terms %s and %s.",
            symbolic::tostring(A).c_str(),
            symbolic::tostring(B).c_str()
        );
    }
    
    S.ki /= B.ki;
    S.kr /= B.kr;
    S.a -= B.a;
    S.c -= B.c;
    
    return S;
}

symbolic::term symbolic::expm (symbolic::rational c)
{
    symbolic::term p;
    
    p.ki = 1.;
    p.kr = 1;
    p.a = 0;
    p.gf = GF_NONE;
    p.b = 0.;
    p.c = c;
    
    return p;
}

symbolic::poly symbolic::operator * (symbolic::poly const & P, symbolic::poly const & Q)
{
    symbolic::poly product_poly;
    
    for (symbolic::term const & p : P) for (symbolic::term const & q : Q)
    {
        symbolic::term product_term;
        product_term.ki = p.ki * q.ki;
        product_term.kr = p.kr * q.kr;
        product_term.a = p.a + q.a;
        product_term.c = p.c + q.c;
        
        if (p.gf == GF_NONE and q.gf == GF_NONE)
        {
            product_term.gf = GF_NONE;
            product_term.b = 0;
            
            product_poly.push_back(product_term);
        }
        else if (p.gf != GF_NONE and q.gf == GF_NONE)
        {
            product_term.gf = p.gf;
            product_term.b = p.b;
            
            product_poly.push_back(product_term);
        }
        else if (p.gf == GF_NONE and q.gf != GF_NONE)
        {
            product_term.gf = q.gf;
            product_term.b = q.b;
            
            product_poly.push_back(product_term);
        }
        else if (p.gf == GF_SIN and q.gf == GF_SIN)
        {
            product_term.kr *= symbolic::onehalf;
            product_term.gf = GF_COS;
            
            product_term.b = p.b - q.b;
            product_poly.push_back(product_term);
            if (product_poly.back().b == 0)
                product_poly.back().gf = GF_NONE;
            
            product_term.kr = -product_term.kr;
            product_term.b = p.b + q.b;
            product_poly.push_back(product_term);
            if (product_poly.back().b == 0)
                product_poly.back().gf = GF_NONE;
        }
        else if (p.gf == GF_SIN and q.gf == GF_COS)
        {
            product_term.kr *= symbolic::onehalf;
            product_term.gf = GF_SIN;
            
            product_term.b = p.b + q.b;
            product_poly.push_back(product_term);
            
            product_term.b = p.b - q.b;
            product_poly.push_back(product_term);
        }
        else if (p.gf == GF_COS and q.gf == GF_SIN)
        {
            product_term.kr *= symbolic::onehalf;
            product_term.gf = GF_SIN;
            
            product_term.b = p.b + q.b;
            product_poly.push_back(product_term);
            
            product_term.kr = -product_term.kr;
            product_term.b = p.b - q.b;
            product_poly.push_back(product_term);
        }
        else if (p.gf == GF_COS and q.gf == GF_COS)
        {
            product_term.kr *= symbolic::onehalf;
            product_term.gf = GF_COS;
            
            product_term.b = p.b + q.b;
            product_poly.push_back(product_term);
            if (product_poly.back().b == 0)
                product_poly.back().gf = GF_NONE;
            
            product_term.b = p.b - q.b;
            product_poly.push_back(product_term);
            if (product_poly.back().b == 0)
                product_poly.back().gf = GF_NONE;
        }
        else
        {
            throw exception ("[operator *] Unknown goniometric function combination.");
        }
    }
    
    product_poly.optimize();
    return product_poly;
}

symbolic::poly symbolic::operator - (symbolic::poly const & P, symbolic::poly const & Q)
{
    // copy P
    symbolic::poly difference = P;
    
    // copy Q
    difference.insert(difference.end(), Q.begin(), Q.end());
    
    // invert signs in the copy of Q
    for (unsigned i = P.size(); i < difference.size(); i++)
        difference[i].kr = -difference[i].kr;
    
    // optimize and return
    difference.optimize();
    return difference;
}

symbolic::poly symbolic::operator + (symbolic::poly const & P, symbolic::poly const & Q)
{
    // copy P
    symbolic::poly sum = P;
    
    // copy Q
    sum.insert(sum.end(), Q.begin(), Q.end());
    
    // optimize and return
    sum.optimize();
    return sum;
}

symbolic::poly symbolic::operator / (symbolic::poly const & P, symbolic::rational const & q)
{
    // copy P
    symbolic::poly result = P;
    
    // divide by rational quotient
    for (symbolic::term & p : result)
        p.kr /= q;
    
    return result;
}

symbolic::poly symbolic::AssociatedLaguerre (int k, int s)
{
    symbolic::poly laguerre;
    
    for (int j = s; j <= k; j++)
    {
        symbolic::term term;
        term.ki = 1.;
        term.kr = cln::factorial(k) * cln::factorial(k) / (cln::factorial(k-j) * cln::factorial(j) * cln::factorial(j-s));
        term.a = j-s;
        term.gf = GF_NONE;
        term.b = 0;
        term.c = 0;
        
        if (j % 2 == 0)
            term.kr = -term.kr;
        
        laguerre.push_back(term);
    }
    
    laguerre.optimize();
    return laguerre;
}

symbolic::poly symbolic::GeneralizedLaguerre (int n, int a)
{
    symbolic::poly asslaguerre;
    
    for (int j = 0; j <= n; j++)
    {
        symbolic::term term;
        term.ki = 1.;
        term.kr = cln::binomial(n + a, n - j) / cln::factorial(j);
        term.a = j;
        term.gf = GF_NONE;
        term.b = 0;
        term.c = 0;
        
        if (j % 2 != 0)
            term.kr = -term.kr;
        
        asslaguerre.push_back(term);
    }
    
    asslaguerre.optimize();
    return asslaguerre;
}

symbolic::term symbolic::HydrogenN (int n, int l)
{
    // compute square of the normalization factor
    symbolic::rational N2 = cln::expt(symbolic::rational(2)/n,3) * cln::factorial(n-l-1) / (2*n*cln::expt(cln::factorial(n+l),3));
    
    // the factor itself
    symbolic::term N;
    
    // try to compute a precise (rational) square root and store it in N.kr
    if (not cln::sqrtp(N2, &N.kr))
    {
        N.kr = 1;
        N.ki = sqrt(cln::double_approx(N2));
    }
    
    return N;
}

symbolic::poly symbolic::HydrogenNP (int n, int l)
{
    poly laguerre = symbolic::AssociatedLaguerre(n+l, 2*l+1);
    for (term & term : laguerre)
        term.kr *= cln::expt(symbolic::rational(2)/n, term.a);
    
    term rest;
    rest.kr = cln::expt(symbolic::rational(2)/n,l);
    rest.a = l+1;
    rest.gf = GF_NONE;
    rest.b = 0;
    rest.c = symbolic::rational(1)/n;
    
    return laguerre * rest;
}

symbolic::poly symbolic::HydrogenP (int n, int l)
{
    return symbolic::HydrogenN(n,l) * symbolic::HydrogenNP(n,l);
}

symbolic::poly symbolic::HydrogenS (int n, int l, symbolic::rational lambda)
{
    symbolic::poly laguerre = symbolic::AssociatedLaguerre(n+l, 2*l+1);
    for (symbolic::term & term : laguerre)
        term.kr *= cln::expt(2*lambda, term.a);
    
    symbolic::term rest;
    rest.kr = n * n * cln::expt(2*lambda,l);
    rest.a = l;
    rest.gf = GF_NONE;
    rest.b = 0;
    rest.c = 1;
    
    return HydrogenN(n,l) * laguerre * rest;
}

symbolic::poly symbolic::LaguerreBasisFunction (int N, int L, symbolic::rational lambda)
{
    poly laguerre = symbolic::GeneralizedLaguerre(N-1, 2*L+2);
    for (term & term : laguerre)
        term.kr *= cln::expt(lambda, term.a);
    
    symbolic::term rest;
    rest.ki = 1.;
    rest.kr = cln::expt(lambda,L+1);
    rest.a = L+1;
    rest.gf = GF_NONE;
    rest.b = 0;
    rest.c = lambda/2;
    
    return laguerre * rest;
}

symbolic::rational symbolic::LaguerreBasisFunctionNsqr (int N, int L, symbolic::rational lambda)
{
    return lambda * cln::factorial(N-1) / cln::factorial(N+2*L+1);
}

/*
 * Riccati-Bessel function
 * BesselF(l,x) := sum(
 *                       sin(x) * (2/x)^(l-2*i) * (l-i)!/i! * combination(-1/2-i,l-2*i) 
 *                     - cos(x) * (2/x)^(l-2*i) * (l-i)!/i! * i * combination(-1/2-i,l-2*i+1),
 *                   i, 0, floor((l+1)/2)
 *                 );
 */
symbolic::poly symbolic::RiccatiBessel (int l, double k)
{
    symbolic::poly ricj;
    
    for (int i = 0; i <= (l+1)/2; i++)
    {
        symbolic::term term;
        
        symbolic::rational comb = symbolic::combination (-symbolic::rational(1)/2-i,l-2*i);
        
        term.ki = pow(2/k, 1-2*i);
        term.kr = cln::factorial(l-i) / cln::factorial(i) * comb;
        term.a = 2*i - l;
        term.gf = GF_SIN;
        term.b = k;
        term.c = 0;
        ricj.push_back(term);
    }
    
    for (int i = 1; i <= (l+1)/2; i++)
    {
        symbolic::term term;
        
        symbolic::rational comb = symbolic::combination (-symbolic::rational(1)/2-i,l-2*i+1);
        
        term.ki = pow(2/k, 1-2*i);
        term.kr = -cln::factorial(l-i) / cln::factorial(i) * i *  comb;
        term.a = 2*i - l;
        term.gf = GF_COS;
        term.b = k;
        term.c = 0;
        ricj.push_back(term);
    }
    
    ricj.optimize();
    return ricj;
}

symbolic::term symbolic::integrate_full (poly const & P)
{
    // the resulting number
    symbolic::term result;
    result.ki = 1;
    result.kr = 0;
    
    // excluded terms
    std::set<int> exclude;
    
    // for all terms of the polynomial
    for (unsigned i = 0; i < P.size(); i++)
    {
        // skip excluded terms (those have been already integrated)
        if (exclude.find(i) != exclude.end())
            continue;
        
        // get term reference
        symbolic::term const & p = P[i];
        
        // check if it has integrable sine-singularity
        if (p.a == -1 and p.gf == GF_SIN)
        {
            result.ki = result.ki * cln::double_approx(result.kr) + p.ki * cln::double_approx(p.kr) * (special::constant::pi_half - atan(cln::double_approx(p.c)/p.b));
            result.kr = 1;
            
            // leave ki positive
            if (result.ki < 0)
            {
                result.ki = -result.ki;
                result.kr = -result.kr;
            }
            
            continue;
        }
        
        // check if it has integrable cosine-singularity
        if (p.a == -1 and p.gf == GF_COS)
        {
            // Integral of this term on its own is divergent. But, when combined with a compatible term that contains
            // no goniometric function and has the opposite sign, the overall integral of those two terms is finite
            // due to compensation near zero: (cos(x) - 1) -> O(x²) -> 0.
            
            // find the compatible term; it must be further due to descending ordering by goniometric function
            unsigned icompat = 0;
            for (unsigned j = i + 1; j < P.size(); j++)
            {
                if (P[j].a == p.a and P[j].c == p.c and P[j].gf == GF_NONE and P[j].ki == p.ki and P[j].kr == -p.kr)
                {
                    icompat = j;
                    break;
                }
            }
            
            // if there is a compatible term
            if (icompat > i)
            {
                // add the other term to already used terms
                exclude.insert(icompat);
                
                // add the integral
                double bc = p.b / cln::double_approx(p.c);
                result.ki = result.ki * cln::double_approx(result.kr) + p.ki * cln::double_approx(p.kr) * 0.5 * log(1. + bc*bc);
                result.kr = 1;
                
                // leave ki positive
                if (result.ki < 0)
                {
                    result.ki = -result.ki;
                    result.kr = -result.kr;
                }
            }
            
            continue;
        }
        
        // check if it is convergent at all
        if (p.a < 0)
        {
            throw exception
            (
                "[symbolic::integrate_full] Integral of \"%s\" is divergent.\nThe full integrand: %s",
                tostring(p).c_str(),
                tostring(P).c_str()
            );
        }
        
        // computed integral of a term
        symbolic::term integ;
        
        // branch according to the goniometric function
        if (p.gf == GF_NONE)
        {
            integ.ki = p.ki;
            integ.kr = p.kr * cln::factorial(p.a) / cln::expt(p.c, p.a + 1);
        }
        else if (p.gf == GF_SIN)
        {
            integ.ki = p.ki * gsl_sf_fact(p.a) * std::pow(Complex(cln::double_approx(p.c),-p.b), -p.a-1).imag();
            integ.kr = p.kr;
        }
        else if (p.gf == GF_COS)
        {
            integ.ki = p.ki * gsl_sf_fact(p.a) * std::pow(Complex(cln::double_approx(p.c),-p.b), -p.a-1).real();
            integ.kr = p.kr;
        }
        else
        {
            throw exception ("[symbolic::integrate_full] Unknown goniom. function.");
        }
        
        // add to result
        result = result + integ;
    }
    
    return result;
}

symbolic::poly symbolic::integrate_inf (symbolic::poly const & P)
{
    symbolic::poly result;
    
    // for all terms in the poly
    for (symbolic::term const & p : P)
    {
        // a!/j! = a(a-1)(a-2)...(j+2)(j+1)
        symbolic::rational afac_over_jfac = 1;
        
        // for all terms in the sum
        for (int j = p.a; j >= 0; j--)
        {
            if (p.gf == GF_NONE)
            {
                symbolic::term term;
                term.ki = p.ki;
                term.kr = p.kr * afac_over_jfac * cln::expt(p.c, j - p.a - 1);
                term.a = j;
                term.gf = GF_NONE;
                term.b = 0.;
                term.c = p.c;
                result.push_back(term);
            }
            else if (p.gf == GF_COS)
            {
                symbolic::term term1, term2;
                term1.a = term2.a = j;
                term1.b = term2.b = p.b;
                term1.c = term2.c = p.c;
                term1.kr = term2.kr = p.kr;
                
                term1.gf = GF_COS;
                term1.ki = p.ki * cln::double_approx(afac_over_jfac) * std::pow(Complex(cln::double_approx(p.c),-p.b), j-p.a-1).real();
                result.push_back(term1);
                
                term2.gf = GF_SIN;
                term2.ki = -p.ki * cln::double_approx(afac_over_jfac) * std::pow(Complex(cln::double_approx(p.c),-p.b), j-p.a-1).imag();
                result.push_back(term2);
            }
            else if (p.gf == GF_SIN)
            {
                symbolic::term term1, term2;
                term1.a = term2.a = j;
                term1.b = term2.b = p.b;
                term1.c = term2.c = p.c;
                term1.kr = term2.kr = p.kr;
                
                term1.gf = GF_COS;
                term1.ki = p.ki * cln::double_approx(afac_over_jfac) * std::pow(Complex(cln::double_approx(p.c),-p.b), j-p.a-1).imag();
                result.push_back(term1);
                
                term2.gf = GF_SIN;
                term2.ki = p.ki * cln::double_approx(afac_over_jfac) * std::pow(Complex(cln::double_approx(p.c),-p.b), j-p.a-1).real();
                result.push_back(term2);
            }
            else
            {
                throw exception ("[symbolic::integrate_inf] Unknown goniom. function.");
            }
            
            afac_over_jfac *= symbolic::rational(j);
        }
    }
    
    result.optimize();
    return result;
}

symbolic::poly symbolic::integrate_low (symbolic::poly const & P)
{
    symbolic::poly result = symbolic::integrate_full(P) - symbolic::integrate_inf(P);
    result.optimize();
    return result;
}

bool symbolic::term::ordering (symbolic::term const & u, symbolic::term const & v)
{
    // Determine if this is the correct strict ordering.
    
    // descending powers
    if (u.a > v.a)  return true;
    if (u.a == v.a) return false;
    
    // descending frequency
    if (u.b > v.b)  return true;
    if (u.b == v.b) return false;
    
    // descending scale factor
    if (u.c > v.c)  return true;
    if (u.c == v.c) return false;
    
    // descending irrational factor
    if (u.ki > v.ki)  return true;
    if (u.ki == v.ki) return false;
    
    // descending id of goniometric function
    if (u.gf > v.gf)  return true;
    if (u.gf == v.gf) return false;
    
    // the terms differ just by rational factor
    return false;
}

void symbolic::poly::optimize ()
{
    // if nothing to optimize, return
    if (terms.empty())
        return;
    
    // sort terms: descending powers
    std::sort (terms.begin(), terms.end(), symbolic::term::ordering);
    
    // merge simillar terms
    std::vector<symbolic::term> new_terms;
    new_terms.reserve(terms.size());
    for (term const & T : terms)
    {
        if (new_terms.empty())
        {
            // insert new term at the beginning
            new_terms.push_back(T);
        }
        else if (T.ki == 0. or T.kr == 0 or (T.gf == GF_SIN and T.b == 0))
        {
            // skip zero terms
            continue;
        }
        else
        {
            symbolic::term & T0 = new_terms.back();
            if (T0.a == T.a and T0.b == T.b and T0.c == T.c and T0.gf == T.gf and T0.ki == T.ki)
            {
                // sum the two terms
                T0.kr += T.kr;
                    
                // if the result is zero, remove it
                if (T0.kr == 0)
                    new_terms.pop_back();
            }
            else
            {
                // add a new term
                new_terms.push_back(T);
            }
        }
    }
    
    // update terms
    terms = new_terms;
}

double symbolic::eval (symbolic::poly const & P, double r)
{
    double result = 0;        // evaluated symbolic polynomial P
    
    int a = P.front().a;    // exponent

    double b = 0;            // wave number
    double Sin = 0;            // evaluates sine
    double Cos = 1;            // evaluated cosine
    
    symbolic::rational c = 0;        // c-coefficient
    double Exp = 1;            // evaluated exponential
    
    for (symbolic::term const & p : P)
    {
        // precompute goniometric functions
        if (p.b != b)
        {
            Sin = sin(p.b*r);
            Cos = cos(p.b*r);
        }
        
        // precompute exponential (if 'c' changed)
        if (p.c != c)
        {
            c = p.c;
            Exp = exp(-cln::double_approx(p.c)*r);
        }
        
        // update the coordinate power by the difference in exponent
        // between the previous and the new term
        for (int i = a; i > p.a; i--)
            result *= r;
        
        // update exponent
        a = p.a;
        
        // add the term without coordinate power
        if (p.gf == GF_NONE)
        {
            result += p.ki * cln::double_approx(p.kr) * Exp;
        }
        else if (p.gf == GF_COS)
        {
            result += p.ki * cln::double_approx(p.kr) * Cos * Exp;
        }
        else if (p.gf == GF_SIN)
        {
            result += p.ki * cln::double_approx(p.kr) * Sin * Exp;
        }
        else
        {
            throw exception ("[Eval] Unknown goniom. function.");
        }
    }
    
    // update coordinate power by the power of the last term
    result *= pow (r, P.back().a);
    
    return result;
}

double symbolic::term::todouble () const
{
    if (a == 0 and gf == GF_NONE and c == 0)
        return ki * symbolic::double_approx(kr);
    
    throw exception ("Cannot convert \"%s\" to double.", symbolic::tostring(*this).c_str());    
}

std::string symbolic::tostring (symbolic::poly const & P)
{
    std::ostringstream os;
    os << P;
    return os.str();
}

std::ostream & symbolic::operator << (std::ostream & os, symbolic::poly const & P)
{
    if (P.size() == 0)
        os << "0";
    
    for (unsigned i = 0; i < P.size(); i++)
    {
        // write zero (expect only one term)
        if (P[i].ki == 0 or P[i].kr == 0)
        {
            os << "0";
            break;
        }
        
        // write numerical factor
        {
            // use sign only for the first term
            double x = P[i].ki * cln::double_approx(P[i].kr);
            if (x < 0)
                os << " - ";
            if (i > 0 and x > 0)
                os << " + ";
            if (std::abs(x) != 1 or (P[i].a == 0 and P[i].gf == GF_NONE and P[i].c == 0))
                os << std::abs(x);
        }
        
        // write power
        if (P[i].a != 0)
        {
            if (P[i].a == 1)
                os << " r";
            else
                os << " r^" << P[i].a;
        }
        
        // write goniometric function
        if (P[i].gf != GF_NONE)
        {
            os << " " << gfname(P[i].gf) << "(";
            if (P[i].b != 1)
                os << P[i].b;
            os << "r)";
        }
        
        // write exponential
        if (P[i].c == 1)
        {
            os << " exp(-r)";
        }
        else if (P[i].c != 0)
        {
            os << " exp(-" << cln::double_approx(P[i].c) << "r)";
        }
    }
    
    return os;
}

void symbolic::collect (symbolic::poly const & P, symbolic::term const & p, symbolic::poly & Q, symbolic::poly & R)
{
    Q.clear();
    R.clear();
    
    for (symbolic::term const & q : P)
    {
        bool exp_match = (p.c == 0 or q.c == p.c);
        bool sin_match = (p.gf != GF_SIN or q.gf == GF_SIN);
        bool cos_match = (p.gf != GF_COS or q.gf == GF_COS);
        
        if (exp_match and sin_match and cos_match)
            Q.push_back(q/p);
        else
            R.push_back(q);
    }
}
