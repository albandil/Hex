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
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <vector>

#include <cln/cln.h>

#include "complex.h"
#include "misc.h"
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
    }
    
    return S;
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
    symbolic::poly difference = P;
    
    // copy Q
    difference.insert(difference.end(), Q.begin(), Q.end());
    
    // optimize and return
    difference.optimize();
    return difference;
}

symbolic::poly symbolic::Laguerre (int k, int s)
{
    symbolic::poly laguerre;
    
    for (int j = s; j <= k; j++)
    {
        term term;
        term.kr = -cln::expt(cln::cl_I(-1),j) * cln::factorial(k) * cln::factorial(k) / (cln::factorial(k-j) * cln::factorial(j) * cln::factorial(j-s));
        term.a = j-s;
        term.gf = GF_NONE;
        term.b = 0;
        term.c = 0;
        laguerre.push_back(term);
    }
    
    laguerre.optimize();
    return laguerre;
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

symbolic::poly symbolic::HydrogenP (int n, int l)
{
    poly laguerre = symbolic::Laguerre(n+l, 2*l+1);
    for (term & term : laguerre)
        term.kr *= cln::expt(symbolic::rational(2)/n, term.a);
    
    term rest;
    rest.kr = cln::expt(symbolic::rational(2)/n,l);
    rest.a = l+1;
    rest.gf = GF_NONE;
    rest.b = 0;
    rest.c = symbolic::rational(1)/n;
    
    return symbolic::HydrogenN(n,l) * laguerre * rest;
}

symbolic::poly symbolic::HydrogenS (int n, int l, symbolic::rational lambda)
{
    symbolic::poly laguerre = symbolic::Laguerre(n+l, 2*l+1);
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
    // TODO
    return symbolic::poly();
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
    
    // for all terms of the polynomial
    for (symbolic::term const & p : P)
    {
        if (p.a < -1 or (p.a == -1 and p.gf != GF_SIN))
            throw exception ("[Integral_xge_full] Divergent integral!");
        
        if (p.a == -1 and p.gf == GF_SIN)
        {
            result.ki = result.ki * cln::double_approx(result.kr) + p.ki * cln::double_approx(p.kr) * (M_PI/2 - atan(cln::double_approx(p.c)/p.b));
            result.kr = 1;
        }
        
        // computed integral of a term
        symbolic::term integ;
        integ.ki = p.ki;
        
        // branch according to the goniometric function
        if (p.gf == GF_NONE)
        {
            integ.kr = p.kr * cln::factorial(p.a) * cln::expt(p.c, -p.a-1);
        }
        else if (p.gf == GF_SIN)
        {
            integ.kr = p.kr * cln::factorial(p.a) * cln::rational(cln::imagpart(cln::expt(cln::complex(p.c,-p.b), -p.a-1)));
        }
        else if (p.gf == GF_COS)
        {
            integ.kr = p.kr * cln::factorial(p.a) * cln::rational(cln::realpart(cln::expt(cln::complex(p.c,-p.b), -p.a-1)));
        }
        else
        {
            throw exception ("[Integral_xge_full] Unknown goniom. function.");
        }
        
        // add to result
        result = result + integ;
    }
    
    return result;
}

symbolic::poly symbolic::integrate_inf (symbolic::poly const & P)
{
    symbolic::poly result;
    
    for (symbolic::term const & p : P)
    {
        for (int j = 0; j <= p.a; j++)
        {
            if (p.gf == GF_NONE)
            {
                symbolic::term term;
                term.ki = p.ki;
                term.kr = p.kr * cln::factorial(p.a) * cln::expt(p.c, j - p.a - 1) / cln::factorial(j);
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
                term1.ki = term2.ki = p.ki;
                
                term1.gf = GF_COS;
                term1.kr = p.kr * cln::factorial(p.a) * cln::rational(cln::realpart(cln::expt(cln::complex(p.c,-p.b), j-p.a-1))) / cln::factorial(j);
                result.push_back(term1);
                
                term2.gf = GF_SIN;
                term2.kr = -p.kr * cln::factorial(p.a) * cln::rational(cln::imagpart(cln::expt(cln::complex(p.c,-p.b), j-p.a-1))) / cln::factorial(j);
                result.push_back(term2);
            }
            else if (p.gf == GF_SIN)
            {
                symbolic::term term1, term2;
                term1.a = term2.a = j;
                term1.b = term2.b = p.b;
                term1.c = term2.c = p.c;
                term1.ki = term2.ki = p.ki;
                
                term1.gf = GF_COS;
                term1.kr = p.kr * cln::factorial(p.a) * cln::rational(cln::imagpart(cln::expt(cln::complex(p.c,-p.b), j-p.a-1))) / cln::factorial(j);
                result.push_back(term1);
                
                term2.gf = GF_SIN;
                term2.kr = p.kr * cln::factorial(p.a) * cln::rational(cln::realpart(cln::expt(cln::complex(p.c,-p.b), j-p.a-1))) / cln::factorial(j);
                result.push_back(term2);
            }
            else
            {
                throw exception ("[Integral_xge_full] Unknown goniom. function.");
            }
        }
    }
    
    result.optimize();
    return result;
}

symbolic::poly symbolic::integrate_low (symbolic::poly const & P)
{
    return symbolic::integrate_full(P) - symbolic::integrate_inf(P);
}

void symbolic::poly::optimize ()
{
    // if nothing to optimize, return
    if (terms.empty())
        return;
    
    // sort terms: descending powers
    std::sort (
        terms.begin(),
        terms.end(),
        [ ] (symbolic::term const & T1, symbolic::term const & T2) -> bool {
            return (T1.a > T2.a) or (T1.b > T2.b) or (T1.c > T2.c) or (T1.gf > T2.gf) or (T1.ki > T2.ki);
        }
    );
    
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

std::ostream & symbolic::poly::operator << (std::ostream & os) const
{
    for (symbolic::term const & p : terms)
    {
        if (p.ki == 0 or p.kr == 0)
            continue;
        if (p.ki != 1 or p.kr != 1)
            os << p.ki * cln::double_approx(p.kr);
        if (p.a != 0)
            os << " r^" << p.a;
        if (p.gf != GF_NONE)
            os << " " << gfname(p.gf) << "(" << p.b << "r)";
        if (p.c == 1)
            os << " exp(-r)";
        else if (p.c != 0)
            os << " exp(-" << cln::double_approx(p.c) << "r)";
        
        if (&p != &terms.back())
            os << " + ";
    }
    os << "\n";
    return os;
}
