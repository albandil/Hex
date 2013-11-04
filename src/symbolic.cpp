/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2013                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <limits>

#include "complex.h"
#include "misc.h"
#include "symbolic.h"

#include <gsl/gsl_sf.h>
#define fac(n) gsl_sf_fact(n)

#include <cln/cln.h>

cln::cl_RA Combination(cln::cl_RA const & alpha, int k)
{
    cln::cl_RA res = 1;
    for (int w = 0; w < k; w++)
        res *= (alpha - w) / (k - w);
    return res;
}

const char* gfname(unsigned i)
{
    static const char * const names[] = {"none", "sin", "cos"};
    return names[i];
}

SymbolicTerm operator + (SymbolicTerm const & A, SymbolicTerm const & B)
{
    if ((A.gf != B.gf) or (A.a != B.a) or (A.b != B.b) or (A.c != B.c))
        throw "Can't add non-simillar symbolic terms.";
    
    SymbolicTerm S;
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

// operators on XGE polynomials
SymbolicPoly operator * (const SymbolicPoly& P, const SymbolicPoly& Q)
{
    SymbolicPoly product_poly;
    
    for (SymbolicTerm const & p : P) for (SymbolicTerm const & q : Q)
    {
        SymbolicTerm product_term;
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
            product_term.kr /= cln::cl_I(2);
            product_term.gf = GF_COS;
            
            product_term.b = p.b - q.b;
            product_poly.push_back(product_term);
            if (product_poly.back().b == 0)
                product_poly.back().gf = GF_NONE;
            
            product_term.kr *= cln::cl_I(-1);
            product_term.b = p.b + q.b;
            product_poly.push_back(product_term);
            if (product_poly.back().b == 0)
                product_poly.back().gf = GF_NONE;
        }
        else if (p.gf == GF_SIN and q.gf == GF_COS)
        {
            product_term.kr /= cln::cl_I(2);
            product_term.gf = GF_SIN;
            
            product_term.b = p.b + q.b;
            product_poly.push_back(product_term);
            
            product_term.b = p.b - q.b;
            product_poly.push_back(product_term);
        }
        else if (p.gf == GF_COS and q.gf == GF_SIN)
        {
            product_term.kr /= cln::cl_I(2);
            product_term.gf = GF_SIN;
            
            product_term.b = p.b + q.b;
            product_poly.push_back(product_term);
            
            product_term.kr *= cln::cl_I(-1);
            product_term.b = p.b - q.b;
            product_poly.push_back(product_term);
        }
        else if (p.gf == GF_COS and q.gf == GF_COS)
        {
            product_term.kr /= cln::cl_I(2);
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

SymbolicPoly operator - (const SymbolicPoly& P, const SymbolicPoly& Q)
{
    // copy P
    SymbolicPoly difference = P;
    
    // copy Q
    difference.insert(difference.end(), Q.begin(), Q.end());
    
    // invert signs in the copy of Q
    for (unsigned i = P.size(); i < difference.size(); i++)
        difference[i].kr *= cln::cl_I(-1);
    
    // optimize and return
    difference.optimize();
    return difference;
}

SymbolicPoly operator + (const SymbolicPoly& P, const SymbolicPoly& Q)
{
    // copy P
    SymbolicPoly difference = P;
    
    // copy Q
    difference.insert(difference.end(), Q.begin(), Q.end());
    
    // optimize and return
    difference.optimize();
    return difference;
}

/*
 * Associated Laguerre polynomial
 * Laguerre(k,s,x) := sum((-1)^j * (k!)^2 * x^(j-s) / ( (k-j)! * j! * (j-s)! ), j, s, k);
 */
SymbolicPoly Laguerre(int k, int s)
{
    SymbolicPoly laguerre;
    
    for (int j = s; j <= k; j++)
    {
        SymbolicTerm term;
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

/*
 * Hydrogen radial function normalization factor
 * sqrt((2/n)^3 * (n-l-1)! / (2*n*((n+l)!)^3));
 */
SymbolicTerm HydrogenN(int n, int l)
{
    // compute square of the normalization factor
    cln::cl_RA N2 = cln::expt(cln::cl_RA(2)/n,3) * cln::factorial(n-l-1) / (2*n*cln::expt(cln::factorial(n+l),3));
    
    // the factor itself
    SymbolicTerm N;
    
    // try to compute a precise (rational) square root and store it in N.kr
    if (not cln::sqrtp(N2, &N.kr))
    {
        N.kr = 1;
        N.ki = sqrt(cln::double_approx(N2));
    }
    
    return N;
}

/*
 * Hydrogen radial function.
 * HydrogenP(n,l,r) := r * HydrogenN(n,l) * (2*r/n)^l * Laguerre(n+l,2*l+1,2*r/n) * exp(-r/n);
 */
SymbolicPoly HydrogenP(int n, int l)
{
    SymbolicPoly laguerre = Laguerre(n+l, 2*l+1);
    for (SymbolicTerm & term : laguerre)
        term.kr *= cln::expt(cln::cl_RA(2)/n, term.a);
    
    SymbolicTerm rest;
    rest.kr = cln::expt(cln::cl_RA(2)/n,l);
    rest.a = l+1;
    rest.gf = GF_NONE;
    rest.b = 0;
    rest.c = cln::cl_RA(1)/n;
    
    return HydrogenN(n,l) * laguerre * rest;
}

/*
 * Hydrogen sturmian function.
 * HydrogenS(n,l,λ,r) := n * n * HydrogenN(n,l) * (2*λ*r)^l * Laguerre(n+l,2*l+1,2*λ*r) * exp(-λ*r);
 */
SymbolicPoly HydrogenS(int n, int l, cln::cl_RA lambda)
{
    SymbolicPoly laguerre = Laguerre(n+l, 2*l+1);
    for (SymbolicTerm & term : laguerre)
        term.kr *= cln::expt(2*lambda, term.a);
    
    SymbolicTerm rest;
    rest.kr = n * n * cln::expt(2*lambda,l);
    rest.a = l;
    rest.gf = GF_NONE;
    rest.b = 0;
    rest.c = 1;
    
    return HydrogenN(n,l) * laguerre * rest;
}

/*
 * Riccati-Bessel function
 * BesselF(l,x) := sum(
 *                       sin(x) * (2/x)^(l-2*i) * (l-i)!/i! * combination(-1/2-i,l-2*i) 
 *                     - cos(x) * (2/x)^(l-2*i) * (l-i)!/i! * i * combination(-1/2-i,l-2*i+1),
 *                   i, 0, floor((l+1)/2)
 *                 );
 */
SymbolicPoly RiccatiBessel(int l, double k)
{
    SymbolicPoly ricj;
    
    for (int i = 0; i <= (l+1)/2; i++)
    {
        SymbolicTerm term;
        
        cln::cl_RA comb = Combination(-cln::cl_RA(1)/2-i,l-2*i);
        
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
        SymbolicTerm term;
        
        cln::cl_RA comb = Combination(-cln::cl_RA(1)/2-i,l-2*i+1);
        
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

/*
 * Integrals
 */
SymbolicTerm integrate_full (SymbolicPoly const & P)
{
    // the resulting number
    SymbolicTerm result;
    result.ki = 1;
    result.kr = 0;
    
    // for all terms of the polynomial
    for (SymbolicTerm const & p : P)
    {
        if (p.a < -1 or (p.a == -1 and p.gf != GF_SIN))
            throw exception ("[Integral_xge_full] Divergent integral!");
        
        if (p.a == -1 and p.gf == GF_SIN)
        {
            result.ki = result.ki * cln::double_approx(result.kr) + p.ki * cln::double_approx(p.kr) * (M_PI/2 - atan(cln::double_approx(p.c)/p.b));
            result.kr = 1;
        }
        
        // computed integral of a term
        SymbolicTerm integ;
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

SymbolicPoly integrate_inf(SymbolicPoly const & P)
{
    SymbolicPoly result;
    
    for (SymbolicTerm const & p : P)
    {
        for (int j = 0; j <= p.a; j++)
        {
            if (p.gf == GF_NONE)
            {
                SymbolicTerm term;
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
                SymbolicTerm term1, term2;
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
                SymbolicTerm term1, term2;
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

SymbolicPoly integrate_low(SymbolicPoly const & P)
{
    return integrate_full(P) - integrate_inf(P);
}

void SymbolicPoly::optimize()
{
    // if nothing to optimize, return
    if (terms.empty())
        return;
    
    // sort terms: descending powers
    std::sort (
        terms.begin(),
        terms.end(),
        [ ] (SymbolicTerm const & T1, SymbolicTerm const & T2) -> bool {
            return (T1.a > T2.a) or (T1.b > T2.b) or (T1.c > T2.c) or (T1.gf > T2.gf) or (T1.ki > T2.ki);
        }
    );
    
    // merge simillar terms
    std::vector<SymbolicTerm> new_terms;
    new_terms.reserve(terms.size());
    for (SymbolicTerm const & T : terms)
    {
        if (new_terms.empty())
        {
            // insert new term at the beginning
            new_terms.push_back(T);
        }
        else
        {
            SymbolicTerm & T0 = new_terms.back();
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

double eval (SymbolicPoly const & P, double r)
{
    double result = 0;        // evaluated symbolic polynomial P
    
    int a = P.front().a;    // exponent

    double b = 0;            // wave number
    double Sin = 0;            // evaluates sine
    double Cos = 1;            // evaluated cosine
    
    cln::cl_RA c = 0;        // c-coefficient
    double Exp = 1;            // evaluated exponential
    
    for (SymbolicTerm const & p : P)
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

void write (SymbolicPoly const & P)
{
    for (SymbolicTerm const & p : P)
    {
        if (p.ki == 0 or p.kr == 0)
            continue;
        if (p.ki != 1 or p.kr != 1)
            printf("%g", p.ki * cln::double_approx(p.kr));
        if (p.a != 0)
            printf(" r^%d", p.a);
        if (p.gf != GF_NONE)
            printf(" %s(%gr)", gfname(p.gf), p.b);
        if (p.c == 1)
            printf(" exp(-r)");
        else if (p.c != 0)
            printf(" exp(-%gr)", cln::double_approx(p.c));
        
        if (&p != &P.back())
            printf(" + ");
    }
    printf("\n");
}
