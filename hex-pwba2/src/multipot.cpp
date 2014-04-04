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

#include <cmath>

#include <gsl/gsl_sf.h>

#include "misc.h"
#include "multipot.h"

MultipolePotential::MultipolePotential (int lambda, int Na, int La, int Nb, int Lb)
    : type_(bound_bound), lambda_(lambda), Na_(Na), La_(La), Nb_(Nb), Lb_(Lb)
{
    // combined exponential scale
    c_ = 1./Na_ + 1./Nb_;
    
    // normalization factors
    norm_ = std::sqrt(std::pow(2./Na_,3) * gsl_sf_fact(Na_-La_-1)/(2.*Na_*gsl_sf_fact(Na_+La_)))
          * std::sqrt(std::pow(2./Nb_,3) * gsl_sf_fact(Nb_-Lb_-1)/(2.*Nb_*gsl_sf_fact(Nb_+Lb_)));
    
    // polynomial coefficients for the bound states
    double acoef = std::pow(2./Na_,La_); rArray acoefs(Na_-La_);
    for (int a = 0; a <= Na_-La_-1; a++)
    {
        acoefs[a] = acoef;
        acoef *= -2 * (Na_ - La_ - 1 - a);
        acoef /= Na_ * (a + 1) * (2*La_ + 2 + a);
    }
    double bcoef = std::pow(2./Nb_,Lb_); rArray bcoefs(Nb_-Lb_);
    for (int b = 0; b <= Nb_-Lb_-1; b++)
    {
        bcoefs[b] = bcoef;
        bcoef *= -2 * (Nb_ - Lb_ - 1 - b);
        bcoef /= Nb_ * (b + 1) * (2*Lb_ + 2 + b);
    }
    
    // polynomial coefficients of the product
    coefs_.resize(acoefs.size() * bcoefs.size());
    for (int a = 0; a <= Na_-La_-1; a++)
    for (int b = 0; b <= Nb_-Lb_-1; b++)
    {
        coefs_[a+b] = acoefs[a] * bcoefs[b];
    }
    
    std::cout << "Laguerre a : " << acoefs << std::endl;
    std::cout << "Laguerre b : " << bcoefs << std::endl;
}

MultipolePotential::MultipolePotential (int lambda, double Ka, int La, int Nb, int Lb)
    : type_(bound_free), lambda_(lambda), La_(La), Nb_(Nb), Lb_(Lb), Ka_(Ka)
{
    // initialize
}

MultipolePotential::MultipolePotential (int lambda, int Na, int La, double Kb, int Lb)
    : type_(free_bound), lambda_(lambda), Na_(Na), La_(La), Lb_(Lb), Kb_(Kb)
{
    // initialize
}

MultipolePotential::MultipolePotential (int lambda, double Ka, int La, double Kb, int Lb)
    : type_(free_free), lambda_(lambda), La_(La), Lb_(Lb), Ka_(Ka), Kb_(Kb)
{
    throw exception
    (
        "Multipole free-free potential not implemented!"
    );
}

double MultipolePotential::operator() (double x)
{
    switch (type_)
    {
        case bound_bound:
        {
            // compute integrals for all powers
            double result = 0.;
            double cx = c_ * x;
            for (unsigned i = 0; i < coefs_.size(); i++)
            {
                // exponent
                int n = La_ + Lb_ + 2 + i;
                
                // term
                double term = 0;
                
                // compute integral
                if (lambda_ == 0)
                {
                    // add integral [x,inf]
                    term += gsl_sf_gamma(n) * gsl_sf_gamma_inc_Q(n,cx) - gsl_sf_gamma(n+1) * gsl_sf_gamma_inc_Q(n+1,cx);
                }
                else
                {
                    // add integral [0,x]
                    term += std::pow(cx,-lambda_-1) * gsl_sf_gamma(n+lambda_+1) * gsl_sf_gamma_inc_P(n+lambda_+1,cx);
                    
                    // add integral [x,inf]
                    term += std::pow(cx,lambda_) * gsl_sf_gamma(n-lambda_) * gsl_sf_gamma_inc_Q(n-lambda_,cx);
                }
                
                // add to sum
                result += coefs_[i] * term;
            }
            
            return norm_ * result * std::pow(c_,-(La_+Lb_+2));
        }
        
        case bound_free:
        {
            
        }
        
        case free_bound:
        {
            
        }
        
        case free_free:
        {
            
        }
        
        default:
        {
            throw exception ("Runtime error in %s, line %d.", __FILE__, __LINE__);
        }
    };
}
