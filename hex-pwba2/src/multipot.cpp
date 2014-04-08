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

#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf.h>

#include "diophantine.h"
#include "gausskronrod.h"
#include "hydrogen.h"
#include "misc.h"
#include "multipot.h"

MultipolePotential::MultipolePotential (int lambda, int Na, int La, int Nb, int Lb)
    : type_(bound_bound), lambda_(lambda), Na_(Na), La_(La), Nb_(Nb), Lb_(Lb)
{
    // combined exponential scales
    c_ = 1./Na_ + 1./Nb_;
    
    // normalization factors
    norm_ = std::sqrt(std::pow(2./Na_,3) * gsl_sf_fact(Na_-La_-1)/(2.*Na_*gsl_sf_fact(Na_+La_)))
          * std::sqrt(std::pow(2./Nb_,3) * gsl_sf_fact(Nb_-Lb_-1)/(2.*Nb_*gsl_sf_fact(Nb_+Lb_)));
    
    // polynomial (Laguerre) coefficients for the bound states
    double acoef = std::pow(2./Na_,La_) * gsl_sf_choose(Na_+La_,Na_-La_-1);
    rArray acoefs(Na_-La_);
    for (int a = 0; a <= Na_-La_-1; a++)
    {
        acoefs[a] = acoef;
        acoef *= -2 * (Na_ - La_ - 1 - a);
        acoef /= Na_ * (a + 1) * (2*La_ + 2 + a);
    }
    double bcoef = std::pow(2./Nb_,Lb_) * gsl_sf_choose(Nb_+Lb_,Nb_-Lb_-1);
    rArray bcoefs(Nb_-Lb_);
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
}

MultipolePotential::MultipolePotential (int lambda, double Ka, int La, int Nb, int Lb)
    : type_(bound_free), lambda_(lambda), La_(La), Nb_(Nb), Lb_(Lb), Ka_(Ka)
{
    // determine sufficiently far radius where the radial function is small
    double Rfar = Hydrogen::getBoundFar(Nb_, Lb_, 1e-15);
    
    // store at 100 samples per wave length
    double wavelength = 2 * special::constant::pi / Ka_;
    int samples = 100 * Rfar / wavelength;
    
    // the multipole integral evaluated at r = y
    auto integral = [&](double y) -> double
    {
        if (lambda_ == 0)
        {
            auto integrand = [&](double x) -> double
            {
                return Hydrogen::P(Nb_, Lb_, x) * (1./x - 1./y) * Hydrogen::F(Ka_, La_, x);
            };
            BesselNodeIntegrator<decltype(integrand),GaussKronrod<decltype(integrand)>> Q(integrand,Ka_,La_);
            Q.integrate(y, special::constant::Inf);
            
            return Q.result();
        }
        else
        {
            auto integrand1 = [&](double x) -> double
            {
                return Hydrogen::P(Nb_, Lb_, x) * std::pow(x/y,lambda_) * Hydrogen::F(Ka_, La_, x);
            };
            BesselNodeIntegrator<decltype(integrand1),GaussKronrod<decltype(integrand1)>> Q1(integrand1,Ka_,La_);
            Q1.integrate(0, y);
            
            auto integrand2 = [&](double x) -> double
            {
                return Hydrogen::P(Nb_, Lb_, x) * std::pow(y/x,lambda_+1) * Hydrogen::F(Ka_, La_, x);
            };
            BesselNodeIntegrator<decltype(integrand2),GaussKronrod<decltype(integrand2)>> Q2(integrand2,Ka_,La_);
            Q2.integrate(y, special::constant::Inf);
            
            return (Q1.result() + Q2.result()) / y;
        }
    };
    
    // evaluate on grid
    r_ = linspace(0., Rfar, samples);
    V_ = r_.transform(integral);
    
    // dump to file
    write_array(r_, V_, format("Vbf[%d]-%d-%d-%g-%d.dat",lambda_,Na_,La_,Kb_,Lb_));
    
    // setup the interpolator
    interp_ = gsl_interp_alloc(gsl_interp_cspline, samples);
    gsl_interp_init(interp_, r_.data(), V_.data(), samples);
    accel_ = gsl_interp_accel_alloc();
}

MultipolePotential::MultipolePotential (int lambda, int Na, int La, double Kb, int Lb)
    : type_(free_bound), lambda_(lambda), Na_(Na), La_(La), Lb_(Lb), Kb_(Kb)
{
    // determine sufficiently far radius where the radial function is small
    double Rfar = Hydrogen::getBoundFar(Na_, La_, 1e-15);
    
    // store at 100 samples per wave length
    double wavelength = 2 * special::constant::pi / Kb_;
    int samples = 100 * Rfar / wavelength;
    
    // the multipole integral evaluated at r = y
    auto integral = [&](double y) -> double
    {
        if (lambda_ == 0)
        {
            auto integrand = [&](double x) -> double
            {
                return Hydrogen::P(Na_, La_, x) * (1./x - 1./y) * Hydrogen::F(Kb_, Lb_, x);
            };
            BesselNodeIntegrator<decltype(integrand),GaussKronrod<decltype(integrand)>> Q(integrand,Kb_,Lb_);
            Q.integrate(y, special::constant::Inf);
            
            return Q.result();
        }
        else
        {
            auto integrand1 = [&](double x) -> double
            {
                return Hydrogen::P(Na_, La_, x) * std::pow(x/y,lambda_) * Hydrogen::F(Kb_, Lb_, x);
            };
            BesselNodeIntegrator<decltype(integrand1),GaussKronrod<decltype(integrand1)>> Q1(integrand1,Kb_,Lb_);
            Q1.integrate(0, y);
            
            auto integrand2 = [&](double x) -> double
            {
                return Hydrogen::P(Na_, La_, x) * std::pow(y/x,lambda_+1) * Hydrogen::F(Kb_, Lb_, x);
            };
            BesselNodeIntegrator<decltype(integrand2),GaussKronrod<decltype(integrand2)>> Q2(integrand2,Kb_,Lb_);
            Q2.integrate(y, special::constant::Inf);
            
            return (Q1.result() + Q2.result()) / y;
        }
    };
    
    // evaluate on grid
    r_ = linspace(0., Rfar, samples);
    V_ = r_.transform(integral);
    
    // dump to file
    write_array(r_, V_, format("Vfb[%d]-%d-%d-%g-%d.dat",lambda_,Ka_,La_,Nb_,Lb_));
    
    // setup the interpolator
    interp_ = gsl_interp_alloc(gsl_interp_cspline, samples);
    gsl_interp_init(interp_, r_.data(), V_.data(), samples);
    accel_ = gsl_interp_accel_alloc();
}

MultipolePotential::MultipolePotential (int lambda, double Ka, int La, double Kb, int Lb)
    : type_(free_free), lambda_(lambda), La_(La), Lb_(Lb), Ka_(Ka), Kb_(Kb)
{
    throw exception
    (
        "Multipole free-free potential not implemented!"
    );
}

double MultipolePotential::operator() (double x) const
{
    if (x == 0.)
        return 0;
    
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
                    term += gsl_sf_gamma(n) * gsl_sf_gamma_inc_Q(n,cx) - gsl_sf_gamma(n+1) * gsl_sf_gamma_inc_Q(n+1,cx) / cx;
                }
                else
                {
                    // add integral [0,x]
                    term += std::pow(cx,-lambda_-1) * gsl_sf_gamma(n+lambda_+1) * gsl_sf_gamma_inc_P(n+lambda_+1,cx);
                    
                    // add integral [x,inf]
                    term += std::pow(cx,lambda_) * gsl_sf_gamma(n-lambda_) * gsl_sf_gamma_inc_Q(n-lambda_,cx);
                }
                
                // add to sum
                result += coefs_[i] * term * std::pow(c_,-n);
            }
            
            return norm_ * result;
        }
        
        case bound_free:
        case free_bound:
        case free_free:
        {
            return gsl_interp_eval(interp_, r_.data(), V_.data(), x, accel_);
        }
        
        default:
        {
            throw exception ("Runtime error in %s, line %d.", __FILE__, __LINE__);
        }
    };
}
