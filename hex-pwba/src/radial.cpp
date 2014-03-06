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

#include <iostream>

#include "arrays.h"
#include "gausskronrod.h"
#include "radial.h"
#include "romberg.h"
#include "specf.h"

template <class Functor, class Integrator>
class BesselNodeIntegrator
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
        
        BesselNodeIntegrator (Functor f, double k, int l)
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

double compute_Idir (int li, int lf, int lambda, int Ni, int Li, double ki, int Nf, int Lf, double kf)
{
    std::cout << format
    (
        "Precompute Idir\n"
        "\tlambda = %d\n"
        "\tNi = %d, Li = %d, ki = %g, li = %d\n"
        "\tNf = %d, Lf = %d, kf = %g, lf = %d\n",
        lambda, Ni, Li, ki, li, Nf, Lf, kf, lf
    );
    
    if (lambda == 0)
    {
        //
        // r1 > r2
        //
        
        auto integrand = [Ni,Li,Nf,Lf,li,ki,lf,kf](double r2) -> double
        {
            // inner integrand
            auto iintegrand = [Ni,Li,Nf,Lf,r2](double r1) -> double { return hydro_P(Ni,Li,r1) * hydro_P(Nf,Lf,r1) * (1./r1 - 1./r2); };
            
            // inner integrator
            GaussKronrod<decltype(iintegrand)> Qi(iintegrand);
            
            // integrate and check success
            if (not Qi.integrate(r2,Inf))
            {
                throw exception
                (
                    "compute_Idir (inner) failed for λ=0, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d, r2=%g (\"%s\").\n\tresult = %g\n",
                    Ni, Li, ki, li, Nf, Lf, kf, lf, r2, Qi.status().c_str(), Qi.result()
                );
            }
            
            return Qi.result() * ric_j(li,ki*r2) * ric_j(lf,kf*r2);
        };
        
        // which Bessel function oscillates slowlier ?
        int    l = (ki < kf ? li : lf);
        double k = (ki < kf ? ki : kf);
        
        // outer integrator
        BesselNodeIntegrator<decltype(integrand),GaussKronrod<decltype(integrand)>> R(integrand, k, l);
        R.integrate(0, Inf);
        
        return R.result();
    }
    else
    {
        auto integrand = [Ni,Li,Nf,Lf,li,ki,lf,kf,lambda](double r2) -> double
        {
            // inner integrand
            auto iintegrand1 = [Ni,Li,Nf,Lf,r2,lambda](double r1) -> double { return hydro_P(Ni,Li,r1) * hydro_P(Nf,Lf,r1) * std::pow(r1/r2,lambda); };
            auto iintegrand2 = [Ni,Li,Nf,Lf,r2,lambda](double r1) -> double { return hydro_P(Ni,Li,r1) * hydro_P(Nf,Lf,r1) * std::pow(r2/r1,lambda+1); };
            
            // inner integrator
            GaussKronrod<decltype(iintegrand1)> Q1(iintegrand1);
            GaussKronrod<decltype(iintegrand2)> Q2(iintegrand2);
            
            // integrate and check success
            if (not Q1.integrate(0.,r2))
            {
                throw exception
                (
                    "compute_Idir (inner1) failed for λ=lambda, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d, r2=%g (\"%s\").\n\tresult = %g\n",
                    lambda, Ni, Li, ki, li, Nf, Lf, kf, lf, r2, Q1.status().c_str(), Q1.result()
                );
            }
            if (not Q2.integrate(r2,Inf))
            {
                throw exception
                (
                    "compute_Idir (inner2) failed for λ=lambda, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d, r2=%g (\"%s\").\n\tresult = %g\n",
                    lambda, Ni, Li, ki, li, Nf, Lf, kf, lf, r2, Q2.status().c_str(), Q2.result()
                );
            }
            
            return (Q1.result() + Q2.result()) / r2 * ric_j(li,ki*r2) * ric_j(lf,kf*r2);
        };
        
        // outer integral
        double integral = 0;
        
        // outer integrator
        GaussKronrod<decltype(integrand)> Q(integrand);
        
        // which Bessel function oscillates slowlier ?
        int    l = (ki < kf ? li : lf);
        double k = (ki < kf ? ki : kf);
        
        // end of previous integration parcel
        double prevR = 0;
        
        // for all integration parcels (nodes of one of the Bessel function)
        const double eps = 1e-5;
        for (int inode = 1; inode == 1 or std::abs(Q.result()) > eps * std::abs(integral); inode++)
        {
            // get integration bounds
            gsl_sf_result res;
            int err = gsl_sf_bessel_zero_Jnu_e(l+0.5,inode,&res);
            if (err != GSL_SUCCESS)
            {
                throw exception
                (
                    "Cannot find %d-th root of the Bessel function j[%d](%g r) -- %s.",
                    inode, l, k, gsl_strerror(err)
                );
            }
            double rmin = prevR;
            double rmax = res.val / k;
            
            // integrate and check success
            if (not Q.integrate(rmin,rmax))
            {
                throw exception
                (
                    "compute_Idir (outer) failed on inode %d for λ=%d, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d (\"%s\").\n\tresult = %g",
                    inode, lambda, Ni, Li, ki, li, Nf, Lf, kf, lf, Q.status().c_str(), Q.result()
                );
            }
            
            integral += Q.result();
            prevR = rmax;
        }
        
        return integral;
    }
}

double compute_Iexc (int li, int lf, int lambda, int Ni, int Li, double ki, int Nf, int Lf, double kf)
{
    std::cout << format
    (
        "Precompute Iexc\n"
        "\tlambda = %d\n"
        "\tNi = %d, Li = %d, ki = %g, li = %d\n"
        "\tNf = %d, Lf = %d, kf = %g, lf = %d\n",
        lambda, Ni, Li, ki, li, Nf, Lf, kf, lf
    );
    
    if (lambda == 0)
    {
        //
        // r1 > r2
        //
        
        auto integrand = [Ni,Li,Nf,Lf,li,ki,lf,kf](double r2) -> double
        {
            // inner integrand
            auto iintegrand = [Ni,Li,kf,lf,r2](double r1) -> double { return hydro_P(Ni,Li,r1) * ric_j(lf,kf*r1) * (1./r1 - 1./r2); };
            
            // inner integrator
            GaussKronrod<decltype(iintegrand)> Qi(iintegrand);
            
            // integrate and check success
            if (not Qi.integrate(r2,Inf))
            {
                throw exception
                (
                    "compute_Iexc (inner) failed for λ=0, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d, r2=%g (\"%s\").\n\tresult = %g\n",
                    Ni, Li, ki, li, Nf, Lf, kf, lf, r2, Qi.status().c_str(), Qi.result()
                );
            }
            
            return Qi.result() * ric_j(li,ki*r2) * hydro_P(Nf,Lf,r2);
        };
        
        // outer integrator
        GaussKronrod<decltype(integrand)> Q(integrand);
        
        // integrate and check success
        if (not Q.integrate(0.,Inf))
        {
            throw exception
            (
                "compute_Iexc (outer) failed for λ=0, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d (\"%s\").\n\tresult = %g",
                Ni, Li, ki, li, Nf, Lf, kf, lf, Q.status().c_str(), Q.result()
            );
        }
        
        return Q.result();
    }
    else
    {
        auto integrand = [Ni,Li,Nf,Lf,li,ki,lf,kf,lambda](double r2) -> double
        {
            // inner integrand
            auto iintegrand1 = [Ni,Li,kf,lf,r2,lambda](double r1) -> double { return hydro_P(Ni,Li,r1) * ric_j(lf,kf*r1) * std::pow(r1/r2,lambda); };
            auto iintegrand2 = [Ni,Li,kf,lf,r2,lambda](double r1) -> double { return hydro_P(Ni,Li,r1) * ric_j(lf,kf*r1) * std::pow(r2/r1,lambda+1); };
            
            // inner integrator
            BesselNodeIntegrator<decltype(iintegrand1),GaussKronrod<decltype(iintegrand1)>> Q1 (iintegrand1, kf, lf);
            BesselNodeIntegrator<decltype(iintegrand2),GaussKronrod<decltype(iintegrand2)>> Q2 (iintegrand2, kf, lf);
            
            // integrate and check success
            if (not Q1.integrate(0.,r2))
            {
                throw exception
                (
                    "compute_Iexc (inner1) failed for λ=%d, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d, r2=%g (\"%s\").\n\tresult = %g\n",
                    lambda, Ni, Li, ki, li, Nf, Lf, kf, lf, r2, Q1.status().c_str(), Q1.result()
                );
            }
            if (not Q2.integrate(r2,Inf))
            {
                throw exception
                (
                    "compute_Iexc (inner2) failed for λ=%d, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d, r2=%g (\"%s\").\n\tresult = %g\n",
                    lambda, Ni, Li, ki, li, Nf, Lf, kf, lf, r2, Q2.status().c_str(), Q2.result()
                );
            }
            
            return (Q1.result() + Q2.result()) / r2 * ric_j(li,ki*r2) * hydro_P(Nf,Lf,r2);
        };
        
        // outer integrator
        GaussKronrod<decltype(integrand)> Q(integrand);
        
        // integrate and check success
        if (not Q.integrate(0.,Inf))
        {
            throw exception
            (
                "compute_Iexc (outer) failed for λ=%d, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d (\"%s\").\n\tresult = %g",
                lambda, Ni, Li, ki, li, Nf, Lf, kf, lf, Q.status().c_str(), Q.result()
            );
        }
        
        return Q.result();
    }
}
