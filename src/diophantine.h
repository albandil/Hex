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

#ifndef HEX_DIOPHANTINE
#define HEX_DIOPHANTINE

#include <iostream>
#include <string>

#include <gsl/gsl_sf.h>

#include "misc.h"

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
                
//                 std::cout << "integrating node : " << rmin << " to " << rmax << std::flush;
                
                // integrate and check success
                if (not Q_.integrate(rmin,rmax))
                {
                    ok_ = false;
                    status_ = Q_.status();
                    result_ = integral;
                    return ok_;
                }
                
//                 std::cout << " : " << Q_.result() << std::endl;
                
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

#endif
