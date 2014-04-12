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
#include "specf.h"

template <class Functor, class Integrator>
class NodeIntegrator
{
    private:
        
        Integrator Q_;
        Functor f_;
        
        double result_;
        bool ok_;
        std::string status_;
        
        double epsabs_;
        double epsrel_;
        int limit_;
        
    
    public:
        
        NodeIntegrator (Functor f)
            : Q_(f),  f_(f), result_(0), ok_(true), status_(),
              epsabs_(1e-8), epsrel_(1e-5), limit_(1000)
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
        
        virtual double nthNode (int n) const = 0;
        virtual double turningPoint () const { return 0; }
        
        bool integrate (double a, double b)
        {
            // overall integral
            double integral = 0;
            
            // end of previous integration parcel
            double prevR = 0;
            
            // start by integrating in the classically forbidden region (if any)
            double rt = turningPoint();
            if (rt != 0 and a < rt and b > rt)
            {
                double quotient = 0.75; // FIXME variable ?
                for (int samples = 8; ; samples *= 2)
                {
                    rArray grid = geomspace(a,rt,samples,quotient);
                    rArray eval = grid.transform(f_);
                    double estimate = special::integral::trapz(grid, eval);
                    
                    // check iteration limit
                    if ((int)log2(samples) == limit_)
                    {
                        throw exception
                        (
                            "Cannot integrate classically forbidden region - reached maximal number of iterations."
                        );
                    }
                    
                    // check convergence
                    if (std::abs(estimate - integral) < epsrel_ * std::abs(estimate))
                    {
                        integral = estimate;
                        break;
                    }
                    
                    // update estimate and go to the next iteration
                    integral = estimate;
                }
                prevR = rt;
            }
            
            // for all integration parcels (nodes of the Bessel function)
            for (int inode = 1; inode < limit_; inode++)
            {
                // get integration bounds
                double rmin = prevR;
                double rmax = nthNode(inode);
                
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
                    status_ = format("Node integrator failed on <%g,%g>: %s", rmin, rmax, Q_.status().c_str());
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

template <class Functor, class Integrator>
class BesselNodeIntegrator : public NodeIntegrator<Functor,Integrator>
{
    public:
        
        BesselNodeIntegrator (Functor f, double k, int l)
            : NodeIntegrator<Functor,Integrator>(f), k_(k), l_(l) {}
        
        virtual double nthNode (int n) const
        {
            gsl_sf_result res;
            int err = gsl_sf_bessel_zero_Jnu_e(l_+0.5,n,&res);
            if (err != GSL_SUCCESS)
            {
                throw exception
                (
                    "Cannot find %d-th root of the Bessel function j[%d](%g r) -- %s.",
                    n, l_, k_, gsl_strerror(err)
                );
            }
            return res.val / k_;
        }
        
        virtual double turningPoint () const
        {
            return std::sqrt(l_*(l_+1)) / k_;
        };
        
    private:
        
        double k_;
        int l_;
};

template <class Functor, class Integrator>
class CoulombNodeIntegrator : public NodeIntegrator<Functor,Integrator>
{
    public:
        
        CoulombNodeIntegrator (Functor f, double k, int l)
            : NodeIntegrator<Functor,Integrator>(f), k_(k), l_(l), nzeros_(0), zeros_(nullptr) {}
        
        ~CoulombNodeIntegrator()
        {
            if (zeros_ != nullptr)
                delete [] zeros_;
        }
        
        virtual double nthNode (int n) const
        {
            if (n <= 0)
                return 0;
            
            // compute more zeros if necessary
            if (n > nzeros_)
            {
                // release old data
                if (zeros_ != nullptr)
                    delete [] zeros_;
                
                // compute new number of nodes
                nzeros_ = 2 * std::max(n, nzeros_);
                
                // compute the zeros
                zeros_ = new double [nzeros_];
                special::coulomb_zeros(-1/k_, l_, nzeros_, zeros_);
            }
            
            // return the requested zero
            return zeros_[n-1] / k_;
        }
        
        virtual double turningPoint () const
        {
            return (std::sqrt(1 + l_*(l_+1)*k_*k_) - 1) / (k_*k_);
        };
        
    private:
        
        double k_;
        int l_;
        
        mutable int nzeros_;
        mutable double * zeros_;
};

template <class Functor, class Integrator>
class FixedNodeIntegrator : public NodeIntegrator<Functor,Integrator>
{
    public:
        
        FixedNodeIntegrator (Functor f, const rArrayView zeros, double rt = 0)
            : NodeIntegrator<Functor,Integrator>(f), zeros_(zeros), rt_(rt) {}
        
        virtual double nthNode (int n) const
        {
            if (n <= 0)
                return 0;
            
            if (n <= (int)zeros_.size())
                return zeros_[n-1];
                
            return special::constant::Inf;
        }
        
        virtual double turningPoint () const
        {
            return rt_;
        };
        
    private:
        
        double k_;
        int l_;
        
        const rArrayView zeros_;
        double rt_;
};

#endif
