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

#ifndef HEX_DIOPHANTINE
#define HEX_DIOPHANTINE

#include <cmath>
#include <iostream>
#include <string>

#include <gsl/gsl_sf.h>

#include "misc.h"
#include "special.h"

/**
 * @brief Abstract base for node integrators.
 * 
 * This class is the base class for so called "node integrators", which compute
 * integral as a sum of sub-integrals from specific sub-intervals. For example the
 * sub-intervals can be given by zeros of some special function occuring in the
 * integrand. This particular class defines the integration algorithm; however,
 * it is general and not intended for direct use -- the member function that returns
 * the n-th node is pure virtual and expected to be defined in some derived class.
 * See e.g. @ref FixedNodeIntegrator for example of such class.
 */
template <class Functor, class Integrator>
class NodeIntegrator
{
    private:
        
        /// Integration rule to use on sub-intervals.
        Integrator Q_;
        
        /// Function to integrate (double -> double).
        Functor f_;
        
        /// Result of integration.
        double result_;
        
        /// Whether the integration went all right.
        bool ok_;
        
        /// Error text if the integration failed.
        std::string status_;
        
        /// Absolute tolerance.
        double epsabs_;
        
        /// Relative tolerance.
        double epsrel_;
        
        /// Maximum number of nodes to use.
        int limit_;
        
    
    public:
        
        //
        // constructor
        //
        
        NodeIntegrator (Functor f)
            : Q_(f),  f_(f), result_(0), ok_(true), status_(),
              epsabs_(1e-8), epsrel_(1e-6), limit_(1000) {}
        
        //
        // getters and setters
        //
        
        double result () const { return result_; }
        bool ok () const { return ok_; }
        std::string const & status () const { return status_; }
        
        double epsabs () const { return epsabs_; }
        void setEpsAbs (double eps) { epsabs_ = eps; }
        
        double epsrel () const { return epsrel_; }
        void setEpsRel (double eps) { epsrel_ = eps; }
        
        int limit () const { return limit_; }
        void setLimit (double n) { limit_ = n; }
        
        //
        // redefinable node getters
        //
        
        /// Get n-th node.
        virtual double nthNode (int n) const = 0;
        
        /// Get end of classically forbidden region (= exponential rise, no nodes).
        virtual double turningPoint () const { return 0; }
        
        //
        // integration routine
        //
        
        /**
         * @brief Integration routine.
         * 
         * The function will compute integral of function over the specified finite (!) interval.
         * The integration in the potential classically forbidden region is done by simple
         * trapezoidal rule which is refined up to the convergence. The rest is integrated
         * node-to-node by the supplied integrator class (see e.g. @ref GaussKronrod or
         * @ref ClenshawCurtis).
         * 
         * @todo Use Romberg integration in the classically forbidden region instead of the plain
         * trapezoidal integration.
         */
        bool integrate (double a, double b)
        {
            // initialize
            ok_ = true;
            status_ = "";
            result_ = 0;
            
            // end of previous integration parcel
            double prevR = 0;
            
            // start by integrating in the classically forbidden region (if any)
            double rt = this->turningPoint();
            if (rt != 0 and a < rt)
            {
                double r0 = a;
                double r1 = std::min(rt,b);
                
                double quotient = 0.75; // FIXME variable ?
                for (int samples = 8; ; samples *= 2)
                {
                    rArray grid = geomspace(r0,r1,samples,quotient);
                    rArray eval = grid.transform(f_);
                    double estimate = special::integral::trapz(grid, eval);
                    
                    // check iteration limit
                    if ((int)std::log2(samples) == limit_)
                    {
                        throw exception
                        (
                            "Cannot integrate classically forbidden region - reached maximal number of iterations."
                        );
                    }
                    
                    // check convergence
                    if (std::abs(estimate - result_) < epsrel_ * std::abs(estimate))
                    {
                        result_ = estimate;
                        break;
                    }
                    
                    // check hopelessness
                    if (estimate == 0 and result_ == 0)
                    {
                        // this is probably hopeless...
                        break;
                    }
                    
                    // update estimate and go to the next iteration
                    result_ = estimate;
                }
                
                prevR = r1;
            }
            
            // check if the classically forbidden region was the only one to integrate
            if (prevR == b)
                return ok_;
            
            // for all integration parcels (nodes)
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
                
                // integrate and check success
                if (not Q_.integrate(rmin,rmax))
                {
                    ok_ = false;
                    status_ = format("Node integrator failed on <%g,%g>: %s", rmin, rmax, Q_.status().c_str());
                    return ok_;
                }
                
                // update result
                result_ += Q_.result();
                prevR = rmax;
                
                // check convergence
                if (std::abs(Q_.result()) < epsabs_ * std::abs(rmax - rmin) or std::abs(Q_.result()) < epsrel_ * std::abs(result_))
                    break;
            }
            
            return ok_;
        }
};

/**
 * @brief Bessel node integrator.
 * 
 * This class derives from @ref NodeIntegrator and uses the Bessel function
 * zeros as the integration nodes.
 */
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

/**
 * @brief Coulomb node integrator.
 * 
 * This class derives from @ref NodeIntegrator and uses the Coulomb function
 * zeros as the integration nodes.
 */
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

/**
 * @brief Fixed node integrator.
 * 
 * This class derives from @ref NodeIntegrator and uses the node array supplied
 * by the user.
 */
template <class Functor, class Integrator>
class FixedNodeIntegrator : public NodeIntegrator<Functor,Integrator>
{
    public:
        
        FixedNodeIntegrator (Functor f, const rArrayView zeros, double rt = 0)
            : NodeIntegrator<Functor,Integrator>(f), zeros_(zeros), rt_(rt)
        {
            NodeIntegrator<Functor,Integrator>::setLimit(zeros_.size());
        }
        
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
        
        const rArray zeros_;
        double rt_;
};

#endif
