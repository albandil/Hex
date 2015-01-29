//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2015, Jakub Benda, Charles University in Prague                    //
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

#ifndef HEX_ROMBERG
#define HEX_ROMBERG

#include "arrays.h"
#include "special.h"

/**
 * @brief One-dimensional romberg integration.
 */
template <class T>
class Romberg
{
    private:
        
        // Whether the last computation has been successful.
        bool ok_;
        
        /// Reason for unsuccessful integration.
        std::string status_;
        
        /// Result of the last integration.
        T result_;
        
        /// Minimal level of subdivision.
        unsigned minlevel_;
        
        /// Maximal level of subdivision.
        unsigned maxlevel_;
        
        /// Whether to write verbose output.
        bool verbose_;
        
        /// String prepended to verbose output.
        std::string verbose_pre_;
        
        /// Absolute integration tolerance.
        double epsabs_;
        
        /// Relative integration tolerance.
        double epsrel_;
        
        /// Maximal Romberg order.
        unsigned maxromblevel_;
        
    public:
        
        //
        // constructors
        //
        
        Romberg ()
            : ok_(false), status_(""), result_(0), minlevel_(0), maxlevel_(10),
              verbose_(false), verbose_pre_(""), epsabs_(1e-8), epsrel_(1e-6),
              maxromblevel_(10)
        {
            // nothing to do
        }
        
        //
        // integration routine
        //
        
        /**
         * @brief Integrates the function.
         * 
         * @param f Function to integrate (double -> double).
         * @param a Lower bound.
         * @param b Upper bound.
         */
        template <class Functor> bool integrate (Functor f, double a, double b)
        {
            // reset result
            result_ = 0;
            
            // number of intervals for this subdivision level
            unsigned Nintervals = 1;
                
            // interval size for this subdivision level
            double h = b - a;
            
            // trapezoidal estimates
            std::vector<T> integrals = { 0.5 * (f(b) + f(a)) * h };
            
            // set up Romberg table
            std::vector<std::vector<T>> romberg(1);
            romberg[0].push_back(integrals[0]);
            if (verbose_)
            {
                std::cout << verbose_pre_
                          << std::setw(13) << std::left << h << integrals.back() << " "
                          << std::setw(13) << std::left << romberg[0][0] << std::endl;
            }
            
            // for all subdivision levels
            for (unsigned level = 1; level < maxlevel_; level++)
            {
                // update number and size of intervals
                Nintervals *= 2;
                h *= 0.5;
                
                // sum of evaluations for this subdivision level
                T suma = integrals.back();
                
                // extend storage
                romberg.push_back(std::vector<T>());
                
                // evaluate where not previously evaluated
                for (unsigned i = 1; i < Nintervals; i += 2)
                    suma += f(i * h);
                
                // add the new sum
                integrals.push_back(suma * h);
                
                // update Romberg table
                romberg[level].resize(std::min(level,maxromblevel_) + 1);
                romberg[level][0] = integrals.back();
                if (verbose_) std::cout << verbose_pre_ << std::setw(13) << std::left << h << integrals.back() << " ";
                for (unsigned icol = 1; icol < romberg[level].size(); icol++)
                {
                    T scale = gsl_sf_pow_int(4,icol);
                    romberg[level][icol] = (scale * romberg[level][icol-1] - romberg[level-1][icol-1]) / (scale - 1);
                    if (verbose_) std::cout << std::setw(13) << std::left << romberg[level][icol] << " ";
                }
                if (verbose_) std::cout << std::endl;
                
                // compare estimates
                double Delta = std::abs(romberg[level].back() - romberg[level-1].back());
                if ((Delta < epsabs_ or Delta < epsrel_ * std::abs(romberg[level].back())) and level >= minlevel_)
                {
                    ok_ = true;
                    status_ = "";
                    result_ = romberg[level].back();
                    return ok_;
                }
            }
            
            ok_ = false;
            status_ = format("Subdivision limit (%d) reached.", maxlevel_);
            result_ = romberg.back().back();
            return ok_;
        }
        
        //
        // getters and setters
        //
        
        double epsrel () const { return epsrel_; }
        void setEpsRel (double eps) { epsrel_ = eps; }
        
        double epsabs () const { return epsabs_; }
        void setEpsAbs (double eps) { epsabs_ = eps; }
        
        unsigned minLevel () const { return minlevel_; }
        void setMinLevel (unsigned level) { minlevel_ = level; }
        
        unsigned maxLevel () const { return maxlevel_; }
        void setMaxLevel (unsigned level) { maxlevel_ = level; }
        
        unsigned maxRombLevel () const { return maxromblevel_; }
        void setMaxRombLevel (unsigned level) { maxromblevel_ = level; }
        
        bool verbose () const { return verbose_; }
        void setVerbose (double v, std::string pre = "") { verbose_ = v; verbose_pre_ = pre; }
        
        bool ok () const { return ok_; }
        T result () const { return result_; }
        std::string const & status () const { return status_; }
};

/**
 * @brief Two-dimensional Romberg integration.
 * 
 * This class will integrate a two-dimensional function over a unit square
 * (0,1) × (0,1). It expects the function to be smooth everywhere except for
 * some derivative discontinuities on the diagonal of the integration domain
 * Also, it expects the integrated function to be identically zero on the
 * boundaries, so that all evaluations can be done inside the domain (with
 * a common weight factor).
 *
 * The method used is the Romberg extrapolation applied on simple 2D quadrature
 * rule
 * @f[
 *      \int\!\int_A f(x,y) \mathrm{d}A \approx \frac{1}{3}|A|(f_1 + f_2 + f_3) \ ,
 * @f]
 * where @f$ A @f$ is a triangular integration subdomain and @f$ f_1 @f$,
 * @f$ f_2 @f$ and @f$ f_3 @f$ are function evaluations in its vertices. The error
 * of this rule is @f$ O(h^2) @f$ where @f$ h @f$ is the subgrid step, so the
 * Romberg iterations can be used to reduce the error.
 *
 * The class is a template with the parameter T, which is the arbitrary return
 * type of the integrated function.
 */
template <class T, class Functor>
class UnitSquareRomberg
{
    private:
        
        // function to integrate
        Functor f_;
        
        // relative tolerance
        double epsrel_;
        
        // absolute tolerance
        double epsabs_;
        
        // minimal subdivision limit
        unsigned minlevel_;
        
        // maximal subdivision limit (0 = no evaluations; overrides minlevel_)
        unsigned maxlevel_;
        
        // maximal Romberg aggregation limit (0 = no Romberg)
        unsigned maxromblevel_;
        
        // whether to print diagnostic information
        bool verbose_;
        
        // result of the last integration
        T result_;
        
        // status flag
        bool ok_;
        
        // status message
        std::string status_;
    
    public:
        
        //
        // constructor
        //
        
        UnitSquareRomberg
        (
            Functor f,
            double epsrel = 1e-5,
            double epsabs = 1e-10,
            unsigned minlevel = 0,
            unsigned maxlevel = 10,
            unsigned maxromblevel = 1,
            bool verbose = false
        ) : f_(f), epsrel_(epsrel), epsabs_(epsabs), maxlevel_(maxlevel),
            maxromblevel_(maxromblevel), verbose_(verbose), result_(0),
            ok_(true), status_()
        {
            // do nothing
        }
        
        //
        // getters and setters
        //
        
        double epsrel () const { return epsrel_; }
        void setEpsRel (double eps) { epsrel_ = eps; }
        
        double epsabs () const { return epsabs_; }
        void setEpsAbs (double eps) { epsabs_ = eps; }
        
        unsigned minLevel () const { return minlevel_; }
        void setMinLevel (unsigned level) { minlevel_ = level; }
        
        unsigned maxLevel () const { return maxlevel_; }
        void setMaxLevel (unsigned level) { maxlevel_ = level; }
        
        unsigned maxRombLevel () const { return maxromblevel_; }
        void setMaxRombLevel (unsigned level) { maxromblevel_ = level; }
        
        bool verbose () const { return verbose_; }
        void setVerbose (double v) { verbose_ = v; }
        
        bool ok () const { return ok_; }
        T result () const { return result_; }
        std::string const & status () const { return status_; }
        
        //
        // integration
        //
        
        /**
         * @brief Integrate function.
         * 
         * This method expects the integrated function to accept two coordinates
         * @f$ x @f$ and @f$ y @f$. It will be called for all coordinate pairs.
         */
        bool integrate ()
        {
            // successive estimates
            std::vector<T> integrals(1);
            integrals.push_back(0.);
            
            // edge
            T h = 1.;
            
            // number of cells per dimension
            unsigned n = 1;
            
            // setup Romberg table
            std::vector<std::vector<T>> romberg(1);
            romberg[0].push_back(0.);
            
            // initialize output
            if (verbose_) std::cout << std::setw(13) << std::left << 1. << romberg[0][0] << std::endl;
            
            // for all subdivisions
            for (unsigned level = 1; maxlevel_ == 0 or level <= maxlevel_; level++)
            {
                // extend storage
                integrals.push_back(0);
                romberg.push_back(std::vector<T>());
                
                // update geometry
                h /= 2;
                n *= 2;
                
                // evaluate function at all internal points
                T suma = 0.;
                for (unsigned ix = 1; ix < n; ix++)
                for (unsigned iy = 1; iy < n; iy++)
                    suma += f_ (h*ix,h*iy);
                
                // store integral estimate
                integrals.push_back(suma * h * h);
                
                // update Romberg table
                romberg[level].resize(std::min(level,maxromblevel_) + 1);
                romberg[level][0] = integrals.back();
                if (verbose_) std::cout << std::setw(13) << std::left << h << integrals.back() << " ";
                for (unsigned icol = 1; icol < romberg[level].size(); icol++)
                {
                    T scale = (1 << (2*icol));
                    romberg[level][icol] = (scale * romberg[level][icol-1] - romberg[level-1][icol-1]) / (scale - 1);
                    if (verbose_) std::cout << std::setw(13) << std::left << romberg[level][icol] << " ";
                }
                if (verbose_) std::cout << std::endl;
                
                // compare estimates
                double Delta = std::abs(romberg[level].back() - romberg[level-1].back());
                if ((Delta < epsabs_ or Delta < epsrel_ * std::abs(romberg[level].back())) and level >= minlevel_)
                {
                    ok_ = true;
                    status_ = "";
                    result_ = romberg[level].back();
                    return ok_;
                }
            }
            
            ok_ = false;
            status_ = format("Subdivision limit (%d) reached.", maxlevel_);
            result_ = romberg.back().back();
            return ok_;
        }
        
        /**
         * @brief Integrate function.
         * 
         * This member expects the integrated function to accept list of equidistant coordinates
         * @f$ X = \{x_i \mid i = 1, \dots, n\} @f$, to evaluate the integrand on cartesian product
         * of these sets and to sum the evaluations. I.e.
         * @f[
         *     f(X) \equiv \sum_{i,j} f(x_i,x_j) \,.
         * @f]
         */
        bool integrate_extern ()
        {
            // successive estimates
            std::vector<T> integrals(1);
            integrals.push_back(0.);
            
            // edge
            T h = 1.;
            
            // number of cells per dimension
            unsigned n = 1;
            
            // setup Romberg table
            std::vector<std::vector<T>> romberg(1);
            romberg[0].push_back(0.);
            
            // initialize output
            if (verbose_) std::cout << std::setw(13) << std::left << 1. << romberg[0][0] << std::endl;
            
            // for all subdivisions
            for (unsigned level = 1; maxlevel_ == 0 or level <= maxlevel_; level++)
            {
                // extend storage
                integrals.push_back(0);
                romberg.push_back(std::vector<T>());
                
                // update geometry
                h /= 2;
                n *= 2;
                
                // evaluate function at all internal points (use supplied function)
                T suma = f_ (h * linspace(1u, n-1, n-1));
                
                // store integral estimate
                integrals.push_back(suma * h * h);
                
                // update Romberg table
                romberg[level].resize(std::min(level,maxromblevel_) + 1);
                romberg[level][0] = integrals.back();
                if (verbose_) std::cout << std::setw(13) << std::left << h << integrals.back() << " ";
                for (unsigned icol = 1; icol < romberg[level].size(); icol++)
                {
                    T scale = (1 << (2*icol));
                    romberg[level][icol] = (scale * romberg[level][icol-1] - romberg[level-1][icol-1]) / (scale - 1);
                    if (verbose_) std::cout << std::setw(13) << std::left << romberg[level][icol] << " ";
                }
                if (verbose_) std::cout << std::endl;
                
                // compare estimates
                double Delta = std::abs(romberg[level].back() - romberg[level-1].back());
                if ((Delta < epsabs_ or Delta < epsrel_ * std::abs(romberg[level].back())) and level >= minlevel_)
                {
                    ok_ = true;
                    status_ = "";
                    result_ = romberg[level].back();
                    return ok_;
                }
            }
            
            ok_ = false;
            status_ = format("Subdivision limit (%d) reached.", maxlevel_);
            result_ = romberg.back().back();
            return ok_;
        }
};

#endif /* HEX_ROMBERG */
