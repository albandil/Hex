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

#ifndef HEX_GAUSSKRONROD
#define HEX_GAUSSKRONROD

#include <cmath>

#include <gsl/gsl_integration.h>

#include "misc.h"
#include "special.h"

/** @brief Numerical integrator
 *
 * The class is initialized by a lambda-function (= the integrand) and
 * serves as a QuadPack wrapper. The member function @ref integrate performs
 * the actual integration. The getters @ref result and @ref abserr return
 * the computed numbers.
 * 
 * For the initialization you will mostly want to use the structure
 * @code
 *     GaussKronrod<decltype(integrand)> Q(integrand);
 * @endcode
 */
template <typename Functor> class GaussKronrod : public RadialFunction<double>
{
    private:
        
        /// Integrated function.
        Functor Integrand;
        
        /// Result of the last integration.
        double Result;
        
        /// Estimation of the absolute error of the last integration.
        double AbsErr;
        
        /// Whether the last integration succeeded.
        bool Ok;
        
        /// Absolute tolerance of the integrator.
        double EpsAbs;
        
        /// Relative tolerance of the integrator.
        double EpsRel;
        
        /// Limit on subdivision.
        size_t Limit;
        
        /// Integration workspace pointer.
        gsl_integration_workspace * Workspace;
        
        /// Status information.
        std::string Status;
        
    public:
        
        // constructor
        GaussKronrod (Functor f, size_t limit = 1000)
            : Integrand(f), Result(special::constant::Nan), AbsErr(special::constant::Nan),
              Ok(false), EpsAbs(0.), EpsRel(1e-5), Limit(limit),
              Workspace(gsl_integration_workspace_alloc(Limit)) {}
        
        // destructor
        ~GaussKronrod()
        {
            gsl_integration_workspace_free(Workspace);
        }
        
        // evaluate the integrand
        inline double operator() (double x) const
        {
            return Integrand(x);
        }
        
        /**
         * @brief Evaluates radial function supplied in a pointer.
         * 
         * This function serves as a workaround for impossibility of passing pointer-to-member
         * to the GSL routines.
         */
        static double eval (double x, void * ptr_radf)
        {
            RadialFunction<double> const & f = (* static_cast<RadialFunction<double> const *>(ptr_radf));
            return f(x);
        }
        
        /** @brief Compute the integral.
         *
         * Performs the integration on a given interval. Uses specialized routines
         * from GSL (that is, in this case, nothing than C-port of QuadPack).
         * 
         * You can compute improper integrals. For specifying "infinity" as
         * one or both bounds use either
         * 
          @code
              std::numeric_limits<double>::infinity()
          @endcode
         * 
         * for positive infinity or
         * 
          @code
              -std::numeric_limits<double>::infinity()
          @endcode
         * 
         * for negative infinity.
         * 
         * @return The value of "Ok" (i.e. whether the last integration has
         * been successful according to the library).
         */
        bool integrate (double a, double b)
        {
            // reset
            Result = AbsErr = special::constant::Nan;
            Ok = true;
            
            // skip empty intervals
            if (a == b)
            {
                Ok = true;
                Result = 0;
                Status = "";
                return Ok;
            }
            
            // 
            if (std::isnan(a) or std::isnan(b))
            {
                Ok = false;
                Result = special::constant::Nan;
                Status = "Some of the bounds is not finite.";
                return false;
            }
            
            // check order of bounds
            if (a > b)
            {
                integrate(b, a);
                Result = -Result;
                return Ok;
            }
            
            // setup integrand
            gsl_function F;
            F.function = &GaussKronrod::eval;
            F.params = this;
            
            // setup integrator
            int err = GSL_SUCCESS, err1 = GSL_SUCCESS, err2 = GSL_SUCCESS;
            double Result1, Result2, AbsErr1, AbsErr2;
            
            // use correct integrator
            if (std::isfinite(a) and std::isfinite(b))          /* -∞ < a < b < +∞ */
            {
                err = gsl_integration_qag
                (
                    &F, a, b,
                    EpsAbs, EpsRel, Limit,
                    GSL_INTEG_GAUSS51,
                    Workspace,
                    &Result, &AbsErr
                );
                Status = (err != GSL_SUCCESS ? gsl_strerror(err) : "");
            }
            else if (std::isfinite(a) and not std::isfinite(b))    /* -∞ < a < b = +∞ */
            {
                err = gsl_integration_qagiu
                (
                    &F, a,
                    EpsAbs, EpsRel, Limit,
                    Workspace,
                    &Result, &AbsErr
                );
                Status = (err != GSL_SUCCESS ? gsl_strerror(err) : "");
            }
            else if (not std::isfinite(a) and std::isfinite(b))    /* -∞ = a < b < +∞ */
            {
                err = gsl_integration_qagil
                (
                    &F, b,
                    EpsAbs, EpsRel, Limit,
                    Workspace,
                    &Result, &AbsErr
                );
                Status = (err != GSL_SUCCESS ? gsl_strerror(err) : "");
            }
            else                                 /* -∞ = a < b = +∞ */
            {
                err1 = gsl_integration_qagiu
                (
                    &F, a,
                    EpsAbs, EpsRel, Limit,
                    Workspace,
                    &Result1, &AbsErr1
                );
                
                err2 = gsl_integration_qagil
                (
                    &F, b,
                    EpsAbs, EpsRel, Limit,
                    Workspace,
                    &Result2, &AbsErr2
                );
                
                Result = Result1 + Result2;
                AbsErr = AbsErr1 + AbsErr2;
                
                Status = (err1 != GSL_SUCCESS or err2 != GSL_SUCCESS ? std::string(gsl_strerror(err1)) + std::string(" & ") + std::string(gsl_strerror(err2)) : std::string());
            }
            
            Ok = (err == GSL_SUCCESS and err1 == GSL_SUCCESS and err2 == GSL_SUCCESS);
            return Ok;
        }
        
        //
        // getters & setters
        //
        
        bool ok() const { return Ok; }
        double result() const { return Result; }
        double error() const { return AbsErr; }
        std::string const & status() const { return Status; }
        
        void setEpsAbs (double epsabs) { EpsAbs = epsabs; }
        void setEpsRel (double epsrel) { EpsRel = epsrel; }
        
        double epsabs () const { return EpsAbs; }
        double epsrel () const { return EpsRel; }
};

#endif
