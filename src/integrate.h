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

#ifndef HEX_INTEGRATE
#define HEX_INTEGRATE

#include <cmath>
#include <limits>
#include <string>

#include <o2scl/exception.h>
#include <o2scl/gsl_inte_qag.h>
#include <o2scl/gsl_inte_qagi.h>
#include <o2scl/gsl_inte_qagil.h>
#include <o2scl/gsl_inte_qagiu.h>

#include "misc.h"

/** \brief Numerical integrator
 *
 * The class is initialized by a lambda-function (= the integrand) and
 * serves as a QuadPack wrapper. The member function \ref integrate performs
 * the actual integration. The getters \ref result and \ref abserr return
 * the computed numbers.
 * 
 * For the initialization you will mostly want to use the structure
 * \code
 *     Integrator<decltype(integrand)> Q(integrand);
 * \endcode
 */
template <typename Functor> class Integrator
{
	private:
		Functor F;
		double Result, AbsErr;
		bool Ok;
		std::string Status;
		double eval(double x) { return F(x); }

	public:
		Integrator(Functor f) : F(f) { Result = AbsErr = Nan; Ok = false; }
		~Integrator() {}
		
		/** \brief Compute the integral.
		 *
		 * Performs the integration on a given interval. Uses specialized routines
		 * from the O₂scl library, which itself just wraps GSL (and that is, in this
		 * case, nothing than C-port of QuadPack).
		 * 
		 * You can compute improper integrals. For specifying "infinity" as
		 * one or both bounds use either
		 * 
		 * \code
		 *     std::numeric_limits<double>::infinity()
		 * \endcode
		 * 
		 * for positive infinity or
		 * 
		 * \code
		 *     -std::numeric_limits<double>::infinity()
		 * \endcode
		 * 
		 * for negative infinity.
		 * 
		 * \return The value of \ref Ok (i.e. whether the last integration has
		 * been successful according to the library).
		 */
		bool integrate(double a, double b)
		{
			// reset
			Result = AbsErr = Nan;
			Ok = true;
			Status.clear();
			
			if (a == b) { return 0; }
			
			if (isnan(a) || isnan(b))
			{
				Ok = false;
				Status = std::string("ERROR: Some of the integration bounds is not defined!");
				return false;
			}
			
			// otočíme znaménko
			if (a > b)
			{
				integrate(b, a);
				Result = -Result;
				return Ok;
			}
			
			// EvalPtr je ukazatel na členskou funkci Integrator::eval, kterou budeme posílat
			// do knihovních funkcí.
			o2scl::funct_mfptr<Integrator<Functor>> EvalPtr(const_cast<Integrator<Functor>*>(this), &Integrator<Functor>::eval);
			o2scl::inte<o2scl::funct_mfptr<Integrator<Functor>>> *R = 0;
			
			// podle konečnosti mezí vybere správný integrátor
			if (finite(a) && finite(b))          /* -∞ < a < b < +∞ */
				R = new o2scl::gsl_inte_qag<o2scl::funct_mfptr<Integrator<Functor>>>;
			else if (finite(a) && !finite(b))    /* -∞ < a < b = +∞ */
				R = new o2scl::gsl_inte_qagiu<o2scl::funct_mfptr<Integrator<Functor>>>;
			else if (!finite(a) && finite(b))    /* -∞ = a < b < +∞ */
				R = new o2scl::gsl_inte_qagil<o2scl::funct_mfptr<Integrator<Functor>>>;
			else                                 /* -∞ = a < b = +∞ */
				R = new o2scl::gsl_inte_qagi<o2scl::funct_mfptr<Integrator<Functor>>>;
			
			// pozor! o2scl kope!
			try
			{
				// provede integraci
				R->integ_err(EvalPtr, a, b, Result, AbsErr);
			}
			catch (o2scl::exc_runtime_error e)
			{
				// Integrace se úplně nepovedla, ale pořád lepší než nic. Text výjimy
				// uložíme a necháme uživatele, aby o adlším pokračování rozhodl sám.
				Status = std::string(e.what());
				Ok = false;
			}
			catch (o2scl::exc_exception e)
			{
				// totéž
				Status = std::string(e.what());
				Ok = false;
			}
			
			delete R;
			return Ok;
		}
		
		const std::string& status() const { return Status; }	// text chybového hlášení
		bool ok() const { return Ok; }					// proběhla poslední integrace v pořádku?
		double result() const { return Result; }		// výsledek poslední integrace
		double abserr() const { return AbsErr; }		// a její absolutní chyba
};

#endif
