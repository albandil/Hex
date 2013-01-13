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
		Integrator(Functor f) : F(f) { Result = AbsErr = std::numeric_limits<double>::quiet_NaN(); Ok = false; }
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
			Result = AbsErr = std::numeric_limits<double>::quiet_NaN();
			Ok = true;
			Status.clear();
			
			if (a == b) { return 0; }
			
			if (isnan(a) || isnan(b))
			{
				Ok = false;
				Status = std::string("CHYBA: Některá z integračních mezí není definovaná!");
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

// /** \brief Dvojrozměný numerický integrátor
//  *
//  * Třída se inicializuje požadovaným integrandem (ukazatelem na funkci, funktorem nebo Λ-funkcí)
//  * a pak se zavolá členská funkce \ref integrate. Integrand je dvojdílný - první část je závislá
//  * jen na proměnné ('x'), podle které probíhá vnější integrál, a jde tedy vytknout (což ušetří nějaký čas).
//  * Druhá část závisí na obou proměnných ('x,y'). Výsledky jsou přístupné přes gettry \ref result
//  * a \ref abserr.
//  * \param Funktor1 Datový typ integrandu závislém jen na 'x'.
//  * \param Funktor2 Datový typ integrandu závislém na 'x' i na 'y'.
// **/
// template <typename Funktor1, typename Funktor2> class Integrator2D
// {
// 	private:
// 		Funktor1 F1;	// část integrandu závislá jen na 'x'
// 		Funktor2 F2;	// část integrandu závislá na 'x' i 'y'
// 		
// 		double Result, AbsErr;
// 		bool Ok;
// 		std::string Status;
// 
// 	public:
// 		/** \brief Konstruktor
// 		 *
// 		 * Konstruktor
// 		 * \param f1 Vytknutá x-závislá část integrandu jako objekt s rozhraním 'double operator()(double)'.
// 		 * \param f2 Vnitřní xy-závislá část integrandu jako objekt s rozhraním 'double operator()(double, double)'.
// 		**/
// 		Integrator2D(Funktor1 f1, Funktor2 f2) : F1(f1), F2(f2) { Result = AbsErr = nan; Ok = false; }
// 		~Integrator2D() {}
// 
// 		/** \brief Spočítá dvojrozměrný integrál
// 		 *
// 		 * Funkce provede požadovanou integraci na oblasti vymezené křivkami x = a, x = b, y = y1(x), y = y2(x).
// 		 * Meze y1, y2 musí jít implicitně konvertovat na funkci typu double(*)(double), tedy pokud chce
// 		 * uživatel použít Λ-funkce, tak musí mít obě prázdnou klozuru.
// 		 * \param a Dolní mez pro vnější integrál (přes první proměnnou 'x').
// 		 * \param b Horní mez pro vnější integrál (přes první proměnnou 'x').
// 		 * \param y1 Dolní (x-závislá) mez pro vnitřní integrál (přes druhou proměnnou 'y').
// 		 * \param y2 Horní (x-závislá) mez pro vnitřní integrál (přes druhou proměnnou 'y').
// 		 * \return Vrací hodnotu proměnné Ok (tedy zpravidla zda poslední integrace proběhla bez komentářů ze strany knihovny o2scl).
// 		**/
// 		bool integrate(double a, double b, double (*y1)(double), double (*y2)(double))
// 		{
// 			// reset
// 			Result = AbsErr = nan;
// 			Ok = true;
// 			Status.clear();
// 			
// 			if (a == b) { return 0; }
// 			
// 			if (isnan(a) || isnan(b))
// 			{
// 				Ok = false;
// 				Status = std::string("CHYBA: Některá z integračních mezí není definovaná!");
// 				return false;
// 			}
// 			
// 			// otočíme znaménko
// 			if (a > b)
// 			{
// 				integrate(b, a, y1, y2);
// 				Result = -Result;
// 				return Ok;
// 			}
// 			
// 			// ReducedIntegrand je integrand zintegrovaný podle 'y' od y1(x) do y2(x)
// 			// (je to funkce proměnné 'x')
// 			auto ReducedIntegrand = [this, y1, y2](double x)->double{
// 				// BoundInnerIntegrand je vnitřní část integrandu (ta, co závisí na 'x' i 'y')
// 				// vyčíslená v bodě daném PARAMETREM 'x' a PROMĚNNOU 'ÿ́' (je to tedy funkce
// 				// proměnné 'y')
// 				auto BoundInnerIntegrand = [this, x](double y)->double{ return this->F2(x,y); };
// 				// zintegruje BoundInnerIntegrand od y1(x) do y2(x)
// 				Integrator<decltype(BoundInnerIntegrand)> R(BoundInnerIntegrand);
// 				R.integrate(y1(x), y2(x));
// 				if (!R.ok())
// 					throw o2scl::exc_runtime_error(R.status());
// 				return this->F1(x) * R.result();
// 			};
// 			
// 			// vnější jednorozměrný integrál
// 			Integrator<decltype(ReducedIntegrand)> R(ReducedIntegrand);
// 			R.integrate(a, b);
// 			if (!R.ok())
// 			{
// 				Status = R.status();
// 				Ok = false;
// 			}
// 			
// 			Result = R.result();
// 			AbsErr = R.abserr();
// 			
// 			return Ok;
// 		}
// 		
// 		const std::string & status() const { return Status; }	// text chybového hlášení
// 		bool ok() const { return Ok; }				// proběhla poslední integrace v pořádku?
// 		double result() const { return Result; }		// výsledek poslední integrace
// 		double abserr() const { return AbsErr; }		// a její absolutní chyba
// };

#endif
