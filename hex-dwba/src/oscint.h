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

#ifndef HEX_OSCINT
#define HEX_OSCINT

#include <cmath>
#include <string>

#include "arrays.h"
#include "clenshawcurtis.h"

/**
 * \brief Integrator for oscillating functions.
 * 
 * The method works in the following way: The integration domain is split
 * into finite sections using the callback function WaveCallback. It is supposed
 * that this function returns bounding coordinates of these intervals for every
 * integer \f$ n \ge 1 \f$ supplied. Usually, one wavelength is chosen as the
 * measure and WaveCallback returns \f$ n \f$-multiple of the wavelength.
 * The integrand is integrated using Clenshaw-Curtis quadrature on every single
 * interval starting from the left-most. A sequence of partial sums of the
 * per-interval integrals is formed, accelerated using Aitken \f$ \Delta^2 \f$-process
 * and extended up to convergence.
 */
template <class Functor, class WaveCallback, class FType>
class OscillatingIntegral
{
	public:
		
		OscillatingIntegral(Functor f, WaveCallback w)
		    : F(f), W(w), Verbose(false), vName("OscillatingIntegral"),
		      nIntervals(1000), EpsRel(1e-8), EpsAbs(1e-8) {}
		
		bool verbose() const
			{ return Verbose; }
		void setVerbose(bool flag, std::string name = "OscillatingIntegral")
			{ Verbose = flag; vName = name; }
		
		int intervals() const
			{ return nIntervals; }
		void setIntervals(int inter)
			{ nIntervals = inter; }
		
		double epsilon() const
			{ return EpsRel; }
		void setEps(double eps)
			{ EpsRel = eps; }
		
		double tolerance() const
			{ return EpsAbs; }
		void setTol(double tol)
			{ EpsAbs = tol; }
		
		/**
		 * \brief The integration routine.
		 * 
		 * \param a Lower bound.
		 * \param b Upper bound.
		 * \param n On return, the evaluation count.
		 * 
		 * \return Value of the definite integral.
		 */
		FType integrate(double a, double b, int *n = nullptr) const
		{
			// skip zero-length interval
			if (a == b)
			{
				if (n != nullptr)
					*n = 0;
				return 0.;
			}
			
			// reverse bounds
			if (b < a)
				return -integrate(b, a, n);
			
			// partial sums
			Array<FType> psums;
			
			// accelerated partial sums
			Array<FType> accel;
			
			// Clenshaw-Curtis integrator
			ClenshawCurtis<decltype(F),FType> Q(F);
			Q.setEps(EpsRel);
			
			// loop over intervals
			double left = a, right;
			for (int n = 1; n <= nIntervals; n++)
			{
				// get right bound
				right = W(n);
				if (not std::isfinite(right) or right > b)
				{
					// forced termination by the callback control
					break;
				}
				
				// copy absolute tolerance to the subordinate integrator
				if (std::isfinite(a) and std::isfinite(b))
					Q.setTol((right-left)*EpsAbs/(b-a));
				else
					Q.setTol(EpsAbs);
				
				// integrate using the Clenshaw-Curtis quadrature
				FType integ = Q.integrate(left,right);
				
				// add to the parial sums
				if (psums.size() > 0)
					psums.push_back(psums.back() + integ);
				else
					psums.push_back(integ);
				
				// accelerate the series
				if (psums.size() > 2)
				{
					FType f   = psums.back(2);
					FType fp1 = psums.back(1);
					FType fp2 = psums.back(0);
					accel.push_back(fp2 - (fp2-fp1)*(fp2-fp1)/(fp2-2.*fp1+f));
					
					// check (relative) convergence on the accelerated series
					if (std::abs(accel.back(0)-accel.back(1)) < EpsRel * std::abs(accel.back()))
						return accel.back();
					
					// check (absolute) convergence on the accelerated series
					if (std::isfinite(a) and std::isfinite(b) and std::abs(accel.back()) < EpsAbs * (b - a))
						return accel.back();
				}
				
				// advance left bound
				left = right;
			}
			return 0;
		}
		
	private:
		
		/// Function to integrate.
		Functor F;
		
		/// Callback function for determining the oscillation period.
		WaveCallback W;
		
		/// Verbosity control.
		bool Verbose;
		
		/// Verbosity identificator.
		std::string vName;
		
		/// Maximum number of sub-intervals.
		int nIntervals;
		
		/// Relative tolerance.
		double EpsRel;
		
		/// Absolute tolerance.
		double EpsAbs;
};

#endif
