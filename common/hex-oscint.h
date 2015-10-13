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

#ifndef HEX_OSCINT
#define HEX_OSCINT

#include <cmath>
#include <string>

#include "hex-arrays.h"
#include "hex-clenshawcurtis.h"

/**
 * \brief Integrator for oscillating functions.
 * 
 * The method works in the following way: The integration domain is split
 * into finite sections using the callback function WaveCallback. It is supposed
 * that this function returns bounding coordinates of these intervals for every
 * integer \f$ n \ge 1 \f$ supplied. Usually, one wavelength is chosen as the
 * measure and WaveCallback returns \f$ 1.5n \f$-multiple of the wavelength.
 * Creating intervals of length equal to 1.5×λ makes the sequence of contributions
 * alternate and <b>dramatically</b> improves the convergence (e.g. 70 times
 * less intervals have to be integrated for sin(x)/x when Aitken Δ²-process is
 * used twice). This reducion also lowers possible rounding errors.
 * 
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
		      nIntervals(1000), EpsRel(1e-8), EpsAbs(1e-8), ThrowAll(true) {}
		
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
		
		bool throwall() const
			{ return ThrowAll; }
		void setThrowAll(bool flag)
			{ ThrowAll = flag; }
		
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
			
			if (Verbose)
				std::cout << "[" << vName << "] integrate " << a << " to " << b << "\n";
			
			// reverse bounds
			if (b < a)
				return -integrate(b, a, n);
			
			// partial sums
			Array<FType> psums;
			
			// accelerated partial sums
			Array<FType> accel1, accel2;
			
			// Clenshaw-Curtis integrator
			ClenshawCurtis<decltype(F),FType> Q(F);
			Q.setEps(EpsRel);
			Q.setVerbose(false, "OscillatingIntegral");
			Q.setRec(false);
			Q.setSubdiv(15);
			Q.setThrowAll(ThrowAll);
			
			// loop over intervals
			double left = a, right;
			int k;
			for (k = 1; k <= nIntervals and left < b; k++)
			{
				// get right bound
				right = W(a,k);
				if (not std::isfinite(right))
				{
					// forced termination by the callback control
					break;
				}
				if (right < left)
				{
					// wrong supplied period
					continue;
				}
				if (right > b)
				{
					// do not step out of integration domain
					right = b;
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
				
				// accelerate the series for the first time
				if (psums.size() > 2)
				{
					FType f   = psums.back(2);
					FType fp1 = psums.back(1);
					FType fp2 = psums.back(0);
					accel1.push_back(fp2 - (fp2-fp1)*(fp2-fp1)/(fp2-2.*fp1+f));
				}
				
				// accelerate the series for the second time
				if (accel1.size() > 2)
				{
					FType f   = accel1.back(2);
					FType fp1 = accel1.back(1);
					FType fp2 = accel1.back(0);
					accel2.push_back(fp2 - (fp2-fp1)*(fp2-fp1)/(fp2-2.*fp1+f));
					
					// check (relative) convergence on the accelerated series
					if (std::abs(accel2.back(0)-accel2.back(1)) < EpsRel * std::abs(accel2.back()))
					{
						if (n != nullptr) *n = k;
						return accel2.back();
					}
					
					// check (absolute) convergence on the accelerated series
					if (std::isfinite(a) and std::isfinite(b) and std::abs(accel2.back()) < EpsAbs * (b - a))
					{
						if (n != nullptr) *n = k;
						return accel2.back();
					}
				}
				
				if (Verbose)
				{
					std::cout << "[" << vName << "] " << k << " " 
						<< left << " " << right << " " 
						<< psums.back();
						
					if (accel1.size() > 0)
						std::cout << " " << accel1.back();
					if (accel2.size() > 0)
						std::cout << " " << accel2.back();
					
					std::cout << "\n";
				}
				
				// advance left bound
				left = right;
			}
			
			// store iterations
			if (n != nullptr) *n = k;
			
			// out of integration domain: return sum
			return psums.back();
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
		
		/// Whether to be choleric about errors (and throw every time).
		bool ThrowAll;
};

#endif
