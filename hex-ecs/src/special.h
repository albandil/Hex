/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2012                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef _SPECIAL_H_
#define _SPECIAL_H_

#include <algorithm>
#include <complex>
#include <vector>

#include <gsl/gsl_sf.h>

#include "arrays.h"
#include "bspline.h"
#include "gauss.h"

/** Riccati-Bessel function
 * 
 * Evaluate Riccati-Bessel function for complex argument. Function is not suitable for
 * large degrees, it uses the most naïve (and least stable) evaluation method.
 * Starting from the expressions for zeroth and first Riccati-Bessel function
 * \f[
 *      j_0(z) = \sin z, \qquad j_1(z) = \frac{\sin z}{z} - \cos z
 * \f]
 * the function employs the forward(!) recurrence relation
 * \f[
 *      j_{n+1}(z) = \frac{2n+1}{z} j_n(z) - j_{n-1}(z) .
 * \f]
 * 
 * \param n Degree of the Riccati-Bessel function.
 * \param z Complex argument.
 */
LComplex ric_j(int n, LComplex z);

/** Derivative of Riccati-Bessel function
 * 
 * \param n Degree of the function.
 * \param z Complex argument.
 */
LComplex dric_j(int n, LComplex z);

/** Hydrogen radial function (radius-multiplied)
 * Evaluate hydrogen radial function (radius-multiplied) for complex argument
 * \param n Principal quantum number.
 * \param l Orbital quantum number.
 * \param z Radius: real or complex argument.
 */
Complex hydro_P(unsigned n, unsigned l, Complex z);

/** Derivative of Pnl.
 * \param n
 * \param l
 * \param z
 */
Complex dhydro_P(unsigned n, unsigned l, Complex z);

/** Compute P-overlaps
 * Compute overlap vector of B-splines vs. hydrogen Pnl function.
 * \param n Principal quantum number.
 * \param l Orbital quantum number.
 * \param weightf Weight function to multiply every value of the hydrogenic function.
 *                It is expected to have the "double operator() (Complex z)" interface,
 *                where the sent value is the complex coordinate.
 */
template <class Functor> cArray overlapP(int n, int l, Functor weightf)
{
	cArray res(Nspline);
	
	// per interval
	int points = 20;
	
	// evaluated B-spline and hydrogenic functions (auxiliary variables)
	cArray evalB(points);
	cArray evalP(points);
	
	// for all knots
	for (int iknot = 0; iknot < Nknot - 1; iknot++)
	{
		// skip zero length intervals
		if (t[iknot] == t[iknot+1])
			continue;
		
		// which points are to be used here?
		cArray xs = p_points(points, t[iknot], t[iknot+1]);
		cArray ws = p_weights(points, t[iknot], t[iknot+1]);
		
		// evaluate the hydrogenic function
		std::transform(
			xs.begin(), xs.end(), evalP.begin(),
			[ = ](Complex x) -> Complex {
				gsl_sf_result R;
				if (gsl_sf_hydrogenicR_e(n, l, 1., x.real(), &R) == GSL_EUNDRFLW)
					return 0.;
				else
					return weightf(x) * x * R.val;
			}
		);
		
		// for all relevant B-splines
		for (int ispline = std::max((int)iknot-(int)order,0); ispline < Nspline and ispline <= iknot; ispline++)
		{
			// evaluate the B-spline
			B(ispline, iknot, points, xs.data(), evalB.data());
			
			// sum with weights
			Complex sum = 0.;
			for (int ipoint = 0; ipoint < points; ipoint++)
				sum += ws[ipoint] * evalP[ipoint] * evalB[ipoint];
			
			// store the overlap
			res[ispline] += sum;
		}
	}
	
	return res;
}

/** Compute j-overlaps
 * Compute overlap integrals for Riccati-Bessel function.
 * \param maxL2 Maximal degree of the Riccati-Bessel function.
 * \param vk Array containing linear momenta.
 * \param weightf Weight function to multiply every value of the Bessel function.
 *                It is expected to have the "double operator() (Complex z)" interface,
 *                where the sent value is the complex coordinate.
 * \return Array of shape [vk.size() × (maxL2+1) × Nspline] in column-major format.
 */
template <class Functor> cArray overlapj(int maxL2, std::vector<double> vk, Functor weightf)
{
	// shorthand for energy count
	int Nenergy = vk.size();
	
	// reserve space for the output array
	size_t size = Nspline * Nenergy * (maxL2 + 1);
	cArray res(size);
	
	// per interval
	int points = 20;
		
	// for all knots
	# pragma omp parallel for
	for (int iknot = 0; iknot < Nknot - 1; iknot++)
	{
		// skip zero length intervals
		if (t[iknot] == t[iknot+1])
			continue;
		
		// which points are to be used here?
		cArray xs = p_points(points, t[iknot], t[iknot+1]);
		cArray ws = p_weights(points, t[iknot], t[iknot+1]);
		
		// for all linear momenta (= energies)
		for (int ie = 0; ie < Nenergy; ie++)
		{
			// for all angular momenta
			for (int l = 0; l <= maxL2; l++)
			{
				// evaluate the Riccati-Bessel function
				std::vector<LComplex> evalj(points);
				std::transform(
					xs.begin(), xs.end(), evalj.begin(),
					[ = ](Complex x) -> LComplex {
						return LComplex(weightf(x)) * ric_j(l, vk[ie] * x);
					}
				);
				
				// for all relevant B-splines
				for (int ispline = std::max((int)iknot-(int)order,0); ispline < Nspline and ispline <= iknot; ispline++)
				{
					// evaluate the B-spline
					cArray evalB(points);
					B(ispline, iknot, points, xs.data(), evalB.data());
					
					// sum with weights
					LComplex sum = 0.;
					for (int ipoint = 0; ipoint < points; ipoint++)
						sum += LComplex(ws[ipoint]) * evalj[ipoint] * LComplex(evalB[ipoint]);
					
					// store the overlap; keep the shape Nmomenta × Nspline × (maxl+1)
					res[(ie * (maxL2 + 1) + l) * Nspline + ispline] += Complex(sum);
				}
			}
		}
	}
	
	return res;
}

#endif
