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

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <cstdlib>
#include <vector>

#include <omp.h>

#include "angs.h"
#include "arrays.h"
#include "bspline.h"
#include "special.h"
#include "spmatrix.h"

cArray computeLambda (
	rArray const & kf, rArray const & ki,
	int maxell, int L, int Spin,
	int ni, int li, int mi,
	rArray const & Ei, int lf,
	cArray const & Pf_overlaps
) {
	// shorthands
	unsigned Nenergy = kf.size();                       // energy count
	Complex const * const t = &(Bspline::ECS().t(0));   // B-spline knots
	int order   = Bspline::ECS().order();               // B-spline order
	int Nspline = Bspline::ECS().Nspline();             // B-spline count
	int Nknot   = Bspline::ECS().Nknot();               // number of all knots
	int Nreknot = Bspline::ECS().Nreknot();             // number of real knots
	
	// for all energies, compute the radial factors
	cArray rads(Nenergy * (maxell + 1));
	for (unsigned ie = 0; ie < Nenergy; ie++)
	{
		// compose filename of the data file for this solution
		std::ostringstream oss;
		oss << "psi-" << L << "-" << Spin << "-" << ni << "-" << li << "-" << mi << "-" << Ei[ie] << ".hdf";
		
		// load the solution
		cArray solution;
		
		#pragma omp critical
		solution.hdfload(oss.str().c_str());
		
		// The cross section oscillates, so we will do some averaging
		// As recommended by Bartlett, we will compute several amplitudes
		// separated by π/(n*kf[ie]) near the R₀ turning point.
		double wavelength = M_PI / kf[ie];
		int samples = 10;
		double R0 = t[Nreknot - 1].real();
		
		// skip impact energies with undefined outgoing momentum
		if (isnan(kf[ie]))
			continue;
		
		for (int n = 1; n <= samples; n++)
		{
			// this is the evaluation point
			double eval_r = R0 - wavelength * n / samples;
			
			// determine knot
			int eval_knot = std::lower_bound (
				t,
				t + Nknot,
				Complex(eval_r, 0.),
				[](Complex a, Complex b) -> bool {
					return a.real() < b.real();
				}
			) - t;
			
			// evaluate j and dj at far radius
			cArray j_R0(maxell + 1);
			cArray dj_R0(maxell + 1);
			for (int l = 0; l <= maxell; l++)
			{
				//evaluate the functions
				j_R0[l] = ric_j(l, kf[ie] * eval_r);
				dj_R0[l] = kf[ie] * Complex(dric_j(l, kf[ie] * eval_r));
			}
			
			// evaluate B-splines and their derivatives at evaluation radius
			CooMatrix Bspline_R0(Nspline, 1), Dspline_R0(Nspline, 1);
			for (int ispline = 0; ispline < Nspline; ispline++)
			{
				Complex val;
				
				// evaluate B-spline
				val = Bspline::ECS().bspline(ispline, eval_knot-1, order, eval_r);
				if (val != 0.)
					Bspline_R0.add(ispline, 0, val);
				
				// evaluate B-spline derivative
				val = Bspline::ECS().dspline(ispline, eval_knot-1, order, eval_r);
				if (val != 0.)
					Dspline_R0.add(ispline, 0, val);
			}
			
			// evaluate Wronskians
			CooMatrix Wj[maxell + 1];
			for (int l = 0; l <= maxell; l++)
				Wj[l] = dj_R0[l] * Bspline_R0 - j_R0[l] * Dspline_R0;
				
			// we need "P_overlaps" to have a 'dot' method
			CooMatrix Sp(Nspline, 1, Pf_overlaps.begin());
			
			// compute radial factor
			#pragma omp parallel for
			for (int l = 0; l <= maxell; l++)
			{
				// we don't need to compute forbidden transition
				if (l < abs(lf-L) or l > lf + L)
					continue;
				
				// get correct solution (for this ang. mom.)
				cArrayView PsiSc (
					solution, 
					(lf * (maxell + 1) + l) * Nspline * Nspline, 
					Nspline * Nspline
				);
				
				rads[ie * (maxell + 1) + l] += Sp.transpose().dot(PsiSc).dot(Wj[l].todense()).todense()[0] / double(samples);
			}
		}
	}
	
	return rads;
}
