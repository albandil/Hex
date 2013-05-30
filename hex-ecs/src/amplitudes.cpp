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

#include "angs.h"
#include "arrays.h"
#include "bspline.h"
#include "chebyshev.h"
#include "clenshawcurtis.h"
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

cArrays computeXi(int maxell, int L, int Spin, int ni, int li, int mi, rArray const & Ei)
{
	cArrays results;
	
	// shorthands
	Complex const * const t = &(Bspline::ECS().t(0));   // B-spline knots
	int order   = Bspline::ECS().order();               // B-spline order
	int Nspline = Bspline::ECS().Nspline();             // B-spline count
	int Nreknot = Bspline::ECS().Nreknot();             // number of real knots
	
	// determine evaluation radius
	double rho = t[Nreknot-2].real();
	
	// auxiliary variables
	cArray B1(Nspline), dB1(Nspline), B2(Nspline), dB2(Nspline);
	gsl_sf_result F1, F1p, G1, G1p, F2, F2p, G2, G2p;
	double expF1, expG1, expF2, expG2, cos_alpha, sin_alpha;
	
	std::cout << "\n";
	
	// for all energies
	for (size_t ie = 0; ie < Ei.size(); ie++)
	{
		// compose filename of the data file for this solution
		std::ostringstream oss;
		oss << "psi-" << L << "-" << Spin << "-" << ni << "-" << li << "-" << mi << "-" << Ei[ie] << ".hdf";
		
		// load the solution
		cArray solution;
		#pragma omp critical
		solution.hdfload(oss.str().c_str());
		
		// for all angular states (triangle ℓ₂ ≤ ℓ₁)
		for (int l1 = 0; l1 <= maxell; l1++)
		for (int l2 = std::abs(l1-L); l2 <= l1; l2++)
		{
			std::cout << "\tEi[" << ie << "] = " << Ei[ie] << ", l1 = " << l1 << ", l2 = " << l2 << "\n";
			
			// create subset of the solution
			cArrayView PsiSc (
				solution,
				(l1 * (maxell + 1) + l2) * Nspline * Nspline, 
				Nspline * Nspline
			);
			
			// we want to approximate the following function f_{ℓ₁ℓ₂}^{LS}(k₁,k₂)
			auto fLSl1l2k1k2 = [&](double k1) -> Complex {
				
				if (k1 == 0 or k1*k1 >= Ei[ie])
					return 0.;
				
				// compute momentum of the other electron
				double k2 = sqrt(Ei[ie] - 1./(ni*ni) - k1*k1);
				
				// Xi integrand
				auto integrand = [&](double alpha) -> Complex {
					
					// precompute projectors
					cos_alpha = cos(alpha);
					sin_alpha = sin(alpha);
					
					// precompute coordinates
					double r1 = rho * cos_alpha;
					double r2 = rho * sin_alpha;
					
					// evaluate Coulomb wave functions
					gsl_sf_coulomb_wave_FG_e(-1./k1, k1*r1, l1, 0, &F1, &F1p, &G1, &G1p, &expF1, &expG1);
					gsl_sf_coulomb_wave_FG_e(-1./k2, k2*r2, l2, 0, &F2, &F2p, &G2, &G2p, &expF2, &expG2);
					
					double F1F2 = F1.val * F2.val;
					double ddrho_F1F2 = F1p.val*cos_alpha*F2.val + F1.val*F2p.val*sin_alpha;
					
					// get B-spline knots
					int iknot1 = Bspline::ECS().knot(r1);
					int iknot2 = Bspline::ECS().knot(r2);
					
					// evaluate the B-splines
					for (int ispline1 = std::max(0,iknot1-order); ispline1 <= iknot1; ispline1++)
					{
						B1[ispline1]  = Bspline::ECS().bspline(ispline1,iknot1,order,r1);
						dB1[ispline1] = Bspline::ECS().dspline(ispline1,iknot1,order,r1);
					}
					for (int ispline2 = std::max(0,iknot2-order); ispline2 <= iknot2; ispline2++)
					{
						B2[ispline2]  = Bspline::ECS().bspline(ispline2,iknot2,order,r2);
						dB2[ispline2] = Bspline::ECS().dspline(ispline2,iknot2,order,r2);
					}
					
					// evaluate the solution
					Complex Psi = 0., ddr1_Psi = 0., ddr2_Psi = 0., ddrho_Psi = 0.;
					for (int ispline1 = std::max(0,iknot1-order); ispline1 <= iknot1; ispline1++)
					for (int ispline2 = std::max(0,iknot2-order); ispline2 <= iknot2; ispline2++)
					{
						int idx = ispline1 * Nspline + ispline2;
						
						Psi      += PsiSc[idx] *  B1[ispline1] *  B2[ispline2];
						ddr1_Psi += PsiSc[idx] * dB1[ispline1] *  B2[ispline2];
						ddr2_Psi += PsiSc[idx] *  B1[ispline1] * dB2[ispline2];
					}
					ddrho_Psi = ddr1_Psi * cos_alpha + ddr2_Psi * sin_alpha;
					
					// evaluate the integrand
					return F1F2*ddrho_Psi - Psi*ddrho_F1F2;
				};
				
				/// DEBUG
// 				std::ofstream ofs("integrand.dat");
// 				for (int ia = 0; ia < 1000; ia++)
// 				{
// 					double alpha = 0.5 * M_PI * ia / 1000;
// 					Complex v = integrand(alpha);
// 					ofs << alpha << "\t" << v.real() << "\t" << v.imag() << "\n";
// 				}
// 				ofs.close();
// 				exit(0);
				///
				
				// integrator
				ClenshawCurtis<decltype(integrand),Complex> Q(integrand);
				Q.setEps(1e-6);
				Complex res = 2. * k1 * k2 * rho * Q.integrate(0., 0.5 * M_PI) / (sqrt(M_PI));
				
// 				std::cout << k1 << "\t" << res << "\n";
				
				return res;
				
			};
			
			/// DEBUG
// 			std::ofstream ofs("fLSl1l2k1k2.dat");
// 			double kmax = sqrt(Ei[ie] - 1./(ni*ni));
// 			for (int ik = 0; ik < 1000; ik++)
// 			{
// 				double k = kmax * ik / 1000.;
// 				Complex f = fLSl1l2k1k2(k);
// 				ofs << k/kmax << "\t" << k*sqrt(kmax*kmax-k*k) << "\t" << f.real() << "\t" << f.imag() << "\n";
// 			}
// 			ofs.close();
// 			exit(0);
			///
			
			// Chebyshev approximation
			Chebyshev<double,Complex> CB;
			
			// convergence loop
			for (int N = 4; ; N *= 2)
			{
				
// 				std::cout << "N = " << N << "\n";
				
				// build the approximation
				CB.generate(fLSl1l2k1k2, N, 0., sqrt(Ei[ie] - 1./(ni*ni)));
				
				/// DEBUG
// 				std::ostringstream os;
// 				os << "cb" << N << ".hdf";
// 				CB.coeffs().hdfsave(os.str().c_str());
				///
				
				// check tail
				if (CB.tail(1e-5) != N)
					break;
				
// 				std::cout << CB.str() << "\n";
				
				// limit subdivision
				if (N > 32768)
					throw exception("ERROR: Non-convergent Chebyshev expansion.");
			}
			
// 			exit(0);
			
			results.push_back(CB.coeffs());
		}
	}
	
	return results;
}
