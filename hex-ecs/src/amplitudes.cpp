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

#include <fftw3.h>

#include "angs.h"
#include "arrays.h"
#include "bspline.h"
#include "chebyshev.h"
#include "clenshawcurtis.h"
#include "moments.h"
#include "specf.h"
#include "spmatrix.h"

bool debug = false;

cArray computeLambda (
	rArray const & kf, rArray const & ki,
	int maxell, int L, int Spin,
	int ni, int li, int mi,
	rArray const & Ei, int lf,
	cArray const & Pf_overlaps,
	std::vector<std::pair<int,int>> const & coupled_states
) {
	// shorthands
	unsigned Nenergy = kf.size();                       // energy count
	Complex const * const t = &(Bspline::ECS().t(0));   // B-spline knots
	int order   = Bspline::ECS().order();               // B-spline order
	int Nspline = Bspline::ECS().Nspline();             // B-spline count
	int Nknot   = Bspline::ECS().Nknot();               // number of all knots
	int Nreknot = Bspline::ECS().Nreknot();             // number of real knots
	
	cArray rads(Nenergy * (maxell + 1));
	
	// for all energies, compute the radial factors
	for (unsigned ie = 0; ie < Nenergy; ie++)
	{
		// compose filename of the data file for this solution
		std::ostringstream oss;
		oss << "psi-" << L << "-" << Spin << "-" << ni << "-" << li << "-" << mi << "-" << Ei[ie] << ".hdf";
		
		// load the solution
		cArray solution;
		
		#pragma omp critical
		{
			if (not solution.hdfload(oss.str().c_str()))
				throw exception ("Failed to load \"%s\"\n", oss.str().c_str());
		}
		
		// The cross section oscillates, so we will do some averaging
		// As recommended by Bartlett, we will compute several amplitudes
		// separated by π/(n*kf[ie]) near the R₀ turning point.
		double wavelength = M_PI / kf[ie];
		char const * HEX_RHO = getenv("HEX_RHO");
		char const * HEX_SAMPLES = getenv("HEX_SAMPLES");
		int samples = (HEX_SAMPLES == nullptr) ? 10 : atoi(HEX_SAMPLES);
		double R0 = (HEX_RHO == nullptr) ? t[Nreknot - 1].real() : atof(HEX_RHO);
		
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
			for (unsigned ill = 0; ill < coupled_states.size(); ill++)
			{
				// use blocks that result in requested final lf
				int l1 = coupled_states[ill].first;
				int l2 = coupled_states[ill].second;
				if (l1 != lf)
					continue;
				
				// get correct solution (for this ang. mom.)
				cArrayView PsiSc (solution, ill * Nspline * Nspline, Nspline * Nspline);
				
				rads[ie * (maxell + 1) + l2] += Sp.transpose().dot(PsiSc).dot(Wj[l2].todense()).todense()[0] / double(samples);
			}
		}
	}
	
	return rads;
}

cArrays computeXi(int maxell, int L, int Spin, int ni, int li, int mi, rArray const & Ei, rArray & ics, std::vector<std::pair<int,int>> const & coupled_states)
{
	// resize and clear the output storage for integral cross sections
	ics.resize(Ei.size());
	ics.clear();
	
	// array of Chebyshev expansions for every Ei and angular state
	cArrays results;
	
	// shorthands
	Complex const * const t = &(Bspline::ECS().t(0));   // B-spline knots
	int order   = Bspline::ECS().order();               // B-spline order
	int Nspline = Bspline::ECS().Nspline();             // B-spline count
	int Nreknot = Bspline::ECS().Nreknot();             // number of real knots
	
	// determine evaluation radius
	char const * HEX_RHO = getenv("HEX_RHO");
 	double rho = (HEX_RHO == nullptr) ? t[Nreknot-2].real() : atof(HEX_RHO);
	
	// auxiliary variables
	cArray B1(Nspline), dB1(Nspline), B2(Nspline), dB2(Nspline);
	
	// for all energies
	for (size_t ie = 0; ie < Ei.size(); ie++)
	{
		// compose filename of the data file for this solution
		std::ostringstream oss;
		oss << "psi-" << L << "-" << Spin << "-" << ni << "-" << li << "-" << mi << "-" << Ei[ie] << ".hdf";
		
		// load the solution
		cArray solution;
		#pragma omp critical
		if (not solution.hdfload(oss.str().c_str()))
			throw exception ("Can't open the solution file \"%s\"!", oss.str().c_str());
		
		// for all angular states ???: (triangle ℓ₂ ≤ ℓ₁)
		for (unsigned ill = 0; ill < coupled_states.size(); ill++)
		{
			int l1 = coupled_states[ill].first;
			int l2 = coupled_states[ill].second;
			
			std::cout << "\tEi[" << ie << "] = " << Ei[ie] << ", l1 = " << l1 << ", l2 = " << l2 << "\n";
			
			// create subset of the solution
			cArrayView PsiSc (solution, ill * Nspline * Nspline, Nspline * Nspline);
			
			/// DEBUG
// 			write_2D_data (
// 				Nspline, Nspline,
// 				"psi.dat",
// 				[&](int i, int j) -> double {
// 					return PsiSc[i*Nspline+j].real();
// 				}
// 			);
// 			bool debug = true;
			///
			
			// we want to approximate the following function f_{ℓ₁ℓ₂}^{LS}(k₁,k₂)
			auto fLSl1l2k1k2 = [&](double k1) -> Complex {
				
				if (k1 == 0 or k1*k1 >= Ei[ie])
					return 0.;
				
				// compute momentum of the other electron
				double k2 = sqrt(Ei[ie] - 1./(ni*ni) - k1*k1);
				
				// Xi integrand
				auto integrand = [&](double alpha) -> Complex {
					
					// precompute projectors
					double cos_alpha = (alpha == 0.5 * M_PI) ? 0. : cos(alpha);
					double sin_alpha = sin(alpha);
					
					// precompute coordinates
					double r1 = rho * cos_alpha;
					double r2 = rho * sin_alpha;
					
					// evaluate Coulomb wave functions and derivatives
					double F1, F2, F1p, F2p;
					int err1 = coul_F(l1,k1,r1, F1,F1p);
					int err2 = coul_F(l2,k2,r2, F2,F2p);
					
					/// DEBUG
					if (err1 != GSL_SUCCESS or err2 != GSL_SUCCESS)
					{
						std::cerr << "Errors while evaluating Coulomb function:\n";
						std::cerr << "\terr1 = " << err1 << "\n";
						std::cerr << "\terr2 = " << err2 << "\n";
						exit(-1);
					}
					///
					
					/// DEBUG
// 					double expF, expG;
// 					gsl_sf_result f1, g1, f1p, g1p, f2, g2, f2p, g2p;
// 					int gslerr1, gslerr2;
// 					gslerr1 = gsl_sf_coulomb_wave_FG_e(-1/k1, k1*r1, l1, 0, &f1, &f1p, &g1, &g1p, &expF, &expG);
// 					F1 = F.val; F1p = Fp.val;
// 					if (F1 != f1.val)
// 					{
// 						std::cerr << "coul_F: " << F1 << "," << F1p << "\n";
// 						std::cerr << "GSL: " << f1.val << "," << f1p.val << "\n";
// 					}
// 					gslerr2 = gsl_sf_coulomb_wave_FG_e(-1/k2, k2*r2, l2, 0, &f2, &f2p, &g2, &g2p, &expF, &expG);
// 					F2 = F.val; F2p = Fp.val;
// 					if (F2 != f2.val)
// 					{
// 						std::cerr << "coul_F: " << F2 << "," << F2p << "\n";
// 						std::cerr << "GSL: " << f2.val << "," << f2p.val << "\n";
// 					}
					///
					
					double F1F2 = F1 * F2;
					double ddrho_F1F2 = 0.;
					
// 					std::cerr << "F1p[" << l1 << "," << -1/k1 << "](" << k1*r1 << ") = " << F1p << " (" << f1p.val << ", err " << gslerr1 << ")\n";
// 					std::cerr << "F2p[" << l2 << "," << -1/k2 << "](" << k2*r2 << ") = " << F2p << " (" << f2p.val << ", err " << gslerr2 << ")\n";
					
					if (cos_alpha != 0.)
						ddrho_F1F2 += k1*F1p*cos_alpha*F2;
					if (sin_alpha != 0.)
						ddrho_F1F2 += k2*F1*F2p*sin_alpha;
					
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
					
					if (cos_alpha != 0.)
						ddrho_Psi += ddr1_Psi * cos_alpha;
					if (sin_alpha != 0.)
						ddrho_Psi += ddr2_Psi * sin_alpha;
					
					/// DEBUG
					if (not finite(F1F2))
						std::cerr << "F1F2 = " << F1F2 << "\n";
					if (not finite(std::abs(ddrho_Psi)))
						std::cerr << "ddrho_Psi = " << ddrho_Psi << "\n";
					if (not finite(std::abs(Psi)))
						std::cerr << "Psi = " << Psi << "\n";
					if (not finite(ddrho_F1F2))
						std::cerr << "ddrho_F1F2 = " << ddrho_F1F2 << "\n";
					///
					
					// evaluate the integrand
					return F1F2*ddrho_Psi - Psi*ddrho_F1F2;
				};
				
// 				std::cout << integrand(2.45436e-05) << "\n";
// 				exit(0);
				
				// integrator
				ClenshawCurtis<decltype(integrand),Complex> Q(integrand);
				Q.setEps(1e-6);
				Complex res = 2. * rho * Q.integrate(0., 0.5 * M_PI) / sqrt(M_PI);
				
				/// DEBUG
// 				if (debug)
// 				{
// 					std::ofstream ofs("integrand.dat");
// 					for (int ia = 0; ia < 1000; ia++)
// 					{
// 						double alpha = 0.0122718 /* 0.5 * M_PI */ * ia / 1000;
// 						Complex v = integrand(alpha);
// 						ofs << alpha << "\t" << v.real() << "\t" << v.imag() << "\n";
// 					}
// 					ofs.close();
// 					Q.setVerbose(true);
// 					Q.setStack(10);
// 					Complex inte = Q.integrate(0.,0.0122718);
// 					std::cerr << "Integrated value: " << inte << "\n";
// 					exit(0);
// 				}
				///
				
				return res;
				
			};
			
			/// DEBUG
// 			debug = true;
// 			fLSl1l2k1k2(sqrt(0.5*(Ei[ie]-1./(ni*ni))));
			///
			
			/// DEBUG
// 			std::ostringstream name;
// 			name << "f_" << l1 << "_" << l2 << ".dat";
// 			std::ofstream ofs(name.str().c_str());
// 			double kmax = sqrt(Ei[ie] - 1./(ni*ni));
// 			for (int ik = 0; ik < 1000; ik++)
// 			{
// 				double k = kmax * ik / 1000.;
// 				Complex f = fLSl1l2k1k2(k);
// 				
// 				if (std::abs(f) > 1e+5)
// 				{
// 					debug = true;
// 					fLSl1l2k1k2(k);
// 					abort();
// 				}
// 				
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
				// build the approximation
				CB.generate(fLSl1l2k1k2, N, 0., sqrt(Ei[ie] - 1./(ni*ni)));
				
				/// DEBUG
// 				std::ostringstream os;
// 				os << "cb_" << l1 << "_" << l2 << "_" << N << ".hdf";
// 				CB.coeffs().hdfsave(os.str().c_str());
				///
				
				// check tail
				if (CB.tail(1e-5) != N)
					break;
				
				// limit subdivision
				if (N > 32768)
					throw exception("ERROR: Non-convergent Chebyshev expansion.");
			}
			
			results.push_back(CB.coeffs());
			
			//
			// integrate the expansion
			//
			
			// setup FFTW
			int N = CB.coeffs().size();
			cArray mirror_coeffs(4*N+1), evalf(4*N);
			fftw_plan plan = fftw_plan_dft_1d (
				4*N,
				reinterpret_cast<fftw_complex*>(&mirror_coeffs[0]),
				reinterpret_cast<fftw_complex*>(&evalf[0]),
				FFTW_FORWARD,
				0
			);
			
			// mirror oddly around N (3N), evenly around 2N (0,4N)
			for (int k = 0; k < N; k++)
			{
				mirror_coeffs[4*N-k] = mirror_coeffs[k] = CB.coeffs()[k];
				mirror_coeffs[2*N+k] = mirror_coeffs[2*N-k] = -CB.coeffs()[k];
			}
			
			// integrate
			//    ₁                       n/2-1                
			//   ⌠              dx     2π ===  |       2j+1     |²
			// 2 ⎮ |f(|x|)|² ——————— = —— >    | f(cos(———— π)) |  
			//   ⌡           √(1-x²)   n  ===  |        2n      |
			//  ⁰                         j=0
			// where
			//                        N-1
			//        2j+1       c₀   ===         j+½
			//  f(cos(———— π)) = —— + >   ck cos( ——— kπ )
			//         2n        2    ===          n
			//                        k=1
			// can be evaluated by DCT-III (inverse DCT-II) if full precision is used,
			// i.e. n = N; the result will be stored in odd elements (multiplied by 4)
			
			// evaluate the function using the FFT
			fftw_execute(plan);
			fftw_destroy_plan(plan);
			
			// sum contributions
			double cs = 0;
			for (int j = 0; j < N/2; j++)
				cs += sqrabs(evalf[2*j+1]);            // (FFTW magic) odd elements only
			cs *= 0.0625 * M_PI / CB.coeffs().size();  // (FFTW magic) 1/4²
			
			std::cout << "\t\t- contrib to ics: " << cs << "\n";
			ics[ie] += cs;
		}
	}
	
	return results;
}
