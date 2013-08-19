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
#include <regex>
#include <vector>

#include <fftw3.h>

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
	int maxell, int L, int Spin, int Pi,
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
		oss << "psi-" << L << "-" << Spin << "-" << Pi << "-" << ni << "-" << li << "-" << mi << "-" << Ei[ie] << ".hdf";
		
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
			
			// evaluate j and dj at far radius for all angular momenta up to maxell
			cArray j_R0 = ric_jv(maxell, kf[ie] * eval_r);
			cArray dj_R0 = dric_jv(maxell, kf[ie] * eval_r) * kf[ie];
			
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

Chebyshev<double,Complex> fcheb(cArrayView const & PsiSc, double kmax, int l1, int l2)
{
	// shorthands
	Complex const * const t = &(Bspline::ECS().t(0));   // B-spline knots
	int Nspline = Bspline::ECS().Nspline();             // number of real knots
	int Nreknot = Bspline::ECS().Nreknot();             // number of real knots
	int order   = Bspline::ECS().order();               // B-spline order
	
	// determine evaluation radius
	char const * HEX_RHO = getenv("HEX_RHO");
 	double rho = (HEX_RHO == nullptr) ? t[Nreknot-2].real() : atof(HEX_RHO);
	
	// we want to approximate the following function f_{ℓ₁ℓ₂}^{LS}(k₁,k₂)
	auto fLSl1l2k1k2 = [&](double k1) -> Complex {
		
		if (k1 == 0 or k1*k1 >= kmax*kmax)
			return 0.;
		
		// compute momentum of the other electron
		double k2 = sqrt(kmax*kmax - k1*k1);
		
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
			
			double F1F2 = F1 * F2;
			double ddrho_F1F2 = 0.;
			
			if (cos_alpha != 0.)
				ddrho_F1F2 += k1*F1p*cos_alpha*F2;
			if (sin_alpha != 0.)
				ddrho_F1F2 += k2*F1*F2p*sin_alpha;
			
			// get B-spline knots
			int iknot1 = Bspline::ECS().knot(r1);
			int iknot2 = Bspline::ECS().knot(r2);
			
			// auxiliary variables
			cArray B1(Nspline), dB1(Nspline), B2(Nspline), dB2(Nspline);
			
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
		
		// integrator
		ClenshawCurtis<decltype(integrand),Complex> Q(integrand);
		Q.setEps(1e-6);
		Complex res = 2. * rho * Q.integrate(0., 0.5 * M_PI) / sqrt(M_PI);
		
		return res;
		
	};
	
	// Chebyshev approximation
	Chebyshev<double,Complex> CB;
	
	// convergence loop
	for (int N = 4; ; N *= 2)
	{
		// build the approximation
		CB.generate(fLSl1l2k1k2, N, 0., kmax);
		
		// check tail
		if (CB.tail(1e-5) != N)
			break;
		
		// limit subdivision
		if (N > 32768)
			throw exception("ERROR: Non-convergent Chebyshev expansion.");
	}
	
	return CB;
}

cArrays computeXi(int maxell, int L, int Spin, int Pi, int ni, int li, int mi, rArray const & Ei, rArray & ics, std::vector<std::pair<int,int>> const & coupled_states)
{
	// resize and clear the output storage for integral cross sections
	ics.resize(Ei.size());
	ics.clear();
	
	// array of Chebyshev expansions for every Ei and angular state
	cArrays results;
	
	// B-spline count
	int Nspline = Bspline::ECS().Nspline();
	
	// for all energies
	for (size_t ie = 0; ie < Ei.size(); ie++)
	{
		// compose filename of the data file for this solution
		std::ostringstream oss;
		oss << "psi-" << L << "-" << Spin << "-" << Pi << "-" << ni << "-" << li << "-" << mi << "-" << Ei[ie] << ".hdf";
		
		// load the solution
		cArray solution;
		#pragma omp critical
		if (not solution.hdfload(oss.str().c_str()))
			throw exception ("Can't open the solution file \"%s\"!", oss.str().c_str());
		
		// maximal available momentum
		double kmax = sqrt(Ei[ie] - 1./(ni*ni));
		
		// for all angular states ???: (triangle ℓ₂ ≤ ℓ₁)
		for (unsigned ill = 0; ill < coupled_states.size(); ill++)
		{
			int l1 = coupled_states[ill].first;
			int l2 = coupled_states[ill].second;
			
			std::cout << "\tEi[" << ie << "] = " << Ei[ie] << ", l1 = " << l1 << ", l2 = " << l2 << "\n";
			
			// create subset of the solution
			cArrayView PsiSc (solution, ill * Nspline * Nspline, Nspline * Nspline);
			
			// compute new ionization amplitude
			Chebyshev<double,Complex> CB = fcheb(PsiSc, kmax, l1, l2);
			results.push_back(CB.coeffs());
			
			// integrate the expansion
			int tail = CB.tail(1e-10);
			int n;
			auto fsqr = [&](double beta) -> double { return sqrabs(CB.clenshaw(kmax*sin(beta), tail)); };
			ClenshawCurtis<decltype(fsqr),double> integrator(fsqr);
			double cs = integrator.integrate(0, 0.25 * M_PI, &n) / sqrt(Ei[ie]);
			std::cout << "\t\t- contrib to ics: " << cs << " (" << n << " evaluations)\n";
			ics[ie] += cs;
		}
	}
	
	return results;
}

void TDCS (
	std::string solutionfile, std::vector<std::pair<int,int>> const & coupled_states,
	double kmax, int ni, int L, int S
) {
	std::cout << "\nEvaluating ionization amplitudes for use in TDCS.\n";
	std::cout << "\tusing file \"" << solutionfile << "\"\n";
	std::cout << "\tusing total energy " << kmax*kmax << " Ry\n";
	
	// state information... C++11 regex doesn't work in GCC < 4.9
	int li = 0, mi = 0;
	double Ei = 0;
	
#if __GNUC__ > 4 && __GNUC_MINOR__ > 8
	// try to match the solution file to Hex-ecs naming scheme
	std::regex e("psi-(+\\d)-(+\\d)-(+\\d)-(+\\d)-(+\\d)-(+[\\.\\d])\\.hdf");
	
	std::smatch sm;
	std::regex_match(solutionfile, sm, e);
	std::cout << "c\n";
	if (sm.size() == 7)
	{
		std::cout << "Assumed solution data:\n";
		ni = atoi(sm[1].str().c_str()); std::cout << "\tni = " << ni << "\n";
		li = atoi(sm[2].str().c_str()); std::cout << "\tli = " << li << "\n";
		mi = atoi(sm[3].str().c_str()); std::cout << "\tmi = " << mi << "\n";
		L  = atoi(sm[4].str().c_str()); std::cout << "\tL  = " << L  << "\n";
		S  = atoi(sm[5].str().c_str()); std::cout << "\tS  = " << S  << "\n";
		Ei = atof(sm[6].str().c_str()); std::cout << "\tEi = " << Ei << "\n";
	}
#endif
	
	// output file
	std::string fsqlname = solutionfile + ".sql";
	std::ofstream fsql(fsqlname.c_str());
	fsql << "BEGIN TRANSACTION;\n";
	
	// B-spline count
	int Nspline = Bspline::ECS().Nspline();
	
	// load the solution
	cArray solution;
	#pragma omp critical
	if (not solution.hdfload(solutionfile.c_str()))
		throw exception ("Can't open the solution file \"%s\"!", solutionfile.c_str());
	
	// for all angular states ???: (triangle ℓ₂ ≤ ℓ₁)
	for (unsigned ill = 0; ill < coupled_states.size(); ill++)
	{
		int l1 = coupled_states[ill].first;
		int l2 = coupled_states[ill].second;
		
		std::cout << "\tEmax = " << kmax*kmax << ", l1 = " << l1 << ", l2 = " << l2 << "\n";
		
		// create subset of the solution
		cArrayView PsiSc (solution, ill * Nspline * Nspline, Nspline * Nspline);
		
		// compute new ionization amplitude
		cArray coeffs = fcheb(PsiSc, kmax, l1, l2).coeffs();
		
		// save data as BLOBs
		fsql << "INSERT OR REPLACE INTO \"ionf\" VALUES ("
				<< ni << "," << li << "," << mi << ","
				<< L  << "," << S  << ","
				<< Ei << "," << coupled_states[ill].first << ","
				<< coupled_states[ill].second << ","
				<< coeffs.toBlob() << ");\n";
	}
	
	fsql << "COMMIT;\n";
	fsql << std::flush;
	fsql.close();
	
	// print some info
	std::cout << "Triple-differential cross section data written to the file \"" << fsqlname << "\"\n";

#if __GNUC__ > 4 && __GNUC_MINOR__ > 8
	if (sm.size() != 7)
	{
#endif
		std::cout << "Please be aware that it does not contain state data (ni,li,mi,L,S,Ei) and these must be added manually!\n";
#if __GNUC__ > 4 && __GNUC_MINOR__ > 8
	}
#endif
}
