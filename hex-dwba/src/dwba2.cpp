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

#include <cstdio>
#include <cmath>

#include "angs.h"
#include "complex.h"
#include "dwba2.h"
#include "hydrogen.h"
#include "integrate.h"
#include "potential.h"
#include "wave_distort.h"
#include "spmatrix.h"

#define EPS_CONTRIB 1e-8

double DWBA2::computeI (
	int lam_a,
	HydrogenFunction const & psi_a,
	DistortedWave const & chi_a,
	DistortingPotential const & U_a,
	SturmianFunction const & sturm1,
	SturmianFunction const & sturm2
) {
	if (lam_a == 0)
	{
		// setup inner integrals
		PhiFunction phi1(integrate_inf, 0, psi_a, sturm1), phi2(integrate_inf, 1, psi_a, sturm1);
		RadialFunction<double> const *phi1r = &phi1, *phi2r = &phi2;
		double xi = XiIntegral(psi_a, sturm1);
		
		// setup outer integrand
		auto integrand = [ & ] (double x) -> double {
			return chi_a(x) * ((*phi1r)(x) - (*phi2r)(x) - xi * U_a(x)) * sturm2(x);
		};
		
		// integrate
		Integrator<decltype(integrand)> Q(integrand);
		Q.integrate(0., std::numeric_limits<double>::infinity());
		
		/// DEBUG
		/*write_1D_data (
			1000,
			"integrand.out",
			[ & ] (size_t i) -> double { return integrand(i * 0.01); }
		);
		write_1D_data (
			1000,
			"sturm2.out",
			[ & ] (size_t i) -> double { return sturm2(i * 0.01); }
		);
		SymbolicPoly phi1_s = (*(SymbolicPoly const *)(&phi1));
		SymbolicPoly phi2_s = (*(SymbolicPoly const *)(&phi2));
		SymbolicPoly phi12  = phi1_s - phi2_s;
		SymbolicPoly sturm1_s = (*(SymbolicPoly const *)(&sturm1));
		SymbolicPoly sturm2_s = (*(SymbolicPoly const *)(&sturm2));
		SymbolicPoly psia_s = (*(SymbolicPoly const *)(&psi_a));
		printf("phi1: "); write(phi1_s);
		printf("phi2: "); write(phi2_s);
		printf("phi1-phi2: "); write(phi12);
		printf("psi_a: "); write(psia_s);
		printf("sturm1: "); write(sturm1_s);
		printf("sturm2: "); write(sturm2_s);
		printf("psi_a*sturm1: "); write(psi_a*sturm1_s);
		chi_a.toFile("chii.out");
		printf("xi = %g\n", xi);
		printf("res = %g\n", Q.result());
		abort();*/
		
		return Q.result();
	}
	else
	{
		// setup inner integrals
		PhiFunction phi1(integrate_inf, -lam_a, psi_a, sturm1), phi2(integrate_low, lam_a+1, psi_a, sturm1);
		RadialFunction<double> const *phi1r = &phi1, *phi2r = &phi2;
		
		// setup outer integrand
		auto integrand = [ & ] (double x) -> double {
			return chi_a(x) * ((*phi1r)(x) + (*phi2r)(x)) * sturm2(x);
		};
		
		// integrate
		Integrator<decltype(integrand)> Q(integrand);
		Q.integrate(0., std::numeric_limits<double>::infinity());
		return Q.result();
	}
}

double DWBA2::XiIntegral(HydrogenFunction const & psi, SturmianFunction const & S)
{
	SymbolicPoly psi_r = (*(SymbolicPoly const *)(&psi));
	SymbolicPoly S_r   = (*(SymbolicPoly const *)(&S));
	SymbolicTerm r;
	r.ki = 1;
	r.kr = 1;
	r.a = 1;
	
	return eval(integrate_full(r * psi_r * S_r), 0.);
}

double DWBA2::computeA (double E, int Nf1, int Nf2, int Ni1, int Ni2, int L1, int L2)
{
	// overlap integral of r1-dependent states
	double Sigma1 = 0;
	if (Nf1 == Ni1)
		Sigma1 = Nf1;
	if (Nf1 == Ni1 + 1)
		Sigma1 = -0.5 * ALPHA_PLUS(Ni1,L1);
	if (Nf1 == Ni1 - 1)
		Sigma1 = -0.5 * ALPHA_MINUS(Ni1,L1);
	
	// overlap integral of r2-dependent states
	double Sigma2 = 0;
	if (Nf2 == Ni2)
		Sigma2 = Nf2;
	if (Nf2 == Ni2 + 1)
		Sigma2 = -0.5 * ALPHA_PLUS(Ni2,L2);
	if (Nf2 == Ni2 - 1)
		Sigma2 = -0.5 * ALPHA_MINUS(Ni2,L2);
	
	// atomic hamiltonian matrix element
	double Hat1 = (Ni1 - 1) * DELTA(Nf1,Ni1) - 0.5 * Sigma1;
	double Hfree2 = Ni2 * DELTA(Nf2,Ni2) - 0.5 * Sigma2;	// FIXME rozdělení na U2 a U12
	
	// return the result
	return 0.5 * E * Sigma1 * Sigma2 - Hat1 * Sigma2 - Sigma1 * Hfree2 /* -U12 */;
}

void DWBA2::DWBA2_Ln (
	double Ei, int li, int lf, double ki, double kf, 
	int Ni, int Nf, int Li, int Lf, int L1, int L2,
	DistortingPotential const & Ui,
	DistortingPotential const & Uf,
	HydrogenFunction    const & psii,
	HydrogenFunction    const & psif,
	DistortedWave       const & chii,
	DistortedWave       const & chif,
	cArray & DD_lf_li_Ln,
	cArray & DE_lf_li_Ln,
	cArray & ED_lf_li_Ln,
	cArray & EE_lf_li_Ln,
	bool & compute_DD,
	bool & compute_DE,
	bool & compute_ED,
	bool & compute_EE
) {
	// multipole bounds
	//   λi = lami_min ... lami_max
	//   λf = lamf_min ... lamf_max
	
	int lami_min = INT_MAX;
	int lami_max = INT_MIN;
	int lamf_min = INT_MAX;
	int lamf_max = INT_MIN;
	
	int lami_min_D = std::max(abs(L2-li),abs(L1-Li));
	int lami_max_D = std::min(    L2+li ,    L1+Li );
	
	int lamf_min_D = std::max(abs(L2-lf),abs(L1-Lf));
	int lamf_max_D = std::min(    L2+lf ,    L1+Lf );
	
	int lami_min_E = std::max(abs(L2-Li),abs(L1-li));
	int lami_max_E = std::min(    L2+Li ,    L1+li );
	
	int lamf_min_E = std::max(abs(L2-Lf),abs(L1-lf));
	int lamf_max_E = std::min(    L2+Lf ,    L1+lf );
	
	if (compute_DD)
	{
		lami_min = std::min(lami_min, lami_min_D);
		lami_max = std::max(lami_max, lami_max_D);
		lamf_min = std::min(lamf_min, lamf_min_D);
		lamf_max = std::max(lamf_max, lamf_max_D);
	}
	
	if (compute_DE)
	{
		lami_min = std::min(lami_min, lami_min_E);
		lami_max = std::max(lami_max, lami_max_E);
		lamf_min = std::min(lamf_min, lamf_min_D);
		lamf_max = std::max(lamf_max, lamf_max_D);
	}
	
	if (compute_ED)
	{
		lami_min = std::min(lami_min, lami_min_D);
		lami_max = std::max(lami_max, lami_max_D);
		lamf_min = std::min(lamf_min, lamf_min_E);
		lamf_max = std::max(lamf_max, lamf_max_E);
	}
	
	if (compute_EE)
	{
		lami_min = std::min(lami_min, lami_min_E);
		lami_max = std::max(lami_max, lami_max_E);
		lamf_min = std::min(lamf_min, lamf_min_E);
		lamf_max = std::max(lamf_max, lamf_max_E);
	}
	
	// multipole loop
	for (int lami = lami_min; lami <= lami_max; lami++)
	for (int lamf = lamf_min; lamf <= lamf_max; lamf++)
	{
		// radial integral matrices
		CooMatrix Ii, If;
		
		// contributions to the scattering amplitude
		// FIXME: "double" would suffice, but it means to rewrite sparse matrices
		//        to templates
		Complex DD_N_prev = 0, DE_N_prev = 0, ED_N_prev = 0, EE_N_prev = 0;
		
		// Sturmian basis convergence loop
		for (int N = 1; ; N++)
		{
			// N1i = L1+1, ..., L1+N
			// N2i = L2+1, ..., L2+N
			
			/// DEBUG
			char sturm_file[20];
			sprintf(sturm_file, "sturm-%.2d.dat", N);
			SturmianFunction s1(N,0);
			write_1D_data (
				1000,
				sturm_file,
				[&](size_t i) -> double {
					return s1(i*0.05);
				}
			);
			
			// resize matrices
			CooMatrix A(N*N,N*N);
			Ii.resize(N,N);
			If.resize(N,N);
			
			/// DEBUG
			PhiFunction phi1(integrate_inf, 0, psii, SturmianFunction(N,L1));
			PhiFunction phi2(integrate_inf, 1, psii, SturmianFunction(N,L1));
			RadialFunction<double> const *phi1r = &phi1, *phi2r = &phi2;
			char phi1name[20], phi2name[20];
			sprintf(phi1name, "phi1-%.2d.dat", N);
			sprintf(phi2name, "phi2-%.2d.dat", N);
			write_1D_data (
				10000,
				phi1name,
				[&](size_t i) -> double { return (*phi1r)(0.01 * i); }
			);
			write_1D_data (
				10000,
				phi2name,
				[&](size_t i) -> double { return (*phi2r)(0.01 * i); }
			);
			
			// update initial radial integral matrix
			# pragma omp parallel for collapse(2) schedule(dynamic,1)
			for (int Ni1 = L1 + 1; Ni1 <= L1 + N; Ni1++)
			for (int Ni2 = L2 + 1; Ni2 <= L2 + N; Ni2++)
			{
				// if already done, skip
				if ( Ii(Ni1-L1-1,Ni2-L2-1) != 0. )
					continue;
				
				// otherwise compute the value
				double I = computeI (lami, psii, chii, Ui, SturmianFunction(Ni1,L1), SturmianFunction(Ni2,L2));
				
				# pragma omp critical
				Ii.add(Ni1-L1-1, Ni2-L2-1, I);
			}
			
			// update final radial integral matrix
			# pragma omp parallel for collapse(2) schedule(dynamic,1)
			for (int Nf1 = L1 + 1; Nf1 <= L1 + N; Nf1++)
			for (int Nf2 = L2 + 1; Nf2 <= L2 + N; Nf2++)
			{
				// if already done, skip
				if ( If(Nf1-L1-1,Nf2-L2-1) != 0. )
					continue;
				
				// otherwise compute the value
				double I = computeI (lamf, psif, chif, Uf, SturmianFunction(Nf1,L1), SturmianFunction(Nf2,L2));
				
				# pragma omp critical
				If.add(Nf1-L1-1, Nf2-L2-1, I);
			}
			
			// update matrix of the inverse Green operator
			# pragma omp parallel for collapse(4) schedule(dynamic,1)
			for (int Ni1 = L1 + 1; Ni1 <= L1 + N; Ni1++)
			for (int Nf1 = L1 + 1; Nf1 <= L1 + N; Nf1++)
			for (int Ni2 = L2 + 1; Ni2 <= L2 + N; Ni2++)
			for (int Nf2 = L2 + 1; Nf2 <= L2 + N; Nf2++)
			{
				// otherwise compute the value
				double A_val = computeA(ki*ki - 1./(Ni*Ni), Nf1, Nf2, Ni1, Ni2, L1, L2);
				
				// store, if nonzero
				# pragma omp critical
				if (A_val != 0.)
					A.add ((Nf1-L1-1) * N + (Nf2-L2-1), (Ni1-L1-1) * N + (Ni2-L2-1), A_val);
			}
			
			// prepare integral matrices as 1D arrays
			cArray Ii_dense  = Ii.todense();
			cArray Iit_dense = Ii.transpose().todense();
			cArray If_dense  = If.todense();
			cArray Ift_dense = If.transpose().todense();
			
			// factorize the Green function inverse matrix A
			CsrMatrix A_csr = A.tocsr();	// necessary anchor for LUft !!!
			CsrMatrix::LUft A_LU = A_csr.factorize();
			
// 			char A_dat_name[20], A_mat_name[20], A_png_name[20];
// 			sprintf(A_dat_name, "A-%.2d.dat", N);
// 			sprintf(A_mat_name, "A-%.2d.mat", N);
// 			sprintf(A_png_name, "A-%.2d.png", N);
// 			A.write(A_dat_name);
// 			A_csr.plot(A_png_name);
// 			write_2D_data (
// 				A.rows(), A.cols(),
// 				A_mat_name,
// 				[&](size_t i, size_t j) -> double {
// 					return A(i,j).real();
// 				}
//  			);
// 			char I_dat_name[20], I_mat_name[20];
// 			sprintf(I_dat_name, "I-%.2d.dat", N);
// 			sprintf(I_mat_name, "I-%.2d.mat", N);
// 			Ii.write(I_dat_name);
// 			write_2D_data (
// 				Ii.rows(), Ii.cols(),
// 				I_mat_name,
// 				[&](size_t i, size_t j) -> double {
// 					return Ii(i,j).real();
// 				}
// 			);
			
			// convergence indicator
			bool converged = true;
			
			// solve, project and test convergence
			printf("%d\t", N); fflush(stdout);
			if (compute_DD)
			{
				cArray sol = A_LU.solve(Ii_dense);
				Complex DD_N = (If_dense  | sol);
				
				if (DD_N == 0. or abs((DD_N - DD_N_prev) / DD_N ) > EPS_CONTRIB)
				{
					converged = false;
					printf("%g\t", abs(DD_N)); fflush(stdout);
				}
				DD_N_prev = DD_N;
			}
			if (compute_DE)
			{
				Complex DE_N = (If_dense  | A_LU.solve(Iit_dense));
				if (DE_N == 0. or abs((DE_N - DE_N_prev) / DE_N) > EPS_CONTRIB)
				{
					converged = false;
					printf("%g\t", abs(DE_N)); fflush(stdout);
				}
				DE_N_prev = DE_N;
			}
			if (compute_ED)
			{
				Complex ED_N = (Ift_dense | A_LU.solve(Ii_dense));
				if (ED_N == 0. or abs((ED_N - ED_N_prev) / ED_N) > EPS_CONTRIB)
				{
					converged = false;
					printf("%g\t", abs(ED_N)); fflush(stdout);
				}
				ED_N_prev = ED_N;
			}
			if (compute_EE)
			{
				Complex EE_N = (Ift_dense | A_LU.solve(Iit_dense));
				if (EE_N == 0. or abs((EE_N - EE_N_prev) / EE_N) > EPS_CONTRIB)
				{
					converged = false;
					printf("%g\t", abs(EE_N)); fflush(stdout);
				}
				EE_N_prev = EE_N;
			}
			printf("\n"); fflush(stdout);
			
			// if converged, exit the Sturmian convergence loop
			if (converged)
				break;
			
		} // end For N
		
		// add all lambda- and M- dependent angular factors
		for (int Mi = -Li; Mi <= Li; Mi++)
		for (int Mf = -Lf; Mf <= Lf; Mf++)
		for (int mui = -lami; mui <= lami; mui++)
		for (int muf = -lamf; muf <= lamf; muf++)
		{
			int mindex = (Mi + Li) * (2*Lf + 1) + Mf + Lf;
			
			{
				int mi = 0;
				int M1 = Mi - mui;
				int M2 = mui + mi;
				int mf = M2 - muf;
				DD_lf_li_Ln[mindex] += 
					Gaunt(lami,mui,li,mi,L2,M2) * Gaunt(L1,M1,lami,mui,Li,Mi) * // d
					Gaunt(lamf,muf,lf,mf,L2,M2) * Gaunt(L1,M1,lamf,muf,Lf,Mf) * // d
					DD_N_prev;
			}
			
			{
				int mi = 0;
				int M1 = mi - mui;
				int M2 = mui + Mi;
				int mf = M2 - mui;
				DE_lf_li_Ln[mindex] += 
					Gaunt(lami,mui,Li,Mi,L2,M2) * Gaunt(L1,M1,lami,mui,li,mi) * // e
					Gaunt(lamf,muf,lf,mf,L2,M2) * Gaunt(L1,M1,lamf,muf,Lf,Mf) * // d
					DE_N_prev;
			}
			
			{
				int mi = 0;
				int M1 = Mi - mui;
				int M2 = mui + mi;
				int mf = M1 + muf;
				ED_lf_li_Ln[mindex] += 
					Gaunt(lami,mui,li,mi,L2,M2) * Gaunt(L1,M1,lami,mui,Li,Mi) * // d
					Gaunt(lamf,muf,Lf,Mf,L2,M2) * Gaunt(L1,M1,lamf,muf,lf,mf) * // e
					ED_N_prev;
			}
			
			{
				int mi = 0;
				int M1 = mi - mui;
				int M2 = mui + Mi;
				int mf = M1 + muf;
				EE_lf_li_Ln[mindex] += 
					Gaunt(lami,mui,Li,Mi,L2,M2) * Gaunt(L1,M1,lami,mui,li,mi) * // e
					Gaunt(lamf,muf,Lf,Mf,L2,M2) * Gaunt(L1,M1,lamf,muf,lf,mf) * // e
					EE_N_prev;
			}
		}
		
		// slash by lambdas
		double lam_factor = (2*lami + 1) * (2*lamf + 1);
		DD_lf_li_Ln /= lam_factor;
		DE_lf_li_Ln /= lam_factor;
		ED_lf_li_Ln /= lam_factor;
		EE_lf_li_Ln /= lam_factor;
		
	} // end For lambdas
	
	// add the remaining constant complex factors
	Complex factor = pow(4*M_PI, 4) * pow(Complex(0.,1.), li-lf) / (ki * kf);
	DD_lf_li_Ln *= factor;
	DE_lf_li_Ln *= factor;
	ED_lf_li_Ln *= factor;
	EE_lf_li_Ln *= factor;
}
