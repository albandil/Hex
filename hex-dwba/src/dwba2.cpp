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
#include <cstdio>
#include <cmath>

#include "angs.h"
#include "complex.h"
#include "dwba2.h"
#include "hydrogen.h"
#include "integrate.h"
#include "potential.h"
#include "wave_distort.h"
#include "wave_irreg.h"
#include "wave_forbid.h"
#include "wave_hyperb.h"
#include "spmatrix.h"

#define EPS_CONTRIB 1e-8

DWBA2::PhiFunctionDirIntegral::PhiFunctionDirIntegral (
	HydrogenFunction const & psin, 
	int lam, 
	DistortingPotential const & U, 
	HydrogenFunction const & psi
) : Lam(lam), U(U), Psi(psi), Psin(psin) {}

double DWBA2::PhiFunctionDirIntegral::operator()(double x2) const
{
	if (not finite(x2))
		return 0.;
	
	//
	// finite integrand
	//
		
	auto integrand1 = [&](double x1) -> double {
		if (x2 == 0.)
			return 0.;
		return Psin(x1) * pow(x1/x2, Lam) * Psi(x1);
	};
	
	// compactification
	CompactIntegrand<decltype(integrand1),double> R1(integrand1, 0., Inf, 1.0);
	
	// integration system
	ClenshawCurtis<decltype(R1),double> Q1(R1);
	Q1.setEps(1e-10);	// be precise
	
	// integrate
	double i1 = (x2 == 0) ? 0. : Q1.integrate(R1.scale(0.), R1.scale(x2)) / x2;
	
	//
	// infinite integrand:
	//
	
	auto integrand2 = [&](double x1) -> double {
		if (x1 == 0.)
			return 0.;
		return Psin(x1) * pow(x2/x1, Lam) * Psi(x1) / x1;
	};
	
	// compactification
	CompactIntegrand<decltype(integrand2),double> R2(integrand2, 0., Inf, 1.0);
	
	// integration system
	ClenshawCurtis<decltype(R2),double> Q2(R2);
	Q2.setEps(1e-10);
	
	// integrate
	double i2 = Q2.integrate(R2.scale(x2), R2.scale(Inf));
	
	// sum the two integrals
	return i1 + i2;
};

DWBA2::PhiFunctionDir::PhiFunctionDir (
	HydrogenFunction const & psin, 
	int lam, 
	DistortingPotential const & U, 
	HydrogenFunction const & psi
) : Cb_inf(0.), Lam(lam), U(U), Diag(psin == psi), 
	Zero(Diag and (DistortingPotential(psi) == U)),
	Integrand(psin,lam,U,psi),
	CompactIntegral(Integrand,0.,false,1.0)
{
	// the case λ = 0 is easy:
	if (lam == 0)
	{
		if (Zero)
			return;
		
		if (Diag)
		{
			// integrand
			auto integrand = [&](double x1) -> double {
				return (x1 == 0.) ? 0. : gsl_sf_pow_int(psi(x1), 2) / x1;
			};
			
			// integration system
			ClenshawCurtis<decltype(integrand),double> Q(integrand);
			
			// integrate
			Cb_inf = Q.integrate(0., Inf);
			
			return;
		}
		
		/* otherwise continue */
	}
	
	// convergence loop
	for (int N = 16; ; N *= 2)
	{
		// construct a Chebyshev approximation of the compactified function
		CompactIntegralCb = Chebyshev<double,double> (CompactIntegral, N, -1., 1.);

		// check convergence
		if (CompactIntegralCb.tail(1e-10) != N)
			break;
	}
	
	// get optimal truncation index
	Tail = CompactIntegralCb.tail(1e-10);
	
	// evaluate Phi at positive infinity
	Cb_inf = CompactIntegralCb.clenshaw(1., Tail);
}

double DWBA2::PhiFunctionDir::operator() (double x) const
{
	if (Lam == 0)
	{
		if (Zero) return 0.;
		if (Diag) return Cb_inf - U.plusMonopole(x);
	}
	
	return Cb_inf - CompactIntegralCb.clenshaw(CompactIntegral.scale(x), Tail);
}

void DWBA2::DWBA2_Ln (
	double Ei, int li, int lf, double ki, double kf, 
	int Ni, int Nf, int Li, int Lf, int Ln,
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
	
	int lami_min_D = std::max(abs(Ln-li),abs(Ln-Li));
	int lami_max_D = std::min(    Ln+li ,    Ln+Li );
	
	int lamf_min_D = std::max(abs(Ln-lf),abs(Ln-Lf));
	int lamf_max_D = std::min(    Ln+lf ,    Ln+Lf );
	
	lami_min = std::min(lami_min, lami_min_D);
	lami_max = std::max(lami_max, lami_max_D);
	lamf_min = std::min(lamf_min, lamf_min_D);
	lamf_max = std::max(lamf_max, lamf_max_D);
	
	// multipole loop
	for (int lami = lami_min; lami <= lami_max; lami++)
	for (int lamf = lamf_min; lamf <= lamf_max; lamf++)
	{
		// contributions to the scattering amplitude
		// FIXME: "double" would suffice, but it means to rewrite sparse matrices
		//        to templates
		Complex DD_N_prev = 0;
		
		// compute energy sum/integral
		DWBA2::DWBA2_energy_driver (
			Ei, li, lf, ki, kf, 
			Ni, Nf, Li, Lf, Ln,
			lami, lamf, 
			psii, psif, Ui, Uf, chii, chif, DD_N_prev
		);
		
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
					Gaunt(lami,mui,li,mi,Ln,M2) * Gaunt(Ln,M1,lami,mui,Li,Mi) * // d
					Gaunt(lamf,muf,lf,mf,Ln,M2) * Gaunt(Ln,M1,lamf,muf,Lf,Mf) * // d
					DD_N_prev;
			}
		}
		
		// slash by lambdas
		double lam_factor = (2*lami + 1) * (2*lamf + 1);
		DD_lf_li_Ln /= lam_factor;
		
	} // end For lambdas
	
	// add the remaining constant complex factors
	Complex factor = pow(4*M_PI, 4) * pow(Complex(0.,1.), li-lf) / (ki * kf);
	DD_lf_li_Ln *= factor;
}

void DWBA2::DWBA2_energy_driver (
	double Ei, int li, int lf, double ki, double kf, 
	int Ni, int Nf, int Li, int Lf, int Ln,
	int lami, int lamf, 
	HydrogenFunction const & psii, 
	HydrogenFunction const & psif, 
	DistortingPotential const & Ui,
	DistortingPotential const & Uf,
	DistortedWave const & chii,
	DistortedWave const & chif, Complex & cDD
) {
	
	static const double EPS = 1e-3;
	
	// set Green's function distorting potential
	DistortingPotential Ug(1);	// = U(1s)
	
	// sum over discrete intermadiate states
	for (int Nn = 1; ; Nn++)
	{
		// info
		std::cout << "[DWBA2_energy_part] Nn = " << Nn << std::endl;
		
		// get intermediate state
		HydrogenFunction psin(Nn,Ln);
		
		// compute contribution from this discrete intermediate state
		Complex cDD_Nn;
		DWBA2::DWBA2_En (
			Ei, Ni, -1./(Nn*Nn), Ln, lami, lamf,
			psii, psif, psin,
			Ui, Uf, Ug,
			chii, chif,
			cDD_Nn
		);
		
		// update sums
		cDD += cDD_Nn;
		
		std::cout << "[DWBA2_energy_part] cDD_Nn = " << cDD_Nn << ", cDD = " << cDD << " [" << std::abs(cDD_Nn)/std::abs(cDD) << "]" << std::endl;
		
		// check convergence
		if (std::abs(cDD_Nn) < EPS * std::abs(cDD))
			break;
	}
	
	// integrate over intermediate continuum using Clenshaw-Curtis quadrature
	//  - only integrate over allowed energies
	//  - disable recurrence for maximal reuse of evaluations
	
	double min_Kn = 0;						// just after ionization
	double max_Kn = sqrt(Ei - 1./(Ni*Ni));	// all energy of the projectile
	
	// integrand
	auto DD_integrand = [&](double Kn) -> Complex {
		
		// get intermediate state
		HydrogenFunction psin(Kn,Ln);
		
		// compute amplitude
		Complex dd;
		DWBA2::DWBA2_En (
			Ei, Ni, Kn*Kn, Ln, lami, lamf,
			psii, psif, psin,
			Ui, Uf, Ug,
			chii, chif,
			dd
		);
		
		return dd;
		
	};
	
	// integration system
	ClenshawCurtis<decltype(DD_integrand),Complex> QDD(DD_integrand);
	QDD.setLim(false);
	QDD.setRec(false);
	QDD.setEps(1e-5);
	cDD += QDD.integrate(min_Kn, max_Kn);
}

void DWBA2::DWBA2_En (
	double Ei, int Ni,
	double Eatn, int Ln,
	int lami, int lamf,
	HydrogenFunction const & psii,
	HydrogenFunction const & psif,
	HydrogenFunction const & psin,
	DistortingPotential const & Ui,
	DistortingPotential const & Uf,
	DistortingPotential const & Ug,
	DistortedWave const & chii,
	DistortedWave const & chif,
	Complex & DD
) {
	// get projectile energy
	double Kn = sqrt(Ei - 1./(Ni*Ni) - Eatn);
	
	std::cout << "[DWBA2_En] Kn = " << Kn << std::endl;
	
	// construct Green's function parts
	DistortedWave gphi = Ug.getDistortedWave(Kn,Ln);
	std::cout << "gphi OK\n";
	IrregularWave geta = Ug.getIrregularWave(Kn,Ln);
	std::cout << "geta OK\n";
	
	
	// construct direct inner integrals
	PhiFunctionDir phii(psin, lami, Ui, psii);
	PhiFunctionDir phif(psin, lamf, Uf, psif);
	
	// get integration upper bound
	double fari = 5. * psii.far(1e-8);
	double farf = 5. * psif.far(1e-8);
	
	// Green's function integrand
	auto integrand = [&](double r1) -> Complex {
		
		// inner integrand variations
		
		auto integrand1 = [&](double r2) -> Complex {
			return gphi(r2) * phii(r2) * chii(r2);
		};
		auto integrand2 = [&](double r2) -> Complex {
			if (not finite(r2))
				return 0.;
			return geta(r2) * phii(r2) * chii(r2);
		};
		
		// integration systems
		
		ClenshawCurtis<decltype(integrand1),Complex> Q1(integrand1);
		Q1.setLim(false);
		Q1.setRec(false);
		Q1.setEps(1e-7);
		
		ClenshawCurtis<decltype(integrand2),Complex> Q2(integrand2);
		Q2.setLim(false);
		Q2.setRec(false);
		Q2.setEps(1e-7);
		
		// integrate
		
		Complex q1 = Q1.integrate(0., std::min(r1,fari));
		Complex q2 = Q2.integrate(std::min(r1,fari), fari);
		
		// evaulate outer integrand
		
		return phif(r1) * (geta(r1) * q1 + gphi(r1) * q2);
	};
	
	// outer integration system
	ClenshawCurtis<decltype(integrand),Complex> Q(integrand);
	Q.setLim(false);
	Q.setRec(false);
	Q.setEps(1e-5);
	
	std::cout << "[DWBA2_En] Green's integral\n";
	
	// integrate
	DD = Q.integrate(0., farf);
}
