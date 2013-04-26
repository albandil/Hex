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
#include "multi.h"
#include "potential.h"
#include "wave_distort.h"
#include "wave_irreg.h"
#include "wave_forbid.h"
#include "wave_hyperb.h"

#define EPS_CONTRIB 1e-8

void DWBA2_Ln (
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
		Complex DD_N_prev = 0, DD_N_prev_1 = 0;
		
		// for all Green's function angular momenta
		for (int ln = 0; ; ln++)
		{
			// skip too low ln-s
			if (Ln + ln < std::max(std::abs(li - Li), std::abs(lf - Lf)) or 
				Ln - ln > std::min(         li + Li,           lf + Lf))
				continue;
			
			// exit loop if ln too high
			if (ln - Ln > std::min(Li + li, Lf + lf))
				break;
			
			std::cout << "---------------------------------------\n";
			std::cout << "ln = " << ln << "\n";
			std::cout << "---------------------------------------\n";
			
			// compute energy sum/integral
			DWBA2_energy_driver (
				Ei, li, lf, ki, kf, 
				Ni, Nf, Li, Lf, Ln,
				lami, lamf, ln,
				psii, psif, Ui, Uf, chii, chif, DD_N_prev_1
			);
			
			// update T-matrix contribution
			DD_N_prev += DD_N_prev_1;
		}
		
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

void DWBA2_energy_driver (
	double Ei, int li, int lf, double ki, double kf, 
	int Ni, int Nf, int Li, int Lf, int Ln,
	int lami, int lamf, int ln,
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
	
	// computes contribution from a discrete state
	auto DiscreteContribution = [ & ](int Nn) -> Complex {
		
		// get intermediate state
		HydrogenFunction psin(Nn,Ln);
		
		// compute contribution from this discrete intermediate state
		Complex contrib;
		DWBA2_En (
			Ei, Ni, -1./(Nn*Nn), Ln, lami, lamf, ln,
			psii, psif, psin,
			Ui, Uf, Ug,
			chii, chif,
			contrib
		);
		
		std::cerr << -1./(Nn*Nn) << "\t" << contrib.real() << "\t" << contrib.imag() << "\n";
		
		return contrib;
	};
	
	// computes contribution density from a continuum state
	auto ContinuumContribution = [ & ](double Kn) -> Complex {
		
		if (Kn == 0. or not finite(Kn) or std::abs(Ei - 1./(Ni*Ni) - Kn*Kn) < 1e-8)
		{
			std::cerr << Kn*Kn << "\t" << 0. << "\t" << 0. << "\n";
			return Complex(0.);
		}
		
		// get intermediate state
		HydrogenFunction psin(Kn,Ln);
		
		// compute amplitude
		Complex contrib;
		DWBA2_En (
			Ei, Ni, Kn*Kn, Ln, lami, lamf, ln,
			psii, psif, psin,
			Ui, Uf, Ug,
			chii, chif,
			contrib
		);
		
		std::cerr << Kn*Kn << "\t" << Kn*Kn*contrib.real() << "\t" << Kn*Kn*contrib.imag() << "\n";
		
		return contrib * Kn * Kn;
	};
	

#if 0
	
	// compactification of the previous function
// 	CompactificationR<decltype(ContinuumContribution),Complex>
// 		compact(ContinuumContribution, 0., false, 1.0);
	
// 	for (int Nn = 1; Nn <= 20; Nn++)
// 		DiscreteContribution(Nn);
		
	
// 	for (double Kn = 0; Kn < 5; Kn += 0.01)
// 		compact(compact.scale(Kn));
	
	ContinuumContribution(27.0062);
		
	exit(0);
	
#endif
		
	
	// sum over discrete intermadiate states
	for (int Nn = Ln + 1; ; Nn++)
	{
		// update sums
		Complex contrib = DiscreteContribution(Nn);
		cDD += contrib;
		
		if (std::abs(contrib) < EPS * std::abs(cDD))
			break;
	}
	
	std::cout << "Discrete DD: " << cDD << "\n";
	
	// integrate over intermediate continuum using Clenshaw-Curtis quadrature
	//  - disable recurrence for maximal reuse of evaluations
	
	//
	// allowed regime
	//
	
	ClenshawCurtis<decltype(ContinuumContribution),Complex> QDD(ContinuumContribution);
	QDD.setLim(false);
	QDD.setRec(false);
	QDD.setSubdiv(15);
	QDD.setTol(1e-4);
	QDD.setEps(1e-4);
	QDD.setVerbose(false, "Allowed regime integral");
	
	// integrate over allowed region
	Complex DD_allowed = QDD.integrate (
		0.,
		sqrt(Ei - 1./(Ni*Ni))
	);
	cDD += DD_allowed;
	
	std::cout << "Allowed continuum DD: " << DD_allowed << "\n";
	
	//
	// forbidden regime
	//
	
	CompactIntegrand<decltype(ContinuumContribution),Complex> cQDD (
		ContinuumContribution,
		sqrt(Ei - 1./(Ni*Ni)),
		Inf,
		false,
		1.0
	);
	ClenshawCurtis<decltype(cQDD),Complex> iQDD(cQDD);
	iQDD.setLim(false);
	iQDD.setRec(false);
	iQDD.setSubdiv(15);
	iQDD.setTol(1e-4);
	iQDD.setEps(1e-4);
	iQDD.setVerbose(true, "Forbidden regime integral");
	
	/// DEBUG FIXME
	double momentum_cutoff = Inf;
	
	// integrate over forbidden region
	Complex DD_forbidden = iQDD.integrate (
		cQDD.scale(sqrt(Ei - 1./(Ni*Ni))),
		cQDD.scale(momentum_cutoff)
	);
	cDD += DD_forbidden;
	
	std::cout << "Forbidden continuum DD: " << DD_forbidden << "\n";
	std::cerr << std::endl << std::endl;
}

void DWBA2_En (
	double Ei, int Ni,
	double Eatn, int Ln,
	int lami, int lamf, int ln,
	HydrogenFunction const & psii,
	HydrogenFunction const & psif,
	HydrogenFunction const & psin,
	DistortingPotential const & Ui,
	DistortingPotential const & Uf,
	DistortingPotential const & Ug,
	DistortedWave const & chii,
	DistortedWave const & chif,
	Complex & tmat
) {
	std::cout << "\nDWBA_En\n-------\n\tEatn = " << Eatn << "\n";
	
	// get projectile wave number
	double kn_sqr = Ei - 1./(Ni*Ni) - Eatn;
	double kn = sqrt(std::abs(kn_sqr));
	
	std::cout << "\tpsin.far() = " << psin.getFar() << std::endl;
	
	// construct direct inner integrals
	std::cout << "\tcomputing phii... " << std::flush;
	PhiFunctionDir phii(psin, lami, Ui, psii);
	std::cout << "ok\n\tcomputing phif... " << std::flush;
	PhiFunctionDir phif(psin, lamf, Uf, psif);
	std::cout << "on\n";
	
	// allowed/forbidden regime
	if (kn_sqr > 0)
	{
		std::cout << "\tAllowed regime, kn = " << kn << "\n";
		
		// construct Green's function parts
		std::cout << "\tcomputing gphi... " << std::flush;
		DistortedWave gphi = Ug.getDistortedWave(kn,ln);
		std::cout << "ok\n\tcomputing geta... " << std::flush;
		IrregularWave geta = Ug.getIrregularWave(kn,ln);
		std::cout << "ok\n";
		
		// integrate
		tmat = GreensFunctionIntegral(chif, phif, psif, gphi, geta, psii, phii, chii, false);
	}
	else
	{
		std::cout << "\tForbidden regime, kn = " << kn << "i\n";
		
		// construct Green's function parts
		std::cout << "\tcomputing gtheta... " << std::flush;
		HyperbolicWave gtheta = Ug.getHyperbolicWave(kn,ln);
		std::cout << "ok\n\tcomputing gzeta... " << std::flush;
		ForbiddenWave gzeta = Ug.getForbiddenWave(kn,ln);
		std::cout << "ok\n";
		
		// return only scaled values
		gtheta.scale(true);
		gzeta.scale(true);
		
		// integrate
		tmat = GreensFunctionIntegral(chif, phif, psif, gtheta, gzeta, psii, phii, chii, true);
	}
}
