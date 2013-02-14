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
	//
	
	double min_Kn = 0;						// just after ionization
	double max_Kn = sqrt(Ei - 1./(Ni*Ni));	// all energy of the projectile
	
	rArray eKn;				// evaluation energies
	cArray eDD;				// evaluated values
	Complex pDD = 0;		// previous values
	
	for (int points = 4; ; points *= 2)
	{
		std::cout << "points = " << points << std::endl;
		
		// evaluate T-matrices at Chebyshev roots
		for (int point = 0; point < points; point++)
		{
			// get evaluation point
			double Kn = Chebyshev<double,Complex>::root(points,point,min_Kn,max_Kn);
			
			// did we already evaluate the T-matrices for this point?
			if (std::find(eKn.begin(),eKn.end(),Kn) == eKn.end())
			{
				// no, add element to the arrays
				eKn.push_back(Kn);
				eDD.push_back(Complex(0));
				
				// get intermediate state
				HydrogenFunction psin(Kn,Ln);
				
				// and evaluate the T-matrices
				DWBA2::DWBA2_En (
					Ei, Ni, Kn*Kn, Ln, lami, lamf,
					psii, psif, psin,
					Ui, Uf, Ug,
					chii, chif,
					eDD.back()
				);
			}
		}
		
		// construct Chebyshev approximation and integral of the DD integrand
		auto DDre = [&](double x) -> double { return std::find(eDD.begin(),eDD.end(),x)->real(); };
		auto DDim = [&](double x) -> double { return std::find(eDD.begin(),eDD.end(),x)->imag(); };
		Chebyshev<double,double> cb_DDre(DDre, points, min_Kn, max_Kn);
		Chebyshev<double,double> cb_DDim(DDim, points, min_Kn, max_Kn);
		Chebyshev<double,double> cb_DDre_int = cb_DDre.integrate();
		Chebyshev<double,double> cb_DDim_int = cb_DDim.integrate();
		Complex DD_val (
			cb_DDre_int.clenshaw(max_Kn,1e-8) - cb_DDre_int.clenshaw(min_Kn,1e-8),
			cb_DDim_int.clenshaw(max_Kn,1e-8) - cb_DDim_int.clenshaw(min_Kn,1e-8)
		);
		
		std::cout << "DD_val = " << DD_val << "\n";
		
		// check convergence
		if (abs(DD_val - pDD) < EPS * abs(DD_val))
		{
			// we converged! update sums and return
			cDD += DD_val;
			return;
		}
		else
		{
			// move current values to "previous" values and loop
			pDD = DD_val;
		}
	}
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
	IrregularWave geta = Ug.getIrregularWave(Kn,Ln);
	
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
		
		int n1 = -1;	// disable bisection
		Complex Q1 = ClenshawCurtis<decltype(integrand1),Complex>(integrand1, 0., std::min(r1,fari), false, 1.0, 1e-7, &n1);
		
		int n2 = -1;	// disable bisection
		Complex Q2 = ClenshawCurtis<decltype(integrand2),Complex>(integrand2, std::min(r1,fari), fari, false, 1.0, 1e-7, &n2);
		
		return phif(r1) * (geta(r1) * Q1 + gphi(r1) * Q2);
	};
	
	// integrate
	int n = -1;	// disable bisection
	DD = ClenshawCurtis<decltype(integrand), Complex>(integrand, 0., farf, false, 1.0, 1e-5, &n);
}
