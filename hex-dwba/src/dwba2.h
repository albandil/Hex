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

#ifndef HEX_DWBA_DWBA2
#define HEX_DWBA_DWBA2

#include "arrays.h"
#include "chebyshev.h"
#include "potential.h"
#include "specf.h"

/**
 * Main DWBA2-part function that will compute a scattering amplitude
 * contribution for given angular numbers.
 */
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
);
	
void DWBA2_energy_driver (
	double Ei, int li, int lf, double ki, double kf, 
	int Ni, int Nf, int Li, int Lf, int Ln,
	int lami, int lamf, 
	HydrogenFunction const & psii, 
	HydrogenFunction const & psif, 
	DistortingPotential const & Ui,
	DistortingPotential const & Uf,
	DistortedWave const & chii,
	DistortedWave const & chif, Complex & cDD
);
	
void DWBA2_En (
	double Ei, int Ni,
	double Kn, int Ln,
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
);

template <
	class Chif, class Phif, class Psif,
	class GPhi, class GEta,
	class Psii, class Phii, class Chii
> Complex GreensFunctionIntegral (
	Chif const & chif, Phif const & phif, Psif const & psif,
	GPhi const & gphi, GEta const & geta,
	Psii const & psii, Phii const & phii, Chii const & chii,
	bool scaling
) {
	std::cout << "\tGreen's integral\n";
	
	// get integration upper bound
	double fari = 5. * psii.far(1e-8);
	double farf = 5. * psif.far(1e-8);
	
	// Green's function integrand
	auto integrand = [&](double r1) -> Complex {
		
		// inner integrand variations
		
		auto integrand1 = [&](double r2) -> Complex {
			if (scaling)
				return gphi(r2) * phii(r2) * chii(r2) * exp(r2 - r1);
			else
				return gphi(r2) * phii(r2) * chii(r2);
		};
		auto integrand2 = [&](double r2) -> Complex {
			if (not finite(r2))
				return 0.;
			if (scaling)
				return geta(r2) * phii(r2) * chii(r2) * exp(r1 - r2);
			else
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
	
	// integrate
	return Q.integrate(0., farf);
}

#endif

