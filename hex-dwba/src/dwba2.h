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
#include "potential.h"
#include "specf.h"

/**
 * Class members compute contributions to the second order of the distorted
 * wave Born approximation.
 */
class DWBA2
{
public:

    /**
     * Main DWBA2-part function that will compute a scattering amplitude
     * contribution for given angular numbers.
     */
    static void DWBA2_Ln (
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
	
	static void DWBA2_energy_driver (
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
	
	static void DWBA2_En (
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
};

#endif

