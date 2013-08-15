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

#include <map>
#include <string>
#include <vector>

#include "../arrays.h"
#include "../interpolate.h"
#include "../variables.h"

//
// forward declaration of recycled subroutines from other variables
//

extern cArray scattering_amplitude (
	sqlitepp::session & db,
	int ni, int li, int mi,
	int nf, int lf, int mf,
	int S, double E, rArray const & angles
);

extern rArray differential_cross_section (
	sqlitepp::session & db,
	int ni, int li, int mi,
	int nf, int lf, int mf,
	int S, double E, rArray const & angles
);

//
// StokesParameters members
//

const std::string StokesParameters::Id = "stokes";
const std::string StokesParameters::Description = "Reduced Stokes parameters for the ns->n'p transition.";
const std::vector<std::string> StokesParameters::Dependencies = {
	"ni", /* li = 0, mi = 0 */
	"nf", /* lf = 1, |mf| ≤ 1 */
	"Ei", "theta"
};
const std::vector<std::string> StokesParameters::VecDependencies = { "theta" };

bool StokesParameters::initialize(sqlitepp::session & db) const
{
	return true;
}

std::vector<std::string> const & StokesParameters::SQL_CreateTable() const
{
	static const std::vector<std::string> cmd;
	return cmd;
}

std::vector<std::string> const & StokesParameters::SQL_Update() const
{
	static const std::vector<std::string> cmd;
	return cmd;
}

bool StokesParameters::run (
	sqlitepp::session & db,
	std::map<std::string,std::string> const & sdata
) const {
	
	// manage units
	double efactor = change_units(Eunits, eUnit_Ry);
	double afactor = change_units(Aunits, aUnit_rad);
	
	// atomic and projectile data
	int ni = As<int>(sdata, "ni", Id);
	int nf = As<int>(sdata, "nf", Id);
	int Ei = As<double>(sdata, "Ei", Id) * efactor;
	double ki = sqrt(Ei);
	double kf = sqrt(Ei - 1./(ni*ni) + 1./(nf*nf));
	
	// angles
	rArray angles;
	
	// get angle / angles
	try {
		
		// is there a single angle specified using command line ?
		angles.push_back(As<double>(sdata, "theta", Id));
		
	} catch (std::exception e) {
		
		// are there more angles specified using the STDIN ?
		angles = readStandardInput<double>();
	}
	
	// Should return the following:
	//    P₁ = 2λ - 1
	//    P₂ = -2√2 R
	//    P₃ = 2√2 I
	// and
	//    Pl = √(P₁²+P₂²)
	//    γ  = arg (P₁ + iP₂) / 2
	//    P⁺ = √(P₁²+P₂² + P₃²)
	// where
	//    λ = <|f₀²|> / [ dσ/dΩ ]
	//    R = Re { <f₀* f₁> } / [ dσ/dΩ ]
	//    I = Im { <f₀* f₁> } / [ dσ/dΩ ]
	// where
	//    f₀ ... is the scattering amplitude to mf = 0
	//    f₁ ... is the scattering amplitude to mf = 1
	//    dσ/dΩ ... is the DCS summed over mf
	//    <.> ... stands for averaging over spin states, i.e.
	//        <a> = [ a(S=0) + 3a(S=1) ] / 4
	
	// compute scattering amplitudes
	cArray f0_singlet = scattering_amplitude(db, ni,0,0, nf,1,0, 0, Ei, angles * afactor);
	cArray f0_triplet = scattering_amplitude(db, ni,0,0, nf,1,0, 1, Ei, angles * afactor);
	cArray f1_singlet = scattering_amplitude(db, ni,0,0, nf,1,1, 0, Ei, angles * afactor);
	cArray f1_triplet = scattering_amplitude(db, ni,0,0, nf,1,1, 1, Ei, angles * afactor);
	
	// compute differential cross sections
	rArray dcs = differential_cross_section(db, ni,0,0, nf,1, 0, 0, Ei, angles * afactor)
	           + differential_cross_section(db, ni,0,0, nf,1, 1, 0, Ei, angles * afactor)
	           + differential_cross_section(db, ni,0,0, nf,1,-1, 0, Ei, angles * afactor)
	           + differential_cross_section(db, ni,0,0, nf,1, 0, 1, Ei, angles * afactor)
	           + differential_cross_section(db, ni,0,0, nf,1, 1, 1, Ei, angles * afactor)
	           + differential_cross_section(db, ni,0,0, nf,1,-1, 1, Ei, angles * afactor);
	
	// compute basic parameters
	rArray lambda = (0.25 * abs(f0_singlet*f0_singlet) + 0.75 * abs(f0_triplet*f0_triplet)) / dcs * (kf/ki);
	rArray R = realpart(0.25 * f1_singlet*f0_singlet.conj() + 0.75 * f1_triplet*f0_triplet.conj()) / dcs * (kf/ki);
	rArray I = imagpart(0.25 * f1_singlet*f0_singlet.conj() + 0.75 * f1_triplet*f0_triplet.conj()) / dcs * (kf/ki);
	
	// compute derived parameters
	rArray P1 = 2. * lambda - 1.;
	rArray P2 = 2 * M_SQRT2 * R;
	rArray P3 = -2 * M_SQRT2 * I;
	rArray Pl = hypot(P1, P2);
	rArray gamma = 0.5 * atan2(P2, P1);
	rArray Pplus = hypot(Pl, P3);
	
	// write out
	std::cout << this->logo() <<
		"# Stokes parameters for\n"
		"#     ni = " << ni << ", li = " << 0 << ",\n"
	    "#     nf = " << nf << ", lf = " << 1 << ",\n"
	    "#     E = " << Ei/efactor << " " << unit_name(Eunits) << "\n"
	    "# ordered by angle in " << unit_name(Aunits) << "\n"
		"# \n"
	    "# θ\t λ\t R\t I\t P₁\t P₂\t P₃\t Pl\t γ\t P⁺\n";
	for (size_t i = 0; i < angles.size(); i++)
	{
		std::cout << 
			angles[i] << "\t" << 
			lambda[i] << "\t" <<
			R[i] << "\t" <<
			I[i] << "\t" <<
			P1[i] << "\t" <<
			P2[i] << "\t" <<
			P3[i] << "\t" <<
			Pl[i] << "\t" <<
			gamma[i] / afactor << "\t" <<
			Pplus[i] << "\n";
	}
	
	return true;
}
