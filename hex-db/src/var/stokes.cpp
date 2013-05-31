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

const std::string StokesParameters::Id = "stokes";
const std::string StokesParameters::Description = "Reduced Stokes parameters.";
const std::vector<std::string> StokesParameters::Dependencies = {
	"ni", /* li = 0, mi = 0 */
	"nf", /* lf = 1, |mf| ≤ 1 */
	"Ei", "theta"
};

std::string const & StokesParameters::SQL_CreateTable() const
{
	static const std::string cmd = "";
	return cmd;
}

std::string const & StokesParameters::SQL_Update() const
{
	static const std::string cmd = "";
	return cmd;
}

bool StokesParameters::run (
	eUnit Eunits, lUnit Lunits,
	sqlitepp::session & db,
	std::map<std::string,std::string> const & sdata,
	rArray const & angles
) const {
	
	// manage units
	double efactor = change_units(Eunits, eUnit_Ry);
	double lfactor = change_units(lUnit_au, Lunits);
	
	// atomic and projectile data
	int ni = As<int>(sdata, "ni", Id);
	int nf = As<int>(sdata, "nf", Id);
	int Ei = As<double>(sdata, "Ei", Id) * efactor;
	
	// Should return the following:
	//    P₁ = 2λ - 1
	//    P₂ = -2√2 R
	//    P₃ = 2√2 I
	// where
	//    λ = <|f₀²|>
	//    R = Re { <f₀* f₁> } / [ dσ/dΩ ]
	//    I = Im { <f₀* f₁> } / [ dσ/dΩ ]
	// where
	//    f₀ ... is the scattering amplitude to mf = 0
	//    f₁ ... is the scattering amplitude to mf = 1
	//    dσ/dΩ ... is the DCS summed over mf
	//    <.> ... stands for averaging over spin states, i.e.
	//        <a> = [ a(S=0) + 3a(S=1) ] / 4
	
	// TODO
	
	// energy and real and imarinary part of the T-matrix
// 	double E, Re_T_ell, Im_T_ell;
	
	// create query statement
// 	sqlitepp::statement st(db);
// 	st << "SELECT Ei, Re_T_ell, Im_T_ell FROM " + TMatrix::Id + " "
// 	      "WHERE ni = :ni "
// 		  "  AND li = :li "
// 		  "  AND mi = :mi "
// 		  "  AND nf = :nf "
// 		  "  AND lf = :lf "
// 		  "  AND mf = :mf "
// 		  "  AND  L = :L  "
// 		  "  AND  S = :S  "
// 		  "  AND ell=:ell "
// 		  "ORDER BY Ei ASC",
// 	   sqlitepp::into(E), sqlitepp::into(Re_T_ell), sqlitepp::into(Im_T_ell),
// 	   sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
// 	   sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
// 	   sqlitepp::use(L), sqlitepp::use(S), sqlitepp::use(ell);
	
	// get T-matrices
// 	rArray E_arr;
// 	cArray T_arr;
// 	while (st.exec())
// 	{
// 		E_arr.push_back(E);
// 		T_arr.push_back(Complex(Re_T_ell,Im_T_ell));
// 	}
	
	// interpolate
// 	cArray T_out = interpolate(E_arr, T_arr, energies * efactor);
	
	// write out
// 	std::cout << this->logo() <<
// 		"# T-matrices in " << unit_name(Lunits) << " for\n" <<
// 		"#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
// 	    "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n" <<
// 	    "#     L = " << L << ", S = " << S << ", ℓ = " << ell << "\n" <<
// 	    "# ordered by energy in " << unit_name(Eunits) << "\n" <<
// 	    "# \n" <<
// 	    "# E\t Re T\t Im T\n";
// 	for (size_t i = 0; i < energies.size(); i++)
// 	{
// 		std::cout << 
// 			energies[i] << "\t" << 
// 			T_out[i].real()*lfactor << "\t" <<
// 			T_out[i].imag()*lfactor << "\n";
// 	}
	
// 	return true;
	return false;
}
