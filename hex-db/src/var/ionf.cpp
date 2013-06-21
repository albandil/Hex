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

#include <gsl/gsl_sf.h>

#include "../arrays.h"
#include "../interpolate.h"
#include "../chebyshev.h"
#include "../variables.h"
#include "../vec3d.h"

const std::string IonizationF::Id = "ionf";
const std::string IonizationF::Description = "Ionization amplitude radial part.";
const std::vector<std::string> IonizationF::Dependencies = {
	"ni", "li", "mi", 
	"L", "S",
	"Ei", "l1", "l2",
	"Eshare"
};
const std::vector<std::string> IonizationF::VecDependencies = { "Eshare" };

bool IonizationF::initialize(sqlitepp::session & db) const
{
	return true;
}

std::vector<std::string> const & IonizationF::SQL_CreateTable() const
{
	static const std::vector<std::string> cmd = {
		"CREATE TABLE '" + IonizationF::Id + "' ("
			"ni INTEGER, "
			"li INTEGER, "
			"mi INTEGER, "
			"L  INTEGER, "
			"S  INTEGER, "
			"Ei DOUBLE PRECISION, "
			"l1 INTEGER, "
			"l2 INTEGER, "
			"cheb BLOB, "
			"PRIMARY KEY (ni,li,mi,L,S,Ei,l1,l2)"
		")"
	};
	
	return cmd;
}

std::vector<std::string> const & IonizationF::SQL_Update() const
{
	static const std::vector<std::string> cmd;
	return cmd;
}

bool IonizationF::run (
	eUnit Eunits, lUnit Lunits,
	sqlitepp::session & db,
	std::map<std::string,std::string> const & sdata
) const {
	
	// manage units
	double efactor = change_units(Eunits, eUnit_Ry);
	double lfactor = change_units(lUnit_au, Lunits);
	
	// atomic and projectile data
	int ni = As<int>(sdata, "ni", Id);
	int li = As<int>(sdata, "li", Id);
	int mi = As<int>(sdata, "mi", Id);
	int  L = As<int>(sdata,  "L", Id);
	int  S = As<int>(sdata,  "S", Id);
	int l1 = As<int>(sdata, "l1", Id);
	int l2 = As<int>(sdata, "l2", Id);
	double Ei = As<double>(sdata, "Ei", Id) * efactor;
	
	// read energy sharing (in user units)
	rArray Eshare;
	try {
		Eshare.push_back(As<double>(sdata, "Eshare", Id));
	} catch (std::exception e) {
		Eshare = readStandardInput<double>();
	}
	
	// energy and encoded Chebyshev approximation
	double E;
	std::string blob;
	
	// create query statement
	sqlitepp::statement st(db);
	st << "SELECT Ei, QUOTE(cheb) FROM " + IonizationF::Id + " "
	      "WHERE ni = :ni "
		  "  AND li = :li "
		  "  AND mi = :mi "
		  "  AND  L = :L  "
		  "  AND  S = :S  "
		  "  AND l1 = :l1 "
		  "  AND l2 = :l2 "
		  "ORDER BY Ei ASC",
		  
	   sqlitepp::into(E), sqlitepp::into(blob),
	   sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
	   sqlitepp::use(L), sqlitepp::use(S),
	   sqlitepp::use(l1), sqlitepp::use(l2);
	
	// get Chebyshev expansions
	rArray E_arr;
	cArray cb;
	std::vector<Chebyshev<double,Complex>> cheb_arr;
	while (st.exec())
	{
		// save energy
		E_arr.push_back(E);
		
		// decode Chebyshev expansion from hexadecimal format
		cb.fromBlob(blob);
		
		// save Chebyshev expansion
		cheb_arr.push_back (
			Chebyshev<double,Complex> (
				cb,                   // expansion coefficients
				0.,                   // lowest energy
				sqrt(E - 1./(ni*ni))  // highest energy
			)
		);
	}
	
	// for all energy shares
	cArray f_out(Eshare.size());
	for (size_t i = 0; i < Eshare.size(); i++)
	{
		// initialize k₁ and k₂ so that
		//   1) (k₁)² + (k₂)² = Ei
		//   2) (k₁)² / (k₂)² = Eshare / (1 - Eshare)
		rArray k1 = sqrt((E_arr - 1./(ni*ni)) * Eshare[i]);
		rArray k2 = sqrt((E_arr - 1./(ni*ni)) * (1 - Eshare[i]));
		
		// for all impact energies evaluate the radial part
		cArray f0(E_arr.size());
		for (size_t ie = 0; ie < E_arr.size(); ie++)
			f0[ie] = cheb_arr[ie].clenshaw(k1[ie], cheb_arr[ie].tail(1e-8)) / gsl_sf_pow_int(k1[ie] * k2[ie], 2);
		
		// interpolate
		f_out[i] = interpolate(E_arr, f0, { Ei })[0];
	}
	
	// write out
	std::cout << this->logo() <<
		"# Ionization amplitudes (radial part) in " << unit_name(Lunits) << " for\n" <<
		"#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
	    "#     L = " << L << ", S = " << S << ", ℓ₁ = " << l1 << ", ℓ₂ = " << l2 << "\n" <<
		"# and impact energy\n" <<
	    "#     Ei = " << Ei << " in " << unit_name(Eunits) << "\n" <<
	    "# ordered by energy share " << 
	    "# \n" <<
	    "# Eshare\t Re f\t Im f\n";
	for (size_t i = 0; i < Eshare.size(); i++)
	{
		std::cout << 
			Eshare[i] << "\t" << 
			f_out[i].real()*lfactor << "\t" <<
			f_out[i].imag()*lfactor << "\n";
	}
	
	return true;
}
