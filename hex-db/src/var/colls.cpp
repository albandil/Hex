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
#include <map>
#include <string>
#include <vector>

#include "../interpolate.h"
#include "../variables.h"

const std::string CollisionStrength::Id = "colls";
const std::string CollisionStrength::Description = "Collision strength (energy scaled integral cross section).";
const std::vector<std::string> CollisionStrength::Dependencies = {
	"ni", "li", "mi",
	"nf", "lf", "mf",
	"L", "S",
	"Ei"
};

std::string const & CollisionStrength::SQL_CreateTable() const
{
	static std::string cmd = "";
	
	return cmd;
}

std::string const & CollisionStrength::SQL_Update() const
{
	static std::string cmd = "";
	
	return cmd;
}

bool CollisionStrength::run (
	sqlitepp::session & db,
	std::map<std::string,std::string> const & sdata,
	rArray const & energies
) const {
	
	// scattering event parameters
	int ni = As<int>(sdata, "ni", Id);
	int li = As<int>(sdata, "li", Id);
	int mi = As<int>(sdata, "mi", Id);
	int nf = As<int>(sdata, "nf", Id);
	int lf = As<int>(sdata, "lf", Id);
	int mf = As<int>(sdata, "mf", Id);
	int  L = As<int>(sdata, "L", Id);
	int  S = As<int>(sdata, "S", Id);
	
	// energies and cross sections
	double E, sigma;
	rArray E_arr, sigma_arr;
	
	// compose query
	sqlitepp::statement st(db);
	st << "SELECT Ei, sigma FROM " + IntegralCrossSection::Id + " "
			"WHERE ni = :ni "
			"  AND li = :li "
			"  AND mi = :mi "
			"  AND nf = :nf "
			"  AND lf = :lf "
			"  AND mf = :mf "
			"  AND  L = :L  "
			"  AND  S = :S  "
			"ORDER BY Ei ASC",
		sqlitepp::into(E), sqlitepp::into(sigma),
		sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
		sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
		sqlitepp::use(L), sqlitepp::use(S);
	
	// retrieve data
	while (st.exec())
	{
		E_arr.push_back(E);
		sigma_arr.push_back(sigma);
	}
	
	// interpolate
	rArray omegas = energies * (2*L+1) * (2*S+1) * interpolate(E_arr, sigma_arr, energies);
	
	// write out
	std::cout << "# Collision strength for "
		"ni = " << ni << ", li = " << li << ", mi = " << mi << ", " <<
	    "nf = " << nf << ", lf = " << lf << ", mf = " << mf << ", " <<
	    "L = " << L << ", S = " << S << " " <<
	    " ordered by energy in Rydbergs\n" <<
	    "# E\t Ω\n";
	for (size_t i = 0; i < energies.size(); i++)
		std::cout << energies[i] << "\t" << omegas[i] << "\n";
	
	return true;
}
