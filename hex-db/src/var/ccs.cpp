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

const std::string CompleteCrossSection::Id = "ccs";
const std::string CompleteCrossSection::Description = "Complete cross section (L- and S-summed integral cross section).";
const std::vector<std::string> CompleteCrossSection::Dependencies = {
	"ni", "li", "mi",
	"nf", "lf", "mf",
	"Ei"
};

std::string const & CompleteCrossSection::SQL_CreateTable() const
{
	static std::string cmd = "CREATE TABLE '" + CompleteCrossSection::Id + "' ("
		"ni INTEGER, "
		"li INTEGER, "
		"mi INTEGER, "
		"nf INTEGER, "
		"lf INTEGER, "
		"mf INTEGER, "
		"Ei DOUBLE PRECISION, "
		"sigma DOUBLE PRECISION, "
		"PRIMARY KEY (ni,li,mi,nf,lf,mf,Ei)"
	")";
	
	return cmd;
}
	
std::string const & CompleteCrossSection::SQL_Update() const
{
	static std::string cmd = "INSERT OR REPLACE INTO '" + CompleteCrossSection::Id + "' "
		"SELECT ni, li, mi, nf, lf, mf, Ei, sum(sigma) "
		"FROM '" + IntegralCrossSection::Id + "' "
		"GROUP BY ni, li, mi, nf, lf, mf, Ei";
	
	return cmd;
}

bool CompleteCrossSection::run (
	eUnit Eunits, lUnit Lunits,
	sqlitepp::session & db,
	std::map<std::string,std::string> const & sdata,
	rArray const & energies
) const {
	
	// manage units
	double efactor = change_units(Eunits, eUnit_Ry);
	double lfactor = change_units(lUnit_au, Lunits);
	
	// scattering event parameters
	int ni = As<int>(sdata, "ni", Id);
	int li = As<int>(sdata, "li", Id);
	int mi = As<int>(sdata, "mi", Id);
	int nf = As<int>(sdata, "nf", Id);
	int lf = As<int>(sdata, "lf", Id);
	int mf = As<int>(sdata, "mf", Id);
	
	// energies and cross sections
	double E, sigma;
	rArray E_arr, sigma_arr;
	
	// compose query
	sqlitepp::statement st(db);
	st << "SELECT Ei, sigma FROM " + CompleteCrossSection::Id + " "
			"WHERE ni = :ni "
			"  AND li = :li "
			"  AND mi = :mi "
			"  AND nf = :nf "
			"  AND lf = :lf "
			"  AND mf = :mf "
			"ORDER BY Ei ASC",
		sqlitepp::into(E), sqlitepp::into(sigma),
		sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
		sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf);
	
	// retrieve data
	while (st.exec())
	{
		E_arr.push_back(E);
		sigma_arr.push_back(sigma);
	}
	
	// interpolate
	rArray ccs = interpolate(E_arr, sigma_arr, energies * efactor);
	
	// write out
	std::cout << this->logo() <<
		"# Complete cross section in " << unit_name(Lunits) << " for\n" <<
		"#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
	    "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n" <<
	    "# ordered by energy in " << unit_name(Eunits) << "\n" <<
	    "# \n" <<
	    "# E\t σ\n";
	for (size_t i = 0; i < energies.size(); i++)
		std::cout << energies[i] << "\t" << (finite(ccs[i]) ? ccs[i] * lfactor * lfactor : 0.) << "\n";
	
	return true;
}
