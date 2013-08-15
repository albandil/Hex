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
const std::vector<std::string> CompleteCrossSection::VecDependencies = { "Ei" };

bool CompleteCrossSection::initialize(sqlitepp::session & db) const
{
	return true;
}

std::vector<std::string> const & CompleteCrossSection::SQL_CreateTable() const
{
	static std::vector<std::string> cmd = {
		"CREATE TABLE '" + CompleteCrossSection::Id + "' ("
			"ni INTEGER, "
			"li INTEGER, "
			"mi INTEGER, "
			"nf INTEGER, "
			"lf INTEGER, "
			"mf INTEGER, "
			"Ei DOUBLE PRECISION, "
			"sigma DOUBLE PRECISION, "
			"PRIMARY KEY (ni,li,mi,nf,lf,mf,Ei)"
		")"
	};
	
	return cmd;
}
	
std::vector<std::string> const & CompleteCrossSection::SQL_Update() const
{
	static std::vector<std::string> cmd = {
		"INSERT OR REPLACE INTO '" + CompleteCrossSection::Id + "' "
			"SELECT ni, li, mi, nf, lf, mf, Ei, sum(sigma) "
			"FROM '" + IntegralCrossSection::Id + "' "
			"GROUP BY ni, li, mi, nf, lf, mf, Ei"
	};
	
	return cmd;
}

bool CompleteCrossSection::run (
	sqlitepp::session & db,
	std::map<std::string,std::string> const & sdata
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
	rArray energies, E_arr, sigma_arr;
	
	// get energy / energies
	try {
		
		// is there a single energy specified using command line ?
		energies.push_back(As<double>(sdata, "Ei", Id));
		
	} catch (std::exception e) {
		
		// are there more energies specified using the STDIN ?
		energies = readStandardInput<double>();
	}
	
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
	
	// write header
	std::cout << this->logo() <<
		"# Complete cross section in " << unit_name(Lunits) << " for\n" <<
		"#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
	    "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n" <<
	    "# ordered by energy in " << unit_name(Eunits) << "\n" <<
	    "# \n" <<
	    "# E\t σ\n";
	
	if (energies[0] < 0.)
	{
		// negative energy indicates full output
		for (size_t i = 0; i < E_arr.size(); i++)
			std::cout << E_arr[i] / efactor << "\t" << sigma_arr[i] * lfactor * lfactor << "\n";
	}
	else
	{
		// threshold for ionization
		double Eion = 1./(ni*ni);
		
		// interpolate
		rArray ccs = (efactor * energies.front() < Eion) ? 
			interpolate_real(E_arr, sigma_arr, energies * efactor, o2scl::itp_linear) :
			interpolate_real(E_arr, sigma_arr, energies * efactor, o2scl::itp_cspline);
			
		// output
		for (size_t i = 0; i < energies.size(); i++)
			std::cout << energies[i] << "\t" << (finite(ccs[i]) ? ccs[i] * lfactor * lfactor : 0.) << "\n";
	}
	
	return true;
}
