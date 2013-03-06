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

const std::string TotalCrossSection::Id = "tcs";
const std::string TotalCrossSection::Description = "Total cross section.";
const std::vector<std::string> TotalCrossSection::Dependencies = {
	"ni", "li", "mi", 
	"nf", "lf", "mf",
	"Ei"
};

std::string const & TotalCrossSection::SQL_Update() const
{
	static const std::string cmd = "";
	return cmd;
}

std::string const & TotalCrossSection::SQL_CreateTable() const
{
	static const std::string cmd = "";
	return cmd;
}

bool TotalCrossSection::run (
	sqlitepp::session & db,
	std::map<std::string,std::string> const & sdata,
	rArray const & energies
) const {
	
	// scattering event parameters
	int ni = As<int>(sdata, "ni", Id);
	int li = As<int>(sdata, "li", Id);
	int mi = As<int>(sdata, "mi", Id);
	
	// energies and cross sections
	double E, sigma;
	rArray E_arr, sigma_arr;
	
	// compose query
	sqlitepp::statement st(db);
	st << "SELECT Ei, sum(sigma) FROM " + CompleteCrossSection::Id + " "
			"WHERE ni = :ni "
			"  AND li = :li "
			"  AND mi = :mi "
			"GROUP BY Ei "
			"ORDER BY Ei ASC",
		sqlitepp::into(E),   sqlitepp::into(sigma),
		sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi);
	
	// retrieve data
	while (st.exec())
	{
		E_arr.push_back(E);
		sigma_arr.push_back(sigma);
	}
	
	// interpolate
	rArray tcs = interpolate(E_arr, sigma_arr, energies);
	
	// write out
	std::cout << this->logo() <<
		"# Total cross section for\n"
		"#     ni = " << ni << ", li = " << li << ", mi = " << mi << "\n" <<
	    "# ordered by energy in Rydbergs\n" <<
	    "#\n" <<
	    "# E\t σ\n";
	for (size_t i = 0; i < energies.size(); i++)
		std::cout << energies[i] << "\t" << tcs[i] << "\n";
	
	return true;
}