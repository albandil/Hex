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
const std::vector<std::string> TotalCrossSection::VecDependencies = { "Ei" };

std::vector<std::string> const & TotalCrossSection::SQL_Update() const
{
	static const std::vector<std::string> cmd;
	return cmd;
}

std::vector<std::string> const & TotalCrossSection::SQL_CreateTable() const
{
	static const std::vector<std::string> cmd;
	return cmd;
}

bool TotalCrossSection::run (
	eUnit Eunits, lUnit Lunits,
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
	
	// write header
	std::cout << this->logo() <<
		"# Total cross section in " << unit_name(Lunits) << " for\n"
		"#     ni = " << ni << ", li = " << li << ", mi = " << mi << "\n" <<
	    "# ordered by energy in " << unit_name(Eunits) << "\n" <<
	    "#\n" <<
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
		rArray tcs = (efactor * energies.front() < Eion) ? 
			interpolate_real(E_arr, sigma_arr, energies * efactor, o2scl::itp_linear) :
			interpolate_real(E_arr, sigma_arr, energies * efactor, o2scl::itp_cspline);
		
		// output
		for (size_t i = 0; i < energies.size(); i++)
			std::cout << energies[i] << "\t" << tcs[i]*lfactor*lfactor << "\n";
	}
	
	return true;
}
