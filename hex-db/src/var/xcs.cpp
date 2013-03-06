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

#include "../interpolate.h"
#include "../specf.h"
#include "../variables.h"

const std::string ExtrapolatedCrossSection::Id = "xcs";
const std::string ExtrapolatedCrossSection::Description = "Extrapolated cross section (using Aitken Δ²-process).";
const std::vector<std::string> ExtrapolatedCrossSection::Dependencies = {
	"ni", "li", "mi",
	"nf", "lf", "mf",
	"Ei"
};

std::string const & ExtrapolatedCrossSection::SQL_CreateTable() const
{
	static const std::string cmd = "";
	return cmd;
}

std::string const & ExtrapolatedCrossSection::SQL_Update() const
{
	static const std::string cmd = "";
	return cmd;
}

bool ExtrapolatedCrossSection::run (
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
	
	// SQL interface variables
	double E, sigma; int delta;
	
	// last contributions to the complete integral cross section
	// for every energy from E0
	rArray E0, dsigma0;
	
	// last but one contributions to the complete integral cross section
	// for every energy from E0
	rArray E1, dsigma1;
	
	// compose query
	sqlitepp::statement st(db);
	st << "SELECT Dat.Ei AS Ei, Lim.L - Dat.L AS Delta, Dat.sigma                    "
	      "FROM (                                                                    "
	      "  SELECT Ei, L, SUM(sigma) AS sigma FROM " + IntegralCrossSection::Id + " "
	      "  WHERE ni = :ni AND li = :li AND mi = :mi                                "
	      "    AND nf = :nf AND lf = :lf AND mf = :mf                                "
	      "  GROUP BY Ei, L                                                          "
	      "  ORDER BY Ei ASC                                                         "
		  ") AS Dat CROSS JOIN (                                                     "
		  "  SELECT Ei, MAX(L) AS L FROM " + IntegralCrossSection::Id + "            "
		  "  WHERE ni = :Ni AND li = :Li AND mi = :Mi                                "
	      "    AND nf = :Nf AND lf = :Lf AND mf = :Mf                                "
		  "  GROUP BY Ei                                                             "
		  "  ORDER BY Ei ASC                                                         "
		  ") AS Lim ON Dat.Ei = Lim.Ei                                               "
		  "WHERE Dat.L = Lim.L OR Dat.L = Lim.L - 1                                  "
		  "ORDER BY Ei ASC, Delta DESC",
		sqlitepp::into(E), sqlitepp::into(delta), sqlitepp::into(sigma),
		sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
		sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
		sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
		sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf);
	
	// retrieve last and last but one contributions to the cross sections
	while (st.exec())
	{
		if (delta == 0)
		{
			E0.push_back(E);
			dsigma0.push_back(sigma);
		}
		else /* delta == 1 */
		{
			E1.push_back(E);
			dsigma1.push_back(sigma);
		}
	}
	
	// full complete cross section
	rArray Efull, sigmafull;
	
	// compose query
	sqlitepp::statement st1(db);
	st1 << "SELECT Ei, sigma FROM " + CompleteCrossSection::Id + " "
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
	while (st1.exec())
	{
		Efull.push_back(E);
		sigmafull.push_back(sigma);
	}
	
	// compute complete cross section
	rArray ccs = interpolate(Efull, sigmafull, energies);
	
	// reshape arrays for compatible indexing by energies
	rArray dsigma0empty(dsigma0.size());	// zero-filled ghost
	rArray dsigma1empty(dsigma1.size());	// zero-filled ghost
	merge (Efull, sigmafull, E0, dsigma0empty);	// inflate Efull and sigmafull to accomodate E0 and dsigma0
	merge (Efull, sigmafull, E1, dsigma1empty);	// inflate Efull and sigmafull to accomodate E1 and dsigma1
	rArray sigmafullempty(sigmafull.size());	// // zero-filled ghost
	merge (E0, dsigma0, Efull, sigmafullempty);	// inflate E0 and dsigma0
	merge (E1, dsigma1, Efull, sigmafullempty);	// inflate E1 and dsigma1
	
	// extrapolate using Aitken Δ²-method
	sigmafull -= dsigma0 * dsigma0 / (dsigma0 - dsigma1);
	
	// interpolate
	rArray xcs = interpolate(Efull, sigmafull, energies);
	
	// write out
	std::cout << this->logo() <<
		"# Extrapolated and complete cross section for\n" <<
		"#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
	    "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n" <<
	    "# ordered by energy in Rydbergs\n" <<
		"#\n" <<
	    "# E\t σx\t σc\n";
	for (size_t i = 0; i < energies.size(); i++)
		std::cout << energies[i] << "\t" << xcs[i] << "\t" << ccs[i] << "\n";
	
	return true;
}