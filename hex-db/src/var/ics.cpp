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
#include "../variables.h"

const std::string IntegralCrossSection::Id = "ics";
const std::string IntegralCrossSection::Description = "Integral cross section.";
const std::vector<std::string> IntegralCrossSection::Dependencies = {
	"ni", "li", "mi",
	"nf", "lf", "mf",
	"L", "S",
	"Ei"
};

std::string const & IntegralCrossSection::SQL_CreateTable() const
{
	static const std::string cmd = "CREATE TABLE '" + IntegralCrossSection::Id + "' ("
		"ni INTEGER, "
		"li INTEGER, "
		"mi INTEGER, "
		"nf INTEGER, "
		"lf INTEGER, "
		"mf INTEGER, "
		"L  INTEGER, "
		"S  INTEGER, "
		"Ei DOUBLE PRECISION, "
		"sigma DOUBLE PRECISION, "
		"PRIMARY KEY (ni,li,mi,nf,lf,mf,L,S,Ei)"
	")";
	
	return cmd;
}

std::string const & IntegralCrossSection::SQL_Update() const
{
	static const std::string cmd = "INSERT OR REPLACE INTO " + IntegralCrossSection::Id + " "
		"SELECT ni, li, mi, nf, lf, mf, L, S, Ei, "
		    "sqrt(Ei-1./(ni*ni)+1./(nf*nf))/sqrt(Ei)*(2*S+1)*SUM(Re_T_ell*Re_T_ell+Im_T_ell*Im_T_ell)/157.91367 " // 16π²
		"FROM " + TMatrix::Id + " "
		"WHERE Ei > 0 AND Ei - 1./(ni*ni) + 1./(nf*nf) > 0 "
		"GROUP BY ni, li, mi, nf, lf, mf, L, S, Ei";
	
	return cmd;
}

bool IntegralCrossSection::run (
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
	st << "SELECT Ei, sigma FROM " + Id + " "
			"WHERE ni = :ni "
			"  AND li = :li "
			"  AND mi = :mi "
			"  AND nf = :nf "
			"  AND lf = :lf "
			"  AND mf = :mf "
			"  AND  L = :L  "
			"  AND  S = :S  "
			"ORDER BY Ei ASC;",
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
	rArray ics = interpolate(E_arr, sigma_arr, energies);
	
	// write out
	std::cout << "# Integral cross section for "
		"ni = " << ni << ", li = " << li << ", mi = " << mi << ", " <<
	    "nf = " << nf << ", lf = " << lf << ", mf = " << mf << ", " <<
	    "L = " << L << ", S = " << S << " " <<
	    " ordered by energy in Rydbergs\n" <<
	    "# E\t σ\n";
	for (size_t i = 0; i < energies.size(); i++)
		std::cout << energies[i] << "\t" << ics[i] << "\n";
	
	return true;
}
