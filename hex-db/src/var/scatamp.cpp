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

const std::string ScatteringAmplitude::Id = "scatamp";
const std::string ScatteringAmplitude::Description = "Scattering amplitude.";
const std::vector<std::string> ScatteringAmplitude::Dependencies = {
	"ni", "li", "mi", 
	"nf", "lf", "mf",
	"S", "Ei"
};

std::string const & ScatteringAmplitude::SQL_Update() const
{
	static const std::string cmd = "";
	return cmd;
}

std::string const & ScatteringAmplitude::SQL_CreateTable() const
{
	static const std::string cmd = "";
	return cmd;
}

bool ScatteringAmplitude::run (
	sqlitepp::session & db,
	std::map<std::string,std::string> const & sdata,
	rArray const & angles
) const {
	
	// scattering event parameters
	int ni = As<int>(sdata, "ni", Id);
	int li = As<int>(sdata, "li", Id);
	int mi = As<int>(sdata, "mi", Id);
	int nf = As<int>(sdata, "nf", Id);
	int lf = As<int>(sdata, "lf", Id);
	int mf = As<int>(sdata, "mf", Id);
	int  S = As<int>(sdata, "S", Id);
	double E = As<double>(sdata, "E", Id);
	
	// the scattering amplitudes
	cArray amplitudes(angles.size());
	
	// total angular momentum projection (given by axis orientation)
	int M = mi;
	
	// get maximal partial wave angular momentum
	int max_L;
	sqlitepp::statement st3(db);
	st3 << "SELECT MAX(L) FROM " + TMatrix::Id + " "
	       "WHERE ni = :ni "
		   "  AND li = :li "
		   "  AND mi = :mi "
		   "  AND nf = :nf "
		   "  AND lf = :lf "
		   "  AND mf = :mf "
		   "  AND  S = :S  ",
		sqlitepp::into(max_L),
		sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
		sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
		sqlitepp::use(S);
	st3.exec();
	
	// for all total angular momenta
	for (int L = M; L <= max_L; L++)
	{
		// get maximal projectile partial wave angular momentum
		int max_ell;
		sqlitepp::statement st1(db);
		st1 << "SELECT MAX(ell) FROM " + TMatrix::Id + " "
		       "WHERE ni = :ni "
			   "  AND li = :li "
			   "  AND mi = :mi "
			   "  AND nf = :nf "
			   "  AND lf = :lf "
			   "  AND mf = :mf "
			   "  AND  L = :L  "
			   "  AND  S = :S  ",
			sqlitepp::into(max_ell),
			sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
			sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
			sqlitepp::use(L), sqlitepp::use(S);
		st1.exec();
		
		// for all outgoing partial waves
		for (int ell = abs(M - mf); ell <= max_ell; ell++)
		{
			// get all relevant lines from database
			double __Ei, __Re_T_ell, __Im_T_ell;
			rArray db_Ei;
			cArray db_T_ell;
			sqlitepp::statement st2(db);
			st2 << "SELECT Ei, Re_T_ell, Im_T_ell FROM " + TMatrix::Id + " "
				   "WHERE ni = :ni "
				   "  AND li = :li "
				   "  AND mi = :mi "
				   "  AND nf = :nf "
				   "  AND lf = :lf "
				   "  AND mf = :mf "
				   "  AND ell=:ell "
				   "  AND  L = :L  "
				   "  AND  S = :S  "
				   "ORDER BY Ei ASC",
				sqlitepp::into(__Ei), sqlitepp::into(__Re_T_ell), sqlitepp::into(__Im_T_ell),
				sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
				sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
				sqlitepp::use(ell), sqlitepp::use(L), sqlitepp::use(S);
			while ( st2.exec() )
			{
				db_Ei.push_back(__Ei);
				db_T_ell.push_back(Complex(__Re_T_ell, __Im_T_ell));
			}
			
			if (db_Ei.size() == 0)
				continue;
			
			// update value of "f"
			Complex Tmatrix = interpolate(db_Ei, db_T_ell, {E});
			for (size_t i = 0; i < angles.size(); i++)
				amplitudes[i] += -1./(2.*M_PI) * Tmatrix * sphY(ell, abs(M-mf), angles[i], 0.);
		}
	}
	
	// write out
	std::cout << "# Scattering amplitudes for "
		"ni = " << ni << ", li = " << li << ", mi = " << mi << ", " <<
	    "nf = " << nf << ", lf = " << lf << ", mf = " << mf << ", " <<
	    "S = " << S << ", E = " << E << " " <<
	    " ordered by angle in radians\n" <<
	    "# θ\t Re f\t Im f\n";
	for (size_t i = 0; i < angles.size(); i++)
	{
		std::cout << 
			angles[i] << "\t" << 
			amplitudes[i].real() << "\t" <<
			amplitudes[i].imag() << "\n";
	}
	
	return true;
}
