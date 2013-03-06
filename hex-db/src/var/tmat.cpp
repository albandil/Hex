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

const std::string TMatrix::Id = "tmat";
const std::string TMatrix::Description = "T-matrix.";
const std::vector<std::string> TMatrix::Dependencies = {
	"ni", "li", "mi", 
	"nf", "lf", "mf",
	"L", "S",
	"Ei", "ell"
};

std::string const & TMatrix::SQL_CreateTable() const
{
	static const std::string cmd = "CREATE TABLE '" + Id + "' ("
		"ni INTEGER, "
		"li INTEGER, "
		"mi INTEGER, "
		"nf INTEGER, "
		"lf INTEGER, "
		"mf INTEGER, "
		"L  INTEGER, "
		"S  INTEGER, "
		"Ei DOUBLE PRECISION, "
		"ell INTEGER, "
		"Re_T_ell DOUBLE PRECISION, "
		"Im_T_ell DOUBLE PRECISION, "
		"PRIMARY KEY (ni,li,mi,nf,lf,mf,L,S,Ei,ell)"
	")";
	
	return cmd;
}

std::string const & TMatrix::SQL_Update() const
{
	static const std::string cmd = "";
	return cmd;
}

bool TMatrix::run (
	sqlitepp::session & db,
	std::map<std::string,std::string> const & sdata,
	rArray const & energies
) const {
	
	// atomic and projectile data
	int ni = As<int>(sdata, "ni", Id);
	int li = As<int>(sdata, "li", Id);
	int mi = As<int>(sdata, "mi", Id);
	int nf = As<int>(sdata, "nf", Id);
	int lf = As<int>(sdata, "lf", Id);
	int mf = As<int>(sdata, "mf", Id);
	int  L = As<int>(sdata,  "L", Id);
	int  S = As<int>(sdata,  "S", Id);
	int ell= As<int>(sdata, "ell",Id);
	
	// energy and real and imarinary part of the T-matrix
	double E, Re_T_ell, Im_T_ell;
	
	// create query statement
	sqlitepp::statement st(db);
	st << "SELECT Ei, Re_T_ell, Im_T_ell FROM " + TMatrix::Id + " "
	      "WHERE ni = :ni "
		  "  AND li = :li "
		  "  AND mi = :mi "
		  "  AND nf = :nf "
		  "  AND lf = :lf "
		  "  AND mf = :mf "
		  "  AND  L = :L  "
		  "  AND  S = :S  "
		  "  AND ell=:ell "
		  "ORDER BY Ei ASC",
	   sqlitepp::into(E), sqlitepp::into(Re_T_ell), sqlitepp::into(Im_T_ell),
	   sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
	   sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
	   sqlitepp::use(L), sqlitepp::use(S), sqlitepp::use(ell);
	
	// get T-matrices
	rArray E_arr;
	cArray T_arr;
	while (st.exec())
	{
		E_arr.push_back(E);
		T_arr.push_back(Complex(Re_T_ell,Im_T_ell));
	}
	
	// interpolate
	cArray T_out = interpolate(E_arr, T_arr, energies);
	
	// write out
	std::cout << this->logo() <<
		"# T-matrices for\n" <<
		"#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
	    "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n" <<
	    "#     L = " << L << ", S = " << S << ", ℓ = " << ell << "\n" <<
	    "# ordered by energy in Rydbergs\n" <<
	    "# \n" <<
	    "# E\t Re T\t Im T\n";
	for (size_t i = 0; i < energies.size(); i++)
	{
		std::cout << 
			energies[i] << "\t" << 
			T_out[i].real() << "\t" <<
			T_out[i].imag() << "\n";
	}
	
	return true;
}