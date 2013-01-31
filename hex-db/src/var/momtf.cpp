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

const std::string MomentumTransfer::Id = "momtf";
const std::string MomentumTransfer::Description = "Momentum transfer.";
const std::vector<std::string> MomentumTransfer::Dependencies = {
	"ni", "li", "mi", 
	"nf", "lf", "mf",
	"S", "Ei"
};

std::string const & MomentumTransfer::SQL_CreateTable() const
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
		"sigma DOUBLE PRECISION, "
		"PRIMARY KEY (ni,li,mi,nf,lf,mf,L,S,Ei)"
	")";
	
	return cmd;
}

std::string const & MomentumTransfer::SQL_Update() const
{
	static const std::string cmd = "";
	
	return cmd;
}

bool MomentumTransfer::run (
	sqlitepp::session & db,
	std::map<std::string,std::string> const & data1,
	rArray const & data2
) const {
	
	// TODO
	return false;
}

/*
rArray MomentumTransfer::compute(int ni, int li, int mi, int nf, int lf, int mf, int L, int S, rArray Ei)
{
	rArray eta(Ei.size());
	
	// get maximal partial wave angular momentum
	int max_ell;
	sqlitepp::statement st1(db);
	
	st1 << "SELECT MAX(ell) FROM hex WHERE ni = :ni AND li = :li AND mi = :mi AND nf = :nf AND lf = :lf AND mf = :mf AND L = :L AND S = :S",
		sqlitepp::into(max_ell),
		sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
		sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
		sqlitepp::use(L), sqlitepp::use(S);
	
	if (not st1.exec())
	{
		printf ("no data for these quantum numbers\n");
		return eta;
	}
	
	// load all partial wave T-matrices
	cArray Tmatrices(Ei.size() * (max_ell + 1));
	for (int ell = 0; ell <= max_ell; ell++)
	{
		
		// get all relevant lines from database
		double __Ei, __Re_T_ell, __Im_T_ell;
		rArray db_Ei;
		cArray db_T_ell;
		sqlitepp::statement st2(db);
		st2 << "SELECT Ei, Re_T_ell, Im_T_ell FROM hex "
		      "WHERE ni = :ni AND li = :li AND mi = :mi AND nf = :nf AND lf = :lf AND mf = :mf AND ell = :ell AND L = :L AND S = :S "
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
		
		// update value of "sigma"
		for (size_t ei = 0; ei < Ei.size(); ei++)
			Tmatrices[ei*(max_ell+1) + ell] = interpolate(db_Ei, db_T_ell, Ei[ei]);
	}
	
	// compute the momentum transfer
	for (size_t ei = 0; ei < Ei.size(); ei++)
	{
		// check conservation of energy
		if (Ei[ei] <= 0. or Ei[ei] - 1./(ni*ni) + 1./(nf*nf) <= 0.)
		{
			printf("energy not conserved for ni = %d, nf = %d, Ei = %g\n", ni, nf, Ei[ei]);
			continue;
		}
		
		// for all angular momenta
		for (int ell = 0; ell <= max_ell; ell++)
		{
			// cross section contribution
			Complex T = Tmatrices[ei * (max_ell + 1) + ell];
			eta[ei] += (T.real() * T.real() + T.imag() * T.imag());
			
			// shift contribution
			for (int elp = abs(ell - 1); elp <= ell + 1; elp++)
			{
				Complex Tpc = conj(Tmatrices[elp]);
				eta[ei] += sqrt(4 * M_PI / 3) * (T*Tpc).real() * Gaunt(ell,mi-mf,1,0,elp,mi-mf);
			}
		}
	}
	
	// apply prefactor
	for (size_t ei = 0; ei < Ei.size(); ei++)
		eta[ei] /= 2.*M_PI*2.*M_PI;
	
	return eta;
}
*/