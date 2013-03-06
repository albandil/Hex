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

#include "../angs.h"
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
	static const std::string cmd = "";
	
	return cmd;
}

std::string const & MomentumTransfer::SQL_Update() const
{
	static const std::string cmd = "";
	
	return cmd;
}

bool MomentumTransfer::run (
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
	int  S = As<int>(sdata,  "S", Id);
	
	// SQL interface variables
	double E, Re_T_ell, Im_T_ell;
	int L, ell;
	
	// compose the statement
	sqlitepp::statement st(db);
	st << "SELECT Ei, L, ell, Re_T_ell, Im_T_ell FROM " + TMatrix::Id + " "
	      "WHERE ni = :ni "
		  "  AND li = :li "
		  "  AND mi = :mi "
		  "  AND nf = :nf "
		  "  AND lf = :lf "
		  "  AND mf = :mf "
		  "  AND  S = :S  "
		  "ORDER BY L, ell, Ei ASC",
		sqlitepp::into(E), sqlitepp::into(L), sqlitepp::into(ell), 
		sqlitepp::into(Re_T_ell), sqlitepp::into(Im_T_ell),
		sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
		sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
		sqlitepp::use(S);
	
	// auxiliary variables
	rArray *pE = nullptr;	// pointer to current energy array
	cArray *pT = nullptr;	// pointer to current T-matrix array
	std::pair<int,int> key(-1,-1);	// current (L,ell) pair
	int maxL = -1;	// maximal L
	int maxl = -1;	// maximal ell
		
	// load T-matrices
	std::map<std::pair<int,int>, std::pair<rArray*, cArray*>> Tmatrices;
	while (st.exec())
	{
		// create new arrays if necessary
		if (key.first != L or key.second != ell)
		{
			pE = new rArray;
			pT = new cArray;
			key = std::make_pair(L,ell);
			Tmatrices[key] = std::make_pair(pE,pT);
			maxL = (L > maxL) ? L : maxL;
			maxl = (ell > maxl) ? ell : maxl;
		}
		
		// insert new data
		pE->push_back(E);
		pT->push_back(Complex(Re_T_ell,Im_T_ell));
	}
	
	// compute momentum transfer
	rArray eta_energies;
	rArray eta;
	for (int L = 0; L <= maxL; L++)
	for (int Lp = 0; Lp <= maxL; Lp++)
	for (int l = 0; l <= maxl; l++)
	for (int lp = 0; lp <= maxl; lp++)
	{
		std::pair<int,int> key(L,l), keyp(Lp,lp);
		
		// if no data, skip
		if (Tmatrices.find(key) == Tmatrices.end() or Tmatrices.find(keyp) == Tmatrices.end())
			continue;
		
		// get energies and T-matrices pointers
		rArray *pE = Tmatrices[key].first;
		cArray *pT = Tmatrices[key].second;
		rArray *pEp = Tmatrices[keyp].first;
		cArray *pTp = Tmatrices[keyp].second;
		
		// normalize
		cArray emptyp(pEp->size());	// zero-filled ghost
		merge (*pE, *pT, *pEp, emptyp);	// inflate *pE and *pT to accomodate *pEp and emptyp
		cArray empty(pE->size()); // zero-filled ghost
		merge (*pEp, *pTp, *pE, empty);
		
		// compute factor
		double factor = 0;
		for (int m = -l; m <= l; m++)
			factor += ClebschGordan(l,m,1,0,lp,m);
		factor *= -sqrt((2.*l+1.)/(2.*lp+1.)) * ClebschGordan(l,0,1,0,lp,0);
		if (l == lp)
			factor += 1.;
		
		// update momentum transfer
		merge (eta_energies, eta, *pE, factor * abs((*pT) * (*pTp).conj()));
	}
	eta *= 1. / (4. * M_PI * M_PI);
	
	// interpolate
	eta = interpolate(eta_energies, eta, energies);
	
	// write out
	std::cout << this->logo() <<
		"# Momentum transfer for\n" <<
		"#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
	    "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n" <<
	    "#     S = " << S << "\n" <<
	    "# ordered by energy in Rydbergs\n" <<
		"#\n" <<
	    "# E\t η\n";
	for (size_t i = 0; i < energies.size(); i++)
		std::cout << energies[i] << "\t" << eta[i] << "\n";
	
	return true;
}