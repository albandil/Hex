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
	rArray const & nums
) const {
	
	// TODO
	return false;
}

/*
rArray ExtrapolatedCrossSection::aitkenD2(int ni, int li, int mi, int nf, int lf, int mf, rArray Ei)
{
	// vector size
	size_t Nenergy = Ei.size();
	
	// the cross section
	rArray sigma(Nenergy), dsigma_prev(Nenergy), dsigma(Nenergy), relchng(Nenergy);
	
	// get maximal "big" partial wave angular momentum
	int max_L;
	sqlitepp::statement st1(db);
	st1 << "SELECT MAX(L) FROM hex "
			"WHERE ni = :ni AND li = :li AND mi = :mi AND nf = :nf AND lf = :lf AND mf = :mf",
		sqlitepp::into(max_L),
		sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
		sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf);
	if (not st1.exec())
	{
		printf ("no data for these quantum numbers\n");
		return sigma;
	}
	
	// check if enough samples are present
	if (max_L < 2)
	{
		printf("not enough data for Aitken Δ²-proces (max L = %d)\n", max_L);
		return CompleteCrossSection::compute(ni,li,mi,nf,lf,mf,Ei);
	}
	
	// load all cross sections
	rArray singlet, triplet;
	rArrays sigmas;
	for (int L = 0; L <= max_L; L++)
	{
		// get spin-dependent cross sections
		singlet = IntegralCrossSection::compute(ni, li, mi, nf, lf, mf, L, 0, Ei);
		triplet = IntegralCrossSection::compute(ni, li, mi, nf, lf, mf, L, 1, Ei);
		
		// compute L-dependent cross section
		if (L == 0)
			sigmas.push_back(singlet + triplet);
		else
			sigmas.push_back(singlet + triplet + sigmas.back());
	}
	
	// extrapolate, σn' = σn - (Δσn)² / (Δ²σn)
	rArray sig_n   = sigmas[max_L];
	rArray sig_nm1 = sigmas[max_L-1];
	rArray sig_nm2 = sigmas[max_L-2];
	
	sigma = sig_n - pow(sig_n - sig_nm1, 2) / (sig_n - 2*sig_nm1 + sig_nm2);
	
	// return the result
	return sigma;
}
*/

/*
rArray ExtrapolatedCrossSection::extrapolate(int ni, int li, int mi, int nf, int lf, int mf, rArray Ei)
{
	// vector size
	size_t Nenergy = Ei.size();
	
	// the cross section
	rArray sigma(Nenergy), dsigma_prev(Nenergy), dsigma(Nenergy), relchng(Nenergy);
	
	// get maximal "big" partial wave angular momentum
	int max_L;
	sqlitepp::statement st1(db);
	st1 << "SELECT MAX(L) FROM hex "
			"WHERE ni = :ni AND li = :li AND mi = :mi AND nf = :nf AND lf = :lf AND mf = :mf",
		sqlitepp::into(max_L),
		sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
		sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf);
	if (not st1.exec())
	{
		printf ("no data for these quantum numbers\n");
		return sigma;
	}
	
	// for all "big" partial waves
	for (int L = 0; L <= max_L; L++)
	{
		rArray singlet = IntegralCrossSection::compute(ni, li, mi, nf, lf, mf, L, 0, Ei);
		rArray triplet = IntegralCrossSection::compute(ni, li, mi, nf, lf, mf, L, 1, Ei);
		
		for (size_t i = 0; i < Nenergy; i++)
		{
			dsigma[i] = singlet[i] + triplet[i];
			sigma[i] += dsigma[i];
		
			if (L > 0)
				relchng[i] = dsigma[i] / dsigma_prev[i];
			
			dsigma_prev[i] = dsigma[i];
		}
	}
	
	// extrapolate
	if (max_L > 1) for (size_t i = 0; i < Nenergy; i++)
		sigma[i] += dsigma[i] / (1. - relchng[i]);
	
	// return the result
	return sigma;
}
*/

