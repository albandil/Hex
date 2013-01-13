#include <algorithm>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "sqlitepp/sqlitepp.hpp"

#include "angs.h"
#include "arrays.h"
#include "complex.h"
#include "hex-db.h"
#include "specf.h"
#include "interpolate.h"

// #define __HEX_DB_NAME__ "hex.db"
sqlitepp::session db;

void initialize(const char* dbname)
{
	// open database
	db.open(dbname);
}

void create_new_database(const char* dbname)
{
	char cmd[1024];
	sprintf(cmd, "sqlite3 %s 'create table \"hex\" (ni integer, li integer, mi integer, "
	"nf integer, lf integer, mf integer, L integer, S integer, Ei double precision, "
	"ell integer, Re_T_ell double precision, Im_T_ell double precision, primary "
	"key (ni,li,mi,nf,lf,mf,L,S,Ei,ell))'", dbname);
	system(cmd);
}

cArray scattering_amplitude(int ni, int li, int mi, int nf, int lf, int mf, int S, double Ei, rArray theta)
{
	// the amplitudes
	cArray amplitudes(theta.size());
	
	// check conservation of energy
	if (Ei <= 0. or Ei - 1./(nf*nf) + 1./(ni*ni) <= 0.)
	{
		printf("energy not conserved for ni = %d, nf = %d, Ei = %g\n", ni, nf, Ei);
		return amplitudes;
	}
	
	// total angular momentum projection (given by axis orientation)
	int M = mi;
	
	// get maximal partial wave angular momentum
	int max_L;
	sqlitepp::statement st3(db);
	st3 << "SELECT MAX(L) FROM hex WHERE ni = :ni AND li = :li AND mi = :mi AND nf = :nf AND lf = :lf AND mf = :mf AND S = :S",
		sqlitepp::into(max_L),
		sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
		sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
		sqlitepp::use(S);
	st3.exec();
	
	for (int L = M; L <= max_L; L++)
	{
		// get maximal projectile partial wave angular momentum
		int max_ell;
		sqlitepp::statement st1(db);
		st1 << "SELECT MAX(ell) FROM hex WHERE ni = :ni AND li = :li AND mi = :mi AND nf = :nf AND lf = :lf AND mf = :mf AND L = :L AND S = :S",
			sqlitepp::into(max_ell),
			sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
			sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
			sqlitepp::use(L), sqlitepp::use(S);
		st1.exec();
		
		// for all outgoing partial waves
		for (int ell = abs(M-mf); ell <= max_ell; ell++)
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
			
			if (db_Ei.size() == 0)
				continue;
			
			// update value of "f"
			Complex Tmatrix = interpolate(db_Ei, db_T_ell, Ei);
			for (size_t i = 0; i < amplitudes.size(); i++)
				amplitudes[i] += -1./(2.*M_PI) * abs(sphY(ell, abs(M-mf), theta[i], 0.)) * Tmatrix;
		}
	}
	
	// return the result
	return amplitudes;
}

rArray differential_cross_section(int ni, int li, int mi, int nf, int lf, int mf, int S, double Ei, rArray theta)
{
	// the scattering amplitudes sum and the differential cross section
	cArray f = scattering_amplitude(ni, li, mi, nf, lf, mf, S, Ei, theta);
	rArray dcs(theta.size());
	
	// compute projectile momenta
	double ki = sqrt(Ei);
	double kf = sqrt(Ei - 1./(ni*ni) + 1./(nf*nf));
	
	// compute cross section from the amplitudes
	dcs = kf / ki * (2. * S + 1.) / 4. * pow(abs(f), 2);
	
	return dcs;
}

rArray momentum_transfer(int ni, int li, int mi, int nf, int lf, int mf, int L, int S, rArray Ei)
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

rArray integral_cross_section(int ni, int li, int mi, int nf, int lf, int mf, int L, int S, rArray Ei)
{
	// the cross section
	rArray sigma(Ei.size());
	
	// get maximal partial wave angular momentum
	int max_ell;
	sqlitepp::statement st1(db);
	
	st1 << "SELECT MAX(ell) FROM hex "
			"WHERE ni = :ni AND li = :li AND mi = :mi AND nf = :nf AND lf = :lf AND mf = :mf AND L = :L AND S = :S",
		sqlitepp::into(max_ell),
		sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
		sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
		sqlitepp::use(L), sqlitepp::use(S);
	
	if (not st1.exec())
	{
		printf ("no data for these quantum numbers\n");
		return sigma;
	}
	
	// for all outgoing partial waves
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
		
		if (db_Ei.size() == 0)
			continue;
		
		// for all energies
		for (size_t i = 0; i < Ei.size(); i++)
		{
			// check conservation of energy
			if (Ei[i] <= 0. or Ei[i] - 1./(ni*ni) + 1./(nf*nf) <= 0.)
				continue;
			
			// compute projectile momenta
			double ki = sqrt(Ei[i]);
			double kf = sqrt(Ei[i] - 1./(ni*ni) + 1./(nf*nf));
			
			// compute "sigma" from T-matrices
			auto compute_sigma_from_T = [ S, kf, ki ](Complex T) -> double {
				return 1./(2.*M_PI*2.*M_PI) * (2*S + 1) / 4. * (T.real() * T.real() + T.imag() * T.imag()) * kf / ki;
			};
			
			// update values of "sigma"
			sigma[i] += interpolate(db_Ei, db_T_ell, Ei[i], compute_sigma_from_T);
		}
	}
	
	return sigma;
}

rArray complete_cross_section(int ni, int li, int mi, int nf, int lf, int mf, rArray Ei)
{
	// vector size
	size_t Nenergy = Ei.size();
	
	// the cross section
	rArray sigma(Nenergy);
	
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
		rArray singlet = integral_cross_section(ni, li, mi, nf, lf, mf, L, 0, Ei);
		rArray triplet = integral_cross_section(ni, li, mi, nf, lf, mf, L, 1, Ei);
		
		sigma = sigma + singlet + triplet;
	}
	
	// return the result
	return sigma;
}

rArray extrapolate_cross_section(int ni, int li, int mi, int nf, int lf, int mf, rArray Ei)
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
		rArray singlet = integral_cross_section(ni, li, mi, nf, lf, mf, L, 0, Ei);
		rArray triplet = integral_cross_section(ni, li, mi, nf, lf, mf, L, 1, Ei);
		
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

rArray aitkenD2_cross_section(int ni, int li, int mi, int nf, int lf, int mf, rArray Ei)
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
		return complete_cross_section(ni,li,mi,nf,lf,mf,Ei);
	}
	
	// load all cross sections
	rArray singlet, triplet;
	rArrays sigmas;
	for (int L = 0; L <= max_L; L++)
	{
		// get spin-dependent cross sections
		singlet = integral_cross_section(ni, li, mi, nf, lf, mf, L, 0, Ei);
		triplet = integral_cross_section(ni, li, mi, nf, lf, mf, L, 1, Ei);
		
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

rArray collision_strength(int ni, int li, int mi, int nf, int lf, int mf, int L, int S, rArray Ei)
{
	rArray omegas = integral_cross_section(ni, li, mi, nf, lf, mf, L, S, Ei);
	omegas = omegas * Ei * ((2.*L+1.) * (2.*S+1.));
	return omegas;
}

rArray total_cross_section(int ni, int li, int mi, rArray Ei)
{
	// the cross section
	rArray sigma(Ei.size());
	
	// get maximal "nf" for this initial state
	int max_nf;
	sqlitepp::statement st1(db);
	st1 << "SELECT MAX(nf) FROM hex "
			"WHERE ni = :ni AND li = :li AND mi = :mi",
		sqlitepp::into(max_nf),
		sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi);
	if (not st1.exec())
	{
		printf ("no data for these quantum numbers\n");
		return sigma;
	}
	
	// for all allowed excitations
	for (int nf = 1; nf <= max_nf; nf++)
		for (int lf = 0; lf < nf; lf++)
			for (int mf = -lf; mf <= lf; mf++)
			{
				rArray dsigma = complete_cross_section(ni, li, mi, nf, lf, mf, Ei);
				for (size_t i = 0; i < Ei.size(); i++)
					sigma[i] += dsigma[i];
			}
	
	return sigma;
}

// -------------------------------------------------------------------------- //
// One-value interfaces                                                       //
// -------------------------------------------------------------------------- //

Complex scattering_amplitude(int ni, int li, int mi, int nf, int lf, int mf, int S, double Ei, double theta)
{
	rArray thetas = { theta };
	cArray amplitudes = scattering_amplitude(ni,li,mi, nf,lf,mf, S, Ei, thetas);
	return amplitudes[0];
}

double differential_cross_section(int ni, int li, int mi, int nf, int lf, int mf, int S, double Ei, double theta)
{
	// compute the scattering amplitude
	Complex f = scattering_amplitude(ni, li, mi, nf, lf, mf, S, Ei, theta);
	
	// compute projectile momenta
	double ki = sqrt(Ei);
	double kf = sqrt(Ei - 1./(ni*ni) + 1./(nf*nf));
	
	// return differential cross section
	return kf / ki * (2. * S + 1.) / 4. * (f.real() * f.real() + f.imag() * f.imag());
}

// -------------------------------------------------------------------------- //
// Fortran interfaces                                                         //
// -------------------------------------------------------------------------- //

void scattering_amplitude(
	int const * ni, int const * li, int const * mi,
	int const * nf, int const * lf, int const * mf,
	int const * S, double const * Ei, int const * N,
	double const * theta, Complex * amplitudes
) {
	// save angles to a C++ structure
	rArray vec_theta(*N);
	for (int i = 0; i < *N; i++)
		vec_theta[i] = theta[i];
	
	// compute
	cArray vec_amplitudes = scattering_amplitude(
		*ni,*li,*mi, 
		*nf,*lf,*mf, 
		*S, *Ei, vec_theta
	);
	
	// copy amplitudes from C++ structure
	for (int i = 0; i < *N; i++)
		amplitudes[i] = vec_amplitudes[i];
}

void differential_cross_section(
	int const * ni, int const * li, int const * mi, 
	int const * nf, int const * lf, int const * mf, 
	int const * S, double const * Ei, int const * N,
	double const * theta, double * dsigma
) {
	// save angles to a C++ structure
	rArray vec_theta(*N);
	for (int i = 0; i < *N; i++)
		vec_theta[i] = theta[i];
	
	// compute
	rArray vec_amplitudes = differential_cross_section(
		*ni,*li,*mi,
		*nf,*lf,*mf,
		*S, *Ei, vec_theta
	);
	
	// copy amplitudes from C++ structure
	for (int i = 0; i < *N; i++)
		dsigma[i] = vec_amplitudes[i];
}

void momentum_transfer(
	int const * ni, int const * li, int const * mi, 
	int const * nf, int const * lf, int const * mf, 
	int const * L, int const * S,
	int const * N, double const * Ei, double * eta
) {
	// copy energies to a C++ structure
	rArray vec_Ei(*N);
	for (int i = 0; i < *N; i++)
		vec_Ei[i] = Ei[i];
	
	rArray vec_eta = momentum_transfer(*ni, *li, *mi, *nf, *lf, *mf, *L, *S, vec_Ei);
	
	// copy cross sections from a C++ structure
	for (int i = 0; i < *N; i++)
		eta[i] = vec_eta[i];
}

void integral_cross_section(
	int const * ni, int const * li, int const * mi,
	int const * nf, int const * lf, int const * mf,
	int const * L, int const * S, int const * N, double const * Ei,
	double * sigma
) {
	// copy energies to a C++ structure
	rArray vec_Ei(*N);
	for (int i = 0; i < *N; i++)
		vec_Ei[i] = Ei[i];
	
	rArray vec_sigma = integral_cross_section(*ni,*li,*mi, *nf,*lf,*mf, *L,*S, vec_Ei);
	
	// copy cross sections from a C++ structure
	for (int i = 0; i < *N; i++)
		sigma[i] = vec_sigma[i];
}

void complete_cross_section(
	int const * ni, int const * li, int const * mi,
	int const * nf, int const * lf, int const * mf,
	int const * N, double const * Ei, double * sigma
) {
	// copy energies to a C++ structure
	rArray vec_Ei(*N);
	for (int i = 0; i < *N; i++)
		vec_Ei[i] = Ei[i];
	
	rArray vec_sigma = complete_cross_section(*ni,*li,*mi, *nf,*lf,*mf, vec_Ei);
	
	// copy cross sections from a C++ structure
	for (int i = 0; i < *N; i++)
		sigma[i] = vec_sigma[i];
}

void collision_strength(
	int const * ni, int const * li, int const * mi,
	int const * nf, int const * lf, int const * mf,
	int const * L, int const * S,
	int const * N, double const * Ei, double * omega
) {
	// copy energies to a C++ structure
	rArray vec_Ei(*N);
	for (int i = 0; i < *N; i++)
		vec_Ei[i] = Ei[i];
	
	rArray vec_omega = collision_strength(*ni,*li,*mi, *nf,*lf,*mf, *L,*S,vec_Ei);
	
	// copy cross sections from a C++ structure
	for (int i = 0; i < *N; i++)
		omega[i] = vec_omega[i];
}

void total_cross_section(
	int const * ni, int const * li, int const * mi,
	int const * N, double const * Ei, double * sigma
) {
	// copy energies to a C++ structure
	rArray vec_Ei(*N);
	for (int i = 0; i < *N; i++)
		vec_Ei[i] = Ei[i];
	
	rArray vec_sigma = total_cross_section(*ni,*li,*mi, vec_Ei);
	
	// copy cross sections from a C++ structure
	for (int i = 0; i < *N; i++)
		sigma[i] = vec_sigma[i];
}
