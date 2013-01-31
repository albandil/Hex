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
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>

#include <sqlite3.h>
#include <sqlitepp/sqlitepp.hpp>

#include "angs.h"
#include "arrays.h"
#include "complex.h"
#include "hex-db.h"
#include "specf.h"
#include "interpolate.h"
#include "variables.h"

sqlitepp::session db;	// database handle
VariableList vlist;		// list of scattering variables

void db_sqrt(sqlite3_context* pdb, int n, sqlite3_value** val)
{
	sqlite3_result_double(pdb, sqrt(sqlite3_value_double(*val)));
}

void initialize(const char* dbname)
{
	// open database
	db.open(dbname);
	
	// define SQRT function
	sqlite3_create_function (
		db.impl(),
		"sqrt",
		1,
		SQLITE_UTF8,
		nullptr,
		&db_sqrt,
		nullptr,
		nullptr
	);
}

void create_new_database()
{
	// create tables
	for (const Variable* var : vlist)
	{
		sqlitepp::statement st(db);
		std::string cmd = var->SQL_CreateTable();
		
		if (cmd.size() == 0)
			continue;
		
		st << cmd;
		
		try {
			
			st.exec();
			
		} catch (sqlitepp::exception & e) {
			
			std::cerr << "ERROR: Creation of tables failed, code = " << e.code() << " (\"" << e.what() << "\")" << std::endl;
			std::cerr << "       Failed SQL command was: \"" << cmd << "\"" << std::endl;
			exit(-1);
			
		}
	}
}

void import (const char* sqlname)
{
	std::ifstream ifs(sqlname);
	
	if (not ifs.good())
	{
		std::cerr << "ERROR: Cannot open file \"" << sqlname << "\"" << std::endl;
		exit(-1);
	}
	
	// NOTE assuming one statement per line
	do
	{
		sqlitepp::statement st(db);
		
		std::string cmd;
		getline(ifs, cmd);
		
		if (cmd.size() == 0)
			continue;
		
		st << cmd;
		
		try {
		
			st.exec();
			
		} catch (sqlitepp::exception & e) {
		
			std::cerr << "ERROR: Import failed, code = " << e.code() << " (\"" << e.what() << "\")" << std::endl;
			std::cerr << "       Failed SQL command was: \"" << cmd << "\"" << std::endl;
			exit(-1);
			
		}
		
	} while (not ifs.eof());
}

void update ()
{
	for (const Variable* var : vlist)
	{
		sqlitepp::statement st(db);
	
		std::string cmd = var->SQL_Update();
	
		if (cmd.size() == 0)
			continue;
	
		st << cmd;
		
		try {
			
			st.exec();
			
		} catch (sqlitepp::exception & e) {
			
			std::cerr << "ERROR: Update failed, code = " << e.code() << " (\"" << e.what() << "\")" << std::endl;
			std::cerr << "       Failed SQL command was: \"" << cmd << "\"" << std::endl;
			exit(-1);
			
		}
	}
}

int run (
	std::vector<std::string> const & vars,
	std::map<std::string,std::string> const & sdata,
	rArray const & nums)
{
	// for all requested variables (mostly there will be just one)
	for (std::string const & varname : vars)
	{
		// pointer to the correct "Variable" object
		Variable const * var;
		
		// get a pointer from the dictionary
		if ((var = vlist.get(varname)) == nullptr)
		{
			// this should never happen
			fprintf(stderr, "[run] Runtime error.\n");
			return EXIT_FAILURE;
		}
		
		// try to compute the results
		if (not var->run(db, sdata, nums))
		{
			// this can easily happen
			std::cerr << "Computation of \"" << varname << "\" failed." << std::endl;
			return EXIT_FAILURE;
		}
	}
	
	return EXIT_SUCCESS;
}


// -------------------------------------------------------------------------- //
// One-value interfaces                                                       //
// -------------------------------------------------------------------------- //
// 
// Complex scattering_amplitude(int ni, int li, int mi, int nf, int lf, int mf, int S, double Ei, double theta)
// {
// 	rArray thetas = { theta };
// 	cArray amplitudes = ScatteringAmplitude::compute(ni,li,mi, nf,lf,mf, S, Ei, thetas);
// 	return amplitudes[0];
// }
// 
// double differential_cross_section(int ni, int li, int mi, int nf, int lf, int mf, int S, double Ei, double theta)
// {
// 	// compute the scattering amplitude
// 	Complex f = scattering_amplitude(ni, li, mi, nf, lf, mf, S, Ei, theta);
// 	
// 	// compute projectile momenta
// 	double ki = sqrt(Ei);
// 	double kf = sqrt(Ei - 1./(ni*ni) + 1./(nf*nf));
// 	
// 	// return differential cross section
// 	return kf / ki * (2. * S + 1.) / 4. * (f.real() * f.real() + f.imag() * f.imag());
// }

// -------------------------------------------------------------------------- //
// Fortran interfaces                                                         //
// -------------------------------------------------------------------------- //
/*
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
	cArray vec_amplitudes = ScatteringAmplitude::compute (
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
	rArray vec_amplitudes = differential_cross_section (
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
	
	rArray vec_eta = MomentumTransfer::compute(*ni, *li, *mi, *nf, *lf, *mf, *L, *S, vec_Ei);
	
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
	
	rArray vec_sigma = IntegralCrossSection::compute(*ni,*li,*mi, *nf,*lf,*mf, *L,*S, vec_Ei);
	
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
	
	rArray vec_sigma = CompleteCrossSection::compute(*ni,*li,*mi, *nf,*lf,*mf, vec_Ei);
	
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
	
	rArray vec_omega = CollisionStrength::compute(*ni,*li,*mi, *nf,*lf,*mf, *L,*S,vec_Ei);
	
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
	
	rArray vec_sigma = TotalCrossSection::compute(*ni,*li,*mi, vec_Ei);
	
	// copy cross sections from a C++ structure
	for (int i = 0; i < *N; i++)
		sigma[i] = vec_sigma[i];
}*/
