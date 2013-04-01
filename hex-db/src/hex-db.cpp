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

void dump (const char* dumpfile)
{
	sqlitepp::statement st(db);
	std::string dumpline, dumpcmd =
		"SELECT 'INSERT INTO " + TMatrix::Id + " VALUES(' || "
		"quote(ni) || ',' || "
		"quote(li) || ',' || "
		"quote(mi) || ',' || "
		"quote(nf) || ',' || "
		"quote(lf) || ',' || "
		"quote(mf) || ',' || "
		"quote(L) || ',' || "
		"quote(S) || ',' || "
		"quote(Ei) || ',' || "
		"quote(ell) || ',' || "
		"quote(Re_t_ell) || ',' || "
		"quote(Im_T_ell) || ')' FROM " + TMatrix::Id + ";";
	
	try {
		
		// create statement
		st << dumpcmd, sqlitepp::into(dumpline);
		
		// open output file
		std::ofstream outfile;
		if (dumpfile != std::string("-"))
		{
			outfile.open(dumpfile);
			if (outfile.bad())
			{
				std::cerr << "Couldn't open file \"" << dumpfile << "\".\n";
				exit(-1);
			}
		}
		
		// get and write data
		if (dumpfile != std::string("-"))
		{
			outfile << "BEGIN TRANSACTION" << std::endl;
			while (st.exec())
				outfile << dumpline << std::endl;
			outfile << "COMMIT" << std::endl;
		}
		else
		{
			std::cout << "BEGIN TRANSACTION" << std::endl;
			while (st.exec())
				std::cout << dumpline << std::endl;
			std::cout << "COMMIT" << std::endl;
		}
		
	} catch (sqlitepp::exception & e) {
		
		std::cerr << "ERROR: Dump failed, code = " << e.code() << " (\"" << e.what() << "\")" << std::endl;
		std::cerr << "       Failed SQL command was: \"" << dumpcmd << "\"" << std::endl;
		exit(-1);
		
	}
}

int run (
	eUnit Eunits, lUnit Lunits,
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
		if (not var->run(Eunits, Lunits, db, sdata, nums))
		{
			// this can easily happen
			std::cerr << "Computation of \"" << varname << "\" failed." << std::endl;
			return EXIT_FAILURE;
		}
	}
	
	return EXIT_SUCCESS;
}
