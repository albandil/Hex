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

#include <sqlitepp/sqlitepp.hpp>

#include "arrays.h"
#include "complex.h"
#include "hex-db.h"
#include "specf.h"
#include "interpolate.h"
#include "variables.h"

sqlitepp::session db;    // database handle
VariableList vlist;        // list of scattering variables

// units
eUnit Eunits = eUnit_Ry;
lUnit Lunits = lUnit_au;
aUnit Aunits = aUnit_deg;

void initialize(const char* dbname)
{
    // open database
    db.open(dbname);
    
    // disable journaling
    sqlitepp::statement st(db);
    st << "PRAGMA journal_mode = OFF";
    st.exec();
    
    // initialize variables
    for (const Variable* var : vlist)
        var->initialize(db);
}

void create_new_database()
{
    // create tables
    for (const Variable* var : vlist)
    for (std::string const & cmd : var->SQL_CreateTable())
    {
        sqlitepp::statement st(db);
        
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
    // open file, if necessary
    std::ifstream ifs;
    if (std::string(sqlname) != std::string("-"))
        ifs.open(sqlname);
    
    // get input stream
    std::istream& is = (std::string(sqlname) != std::string("-")) ? ifs : std::cin;
    
    // test input stream
    if (not is.good())
    {
        std::cerr << "ERROR: Cannot open file \"" << sqlname << "\"" << std::endl;
        exit(-1);
    }
    
    // input line number
    int line = 0;
    
    do {
        
        // query statement
        sqlitepp::statement st(db);
        
        // read line from input stream
        std::string cmd, cmd1;
        getline(is, cmd);
        
        // skip empty lines
        if (cmd.size() == 0)
            continue;
        
        // for all statements on this line
        do {
            
            // position of semicolon
            int semicolon = std::find(cmd.begin(), cmd.end(), ';') - cmd.begin();
            
            // split string
            cmd1 = cmd.substr(0, semicolon);
            cmd = cmd.substr(semicolon + 1, std::string::npos);
            
            // trim spaces
            cmd1 = trim(cmd1);
            cmd = trim(cmd);
            
            // try to execute the first command
            try {
                
                st << cmd1;
                st.exec();
                
            } catch (sqlitepp::exception & e) {
            
                std::cerr << "ERROR: Import failed, code = " << e.code() << " (\"" << e.what() << "\")" << std::endl;
                std::cerr << "       [Line " << line << "] Failed SQL command was: \"" << cmd1 << "\"" << std::endl;
                exit(-1);
                
            }
            
        } while (cmd.size() > 0);
        
        // move on to the next line
        line++;
        
    } while (not is.eof());
}

void update ()
{
    for (const Variable* var : vlist)
    for (std::string const & cmd : var->SQL_Update())
    {
        sqlitepp::statement st(db);
    
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

void optimize()
{
    sqlitepp::statement st(db);
    st << "VACUUM";
    try {
        
        st.exec();
        
    } catch (sqlitepp::exception & e) {
        
        std::cerr << "ERROR: Database optimizaion failed, code = " << e.code() << " (\"" << e.what() << "\")" << std::endl;
        exit(-1);
        
    }
}

void avail()
{
    std::cout << "\"--avail\" not implemented yet. Use the following code directly in sqlite3:\n";
    std::cout << "\n";
    std::cout << "CREATE TEMP TABLE MaxLs AS\n";
    std::cout << "    SELECT\n";
    std::cout << "        ni, li, mi, Ei, MAX(L) AS MaxL\n";
    std::cout << "    FROM\n";
    std::cout << "        tmat\n";
    std::cout << "    GROUP BY\n";
    std::cout << "        ni,li,mi,Ei\n";
    std::cout << "    ORDER BY\n";
    std::cout << "        Ei ASC;\n";
    std::cout << "\n";
    std::cout << "CREATE TEMP TABLE StackedLs AS\n";
    std::cout << "    SELECT\n";
    std::cout << "        LThis.ni, LThis.li, LThis.mi,\n";
    std::cout << "        LThis.Ei AS ThisE, LThis.MaxL AS ThisL\n";
    std::cout << "    FROM\n";
    std::cout << "        MaxLs AS LThis \n";
    std::cout << "    LEFT OUTER JOIN MaxLs AS LPrev\n";
    std::cout << "        ON LThis.rowid = LPrev.rowid + 1 \n";
    std::cout << "        AND LThis.ni = LPrev.ni \n";
    std::cout << "        AND LThis.li = LPrev.li \n";
    std::cout << "        AND LThis.mi = LPrev.mi\n";
    std::cout << "    LEFT OUTER JOIN MaxLs AS LNext\n";
    std::cout << "        ON LThis.rowid = LNext.rowid - 1 \n";
    std::cout << "        AND LThis.ni = LNext.ni \n";
    std::cout << "        AND LThis.li = LNext.li \n";
    std::cout << "        AND LThis.mi = LNext.mi\n";
    std::cout << "    WHERE\n";
    std::cout << "        LNext.MaxL IS NULL OR LPrev.MaxL IS NULL OR LPrev.MaxL <> LNext.MaxL\n";
    std::cout << "    ORDER BY \n";
    std::cout << "        LThis.ni, LThis.li, LThis.mi, ThisE;\n";
    std::cout << "\n";
    std::cout << "CREATE TEMP TABLE SStackedLs AS\n";
    std::cout << "    SELECT\n";
    std::cout << "        T1.ni, T1.li, T1.mi, T1.ThisE AS MinE, (CASE WHEN T2.ThisE IS NULL THEN T1.ThisE ELSE T2.ThisE END) AS MaxE, T1.ThisL AS MaxL\n";
    std::cout << "    FROM \n";
    std::cout << "        StackedLs AS T1\n";
    std::cout << "    LEFT OUTER JOIN StackedLs AS T2\n";
    std::cout << "        ON T1.rowid = T2.rowid - 1\n";
    std::cout << "        AND T1.ni = T2.ni\n";
    std::cout << "        AND T1.li = T2.li\n";
    std::cout << "        AND T1.mi = T2.mi\n";
    std::cout << "    WHERE T1.ThisL = T2.ThisL OR T2.ThisL IS NULL\n";
    std::cout << "\n";
    std::cout << "SELECT * FROM SStackedLs;\n";
    std::cout << "\n";
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
    std::vector<std::string> const & vars,
    std::map<std::string,std::string> const & sdata
)
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
            throw exception ("Runtime error.");
        }
        
        // try to compute the results
        if (not var->run(db, sdata))
        {
            // this can easily happen
            throw exception ("Computation of \"%s\" failed.", varname.c_str());
        }
    }
    
    return EXIT_SUCCESS;
}
