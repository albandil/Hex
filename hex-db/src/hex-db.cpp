//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2015, Jakub Benda, Charles University in Prague                    //
//                                                                                   //
// MIT License:                                                                      //
//                                                                                   //
//  Permission is hereby granted, free of charge, to any person obtaining a          //
// copy of this software and associated documentation files (the "Software"),        //
// to deal in the Software without restriction, including without limitation         //
// the rights to use, copy, modify, merge, publish, distribute, sublicense,          //
// and/or sell copies of the Software, and to permit persons to whom the             //
// Software is furnished to do so, subject to the following conditions:              //
//                                                                                   //
//  The above copyright notice and this permission notice shall be included          //
// in all copies or substantial portions of the Software.                            //
//                                                                                   //
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS          //
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       //
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE       //
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, //
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF         //
// OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.  //
//                                                                                   //
//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //

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
#include "hex-db.h"
#include "special.h"
#include "interpolate.h"
#include "variables.h"

sqlitepp::session db;    // database handle
VariableList vlist;        // list of scattering variables

// units
eUnit Eunits = eUnit_Ry;
lUnit Lunits = lUnit_au;
aUnit Aunits = aUnit_deg;

void hex_initialize (const char* dbname) { hex_initialize_(dbname); }
void hex_initialize_ (const char* dbname)
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

void hex_new () { hex_new_(); }
void hex_new_ ()
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
            std::exit(EXIT_FAILURE);
            
        }
    }
}

void hex_import (const char* sqlname) { hex_import_(sqlname); }
void hex_import_ (const char* sqlname)
{
    // open file, if necessary
    std::ifstream ifs;
    if (std::string(sqlname) != std::string("-"))
        ifs.open(sqlname);
    
    // get input stream
    std::istream& is = (std::string(sqlname) != std::string("-")) ? ifs : std::cin;
    
    // test input stream
    if (not is.good())
        HexException("Cannot open file \"%s\".", sqlname);
    
    // line numbers (current and total)
    unsigned line = 0, lines = 0;
    if (std::string(sqlname) != std::string("-"))
    {
        std::cout << "Counting lines..." << std::flush;
        lines = std::count
        (
            std::istreambuf_iterator<char>(is),
            std::istreambuf_iterator<char>(),
            '\n'
        );
        
        // reset file to the beginning
        is.clear();
        is.seekg(0);
    }
    
    std::cout << "\rImporting data...  " << std::flush;
    
    do
    {
        // query statement
        sqlitepp::statement st(db);
        
        if (lines > 0)
            std::cout << "\rImporting data...  " << std::fixed << std::setprecision(0) << line * 100. / lines << " % " << std::flush;
        
        // read line from input stream
        std::string cmd, cmd1;
        std::getline (is, cmd);
        
        // skip empty lines
        if (cmd.size() == 0)
            continue;
        
        // for all statements on this line
        do
        {
            // trim spaces
            cmd = trim(cmd);
            
            // skip this line if it starts with a comment
            if (cmd.size() > 1 and cmd[0] == '-' and cmd[1] == '-')
                break;
            
            // position of semicolon
            unsigned semicolon = std::find(cmd.begin(), cmd.end(), ';') - cmd.begin();
            
            // split string
            cmd1 = cmd.substr(0, semicolon);
            cmd = (semicolon == cmd.size() ? "" : cmd.substr(semicolon + 1, std::string::npos));
            
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
                std::exit(EXIT_FAILURE);
                
            }
        }
        while (cmd.size() > 0);
        
        // move on to the next line
        line++;
    }
    while (not is.eof());
    
    std::cout << "\rThe SQL batch file \"" << sqlname << "\" has been successfully imported." << std::endl;
}

void hex_update () { hex_update_(); }
void hex_update_ ()
{
    for (const Variable* var : vlist)
    {
        std::cout << "\rUpdating " << var->id() << "...      " << std::flush;
        
        for (std::string const & cmd : var->SQL_Update())
        {
            sqlitepp::statement st(db);
            
            if (cmd.size() == 0)
                continue;
            
            st << cmd;
            
            try
            {
                st.exec();
            }
            catch (sqlitepp::exception & e)
            {
                std::cerr << "ERROR: Update failed, code = " << e.code() << " (\"" << e.what() << "\")" << std::endl;
                std::cerr << "       Failed SQL command was: \"" << cmd << "\"" << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
    }
    std::cout << "\rThe database has been successfully updated." << std::endl;
}

void hex_optimize () { hex_optimize_(); }
void hex_optimize_ ()
{
    sqlitepp::statement st(db);
    st << "VACUUM";
    try
    {
        st.exec();
    }
    catch (sqlitepp::exception & e)
    {
        HexException("Database optimizaion failed, code = %d (\"%s\").", e.code(), e.what());
    }
}

void hex_dump (const char* dumpfile) { hex_dump_(dumpfile); }
void hex_dump_ (const char* dumpfile)
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
        "quote(Re_T_ell) || ',' || "
        "quote(Im_T_ell) || ',' || "
        "quote(Re_TBorn_ell) || ',' || "
        "quote(Im_TBorn_ell) || ')' FROM " + TMatrix::Id + ";";
    
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
                std::exit(EXIT_FAILURE);
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
        std::exit(EXIT_FAILURE);
        
    }
}

int hex_run
(
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
            HexException("Runtime error.");
        }
        
        // try to compute the results
        if (not var->run(sdata))
        {
            // this can easily happen
            HexException("Computation of \"%s\" failed.", varname.c_str());
        }
    }
    
    return EXIT_SUCCESS;
}
