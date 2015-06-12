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
#include <cstring>
#include <iostream>

#include "cmdline.h"
#include "hex-db.h"
#include "variables.h"
#include "version.h"

typedef std::vector<std::string> const & Args;

const std::string HelpText = 
    "Usage:\n"
    "\thex-db [options]\n"
    "\n"
    "Options:\n"
    "\t--help             Display this information.\n"
    "\t--new              Create new database.\n"
    "\t--database <name>  Use database <name> instead of the default \"hex.db\".\n"
    "\t--import <file>    Import SQL batch file.\n"
    "\t--dump <file>      Export all contents as a SQL batch file.\n"
    "\t--update           Update derived quantities e.g. after manual import using \"sqlite3\".\n"
    "\t--optimize         Execute VACUUM command, i.e. minimize the occupied space.\n"
    "\t--<var>            Compute/retrieve scattering variable <var>.\n"
    "\t--vars             Display all available scattering variables.\n"
    "\t--<param> <val>    Set scattering parameter <param> to the value <val>.\n"
    "\t--params <var>     Display all available scattering parameters for variable <var>.\n"
    "\t--Eunits <Eunits>  Set units for energy: Ry, a.u. or eV (default: Ry).\n"
    "\t--Tunits <Tunits>  Set units for output: a.u. or cgs (default: a.u.).\n"
    "\t--Aunits <Aunits>  Set units for angles: deg or rad (default: deg).\n"
    "\n"
    "Usage examples:\n"
    "\n"
    "   # retrieve scattering amplitude\n"
    "   > seq 0 180             | hex-db --database=\"hex.db\" --scatamp --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --S=0 --E=0.75\n"
    "\n"
    "   # retrieve differential cross section:\n"
    "   > seq 0 180             | hex-db --database=\"hex.db\" --dcs     --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --S=0 --E=0.75\n"
    "\n"
    "   # retrieve momentum transfer:\n"
    "   > seq 0.650 0.001 0.850 | hex-db --database=\"hex.db\" --momtf   --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --L=0 --S=0\n"
    "\n"
    "   # retrieve integral cross section:\n"
    "   > seq 0.650 0.001 0.850 | hex-db --database=\"hex.db\" --ics     --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --L=0 --S=0\n"
    "\n"
    "   # retrieve sum integral cross section:\n"
    "   > seq 0.650 0.001 0.850 | hex-db --database=\"hex.db\" --ccs     --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0\n"
    "\n"
    "   # retrieve collision strength:\n"
    "   > seq 0.650 0.001 0.850 | hex-db --database=\"hex.db\" --colls   --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --L=0 --S=0\n"
    "\n"
    "   #retrieve total cross section:\n"
    "   > seq 0.650 0.001 0.850 | hex-db --database=\"hex.db\" --tcs     --ni=1 --li=0 --mi=0\n"
    "\n"
    "\nOn some locales a decimal comma is used instead of decimal point (e.g. on native Czech systems).\n"
    "This is not compatible with raw C/C++, one has to turn off the locale by executing\n"
    "> setenv LC_NUMERIC en_GB.utf8\n"
    "or\n"
    "> export LC_NUMERIC=en_GB.utf8\n"
    "depending on version of your shell.\n"
    "\n";

void print_help_and_exit()
{
    std::cout << logo(" ") << std::endl;
    std::cout << "=== Database interface program ===" << std::endl;
    std::cout << std::endl << HelpText;
    std::exit(EXIT_SUCCESS);
}

int main (int argc, char* argv[])
{
    // if no parameters are given, print usage and return
    if (argc == 1)
    {
        print_help_and_exit();
    }
    
    // program parameters
    std::string dbname = "hex.db", dumpfile;
    std::vector<std::string> sqlfiles;
    bool create_new = false, doupdate = false, doimport = false;
    bool dooptimize = false;
    
    // scattering variables to compute
    std::vector<std::string> vars;
    
    // used scattering event data
    std::map<std::string,std::string> sdata;
    
    // for all argv-s
    ParseCommandLine
    (
        argc, argv,
        "help", "h",     0, [ & ](Args opts) -> bool { print_help_and_exit(); return true; },
        "new", "n",      0, [ & ](Args opts) -> bool { create_new = true; return true; },
        "database", "D", 1, [ & ](Args opts) -> bool { dbname = opts[0]; return true; },
        "update", "u",   0, [ & ](Args opts) -> bool { doupdate = true; return true; },
        "optimize", "o", 0, [ & ](Args opts) -> bool { dooptimize = true; return true; },
        "dump", "d",     1, [ & ](Args opts) -> bool { dumpfile = opts[0]; return true; },
        "import", "i",  -1, [ & ](Args opts) -> bool
        {
            doimport = true;
            sqlfiles.insert(sqlfiles.end(), opts.begin(), opts.end());
            return true;
        },
        "vars", "v",     0, [ & ](Args opts) -> bool
        {
            OutputTable table;
            table.setAlignment(OutputTable::left, OutputTable::left);
            table.setWidth(20,40);
            table.write("(name)", "(description)");
            for (Variable const * var : vlist)
                table.write(var->id(), var->description());
            
            std::exit(EXIT_SUCCESS);
        },
        "params", "p",   1, [ & ](Args opts) -> bool
        {
            VariableList::const_iterator it = std::find_if
            (
                vlist.begin(),
                vlist.end(),
                [ & ](Variable* const & it) -> bool { return it->id() == opts[0]; }
            );
            
            if (it == vlist.end())
            {
                std::cout << "No such variable \"" << opts[0] << "\"" << std::endl;
            }
            else
            {
                std::cout << std::endl;
                std::cout << "Variable \"" << opts[0] << "\" uses the following parameters:" << std::endl;
                std::cout << std::endl;
                
                OutputTable table;
                table.setAlignment(OutputTable::left, OutputTable::left);
                table.setWidth(20,40);
                for (auto v : (*it)->deps())
                    table.write(v.first, v.second);
                
                std::cout << std::endl;
                
                std::cout << "of which one of the following can be read from the standard input:" << std::endl;;
                std::cout << std::endl;
                
                for (std::string v : (*it)->vdeps())
                    std::cout << v << " ";
                
                std::cout << std::endl << std::endl;
            }
            std::exit(EXIT_SUCCESS);
        },
        "Eunits", "", 1, [ & ](Args opts) -> bool
        { 
            if (opts[0] == std::string("Ry"))
                Eunits = eUnit_Ry;
            else if (opts[0] == std::string("a.u."))
                Eunits = eUnit_au;
            else if (opts[0] == std::string("eV"))
                Eunits = eUnit_eV;
            else
                HexException("Unknown units \"%s\"", opts[0].c_str());
            
            return true;
        },
        "Tunits", "", 1, [ & ](Args opts) -> bool
        {
            if (opts[0] == std::string("a.u."))
                Lunits = lUnit_au;
            else if (opts[0] == std::string("cgs"))
                Lunits = lUnit_cgs;
            else
                HexException("Unknown units \"%s\"", opts[0].c_str());
            
            return true;
        },
        "Aunits", "", 1, [ & ](Args opts) -> bool
        {
            if (opts[0] == std::string("deg"))
                Aunits = aUnit_deg;
            else if (opts[0] == std::string("rad"))
                Aunits = aUnit_rad;
            else
                HexException("Unknown units \"%s\"", opts[0].c_str());
            
            return true;
        },
        /* default*/ [ & ](std::string arg, Args opts) -> bool
        {
            // try to find it in the variable ids
            for (const Variable* var : vlist)
            {
                if (var->id() == arg)
                {
                    // insert this id into ToDo list
                    vars.push_back(arg);
                    return true;
                }
            }
            
            // try to find it in the variable dependencies
            for (const Variable* var : vlist)
            {
                // scan the dependencies for 'arg'
                bool this_arg_is_needed = false;
                for (auto iter = var->deps().begin(); iter != var->deps().end(); iter++)
                {
                    if (iter->first == arg)
                        this_arg_is_needed = true;
                }
                
                // if this parameter has been found, process it furher
                if (this_arg_is_needed)
                {
                    // check parameter
                    if (opts.empty())
                    {
                        std::cerr << "ERROR: The option --" << arg << " requires a parameter!" << std::endl;
                        return false;
                    }
                    
                    // insert into scattering data list
                    sdata[arg] = std::string(opts[0]);
                    return true;
                }
            }
            
            return false;
        }
    );
    
    // is the requested filename is available?
    std::ifstream f (dbname.c_str());
    bool exists = f.good();
    f.close();
    
    if (create_new and exists)
    {
        std::cerr << "The file \"" << dbname << "\" already exists. Please eiher delete it or use another name for the new database." << std::endl;
        std::exit(EXIT_FAILURE);
    }
    
    if (not create_new and not exists)
    {
        std::cerr << "The file \"" << dbname << "\" does not exist. Please create a new database first or use existing one." << std::endl;
        std::exit(EXIT_FAILURE);
    }
    
    // open the database, if it is going to be used
    if (create_new or doimport or doupdate or dooptimize or not dumpfile.empty() or not vars.empty())
        hex_initialize(dbname.c_str());
    
    // create new database if asked to
    if (create_new)
        hex_new();
    
    // import SQL data
    if (doimport)
    {
        for (std::string const & sqlfile : sqlfiles)
            hex_import(sqlfile.c_str());
    }
    
    // update
    if (doupdate)
        hex_update();
    
    // vacuum
    if (dooptimize)
        hex_optimize();
    
    // dump contents of "tmat" table
    if (not dumpfile.empty())
        hex_dump(dumpfile.c_str());
    
    // is there anything more to do?
    if (vars.empty())
        return 0;
    
    // retrieve all requested data
    return hex_run (vars, sdata);
}
