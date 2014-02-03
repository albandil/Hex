/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2014                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <algorithm>
#include <cstring>
#include <iostream>

#include "cmdline.h"
#include "hex-db.h"
#include "variables.h"

const std::string HelpText = 
    "Usage: hex-db [options]\n"
    "Options:\n"
    "  --help            Display this information.\n"
    "  --new             Create new database.\n"
    "  --database <name> Use database <name> instead of the default \"hex.db\".\n"
    "  --import <file>   Import SQL batch file.\n"
    "  --dump <file>     Export all contents as a SQL batch file.\n"
    "  --update          Update derived quantities e.g. after manual import using \"sqlite3\".\n"
    "  --optimize        Execute VACUUM command, i.e. minimize the occupied space.\n"
    "  --avail           Print out available data.\n"
    "  --<var>           Compute/retrieve scattering variable <var>.\n"
    "  --vars            Display all available scattering variables.\n"
    "  --<param> <val>   Set scattering parameter <param> to the value <val>.\n"
    "  --params <var>    Display all available scattering parameters for variable <var>.\n"
    "  --Eunits <Eunits> Set units for energy: Ry, a.u. or eV.\n"
    "  --Tunits <Tunits> Set units for output: a.u. or cgs.\n"
    "  --Aunits <Aunits> Set units for angles: deg or rad.\n"
    "\n"
    "Usage examples:\n"
    "\n"
    "   # retrieve scattering amplitude\n"
    "   > seq 0.01 0.01 3.14    | hex-db --database=\"hex.db\" --scatamp    --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --S=0 --E=0.75\n"
    "\n"
    "   # retrieve differential cross section:\n"
    "   > seq 0.01 0.01 3.14    | hex-db --database=\"hex.db\" --dcs        --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --S=0 --E=0.75\n"
    "\n"
    "   # retrieve momentum transfer:\n"
    "   > seq 0.650 0.001 0.850 | hex-db --database=\"hex.db\" --momtransf  --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --L=0 --S=0\n"
    "\n"
    "   # retrieve integral cross section:\n"
    "   > seq 0.650 0.001 0.850 | hex-db --database=\"hex.db\" --integcs    --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --L=0 --S=0\n"
    "\n"
    "   # retrieve sum integral cross section:\n"
    "   > seq 0.650 0.001 0.850 | hex-db --database=\"hex.db\" --sumintegcs --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0\n"
    "\n"
    "   # retrieve collision strength:\n"
    "   > seq 0.650 0.001 0.850 | hex-db --database=\"hex.db\" --collstr    --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --L=0 --S=0\n"
    "\n"
    "   #retrieve total cross section:\n"
    "   > seq 0.650 0.001 0.850 | hex-db --database=\"hex.db\" --tcs        --ni=1 --li=0 --mi=0\n"
    "\n"
    "\nOn some locales a decimal comma is used instead of decimal point (e.g. on native Czech systems).\n"
    "This is not compatible with raw C/C++, one has to turn off the locale by executing\n"
    "> setenv LC_NUMERIC en_GB.utf8\n"
    "or\n"
    "> export LC_NUMERIC=en_GB.utf8\n"
    "depending on version of your shell.\n"
    "\n";


int main(int argc, char* argv[])
{
    // if no parameters are given, print usage and return
    if (argc == 1)
    {
        std::cout << HelpText;
        return EXIT_SUCCESS;
    }
    
    // program parameters
    std::string dbname = "hex.db", sqlfile, dumpfile;
    bool create_new = false, doupdate = false, doimport = false;
    bool dooptimize = false, doavail = false;
    
    // scattering variables to compute
    std::vector<std::string> vars;
    
    // used scattering event data
    std::map<std::string,std::string> sdata;
    
    // for all argv-s
    ParseCommandLine
    (
        argc, argv,
        "help", "h",     0, [ & ](std::string opt) -> bool { std::cout << HelpText; return true; },
        "new", "n",      0, [ & ](std::string opt) -> bool { create_new = true; return true; },
        "database", "D", 1, [ & ](std::string opt) -> bool { dbname = opt; return true; },
        "import", "i",   1, [ & ](std::string opt) -> bool { doimport = true; sqlfile = opt; return true; },
        "update", "u",   0, [ & ](std::string opt) -> bool { doupdate = true; return true; },
        "optimize", "o", 0, [ & ](std::string opt) -> bool { dooptimize = true; return true; },
        "avail", "a",    0, [ & ](std::string opt) -> bool { doavail = true; return true; },
        "dump", "d",     1, [ & ](std::string opt) -> bool { dumpfile = opt; return true; },
        "vars", "v",     0, [ & ](std::string opt) -> bool {
            std::cout << "(name)\t\t(description)" << std::endl;
            for (Variable const * var : vlist)
                std::cout << var->id() << "\t\t" << var->description() << std::endl;
            return true;
        },
        "params", "p",   1, [ & ](std::string opt) -> bool {
            VariableList::const_iterator it = std::find_if (
                vlist.begin(),
                vlist.end(),
                [ & ](Variable* const & it) -> bool { return it->id() == opt; }
            );
            
            if (it == vlist.end())
            {
                std::cout << "No such variable \"" << opt << "\"\n";
            }
            else
            {
                std::cout << "\nVariable \"" << opt << "\" uses the following parameters:\n\n\t";
                for (std::string v : (*it)->deps())
                    std::cout << v << " ";
                std::cout << "\n\nof that one of the following can be read from the standard input:\n\n\t";
                for (std::string v : (*it)->vdeps())
                    std::cout << v << " ";
                std::cout << "\n\n";
            }
            return true;
        },
        "Eunits", "", 1, [ & ](std::string opt) -> bool { 
            if (opt == std::string("Ry"))
                Eunits = eUnit_Ry;
            else if (opt == std::string("a.u."))
                Eunits = eUnit_au;
            else if (opt == std::string("eV"))
                Eunits = eUnit_eV;
            else
            {
                std::cerr << "Unknown units \"" << opt << "\"" << std::endl;
                abort();
            }
            return true;
        },
        "Tunits", "", 1, [ & ](std::string opt) -> bool {
            if (opt == std::string("a.u."))
                Lunits = lUnit_au;
            else if (opt == std::string("cgs"))
                Lunits = lUnit_cgs;
            else
            {
                std::cerr << "Unknown units \"" << opt << "\"" << std::endl;
                abort();
            }
            return true;
        },
        "Aunits", "", 1, [ & ](std::string opt) -> bool {
            if (opt == std::string("deg"))
                Aunits = aUnit_deg;
            else if (opt == std::string("rad"))
                Aunits = aUnit_rad;
            else
            {
                std::cerr << "Unknown units \"" << opt << "\"" << std::endl;
                abort();
            }
            return true;
        },
        /* default*/ [ & ](std::string arg, std::string par) -> bool {
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
                std::vector<std::string> const & deps = var->deps();
                if (std::find(deps.begin(), deps.end(), arg) != deps.end())
                {
                    // check parameter
                    if (par.empty())
                    {
                        std::cerr << "ERROR: The option --" << arg << " requires a parameter!" << std::endl;
                        return false;
                    }
                    
                    // insert into scattering data list
                    sdata[arg] = std::string(par);
                    return true;
                }
            }

            return false;
        }
    );
    
    // open the database
    initialize(dbname.c_str());
    
    // create new database if asked to
    if (create_new)
        create_new_database();
    
    // import SQL data
    if (doimport)
        import(sqlfile.c_str());
    
    // update
    if (doupdate)
        update();
    
    // vacuum
    if (dooptimize)
        optimize();
    
    // print out available data
    if (doavail)
        avail();
    
    // dump contents of "tmat" table
    if (not dumpfile.empty())
        dump(dumpfile.c_str());
    
    // is there anything more to do?
    if (vars.empty())
        return 0;
    
    // retrieve all requested data
    return run (vars, sdata);
}
