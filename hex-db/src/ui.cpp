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
#include <cstring>
#include <iostream>
#include <getopt.h>

#include "hex-db.h"
#include "variables.h"

std::string HelpText = 
	"Usage: hex-db [options]\n"
	"Options:\n"
	"  --help            Display this information.\n"
	"  --new             Create new database.\n"
	"  --database=<name> Use database <name> instead of default \"hex.db\".\n"
	"  --import=<file>   Import SQL batch file.\n"
	"  --update          Update derived quantities e.g. after manual import using \"sqlite3\".\n"
	"  --<var>           Compute/retrieve scattering variable <var>.\n"
	"  --vars            Display all available scattering variables.\n"
	"  --<param>=<val>   Set scattering parameter <param> to the value <val>.\n"
	"  --params=<var>    Display all available scattering parameters for variable <var>.\n"
	"  --Eunits=<Eunits> Set units for energy: Ry, a.u. or eV.\n"
	"  --Tunits=<Tunits> Set units for output: a.u. or cgs.\n"
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

bool SwitchInput (std::string argnam, std::string argpar)
{
	// we can get here only if there is no default callback
	std::cerr << "ERROR: Uknown option --" << argnam << std::endl;
	return false;
}
	
template <typename DefaultCallback>
bool SwitchInput (std::string argnam, std::string argpar, DefaultCallback defaultcallback)
{
	return defaultcallback(argnam, argpar) or SwitchInput(argnam, argpar);
}

template <typename Callback, typename ...Params>
bool SwitchInput (std::string argnam, std::string argpar, std::string optname, int npar, Callback callback, Params ...p)
{
	// check optname
	if (argnam != optname)
		return SwitchInput(argnam, argpar, p...); // parse the rest
	
	// check parameter
	if (npar > 0 and argpar.empty())
	{
		std::cerr << "ERROR: The option --" << argnam << " requires a parameter!" << std::endl;
		return false;
	}
	if (npar == 0 and not argpar.empty())
	{
		std::cerr << "ERROR: The option --" << argnam << " uses no parameter!" << std::endl;
		return false;
	}
	
	// run the callback
	callback(argpar);
	
	// successfully done
	return true;
}

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
	
	// scattering variables to compute
	std::vector<std::string> vars;
	
	// used scattering event data
	std::map<std::string,std::string> sdata;
	
	// units
	eUnit Eunits = eUnit_Ry;
	lUnit Lunits = lUnit_au;
	
	// for all argv-s
	for (int i = 1; i < argc; i++)
	{
		std::string arg = argv[i];
		
		// remove leading dashes
		while (arg.front() == '-')
			arg.erase(0,1);
		
		// split at equation sign (if present)
		std::string::iterator eq;
		eq = std::find(arg.begin(), arg.end(), '=');
		std::string argnam = arg.substr(0, eq - arg.begin());
		std::string argpar;
		if (eq != arg.end()) argpar = arg.substr(eq + 1 - arg.begin(), arg.end() - eq - 1);
		
		// if there is no equation sign, check if the next argument is parameter
		if (eq == arg.end())
		{
			if (i+1 < argc and (strlen(argv[i+1]) == 1 or argv[i+1][0] != '-'))
			{
				argpar = argv[i+1];
				i++;
			}
		}
		
		// switch option name
		bool OK = SwitchInput ( argnam, argpar,
			"help",     0, [ & ](std::string opt) -> void { std::cout << HelpText; },
			"new",      0, [ & ](std::string opt) -> void { create_new = true; },
			"database", 1, [ & ](std::string opt) -> void { dbname = opt; },
			"import",   1, [ & ](std::string opt) -> void { doimport = true; sqlfile = opt; },
			"update",   0, [ & ](std::string opt) -> void { doupdate = true; },
			"dump",     1, [ & ](std::string opt) -> void { dumpfile = opt; },
			"vars",     0, [ & ](std::string opt) -> void {
				std::cout << "(name)\t\t(description)" << std::endl;
				for (Variable const * var : vlist)
					std::cout << var->id() << "\t\t" << var->description() << std::endl;
			},
			"params",   1, [ & ](std::string opt) -> void {
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
					std::cout << "Variable \"" << opt << "\" uses the following parameters:\n\t";
					bool theta_present = false;
					bool Ei_present = false;
					for (std::string param : (*it)->dependencies())
					{
						std::cout << param << " ";
						if (param == std::string("theta")) theta_present = true;
						if (param == std::string("Ei")) Ei_present = true;
					}
					if (theta_present)
						std::cout << "\n\"theta\" is to be specified using STDIN.\n";
					else if (Ei_present)
						std::cout << "\n\"Ei\" is to be specified using STDIN.\n";
				}
			},
			"Eunits",   1, [ & ](std::string opt) -> void { 
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
			},
			"Tunits",   1, [ & ](std::string opt) -> void {
				if (opt == std::string("a.u."))
					Lunits = lUnit_au;
				else if (opt == std::string("cgs"))
					Lunits = lUnit_cgs;
				else
				{
					std::cerr << "Unknown units \"" << opt << "\"" << std::endl;
					abort();
				}
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
					std::vector<std::string> const & deps = var->dependencies();
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
				
				// unable to match the option
				return false;
			}
		);
		
		// chech result
		if (not OK)
			return EXIT_FAILURE;
	}
	
	// open the database
	initialize(dbname.c_str());
	
	// create new database if asked to
	if (create_new)
		create_new_database();
	
	// import SQL data
	if (doimport)
		import(sqlfile.c_str());
	
	// update
	if (doimport or doupdate)
		update();
	
	// dump contents of "tmat" table
	if (not dumpfile.empty())
		dump(dumpfile.c_str());
	
	// is there anything more to do?
	if (vars.empty())
		return 0;
	
	// read all standard input
	char str[100];
	rArray nums;
	while (gets(str) != 0)
	{
		// manage Czech locale
		for (char& c : str) if (c == ',') c = '.';
		
		// convert to number
		char* tail;
		double num = strtod(str, &tail);
		if (*tail == 0)
			nums.push_back(num);
		else
			fprintf(stderr, "The input \"%s\" is not a valid number.", str);
	}
	
	// retrieve all requested data
	return run (Eunits, Lunits, vars, sdata, nums);
}
