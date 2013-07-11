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

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <string>
#include <tuple>

#include <getopt.h>

#include "arrays.h"
#include "spmatrix.h"
#include "input.h"

/*
 * Load data from input file. Default filename is "hex.inp", but different
 * can be specified using the --input/-i option.
 */

void parse_command_line (
	int argc, char* argv[],
	std::ifstream & inputfile,
	std::string & zipfile, int & zipcount, double & zipmax,
	bool & parallel,
	int & itinerary
){
	// set short options
	const char* const short_options  = "eihznRmabc";
	
	// set long options
	const option long_options[] = {
		{"example",           0,   0, 'e'},
		{"input",             1,   0, 'i'},
		{"help",              0,   0, 'h'},
		{"zipfile",           1,   0, 'z'},
		{"zipcount",          1,   0, 'n'},
		{"zipmax",            1,   0, 'R'},
#ifndef NO_MPI
		{"mpi",               0,   0, 'm'},
#endif
		{"stg-integ",         0,   0, 'a'},
		{"stg-integ-solve",   0,   0, 'b'},
		{"stg-extract",       0,   0, 'c'},
		{0,                   0,   0,   0}
	};
	
	// switch all options
	int next_option = 0;
	do
	{
		next_option = getopt_long(argc, argv, short_options, long_options, 0);
		
		switch (next_option)
		{
			case 'e':
			{
				std::cout << "Writing sample input file to \"example.inp\".\n\n";
				
				// produce sample input file
				std::ofstream out("example.inp");
				if (out.bad())
					throw exception("Error: Cannot write to \"example.inp\"\n");
				
				out <<
					"# B-spline parameters \n"
					"# order      θ\n"
					"      4   0.63\n"
					"\n"
					"# real knot sequences\n"
					" 0.0  0.1   3   -1\n"
					" 0.0  2.0  60\n"
					"   4   20  58\n"
					"\n"
					"# complex knot sequences\n"
					"  60    -1\n"
					" 100\n"
					"  41\n"
					"\n"
					"# initial atomic states\n"
					"# ni\n"
					"  1\n"
					"# angular states (li, mi)\n"
					"  0  -1\n"
					"  0\n"
					"\n"
					"# final atomic states (nf, lf)\n"
					"  1  -1\n"
					"  0\n"
					"\n"
					"# angular momenta\n"
					"# L  S  limit\n"
					"  0  0  0\n"
					"\n"
					"# initial energies in Rydbergs\n"
					" 0.65   -1\n"
					" 0.95\n"
					" 301\n"
					"\n"
					"# magnetic field\n"
					" 0\n"
				;
				
				out.close();
				exit(0);
			}
			case 'i':
			{
				// set custom input file
				inputfile.open(optarg);
				if (inputfile.bad())
					throw exception("Error: Input file \"%s\" not found.\n", optarg);
				break;
			}
			case 'h':
			{
				// print usage information
				std::cout <<
					"------------------------------------------------       \n"
					"Available switches (short forms in parentheses):       \n"
					"------------------------------------------------       \n"
					"\nGeneral:                                             \n"
					"\t--example            (-e)  create sample input file  \n"
					"\t--help               (-h)  display this help         \n"
					"\t--input <filename>   (-i)  use custom input file     \n"
					"\t--zipfile <filename> (-z)  solution file to zip      \n"
					"\t--zipcount <number>  (-n)  zip samples               \n"
					"\t--mpi                (-m)  use MPI                   \n"
					"\t--stg-integ          (-a)  only do radial integrals  \n"
					"\t--stg-integ-solve    (-b)  only do integrals & solve \n"
					"\t--stg-extract        (-c)  only extract amplitudes   \n"
					"                                                       \n"
				;
				exit(0);
			}
			case 'z':
			{
				// zip B-spline expansion file
				zipfile = std::string(basename(optarg));
				break;
			}
			case 'n':
			{
				// zip samples
				zipcount = atol(optarg);
				break;
			}
			case 'R':
			{
				// zip bounding box
				zipmax = atof(optarg);
				break;
			}
			case 'a':
			{
				// run only the first part (computation of radial integrals)
				itinerary = StgRadial;
				break;
			}
			case 'b':
			{
				// run only the first two parts (computation of radial integrals and solution of the equations)
				itinerary = StgRadial & StgSolve;
				break;
			}
			case 'c':
			{
				// run only the third part (extraction of amplitudes)
				itinerary = StgExtract;
				break;
			}
#ifndef NO_MPI
			case 'm':
			{
				// use MPI
				parallel = true;
				break;
			}
#endif
			case -1:
			{
				// end of command line option list
				break;
			}
			default:
			{
				// unknown option
				std::cout << "Unknown option." << std::endl;
				abort();
			}
		}; // switch opetion
		
	} while (next_option != -1);
}

long read_int(std::ifstream& f)
{
	// text buffer
	std::string s;
	
	while (not f.eof())
	{
		// read string
		f >> s;
		
		// check length
		if (s.size() == 0)
			continue;
		
		// check if it is a beginning of a comment
		if (s[0] == '#')
		{
			// get the rest of the line
			std::getline(f, s);
			continue;
		}
		
		break;
	}
	
	// convert to long
	char* tail;
	long val = strtol(s.c_str(), &tail, 10);
	if (*tail != 0)
	{
		if (s == "*")
			throw false;
		else
			throw exception ("Can't read int.\n");
	}
	else
	{
		return val;
	}
}

double read_dbl(std::ifstream& f)
{
	// text buffer
	std::string s;
	
	while (not f.eof())
	{
		// read string
		f >> s;
		
		// check length
		if (s.size() == 0)
			continue;
		
		// check if it is a beginning of a comment
		if (s[0] == '#')
		{
			// get the rest of the line
			std::getline(f, s);
			continue;
		}
		
		break;
	}
	
	// convert to double
	char* tail;
	double val = strtod(s.c_str(), &tail);
	if (*tail != 0)
	{
		if (s == "*")
			throw false;
		else
			throw exception ("Can't read double.\n");
	}
	else
	{
		return val;
	}
}

void parse_input_file (
	std::ifstream& inputfile,
	int& order, double& ecstheta, 
	rArray& rknots, rArray& cknots,
	int& ni,
	std::vector<std::tuple<int,int,int>>& instates,
	std::vector<std::tuple<int,int,int>>& outstates,
	int& L, int& Spin, int& maxell,
	rArray& Ei, double& B
){
	double x;
	
	// load B-spline parameters
	try {
		order = read_int(inputfile);
		ecstheta = read_dbl(inputfile);
	} catch (std::exception e) {
		std::cerr << e.what() << std::endl;
		throw exception("Input error: Check B-spline parameters.\n");
	} catch (bool b) {
		throw exception("Wildcard not accepted here.\n");
	}
	
	std::cout << "\n-----   B-spline environment  -------\n";
	std::cout << "order = " << order << "\n";
	std::cout << "ecsθ = " << ecstheta << "\n";
	
	// load real knot data
	std::vector<double> rknots_begin, rknots_end, rknots_samples;
	try {
		while ((x = read_dbl(inputfile)) != -1.)
			rknots_begin.push_back(x);
		for (size_t i = 0; i < rknots_begin.size(); i++)
			rknots_end.push_back(read_dbl(inputfile));
		for (size_t i = 0; i < rknots_begin.size(); i++)
			rknots_samples.push_back(read_dbl(inputfile));
	} catch (std::exception e) {
		std::cerr << e.what() << std::endl;
		throw exception("Input error: Check real knot data.\n");
	} catch (bool b) {
		throw exception("Wildcard not accepted here.\n");
	}
	
	// construct real knot sequence
	for (unsigned i = 0; i < rknots_begin.size(); i++)
	{
		if (rknots_begin[i] > rknots_end[i])
		{
			std::cout << "\t" << rknots_begin[i] << " > " << rknots_end[i] << "\n";
			throw exception("Inconsistent knot specification!");
		}
		
		auto new_knots = linspace(rknots_begin[i], rknots_end[i], rknots_samples[i]);
		rknots = concatenate(rknots, new_knots);
	}
	
	std::cout << "\n----------   Real knots  ------------\n";
	for (auto knot = rknots.begin(); knot != rknots.end(); knot++)
		std::cout << *knot << " ";
	std::cout << std::endl;
	
	// load complex knot data
	std::vector<double> cknots_begin, cknots_end, cknots_samples;
	try {
		while ((x = read_dbl(inputfile)) != -1.)
			cknots_begin.push_back(x);
		for (size_t i = 0; i < cknots_begin.size(); i++)
			cknots_end.push_back(read_dbl(inputfile));
		for (size_t i = 0; i < cknots_begin.size(); i++)
			cknots_samples.push_back(read_int(inputfile));
	} catch (std::exception e) {
		std::cerr << e.what() << std::endl;
		throw exception("Input error: Check complex knot data.\n");
	} catch (bool b) {
		throw exception("Wildcard not accepted here.\n");
	}
	
	// construct complex(-to-be) knot sequence
	for (unsigned i = 0; i < cknots_begin.size(); i++)
		cknots = concatenate (
			cknots,
			linspace(
				cknots_begin[i],
				cknots_end[i],
				cknots_samples[i]
			)
		);

	std::cout << "\n---------  Complex knots  ------------\n";
	for (auto knot = cknots.begin(); knot != cknots.end(); knot++)
		std::cout << *knot << " ";
	std::cout << std::endl;
	
	// load initial principal quantum number
	try {
		ni = read_int(inputfile);
	} catch (std::exception e) {
		std::cerr << e.what() << std::endl;
		throw exception("Input error: Check \"ni\".\n");
	} catch (bool b) {
		throw exception("Wildcard not accepted here.\n");
	}
	
	// load initial atomic angular states
	std::vector<int> lis, mis;
	int maxli = 0;
	try {
		while ((x = read_int(inputfile)) != -1.)
		{
			lis.push_back(x);
			if (lis.back() >= ni)
				throw exception("Input error: Angular momentum greater than \"ni\".\n");
			if (lis.back() > maxli)
				maxli = lis.back();
		}
	} catch (std::exception e) {
		std::cerr << e.what() << std::endl;
		throw exception("Input error: Check initial atomic state data.\n");
	} catch (bool b) {
		throw exception("Wildcard not accepted here.\n");
	}
	
	for (size_t i = 0; i < lis.size(); i++)
	{
		try {
			
			mis.push_back(read_int(inputfile));
			if (std::abs(mis[i]) > lis[i])
				throw exception("Input error: Magnetic number greater than \"li\".\n");
			
			instates.push_back(std::make_tuple(ni,lis[i],mis[i]));
			
		} catch (std::exception e) {
			
			std::cerr << e.what() << std::endl;
			throw exception("Input error: Check initial atomic state data.\n");
			
		} catch (bool b) {
			
			// wildcard "*" found
			for (int j = -lis[j]; j <= lis[i]; j++)
			{
				mis.push_back(j);
				instates.push_back(std::make_tuple(ni,lis[i],j));
			}
			
		}
	}
	
	// load final atomic quantum numbers
	std::vector<int> nfs, lfs;
	int maxlf = 0;
	try {
		while ((x = read_int(inputfile)) != -1.)
			nfs.push_back(x);
	} catch (std::exception e) {
		std::cerr << e.what() << std::endl;
		throw exception("Input error: Check final atomic state data.\n");
	} catch (bool b) {
		throw exception("Wildcard not accepted here.\n");
	}
	
	for (size_t i = 0; i < nfs.size(); i++)
	{
		try {
			
			lfs.push_back(read_int(inputfile));
			
			if (nfs[i] == 0)
			{
				outstates.push_back(std::make_tuple(0,0,0));
				continue;
			}
			
			if (lfs[i] > nfs[i])
				throw exception("Input error: Angular momentum greater than \"nf\".\n");
			if (lfs[i] > maxlf)
				maxlf = lfs[i];
			
			outstates.push_back(std::make_tuple(nfs[i],lfs[i],0));
			
		} catch (std::exception e) {
			
			std::cerr << e.what() << std::endl;
			throw exception("Input error: Check final atomic state data.\n");
			
		} catch (bool b) {
			
			if (nfs[i] == 0)
			{
				outstates.push_back(std::make_tuple(0,0,0));
				continue;
			}
			
			// wildcard "*" found, add all allowed angular momenta
			for (int j = 0; j < nfs[i]; j++)
			{
				lfs.push_back(j);
				outstates.push_back(std::make_tuple(nfs[i],j,0));
			}
			
		}
	}
	
	
	// load total quantum numbers
	try {
		
		L = read_int(inputfile);
		Spin = read_int(inputfile);
		maxell = read_int(inputfile);
		
		if (L + L%2 + maxell < maxli)
			throw exception("Input error: ℓ is smaller than some initial angular momenta.\n");
		
		if (L + L%2 + maxell < maxlf)
			throw exception("Input error: ℓ is smaller than some final angular momenta.\n");
		
	} catch (std::exception e) {
		
		std::cerr << e.what() << std::endl;
		throw exception("Input error: Check angular momentum data.\n");
		
	} catch (bool b) {
		throw exception("Wildcard not accepted here.\n");
	}
	
	std::cout << "\n----------  Angular momentum limits  -------------\n";
	std::cout << "L = " << L << "\n";
	std::cout << "S = " << Spin << "\n";
	std::cout << "Π = " << L % 2 << "\n";
	std::cout << "ℓ = " << maxell << "\n";
	
	std::cout << "\n----------  Initial atomic states  -------------\n";
	for (auto state : instates)
		std::cout << "[" << std::get<0>(state) << " " << std::get<1>(state) << " " << std::get<2>(state) << "] ";
	std::cout << "\n";
	
	std::cout << "\n----------  Final atomic states  -------------\n";
	for (auto state : outstates)
		std::cout << "[" << std::get<0>(state) << " " << std::get<1>(state) << " " << std::get<2>(state) << "] ";
	std::cout << "\n";
	
	// load initial energies
	std::vector<double> Ei_begin, Ei_end, Ei_samples;
	try {
		while ((x = read_dbl(inputfile)) != -1.)
			Ei_begin.push_back(x);
		for (size_t i = 0; i < Ei_begin.size(); i++)
			Ei_end.push_back(read_dbl(inputfile));
		for (size_t i = 0; i < Ei_begin.size(); i++)
			Ei_samples.push_back(read_int(inputfile));
	} catch (std::exception e) {
		std::cerr << e.what() << std::endl;
		throw exception("Input error: Check energy data.\n");
	} catch (bool b) {
		throw exception("Wildcard not accepted here.\n");
	}
	
	// construct energy sequence
	for (unsigned i = 0; i < Ei_begin.size(); i++)
		Ei = concatenate(Ei, linspace(Ei_begin[i], Ei_end[i], Ei_samples[i]));
	
	std::cout << "\n---  Initial projectile energies  ----\n";
	std::cout << "lowest energy: " << Ei.front() << "\n";
	std::cout << "highest energy: " << Ei.back() << "\n";
	std::cout << "total enegies: " << Ei.size() << "\n";
	
	try {
		B = read_dbl(inputfile);
	} catch (std::exception e) {
		std::cerr << e.what() << std::endl;
		throw exception("Input error: Check magnetic field data.\n");
	} catch (bool b) {
		throw exception("Wildcard not accepted here.\n");
	}
	
	std::cout << "\n---------- Other parameters -----------\n";
	std::cout << "magnetic field: " << B << " a.u.\n\n";
}
