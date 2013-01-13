/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2012                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cstdio>
#include <cstring>
#include <cstdlib>

#include <getopt.h>

#include "arrays.h"
#include "spmatrix.h"

/*
 * Load data from input file. Default filename is "hex.inp", but different
 * can be specified using the --input/-i option.
 */

void parse_command_line(int argc, char* argv[], FILE*& inputfile)
{
	// set short options
	const char* const short_options  = "eih";
	
	// set long options
	const option long_options[] = {
		{"example",           0,   0, 'e'},
		{"input",             1,   0, 'i'},
		{"help",              0,   0, 'h'},
		{0,                   0,   0,   0}
	};
	
	// switch all options
	int next_option = 0;
	do {
		next_option = getopt_long(argc, argv, short_options, long_options, 0);
		switch (next_option)
		{
			case 'e':
				// produce sample input file
				FILE *h;
				h = fopen("example.inp", "w");
				if (h == 0)
				{
					printf("Error: Cannot write to \"example.inp\"\n");
					abort();
				}
				fprintf(h,
					"# B-spline parameters \n"
					"# order     R0    Rmax       θ\n"
					"      4   60.0    100.0    0.63\n"
					"\n"
					"# real knot sequences\n"
					"0.0  0.1   3.0   -1\n"
					"0.0  2.0  60.0\n"
					"  4   20    58\n"
					"\n"
					"# complex knot sequences\n"
					"60.0    -1\n"
					"100.0\n"
					"  41\n"
					"\n"
					"# initial / final atomic state\n"
					"# ni minli maxli maxnf maxlf\n"
					"  1     -1    -1     1    -1\n"
					"\n"
					"# angular momenta\n"
					"# L  ℓ\n"
					"  0  0\n"
					"\n"
					"# initial energies in Rydbergs\n"
					"0.65   -1\n"
					"0.95\n"
					"301\n"
					"\n"
					"# magnetic field\n"
					"0\n"
				);
				fclose(h);
				exit(0);
				break;
				
			case 'i':
				// set custom input file
				FILE *g;
				g = fopen(optarg, "r");
				if (g)
					inputfile = g;
				else
				{
					printf("Error: Input file \"%s\" not found.\n", optarg);
					abort();
				}
				break;
				
			case 'h':
				// print usage information
				printf(
					"                                                       \n"
					"         __  __                                        \n"
					"        / /_/ /  ___   __  __                          \n"
					"       /  _  /  / __\\  \\ \\/ /                       \n"
					"      /_/ /_/   \\__\\   /_/\\_\\                      \n"
					"                                                       \n"
					"    Jakub Benda, MFF UK (c) 2011                       \n"
					"                                                       \n"
					"------------------------------------------------       \n"
					"Available switches (short forms in parentheses):       \n"
					"------------------------------------------------       \n"
					"\nGeneral:                                             \n"
					"\t--example           (-e)   create sample input file  \n"
					"\t--help              (-h)   display this help         \n"
					"\t--input <filename>  (-i)   use custom input file     \n"
					"                                                       \n"
				);
				exit(0);
				
			case -1:
				// end of command line option list
				break;
				
			default:
				// unknown option
				printf("Unknown option.\n");
				abort();
		};
	} while (next_option != -1);
}

long read_int(FILE *f) throw (int)
{
	// text buffer
	char buff[1024];
	
	for/*ever*/(;;)
	{
		// read string
		if (fscanf(f, "%s", buff) == EOF)
			throw 0;
		
		if (strlen(buff) == 0)
			continue;
		
		// check if it is a beginning of a comment
		if (buff[0] == '#')
		{
			// get the rest of the line
			fgets(buff, 1024, f);
			continue;
		}
		
		break;
	}
	
	// convert to long
	char* tail;
	long val = strtol(buff, &tail, 10);
	if (*tail != 0)
		throw 0;
	else
		return val;
}

double read_dbl(FILE *f) throw (int)
{
	// text buffer
	char buff[1024];
	
	for/*ever*/(;;)
	{
		// read string
		if (fscanf(f, "%s", buff) == EOF)
			throw 0;
		
		if (strlen(buff) == 0)
			continue;
		
		// check if it is a beginning of a comment
		if (buff[0] == '#')
		{
			// get the rest of the line
			fgets(buff, 1024, f);
			continue;
		}
		
		break;
	}
	
	// convert to double
	char* tail;
	double val = strtod(buff, &tail);
	if (*tail != 0)
		throw 0;
	else
		return val;
}

void parse_input_file(
	FILE* inputfile,
	int& order, double& R0, double& ecstheta, double& Rmax,
	rArray& rknots, rArray& cknots,
	int& ni, int& maxnf,
	int& minli, int& maxli, int& maxlf,
	int& L, int& maxell,
	rArray& Ei, double& B
){
	double x;
	
	// load B-spline parameters
	try {
		order = read_int(inputfile);
		R0 = read_dbl(inputfile);
		Rmax = read_dbl(inputfile);
		ecstheta = read_dbl(inputfile);
	} catch (...) {
		fprintf(stderr, "Input error: Check B-spline parameters.\n");
		abort();
	}
	
	fprintf(stdout, "\n-----   B-spline environment  -------\n");
	fprintf(stdout, "order = %d\n", order);
	fprintf(stdout, "R0 = %g\n", R0);
	fprintf(stdout, "Rmax = %g\n", Rmax);
	fprintf(stdout, "ecsθ = %g\n", ecstheta);
	
	// load real knot data
	std::vector<double> rknots_begin, rknots_end, rknots_samples;
	try {
		while ((x = read_dbl(inputfile)) != -1.)
			rknots_begin.push_back(x);
		for (size_t i = 0; i < rknots_begin.size(); i++)
			rknots_end.push_back(read_dbl(inputfile));
		for (size_t i = 0; i < rknots_begin.size(); i++)
			rknots_samples.push_back(read_dbl(inputfile));
	} catch (...) {
		fprintf(stderr, "Input error: Check real knot data.\n");
		abort();
	}
	
	// construct real knot sequence
	for (unsigned i = 0; i < rknots_begin.size(); i++)
	{
		auto new_knots = linspace(rknots_begin[i], rknots_end[i], rknots_samples[i]);
		rknots = concatenate(rknots, new_knots);
	}
	
	fprintf(stdout, "\n----------   Real knots  ------------\n");
	for (auto knot = rknots.begin(); knot != rknots.end(); knot++)
		fprintf(stdout, "%g  ", *knot);
	fprintf(stdout, "\n");
	
	// load complex knot data
	std::vector<double> cknots_begin, cknots_end, cknots_samples;
	try {
		while ((x = read_dbl(inputfile)) != -1.)
			cknots_begin.push_back(x);
		for (size_t i = 0; i < cknots_begin.size(); i++)
			cknots_end.push_back(read_dbl(inputfile));
		for (size_t i = 0; i < cknots_begin.size(); i++)
			cknots_samples.push_back(read_int(inputfile));
	} catch (...) {
		fprintf(stderr, "Input error: Check complex knot data.\n");
		abort();
	}
	
	// construct complex(-to-be) knot sequence
	for (unsigned i = 0; i < cknots_begin.size(); i++)
		cknots = concatenate(
			cknots,
			linspace(
				cknots_begin[i],
				cknots_end[i],
				cknots_samples[i]
			)
		);

	fprintf(stdout, "\n---------  Complex knots  ------------\n");
	for (auto knot = cknots.begin(); knot != cknots.end(); knot++)
		fprintf(stdout, "%g  ", *knot);
	fprintf(stdout, "\n");
	
	// load atomic states
	try {
		ni = read_int(inputfile);
		minli = read_int(inputfile);
		maxli = read_int(inputfile);
		maxnf = read_int(inputfile);
		maxlf = read_int(inputfile);
	} catch (...) {
		fprintf(stderr, "Input error: Check atomic state data.\n");
		abort();
	}
	
	// load total quantum numbers
	try {
		L = read_int(inputfile);
		maxell = read_int(inputfile);
	} catch (...) {
		fprintf(stderr, "Input error: Check angular momentum data.\n");
		abort();
	}
	
	fprintf(stdout, "\n----------  Angular momentum limits  -------------\n");
	fprintf(stdout, "L = %d, ℓ = %d\n", L, maxell);
	
	fprintf(stdout, "\n----------  Initial atomic states  -------------\n");
	for (int li = (minli < 0) ? 0 : minli; li <= maxell and li < ni and (maxli < 0 or li <= maxli); li++)
		fprintf(stdout, "[%d %d] ", ni, li);
	fprintf(stdout, "\n");
	
	fprintf(stdout, "\n----------  Final atomic states  -------------\n");
	for (int nf = 0; nf <= maxnf; nf++)
		for (int lf = 0; lf <= maxell and lf < nf and (maxlf < 0 or lf <= maxlf); lf++)
			fprintf(stdout, "[%d %d] ", nf, lf);
	fprintf(stdout, "\n");
	
	// load initial energies
	std::vector<double> Ei_begin, Ei_end, Ei_samples;
	try {
		while ((x = read_dbl(inputfile)) != -1.)
			Ei_begin.push_back(x);
		for (size_t i = 0; i < Ei_begin.size(); i++)
			Ei_end.push_back(read_dbl(inputfile));
		for (size_t i = 0; i < Ei_begin.size(); i++)
			Ei_samples.push_back(read_int(inputfile));
	} catch (...) {
		fprintf(stderr, "Input error: Check energy data.\n");
		abort();
	}
	
	// construct energy sequence
	for (unsigned i = 0; i < Ei_begin.size(); i++)
		Ei = concatenate(Ei, linspace(Ei_begin[i], Ei_end[i], Ei_samples[i]));
	
	fprintf(stdout, "\n---  Initial projectile energies  ----\n");
	fprintf(stdout, "lowest energy: %g\n", Ei.front());
	fprintf(stdout, "highest energy: %g\n", Ei.back());
	fprintf(stdout, "total enegies: %ld\n", Ei.size());
	
	try {
		B = read_dbl(inputfile);
		fprintf(stdout, "\n---------- Other parameters -----------\n");
		fprintf(stdout, "magnetic field: %g a.u.\n", B);
	} catch (...) {
		fprintf(stderr, "Input error: Check magnetic field data.\n");
		abort();
	}
	
	fprintf(stdout, "\n");
}
