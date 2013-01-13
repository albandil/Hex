#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <getopt.h>

#include "hex-db.h"

/*

# scattering amplitude
seq 0.01 0.01 3.14    | hex-db --database="hex.db" --scatamp    --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --S=0

# differential cross section
seq 0.01 0.01 3.14    | hex-db --database="hex.db" --dcs        --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --S=0

# momentum transfer
seq 0.650 0.001 0.850 | hex-db --database="hex.db" --momtransf  --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --L=0 --S=0

# integral cross section
seq 0.650 0.001 0.850 | hex-db --database="hex.db" --integcs    --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --L=0 --S=0

# sum integral cross section
seq 0.650 0.001 0.850 | hex-db --database="hex.db" --sumintegcs --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0

# collision strength
seq 0.650 0.001 0.850 | hex-db --database="hex.db" --collstr    --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --L=0 --S=0

# total cross section
seq 0.650 0.001 0.850 | hex-db --database="hex.db" --tcs        --ni=1 --li=0 --mi=0

*/

struct Integer
{
	Integer() : val(0), set(false) {}
	int val;
	bool set;
};

struct Real
{
	Real() : val(0), set(false) {}
	double val;
	bool set;
};

int main(int argc, char* argv[])
{
	// program parameters
	char dbname[1024] = "hex.db";
	Integer ni, li, mi, nf, lf, mf, L, S;
	Real E;
	bool scatamp = false, dcs = false, momtransf = false, integcs = false,
		sumintegcs = false, collstr = false, tcs = false, extrapolat = false,
		create_new = false;
	
	// available command line arguments
	const char* const short_options = "h1234567abcdefLS";
	const struct option long_options[] = {
		{ "help",       0, 0, 'h' },
		{ "new",        0, 0, 'n' },
		{ "database",   1, 0, 'D' },
		{ "scatamp",    0, 0, '1' },
		{ "dcs",        0, 0, '2' },
		{ "momtransf",  0, 0, '3' },
		{ "integcs",    0, 0, '4' },
		{ "sumintegcs", 0, 0, '5' },
		{ "collstr",    0, 0, '6' },
		{ "tcs",        0, 0, '7' },
		{ "extrapolat", 0, 0, '8' },
		{ "ni",         1, 0, 'a' },
		{ "li",         1, 0, 'b' },
		{ "mi",         1, 0, 'c' },
		{ "nf",         1, 0, 'd' },
		{ "lf",         1, 0, 'e' },
		{ "mf",         1, 0, 'f' },
		{ "L",          1, 0, 'L' },
		{ "S",          1, 0, 'S' },
		{ "E",          1, 0, 'E' },
		{   0,          0, 0,  0  }
	};
	
	// loop over command line arguments
	int next_option;
	char * tail;
	while ((next_option = getopt_long(argc, argv, short_options, long_options, 0)) != -1)
	switch (next_option) {
		default:
			// unknown switch
			printf("Unknown switch %s\n", argv[optind]);
			/*break;*/
		case 'h':
			// print usage information
			printf("\nUsage examples:\n");
			printf("Retrieving scattering amplitude:\n");
			printf("# seq 0.01 0.01 3.14    | hex-db --database=\"hex.db\" --scatamp    --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --S=0 --E=0.75\n\n");
			printf("Retrieving differential cross section:\n");
			printf("# seq 0.01 0.01 3.14    | hex-db --database=\"hex.db\" --dcs        --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --S=0 --E=0.75\n\n");
			printf("Retrieving momentum transfer:\n");
			printf("# seq 0.650 0.001 0.850 | hex-db --database=\"hex.db\" --momtransf  --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --L=0 --S=0\n\n");
			printf("Retrieving integral cross section:\n");
			printf("# seq 0.650 0.001 0.850 | hex-db --database=\"hex.db\" --integcs    --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --L=0 --S=0\n\n");
			printf("Retrieving sum integral cross section:\n");
			printf("# seq 0.650 0.001 0.850 | hex-db --database=\"hex.db\" --sumintegcs --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0\n\n");
			printf("Retrieving collision strength:\n");
			printf("# seq 0.650 0.001 0.850 | hex-db --database=\"hex.db\" --collstr    --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --L=0 --S=0\n\n");
			printf("Retrieving total cross section:\n");
			printf("# seq 0.650 0.001 0.850 | hex-db --database=\"hex.db\" --tcs        --ni=1 --li=0 --mi=0\n\n");
			printf("\nOn some locales a decimal comma is used instead of decimal point (e.g. on native Czech systems).\n");
			printf("This is not compatible with raw C/C++, one has to turn off the locale by executing\n");
			printf("# setenv LC_NUMERIC en_GB.utf8\n");
			printf("or\n");
			printf("# export LC_NUMERIC=en_GB.utf8\n");
			printf("depending on version of your shell.\n\n");
			exit(0);
		case 'D':
			// save database name (default "hex.db")
			strcpy(dbname, optarg);
			break;
		case 'n':
			// create new database
			create_new = true;
			break;
		case '1':
			// queue scattering amplitude for computation
			scatamp = true;
			break;
		case '2':
			// queue dcs for computation
			dcs = true;
			break;
		case '3':
			// queue momtransf for computation
			momtransf = true;
			break;
		case '4':
			// queue integcs for computation
			integcs = true;
			break;
		case '5':
			// queue sumintegcs for computation
			sumintegcs = true;
			break;
		case '6':
			// queue collstr for computation
			collstr = true;
			break;
		case '7':
			// queue tcs for computation
			tcs = true;
			break;
		case '8':
			// queue extrapolat for computation
			extrapolat = true;
			break;
		case 'a':
			// save ni
			ni.val = strtol(optarg, &tail, 10);
			if (*tail == 0) ni.set = true;
			break;
		case 'b':
			// save li
			li.val = strtol(optarg, &tail, 10);
			if (*tail == 0) li.set = true;
			break;
		case 'c':
			// save mi
			mi.val = strtol(optarg, &tail, 10);
			if (*tail == 0) mi.set = true;
			break;
		case 'd':
			// save nf
			nf.val = strtol(optarg, &tail, 10);
			if (*tail == 0) nf.set = true;
			break;
		case 'e':
			// save lf
			lf.val = strtol(optarg, &tail, 10);
			if (*tail == 0) lf.set = true;
			break;
		case 'f':
			// save mf
			mf.val = strtol(optarg, &tail, 10);
			if (*tail == 0) mf.set = true;
			break;
		case 'L':
			// save L
			L.val = strtol(optarg, &tail, 10);
			if (*tail == 0) L.set = true;
			break;
		case 'S':
			// save S
			S.val = strtol(optarg, &tail, 10);
			if (*tail == 0) S.set = true;
			break;
		case 'E':
			// save E
			E.val = strtod(optarg, &tail);
			if (*tail == 0) E.set = true;
			break;
	}
	
	// create new database if asked to
	if (create_new)
		create_new_database(dbname);
	
	if (not scatamp and not dcs and not momtransf and not integcs and
		not sumintegcs and not collstr and not tcs and not extrapolat)
		return 0;	// nothing to do
	
	// read all standard input
	char str[100];
	rArray nums;
	while (gets(str) != 0)
	{
		// manage Czech locale
		for (char& c : str) if (c == ',') c = '.';
		
		// convert to number
		double num = strtod(str, &tail);
		if (*tail == 0)
		{
			nums.push_back(num);
		}
		else
		{
			fprintf(stderr, "The input \"%s\" is not a valid number.", str);
// 			if (strchr(str, ',') != 0)
// 				fprintf(stderr, " It contains \',\' symbol... Did you switch your locale?");
// 			fprintf(stderr, "\n");
		}
	}
	
	// open the database
	initialize(dbname);
	
	// retrieve all requested data
	if (scatamp)
	{
		// get scattering amplitude
		if (ni.set and li.set and mi.set and nf.set and lf.set and mf.set and S.set and E.set)
		{
			cArray fs = scattering_amplitude(ni.val, li.val, mi.val, nf.val, lf.val, mf.val, S.val, E.val, nums);
			printf("# theta\tRe f\tIm f\n");
			for (size_t i = 0; i < fs.size(); i++)
				printf("%g\t%.15g\t%.15g\n", nums[i], fs[i].real(), fs[i].imag());
		}
		else
		{
			printf("Some quantum numbers are missing. Run the program with --help switch to find out.\n");
			exit(0);
		}
	}
	if (dcs)
	{
		// check compulsory parameters
		if (not ni.set)
		{
			printf("'ni' is missing.\n");
			exit(0);
		}
		if (not nf.set)
		{
			printf("'nf' is missing.\n");
			exit(0);
		}
		if (not E.set)
		{
			printf("'E' is missing.\n");
			exit(0);
		}
		
		// summed differential cross section
		rArray sigsum(nums.size());
		
		// for all optional parameters
		for (int __S__ = 0; __S__ <= 1; __S__ ++)
		{
			if (S.set and S.val != __S__) continue;
			for (int __li__ = 0; __li__ < ni.val; __li__++)
			{
				if (li.set and li.val != __li__) continue;
				for (int __lf__ = 0; __lf__ < nf.val; __lf__++)
				{
					if (lf.set and lf.val != __lf__) continue;
					for (int __mi__ = -__li__; __mi__ <= __li__; __mi__++)
					{
						if (mi.set and mi.val != __mi__) continue;
						for (int __mf__ = -__lf__; __mf__ <= __lf__; __mf__++)
						{
							if (mf.set and mf.val != __mf__) continue;
							
							// get partial cross section for this seto of qnumbers
							rArray sig = differential_cross_section(
								ni.val, __li__, __mi__,
								nf.val, __lf__, __mf__,
								__S__, E.val, nums
							);
							
							// add to the sum
							for (size_t i = 0; i < sigsum.size(); i++)
								sigsum[i] += sig[i];
						}
					}
				}
			}
		}
		
		printf("# theta\tdσ/dΩ\n");
		for (size_t i = 0; i < sigsum.size(); i++)
			printf("%g\t%.15g\n", nums[i], sigsum[i]);
	}
	if (momtransf)
	{
		// get momentum transfer
		if (ni.set and li.set and mi.set and nf.set and lf.set and mf.set and L.set and S.set)
		{
			rArray eta = momentum_transfer(ni.val, li.val, mi.val, nf.val, lf.val, mf.val, L.val, S.val, nums);
			printf("# theta\tη\n");
			for (size_t i = 0; i < eta.size(); i++)
				printf("%g\t%.15g\n", nums[i], eta[i]);
		}
		else
		{
			printf("Some quantum numbers are missing. Run the program with --help switch to find out.\n");
			exit(0);
		}
	}
	if (integcs)
	{
		// get integral cross section
		if (ni.set and li.set and mi.set and nf.set and lf.set and mf.set and L.set and S.set)
		{
			rArray sig = integral_cross_section(ni.val, li.val, mi.val, nf.val, lf.val, mf.val, L.val, S.val, nums);
			printf("# theta\tσ_LS\n");
			for (size_t i = 0; i < sig.size(); i++)
				printf("%g\t%.15g\n", nums[i], sig[i]);
		}
		else
		{
			printf("Some quantum numbers are missing. Run the program with --help switch to find out.\n");
			exit(0);
		}
	}
	if (sumintegcs)
	{
		// check compulsory parameters
		if (not ni.set)
		{
			printf("'ni' is missing.\n");
			exit(0);
		}
		if (not nf.set)
		{
			printf("'nf' is missing.\n");
			exit(0);
		}
		
		// summed differential cross section
		rArray sigsum(nums.size());
		
		// for all optional parameters
		for (int __li__ = 0; __li__ < ni.val; __li__++)
		{
			if (li.set and li.val != __li__) continue;
			for (int __lf__ = 0; __lf__ < nf.val; __lf__++)
			{
				if (lf.set and lf.val != __lf__) continue;
				for (int __mi__ = -__li__; __mi__ <= __li__; __mi__++)
				{
					if (mi.set and mi.val != __mi__) continue;
					for (int __mf__ = -__lf__; __mf__ <= __lf__; __mf__++)
					{
						if (mf.set and mf.val != __mf__) continue;
						
						rArray sig = complete_cross_section(
							ni.val, __li__, __mi__,
							nf.val, __lf__, __mf__,
							nums
						);
						
						// add to the sum
						for (size_t i = 0; i < sigsum.size(); i++)
							sigsum[i] += sig[i];
					}
				}
			}
		}
		
		printf("# theta\tσ\n");
		for (size_t i = 0; i < sigsum.size(); i++)
			printf("%g\t%.15g\n", nums[i], sigsum[i]);
	}
	if (collstr)
	{
		// get collision strength
		if (ni.set and li.set and mi.set and nf.set and lf.set and mf.set and L.set and S.set)
		{
			rArray omega = collision_strength(ni.val, li.val, mi.val, nf.val, lf.val, mf.val, L.val, S.val, nums);
			printf("# theta\tΩ\n");
			for (size_t i = 0; i < omega.size(); i++)
				printf("%g\t%.15g\n", nums[i], omega[i]);
		}
		else
		{
			printf("Some quantum numbers are missing. Run the program with --help switch to find out.\n");
			exit(0);
		}
	}
	if (tcs)
	{
		// get total cross section
		if (ni.set and li.set and mi.set)
		{
			rArray sig = total_cross_section(ni.val, li.val, mi.val, nums);
			printf("# theta\tσ*\n");
			for (size_t i = 0; i < sig.size(); i++)
				printf("%g\t%.15g\n", nums[i], sig[i]);
		}
		else
		{
			printf("Some quantum numbers are missing. Run the program with --help switch to find out.\n");
			exit(0);
		}
	}
	if (extrapolat)
	{
		// check compulsory parameters
		if (not ni.set)
		{
			printf("'ni' is missing.\n");
			exit(0);
		}
		if (not nf.set)
		{
			printf("'nf' is missing.\n");
			exit(0);
		}
		
		// summed differential cross section
		rArray sigsum(nums.size());
		
		// for all optional parameters
		for (int __li__ = 0; __li__ < ni.val; __li__++)
		{
			if (li.set and li.val != __li__) continue;
			for (int __lf__ = 0; __lf__ < nf.val; __lf__++)
			{
				if (lf.set and lf.val != __lf__) continue;
				for (int __mi__ = -__li__; __mi__ <= __li__; __mi__++)
				{
					if (mi.set and mi.val != __mi__) continue;
					for (int __mf__ = -__lf__; __mf__ <= __lf__; __mf__++)
					{
						if (mf.set and mf.val != __mf__) continue;
						
						rArray sig = aitkenD2_cross_section(
							ni.val, __li__, __mi__,
							nf.val, __lf__, __mf__,
							nums
						);
						
						// add to the sum
						for (size_t i = 0; i < sigsum.size(); i++)
							sigsum[i] += sig[i];
					}
				}
			}
		}
		
		printf("# theta\tσ\n");
		for (size_t i = 0; i < sigsum.size(); i++)
			printf("%g\t%.15f\n", nums[i], sigsum[i]);
	}
	
	return EXIT_SUCCESS;
}
