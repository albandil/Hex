#include <iostream>
#include <cstdlib>
#include <string>

#include "../../src/complex.h"
#include "../../src/arrays.h"
#include "../../src/arrays.cpp"
#include "../../src/misc.h"
#include "../../src/compact.h"
#include "../../src/chebyshev.h"

template<typename T> void load_and_write(const char* hdf, int samples)
{
	// load coefficients
	NumberArray<T> coefs;
	if (not coefs.hdfload(hdf))
	{
		std::cerr << "Can't read file \"" << hdf << "\"\n";
		exit(-1);
	}
	
	// create the expansion
	Chebyshev<double,T> expansion(coefs);
	
	// evaluate the expansion
	for (int i = 0; i <= samples; i++)
	{
		double x = (2.*i-samples) / samples;
		
		Complex y = expansion.clenshaw(x,coefs.size());
		std::cout << x << "\t" << y.real() << "\t" << y.imag() << "\n";
	}
}

int main(int argc, char *argv[])
{
	bool cpx = false;	// whether to zip complex expansion
	std::string hdf;
	int samples = -1;
	
	for (int iarg = 1; iarg < argc; iarg++)
	{
		if (strcmp(argv[iarg],"--complex") == 0)
		{
			cpx = true;
		}
		else if (hdf.empty())
		{
			hdf = std::string(argv[iarg]);
		}
		else
		{
			char* tail;
			samples = strtol(argv[iarg], &tail, 10);
			if (*tail != 0)
				samples = -1;
		}
	}
	
	if (samples < 0)
	{
		std::cerr << "\nUsage: ./chebeval [--complex] <HDFfile> <samples>\n\n";
		exit(-1);
	}
	
	if (cpx)
		load_and_write<Complex>(hdf.c_str(), samples);
	else
		load_and_write<double>(hdf.c_str(), samples);
	
	return 0;
}
