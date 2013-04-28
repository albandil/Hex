#include <iostream>
#include <cstdlib>

#include "../../src/complex.h"
#include "../../src/arrays.h"
#include "../../src/arrays.cpp"
#include "../../src/misc.h"
#include "../../src/compact.h"
#include "../../src/chebyshev.h"

int main(int argc, char *argv[])
{
	if (argc < 3)
	{
		std::cerr << "\nUsage: ./chebeval <HDFfile> <samples>\n\n";
		exit(-1);
	}
	
	// load coefficients
	rArray coefs;
	if (not coefs.hdfload(argv[1]))
	{
		std::cerr << "Can't read file \"" << argv[1] << "\"\n";
		exit(-1);
	}
	
	// create the expansion
	Chebyshev<double,double> expansion(coefs);
	
	// get sample count
	int samples = atoi(argv[2]);
	
	// evaluate the expansion
	for (int i = 0; i <= samples; i++)
	{
		double x = double(2*i-samples) / double(samples);
		std::cout << x << "\t" << expansion.clenshaw(x,coefs.size()) << "\n";
	}
	
	return 0;
}
