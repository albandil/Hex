#include <iostream>

#include "../../src/arrays.h"
#include "../../src/arrays.cpp"

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		std::cout << "\nUsage:\n\t./hdf2txt <HDFfile>\n\n";
		exit(0);
	}
	
	// head HDF file
	rArray a;
	a.hdfload(argv[1]);
	
	// write raw data
	write_array(a);
	
	return 0;
}
