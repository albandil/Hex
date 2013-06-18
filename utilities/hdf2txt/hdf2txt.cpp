#include <iostream>

#include "../../src/arrays.h"
#include "../../src/arrays.cpp"

int main(int argc, char *argv[])
{
	bool cpx = false;	// whether to write complex expansion
	std::string hdf;
	
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
			break;
	}
	
	if (argc < 2 or hdf.empty())
	{
		std::cout << "\nUsage:\n\t./hdf2txt [--complex] <HDFfile>\n\n";
		exit(0);
	}
	
	if (cpx)
	{
		// head HDF file
		cArray a;
		a.hdfload(hdf);
	
		// write raw data
		write_array(a);
	}
	else
	{
		// head HDF file
		rArray a;
		a.hdfload(hdf);
	
		// write raw data
		write_array(a);
	}
	
	return 0;
}
