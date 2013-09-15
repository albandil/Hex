#include <iostream>

#include "../../src/arrays.h"
#include "../../src/arrays.cpp"
#include "../../src/hdffile.cpp"

int main(int argc, char *argv[])
{
	H5::Exception::dontPrint();
	
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
		for (int i = 0; i < a.size(); i++)
			std::cout << a[i] << "\n";
	}
	else
	{
		// head HDF file
		rArray a;
		a.hdfload(hdf);
	
		// write raw data
		for (int i = 0; i < a.size(); i++)
			std::cout << a[i] << "\n";
	}
	
	return 0;
}
