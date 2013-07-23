#include <iostream>
#include <vector>
#include <string>

#include "../../src/arrays.h"
#include "../../src/arrays.cpp"
#include "../../src/hdffile.cpp"

int main(int argc, char *argv[])
{
	H5::Exception::dontPrint();
	
	std::vector<std::string> files;
	bool compress = false;
	
	for (int iarg = 1; iarg < argc; iarg++)
	{
		if (strcmp(argv[iarg],"--compress") == 0)
			compress = true;
		else
			files.push_back(std::string(argv[iarg]));
	}
	
	if (files.size() != 2)
	{
		std::cout << "\nUsage:\n\thdf2hdf [--compress] <file1> <file2>\n";
		return -1;
	}
	
	rArray a;
	
	if (not a.hdfload(files[0].c_str()))
	{
		std::cout << "Can't open the file \"" << argv[1] << "\"\n";
		return -1;
	}
	
	if (not a.hdfsave(files[1].c_str(), compress))
	{
		std::cout << "Error while writing file \"" << argv[2] << "\"\n";
		return -1;
	}
	
	return 0;
}
