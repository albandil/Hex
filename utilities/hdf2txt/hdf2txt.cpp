#include <iostream>

#include "../../src/arrays.h"
#include "../../src/arrays.cpp"

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		std::cout << "HDF file missing!\n";
		return 0;
	}
	
	rArray a;
	load_array(a, argv[1]);
	write_array(a);
	
	return 0;
}
