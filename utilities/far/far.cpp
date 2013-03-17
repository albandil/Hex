#include <iostream>

#include "../../src/hydrogen.h"
#include "../../src/hydrogen.cpp"
#include "../../src/specf.h"
#include "../../src/specf.cpp"

int main(int argc, char* argv[])
{
	// Expecting n,l,eps triplet.
	if (argc != 4)
	{
		std::cout << "Use: far <n> <l> <eps>\n";
		return 0;
	}
	
	int n = atoi(argv[1]);
	int l = atoi(argv[2]);
	double eps = atof(argv[3]);
	
	std::cout << Hydrogen::getFarRadius(n,l,eps);
	
	return 0;
}
