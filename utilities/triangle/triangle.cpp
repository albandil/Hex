#include <iostream>
#include <cstdlib>

#include "angs.h"

int main (int argc, char *argv[])
{
	if (argc < 3)
	{
		std::cout << "Use: triangle <max L> <max ell>\n";
		return 0;
	}
	
	int MaxL = atoi(argv[1]);
	int MaxEll = atoi(argv[2]);
	
	std::cout << "maxl/L\t";

	for (int L = 0; L <= MaxL; L++)
		std::cout << L << "\t";

	std::cout << std::endl;
	
	for (int maxl = 0; maxl <= MaxEll; maxl++)
	{
		std::cout << maxl << "\t";
		
		for (int L = 0; L <= MaxL; L++)
			std::cout << triangle_count(L,maxl) << "\t";
		
		std::cout << std::endl;
	}

	return 0;
}
