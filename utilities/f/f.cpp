/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2013                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "../../src/angs.h"
#include "../../src/angs.cpp"

#include <iostream>

int main (int argc, char* argv[])
{
	if (argc != 7)
	{
		std::cout << "\nUsage:\n\t./f <lambda> <l1> <l2> <l1p> <l2p> <L>\n\n";
		exit(0);
	}
	
	int lambda = atoi(argv[1]);
	int l1 = atoi(argv[2]);
	int l2 = atoi(argv[3]);
	int l1p = atoi(argv[4]);
	int l2p = atoi(argv[5]);
	int L = atoi(argv[6]);
	
	std::cout << "f[" << lambda << "](" << l1 << "," << l2 << "," << l1p << ","
	          << l2p << ") = " << computef(lambda,l1,l2,l1p,l2p,L) << "\n";
	
	return 0;
}
