#include <iostream>
#include <cmath>

int triangles (int maxl, int L)
{
	int n = 0;
	
	for (int l1 = 0; l1 <= maxl; l1++)
	for (int l2 = 0; l2 <= maxl; l2++)
	{
		if (std::abs(l1-l2) <= L and l1+l2 >= L)
			n++;
	}
	
	return n;
}

int main (int argc, char *argv[])
{
	std::cout << "maxl/L\t";

	for (int L = 0; L <= 10; L++)
		std::cout << L << "\t";

	std::cout << std::endl;
	
	for (int maxl = 0; maxl <= 10; maxl++)
	{
		std::cout << maxl << "\t";

		for (int L = 0; L <= 10; L++)
		{
			std::cout << triangles(maxl,L) << "\t";
		}

		std::cout << std::endl;
	}

	return 0;
}
