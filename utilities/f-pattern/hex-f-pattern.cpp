#include <iostream>
#include <cstdlib>

#include "hex-special.h"

// Usage
//   hex-f-pattern <L> <Pi> <nL>

int main (int argc, char * argv[])
{
    if (argc != 4)
    {
        std::cout << std::endl << "Usage:" << std::endl << std::endl;
        std::cout << "  hex-f-pattern <L> <Pi> <nL>" << std::endl << std::endl;
        return 0;
    }

    int L = std::atoi(argv[1]);
    int Pi = std::atoi(argv[2]);
    int nL = std::atoi(argv[3]);

    int maxell = L + nL;
    int maxlambda = 2 * maxell;
    
    // coupled angular momentum pairs
    std::vector<std::pair<int,int>> coupled_states;
   
    std::cout << "Setting up the coupled angular states..." << std::endl;
 
    // for given L, Π and levels list all available (ℓ₁ℓ₂) pairs
    for (int ell = 0; ell <= nL; ell++)
    {
        std::cout << "\t-> [" << ell << "] ";
        
        // get sum of the angular momenta for this angular level
        int sum = 2 * ell + L + Pi;
        
        // for all angular momentum pairs that do compose L
        for (int l1 = ell; l1 <= sum - ell; l1++)
        {
            std::cout << "(" << l1 << "," << sum - l1 << ") ";
            coupled_states.push_back(std::make_pair(l1, sum - l1));
        }
        std::cout << std::endl;
    }
    
    std::cout << std::endl;
    
    for (int lambda = 0; lambda <= maxlambda; lambda++)
    {
        int n = 0;
        std::cout << "lambda = " << lambda << std::endl;
        for (std::pair<int,int> const & ll : coupled_states)
        {
            for (std::pair<int,int> const & llp : coupled_states)
            {
                double f = special::computef(lambda, ll.first, ll.second, llp.first, llp.second, L);
                if (f != 0)
                {
                    std::cout << "* ";
                    n++;
                }
                else
                    std::cout << ". ";
            }
            std::cout << std::endl;
        }
        
        std::cout << "-> Used in " << n << " matrix super-blocks." << std::endl << std::endl;
    }
    
    return 0;
}
