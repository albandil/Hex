/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2014                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>

#include "complex.h"
#include "radial.h"
#include "version.h"

int main (void)
{
    std::cout << logo() << std::endl;
    
    Complex I = Idir_nBound
    (
        0, 0,
        1, 0, 0, 2.,
        1, 0, 0, 2.,
        1, 0, 0, 2.
    );
    
    std::cout << "I = " << I << std::endl;
    
    return 0;
}
