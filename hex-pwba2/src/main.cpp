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

#include <gsl/gsl_errno.h>

#include "complex.h"
#include "radial.h"
#include "version.h"

int main (int argc, char* argv[])
{
    std::cout << logo() << std::endl;
    
    gsl_set_error_handler_off();
    
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
