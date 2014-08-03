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

#include <cstdlib>
#include <iostream>

#include <gsl/gsl_errno.h>

#include "arrays.h"
#include "basis.h"
#include "cmdline.h"
#include "potential.h"
#include "quadrature.h"
#include "symbolic.h"
#include "version.h"

int main (int argc, char *argv[])
{
    // write the program header
    std::cout << logo(" ");
    
    // turn off GSL and HDF exceptions
    gsl_set_error_handler_off();
    
    // disable buffering of the standard output (-> immediate logging)
    setvbuf(stdout, nullptr, _IONBF, 0);
    
    // get input from command line
    ParseCommandLine
    (
        argc, argv,
        [&](std::string optname, std::string optarg) -> bool { return false; }
    );
    
    //
    //  Initialization
    //
    
    // input data
    int J = 0;      // total angular momentum
    int S = 0;      // total spin
    int Pi = 0;     // total parity
    int E = 3.5;    // total energy in Rydbergs
    int maxell = 2; // single-electron angular momentum limit
    
    // Laguerre basis sizes
    iArray Nl (maxell + 1);
    Nl.fill(30);
    
    // screening constants
    Array<symbolic::rational> lambda (maxell + 1);
    lambda.fill(symbolic::rational(2));
    
    //
    //  The main body of the computation
    //
    
    // create the Laguerre bases
    LaguerreBasis basis (maxell, Nl, lambda);

    // setup the quadrature nodes and weights
    QuadratureRule quadrature (basis, E);
    
    // compute the matrix of the potential
    PotentialMatrix potential (basis, quadrature, J, S, Pi);
    
    // construct the matrix equation
    MatrixEquation eqn (quadrature, potential);
    
    // solve the matrix equation
    rArray K_matrix = eqn.solve();
    
    //
    //  Evaluate T-matrices from the K-matrices
    //
    
    // TODO
    
    return EXIT_SUCCESS;
}
