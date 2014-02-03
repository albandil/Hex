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

#include "arrays.h"
#include "basis.h"
#include "potential.h"
#include "quadrature.h"
#include "version.h"

int main (int argc, char *argv[])
{
    // write the program header
    std::cout << logo_raw();
    
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
    rArray lambda (maxell + 1);
    lambda.fill(2.);
    
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
