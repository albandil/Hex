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

#include <cmath>

#include "multipot.h"

MultipolePotential::MultipolePotential (int lambda, int Na, int La, int Nb, int Lb)
{
    // initialize
}

MultipolePotential::MultipolePotential (int lambda, double Ka, int La, int Nb, int Lb)
{
    // initialize
}

MultipolePotential::MultipolePotential (int lambda, int Na, int La, double Kb, int Lb)
{
    // initialize
}

double MultipolePotential::operator() (double x)
{
    // FIXME : implement also other potentials than 1s-1s & lambda = 0
    return -(1 + 1/x) * std::exp(-2*x);
}
