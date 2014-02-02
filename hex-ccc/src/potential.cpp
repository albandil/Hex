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

#include "basis.h"
#include "potential.h"

PotentialMatrix::PotentialMatrix (
    LaguerreBasis const & basis, QuadratureRule const & quadrature,
    int J, int S, int Pi
){
    // maximal angular momentum of the propagator
    int maxpell = basis.size() - 1;
    
    
}

RowMatrix<double> const & PotentialMatrix::matrix () const
{
    return matrix_;
}

MatrixEquation::MatrixEquation (
    QuadratureRule const & quadrature, PotentialMatrix const & potential
){
    // TODO
}

rArray MatrixEquation::solve() const
{
    return rArray();
}

double compute_Vdir (
    LaguerreBasis const & basis, int lambda,
    double k, int l, rArray const & PNL, int L,
    double kp, int lp, rArray const & PNLp, int Lp
){
    // TODO
    return 0;
}

double compute_Vexc (
    LaguerreBasis const & basis, int lambda,
    double k, int l, rArray const & PNL, int L,
    double kp, int lp, rArray const & PNLp, int Lp
){
    // TODO
    return 0;
}
