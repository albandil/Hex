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

#include "arrays.h"
#include "basis.h"
#include "gausskronrod.h"
#include "matrix.h"
#include "potential.h"
#include "specf.h"

PotentialMatrix::PotentialMatrix (
    LaguerreBasis const & basis,
    QuadratureRule const & quadrature,
    int J, int S, int Pi
) : basis_(basis), quadrature_(quadrature), matrix_(quadrature.nodes().size())
{
    // matrix indices
    size_t irow = 0, icol = 0;
    
    //
    // Precompute the potential matrix.
    //
    
    // row iterations
//     # pragma omp parallel for collapse (8)
    for (int L = 0; L < basis_.size(); L++)
    for (int N = 1; N <= basis_.size(L); N++)
    for (int l = 0; l < basis_.size(); l++) // FIXME : maxpell
    for (double k : quadrature_.nodes(L, N, l))
    {
        // column iterations
        for (int Lp = 0; Lp < basis_.size(); Lp++)
        for (int Np = 1; Np <= basis_.size(Lp); Np++)
        for (int lp = 0; lp < basis_.size(); lp++) // FIXME : maxpell
        for (double kp : quadrature_.nodes(Lp, Np, lp))
        {
            // get lambda limits
            int lambdamin = 0; // TODO
            int lambdamax = 0; // TODO
            
            // for all multipoles
            for (int lambda = lambdamin; lambda <= lambdamax; lambda++)
            {
                double f = 0; // TODO
                matrix_(irow,icol) += f * ComputeVdir (lambda, L, N, l, k, Lp, Np, lp, kp);
            }
            
            // move to next column of the matrix
            icol++;
        }
        
        // move to next row of the matrix
        irow++;
    }
}

RowMatrix<double> const & PotentialMatrix::matrix () const
{
    return matrix_;
}

MatrixEquation::MatrixEquation (
    QuadratureRule const & quadrature,
    PotentialMatrix const & potential
){
    // TODO
}

rArray MatrixEquation::solve() const
{
    // TODO
    return rArray();
}

double PotentialMatrix::ComputeIdir (
    int lambda,
    int L, int i, int l, double k,
    int Lp, int ip, int lp, double kp
) const {
    
    // sanity check
    assert (lambda >= 0 and L >= 0 and Lp >= 0 and l >= 0 and lp >= 0);
    assert (i >= 1 and ip >= 1);
    assert (k > 0. and kp > 0.);
    
    // decide according to lambda where to integrate
    if (lambda == 0)
    {
        //
        // only consider r1 > r2
        //
        
        auto integrand = [&](double r1) -> double {
            
            // evaluate hydrogen states
            double xi  = basis_.basestate (L,  i,  r1);
            double xip = basis_.basestate (Lp, ip, r1);
            
            // evaluate inner integral
            auto iintegrand = [&](double r2) -> double {
                return (1./r1-1./r2) * ric_j(l,k*r2) * ric_j(lp,kp*r2);
            };
            GaussKronrod<decltype(iintegrand)> iintegrator(iintegrand);
            iintegrator.integrate(0., r1);
            
            // return result
            return xi * xip * iintegrator.result();
        };
        GaussKronrod<decltype(integrand)> integrator(integrand);
        integrator.integrate(0., Inf);
        return integrator.result();
    }
    else
    {
        //
        // split integration along the line r1 = r2
        //
        
        // first part, r1 > r2
        auto integrand1 = [&](double r1) -> double {
            
            // evaluate hydrogen states
            double xi  = basis_.basestate (L,  i,  r1);
            double xip = basis_.basestate (Lp, ip, r1);
            
            // evaluate inner integral
            auto iintegrand = [&](double r2) -> double {
                return pow(r2/r1,lambda) * ric_j(l,k*r2) * ric_j(lp,kp*r2);
            };
            GaussKronrod<decltype(iintegrand)> iintegrator(iintegrand);
            iintegrator.integrate(0., r1);
            
            // return result
            return xi * xip * iintegrator.result()/r1;
        };
        GaussKronrod<decltype(integrand1)> integrator1(integrand1);
        integrator1.integrate(0., Inf);
        
        // second part, r1 < r2
        auto integrand2 = [&](double r1) -> double {
            
            // evaluate hydrogen states
            double xi  = basis_.basestate (L,  i,  r1);
            double xip = basis_.basestate (Lp, ip, r1);
            
            // evaluate inner integral
            auto iintegrand = [&](double r2) -> double {
                return pow(r1/r2,lambda) * ric_j(l,k*r2) * ric_j(lp,kp*r2) / r2;
            };
            GaussKronrod<decltype(iintegrand)> iintegrator(iintegrand);
            iintegrator.integrate(r1, Inf);
            
            // return result
            return xi * xip * iintegrator.result();
        };
        GaussKronrod<decltype(integrand2)> integrator2(integrand2);
        integrator2.integrate(0., Inf);
        return integrator1.result() + integrator2.result();
    }
}

double PotentialMatrix::ComputeVdir (
    int lambda,
    int L, int N, int l, double k,
    int Lp, int Np, int lp, double kp
) const {
    
    // sanity check
    assert (L >= 0 and Lp >= 0 and l >= 0 and lp >= 0);
    assert (N >= 1 and Np >= 1);
    assert (k >= 0. and kp >= 0.);
    
    // evaluate basic integrals
    RowMatrix<double> I (basis_.size(L), basis_.size(Lp));
    I.populate (
        [&](int i, int ip) -> double {
            return ComputeIdir (lambda, L, i, l, k, Lp, ip, lp, kp);
        }
    );
    
    // get orbital expansions
    rArrayView P = basis_.orbital(L, N);
    rArrayView Pp = basis_.orbital(Lp, Np);
    
    // sum expansions using bilinear form
    return I (P,Pp);
}
