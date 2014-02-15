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
#include "symbolic.h"
#include "specf.h"

PotentialMatrix::PotentialMatrix (
    LaguerreBasis const & basis,
    QuadratureRule const & quadrature,
    int J, int S, int Pi
) : basis_(basis), quadrature_(quadrature), matrix_(quadrature.nodes().size())
{
    //
    // Setup hydrogen radial orbitals.
    //
    
    Array<Array<symbolic::poly>> P(basis_.size());
    Array<Array<symbolic::term>> Nsqr(basis_.size());
    
    for (int L = 0; L < basis_.size(); L++)
    {
        P[L].resize(basis_.size(L));
        Nsqr[L].resize(basis_.size(L));
        
        for (int N = 1; N <= basis_.size(L); N++)
        {
            // get hydrogen function and its squared normalization factor
            P[L][N-1] = symbolic::LaguerreBasisFunction(N,L,basis_.rat_lambda(L));
            Nsqr[L][N-1] = symbolic::LaguerreBasisFunctionNsqr(N,L,basis_.rat_lambda(L));
        }
    }
    
    //
    // Setup Riccati-Bessel functions.
    //
    
    Array<Array<Array<Array<symbolic::poly>>>> j(basis_.size());
    
    for (int L = 0; L < basis_.size(); L++)
    {
        j[L].resize(basis_.size()); // FIXME : maxpell
        
        for (int l = 0; l < basis_.size(); l++) // FIXME : maxpell
        {
            j[L][l].resize(basis_.size(L));
            
            for (int N = 1; N <= basis_.size(L); N++)
            {
                const rArrayView k = quadrature_.nodes(L, l, N);
                j[L][l][N-1].resize(k.size());
                
                for (unsigned ik = 0; ik < k.size(); ik++)
                    j[L][l][N-1][ik] = symbolic::RiccatiBessel(l,k[ik]);
            }
        }
    }
    
    //
    // Precompute the potential matrix in Laguerre basis.
    //
    
    // matrix indices
    size_t irow = 0, icol = 0;
    
    // row iterations
    for (int L = 0; L < basis_.size(); L++)
    for (int l = 0; l < basis_.size(); l++) // FIXME : maxpell
    for (int N = 1; N <= basis_.size(L); N++)
    {
        const rArrayView k = quadrature_.nodes(L, l, N);
        for (unsigned ik = 0; ik < k.size(); ik++)
        {
            // column iterations
            for (int Lp = 0; Lp < basis_.size(); Lp++)
            for (int lp = 0; lp < basis_.size(); lp++) // FIXME : maxpell
            for (int Np = 1; Np <= basis_.size(Lp); Np++)
            {
                const rArrayView kp = quadrature_.nodes(Lp, lp, Np);
                for (unsigned ikp = 0; ikp < kp.size(); ikp++)
                {
                    // the matrix is symmetrical, so we can skip the lower half
                    if (irow > icol)
                        continue;
                    
                    // get lambda limits
                    int lambdamin = 0; // TODO
                    int lambdamax = 0; // TODO
                    
                    // for all multipoles
                    for (int lambda = lambdamin; lambda <= lambdamax; lambda++)
                    {
                        double f = 0; // TODO
                        
                        std::cout << "Double integral\n"
                                  << "   L = " << L  << " l = " << l  << " N = " << N  << " k = "  << k[ik]   << "\n"
                                  << "   L'= " << Lp << " l'= " << lp << " N'= " << Np << " k' = " << kp[ikp] << "\n"
                                  << "   j : " << j[L][l][N-1][ik] << "\n"
                                  << "   j': " << j[Lp][lp][Np-1][ikp] << "\n"
                                  << "   P : " << P[L][N-1] << "\n"
                                  << "   P': " << P[Lp][Np-1] << "\n";
                                  
                        double integral = ComputeIdir (
                            lambda,
                            j[L][l][N-1][ik],
                            j[Lp][lp][Np-1][ikp],
                            P[L][N-1],
                            P[Lp][Np-1]
                        );
                        
                        std::cout << " -> result " << integral << "\n";
                        
                        matrix_(irow,icol) += f * integral;
                    }
                    
                    // mirror the element across the diagonal
                    matrix_(icol,irow) = matrix_(irow,icol);
                    
                    // move to next column of the matrix
                    icol++;
                }
            }
        
            // move to next row of the matrix
            irow++;
        }
    }
    
    //
    // Transform the matrix to eigenstates.
    //
    
    // TODO
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
    symbolic::poly const & j,
    symbolic::poly const & jp,
    symbolic::poly const & P,
    symbolic::poly const & Pp
){
    if (lambda == 0)
    {
        //          ∞
        //         ⌠             ⎛1   1⎞
        //  φ(x) = |  P(y) P'(y) ⎜- - -⎟ dy
        //         ⌡             ⎝y   x⎠
        //        x
        
        // second term
        symbolic::poly integrand = P * Pp;
        symbolic::poly integral_2 = symbolic::integrate_inf(integrand);
        for (symbolic::term & p : integral_2)
        {
            p.a--;
            p.kr = -p.kr;
        }
        
        // first term
        for (symbolic::term & p : integrand)
        {
            p.a--;
        }
        symbolic::poly integral_1 = symbolic::integrate_inf(integrand);
        symbolic::poly integral = integral_1 + integral_2;
        
        //          ∞
        //         ⌠
        //  Idir = |  j(x) j'(x) φ(x) dx
        //         ⌡
        //        0
        
        integrand = integral * j * jp;
        return symbolic::integrate_full(integrand).todouble();
    }
    else
    {
        //          ∞
        //         ⌠              -λ-1
        //  φ(x) = |  P(y) P'(y) y     dy
        //         ⌡
        //        x
        
        symbolic::poly iintegrand_1 = P * Pp;
        for (symbolic::term & p : iintegrand_1)
        {
            p.a -= lambda + 1;
        }
        symbolic::poly phi = symbolic::integrate_inf(iintegrand_1);
        
        //          x
        //         ⌠              λ
        //  ψ(x) = |  P(y) P'(y) y  dy
        //         ⌡
        //        0
        
        symbolic::poly iintegrand_2 = P * Pp;
        for (symbolic::term & p : iintegrand_2)
        {
            p.a += lambda;
        }
        symbolic::poly psi = symbolic::integrate_low(iintegrand_2);
        
        //          ∞
        //         ⌠             ⎛      λ          -λ-1⎞
        //  Idir = |  j(x) j'(x) ⎜φ(x) x  +  ψ(x) x    ⎟ dx
        //         ⌡             ⎝                     ⎠
        //        0
        
        for (symbolic::term & p : phi)
        {
            p.a += lambda;
        }
        for (symbolic::term & p : psi)
        {
            p.a -= lambda + 1;
        }
        symbolic::poly integrand = j * jp * (phi + psi);
        return symbolic::integrate_full(integrand).todouble();
    }
}

/*
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
            double xxx = xi * xip * iintegrator.result();
            std::cout << xxx << std::endl;
            return xxx;
        };
        GaussKronrod<decltype(integrand)> integrator(integrand);
        integrator.integrate(0., Inf);
        
        if (not integrator.ok())
        {
            throw exception
            (
                "ComputeIdir failed (%g) for L=%d i=%d l=%d k=%g Lp=%d ip=%d lp=%d kp=%g\nStatus: %s",
                 integrator.result(), L, i, l, k, Lp, ip, lp, kp, integrator.status().c_str()
            );
        }
        
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
        [&](int i, int ip) -> double
        {
            int N = i + 1, Np = ip + 1;
            return ComputeIdir (lambda, L, N, l, k, Lp, Np, lp, kp);
        }
    );
    
    // get orbital expansions
    rArrayView P = basis_.orbital(L, N);
    rArrayView Pp = basis_.orbital(Lp, Np);
    
    // sum expansions using bilinear form
    return I (P,Pp);
}
*/