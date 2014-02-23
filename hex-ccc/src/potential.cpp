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
#include "gauss.h"
#include "gausskronrod.h"
#include "matrix.h"
#include "potential.h"
#include "symbolic.h"
#include "specf.h"

PotentialMatrix::PotentialMatrix
(
    LaguerreBasis const & basis,
    QuadratureRule const & quadrature,
    int J, int S, int Pi
) : basis_(basis), quadrature_(quadrature), matrix_(quadrature.nodes().size()), g_()
{
    //
    // Precompute symbolic representation of the Laguerre states.
    //
    
    Array<Array<symbolic::poly>> xi(basis_.size());
    Array<Array<double>> xiNsqr(basis_.size());
    for (int L = 0; L < basis_.size(); L++)
    {
        xi[L].resize(basis_.size(L));
        xiNsqr[L].resize(basis_.size(L));
        
        for (int i = 1; i <= basis_.size(L); i++)
        {
            xi[L][i-1] = symbolic::LaguerreBasisFunction(i, L, basis_.rat_lambda(L));
            xiNsqr[L][i-1] = symbolic::double_approx(symbolic::LaguerreBasisFunctionNsqr(i, L, basis_.rat_lambda(L)));
        }
    }
    
    //
    // Precompute the potential matrix.
    //
    
    std::cout << "Precomputing potential matrix...\n\n";
    
    // matrix indices
    size_t irow = 0, icol = 0;
    
    // row iterations
//     # pragma omp parallel for collapse (3)
    for (int L = 0; L < basis_.size(); L++)
    for (int l = 0; l < basis_.size(); l++) // FIXME : maxpell
    for (int i = 1; i <= basis_.size(L); i++)
    {
        const rArrayView k = quadrature_.nodes(L, l, i);
        for (unsigned ik = 0; ik < k.size(); ik++)
        {
            // column iterations
            for (int Lp = 0; Lp < basis_.size(); Lp++)
            for (int lp = 0; lp < basis_.size(); lp++) // FIXME : maxpell
            for (int ip = 1; ip <= basis_.size(Lp); ip++)
            {
                const rArrayView kp = quadrature_.nodes(Lp, lp, ip);
                for (unsigned ikp = 0; ikp < kp.size(); ikp++)
                {
                    // the matrix is symmetrical, so we can skip the lower half
                    if (irow > icol)
                    {
                        icol++;
                        continue;
                    }
                    
                    // get lambda limits
                    int lambdamin = std::max (std::abs(L-Lp), std::abs(l-lp));
                    int lambdamax = std::max (         L+Lp ,          l+lp );
                    
                    // for all multipoles
                    for (int lambda = lambdamin; lambda <= lambdamax; lambda++)
                    {
                        double f = computef (lambda, L, l, Lp, lp, J);
                        
                        // skip uncontributing multipoles
                        if (f == 0.)
                            continue;
                        
                        // compute the direct double integral
//                         double integral = ComputeIdir (lambda, L, i, l, k[ik], Lp, ip, lp, kp[ikp]);
                        double integral = ComputeJdir (lambda, xi[L][i-1], l, k[ik], xi[Lp][ip-1], lp, kp[ikp]) * sqrt(xiNsqr[L][i-1] * xiNsqr[Lp][ip-1]);
                        
                        // add contribution of this multipole to the matrix
                        matrix_(irow,icol) += f * integral;
                    }
                    
                    /// DEBUG
                    if (lambdamin <= lambdamax)
                    {
                        std::cout << "\tV(" << irow << "," << icol << ") = " << matrix_(irow,icol) << "\n";
                    }
                    
                    // mirror the element across the diagonal
                    matrix_(icol,irow) = matrix_(irow,icol);
                    
                    // move to next column of the matrix
                    icol++;
                }
            }
            
            // move to next row of the matrix
            icol = 0;
            irow++;
        }
    }
}

RowMatrix<double> const & PotentialMatrix::matrix () const
{
    return matrix_;
}

MatrixEquation::MatrixEquation
(
    QuadratureRule const & quadrature,
    PotentialMatrix const & potential
)
{
    // TODO
}

rArray MatrixEquation::solve () const
{
    // TODO
    return rArray();
}

/*
double PotentialMatrix::ComputeIdir
(
    int lambda,
    symbolic::poly const & j,
    symbolic::poly const & jp,
    symbolic::poly const & P,
    symbolic::poly const & Pp
)
{
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
*/

/*
double PotentialMatrix::ComputeIdir
(
    int lambda,
    int L, int i, int l, double k,
    int Lp, int ip, int lp, double kp
) const
{
    // sanity check
    assert (lambda >= 0 and L >= 0 and Lp >= 0 and l >= 0 and lp >= 0);
    assert (i >= 1 and ip >= 1);
    assert (k > 0. and kp > 0.);
    
    // order of quadrature rule
    int n = 50;
    
    // get quadrature points
    rArray xs = g_.nodes(n, 0., 1.);
    rArray ws = g_.weights(n, 0., 1.);
    
    // approximate integral
    double result = 0;
    
    if (lambda == 0)
    {
        auto integrand = [&](double u, double t) -> double
        {
            double v = t * u;
            
            double Jac = u / ((1.-u)*(1.-u)*(1.-v)*(1.-v));
            
            double r1 = u / (1. - u);
            double r2 = v / (1. - v);
            
            double xi = basis_.basestate(L,i,r1);
            double xip = basis_.basestate(Lp,ip,r1);
            double j = ric_j(l,k*r2);
            double jp = ric_j(lp,kp*r2);
            
            return xi * xip * (1./r1 - 1./r2) * j * jp * Jac;
        };
        
        for (int iu = 0; iu < n; iu++)
        for (int it = 0; it < n; it++)
            result += ws[iu] * ws[it] * integrand(xs[iu],xs[it]);
    }
    else
    {
        auto integrand1 = [&](double u, double t) -> double
        {
            double v = t * u;
            
            double Jac = u / ((1.-u)*(1.-u)*(1.-v)*(1.-v));
            
            double r1 = u / (1. - u);
            double r2 = v / (1. - v);
            
            double xi = basis_.basestate(L,i,r1);
            double xip = basis_.basestate(Lp,ip,r1);
            double j = ric_j(l,k*r2);
            double jp = ric_j(lp,kp*r2);
            
            return xi * xip * pow(r2/r1,lambda)/r1 * j * jp * Jac;
        };
        
        for (int iu = 0; iu < n; iu++)
        for (int it = 0; it < n; it++)
            result += ws[iu] * ws[it] * integrand1(xs[iu],xs[it]);
        
        auto integrand2 = [&](double t, double v) -> double
        {
            double u = t * v;
            
            double Jac = v / ((1.-u)*(1.-u)*(1.-v)*(1.-v));
            
            double r1 = u / (1. - u);
            double r2 = v / (1. - v);
            
            double xi = basis_.basestate(L,i,r1);
            double xip = basis_.basestate(Lp,ip,r1);
            double j = ric_j(l,k*r2);
            double jp = ric_j(lp,kp*r2);
            
            return xi * xip * pow(r1/r2,lambda)/r2 * j * jp * Jac;
        };
        
        for (int iu = 0; iu < n; iu++)
        for (int it = 0; it < n; it++)
            result += ws[iu] * ws[it] * integrand2(xs[iu],xs[it]);
    }
    
    return result;
}
*/

double PotentialMatrix::ComputeIdir
(
    int lambda,
    int L, int i, int l, double k,
    int Lp, int ip, int lp, double kp
) const
{
    if (lambda == 0)
    {
        //
        // r1 > r2
        //
        
        auto integrand = [&](double u) -> double
        {
            double r1 = u / (1. - u);
            
            double xi  = basis_.basestate(L, i, r1);
            double xip = basis_.basestate(Lp, ip, r1);
            
            auto iintegrand = [&](double v) -> double
            {
                double r2 = v / (1. - v);
                return (1./r1 - 1./r2) * ric_j(l,k*r2) * ric_j(lp,kp*r2) / ((1.-v)*(1.-v));
            };
            
            GaussKronrod<decltype(iintegrand)> Qi(iintegrand);
            Qi.setEpsAbs(1e-10);
            Qi.setEpsRel(1e-6);
            Qi.integrate(0., u);
            
            return xi * xip * Qi.result() / ((1.-u)*(1.-u));
        };
        
        GaussKronrod<decltype(integrand)> Q(integrand);
        Q.setEpsAbs(1e-10);
        Q.setEpsRel(1e-5);
        Q.integrate(0., 1.);
        
        return Q.result();
    }
    else
    {
        //
        // r1 > r2
        //
        
        auto integrand1 = [&](double u) -> double
        {
            double r1 = u / (1. - u);
            
            double xi  = basis_.basestate(L, i, r1);
            double xip = basis_.basestate(Lp, ip, r1);
            
            auto iintegrand1 = [&](double v) -> double
            {
                double r2 = v / (1. - v);
                return  pow(r2/r1,lambda)/r1 * ric_j(l,k*r2) * ric_j(lp,kp*r2) / ((1.-v)*(1.-v));
            };
            
            GaussKronrod<decltype(iintegrand1)> Qi(iintegrand1);
            Qi.integrate(0., u);
            
            return xi * xip * Qi.result() / ((1.-u)*(1.-u));
        };
        
        GaussKronrod<decltype(integrand1)> Q1(integrand1);
        Q1.integrate(0., 1.);
        
        //
        // r1 < r2
        //
        
        auto integrand2 = [&](double u) -> double
        {
            double r1 = u / (1. - u);
            
            double xi  = basis_.basestate(L, i, r1);
            double xip = basis_.basestate(Lp, ip, r1);
            
            auto iintegrand2 = [&](double v) -> double
            {
                double r2 = v / (1. - v);
                return  pow(r1/r2,lambda)/r2 * ric_j(l,k*r2) * ric_j(lp,kp*r2) / ((1.-v)*(1.-v));
            };
            
            GaussKronrod<decltype(iintegrand2)> Qi(iintegrand2);
            Qi.integrate(u, 1.);
            
            return xi * xip * Qi.result() / ((1.-u)*(1.-u));
        };
        
        GaussKronrod<decltype(integrand2)> Q2(integrand2);
        Q2.integrate(0., 1.);
        
        return Q1.result() + Q2.result();
    }
}

double PotentialMatrix::ComputeJdir
(
    int lambda,
    symbolic::poly const & xi, int l, double k,
    symbolic::poly const & xip, int lp, double kp
) const
{
    symbolic::poly xixip = xi * xip;
    
    if (lambda == 0)
    {
        //          ∞
        //         ⌠             ⎛1   1⎞
        //  φ(x) = |  P(y) P'(y) ⎜- - -⎟ dy
        //         ⌡             ⎝y   x⎠
        //        x
        
        // second term
        symbolic::poly iintegrand = xixip;
        symbolic::poly integral_2 = symbolic::integrate_inf(iintegrand);
        for (symbolic::term & p : integral_2)
        {
            p.a--;
            p.kr = -p.kr;
        }
        
        // first term
        for (symbolic::term & p : iintegrand)
        {
            p.a--;
        }
        symbolic::poly integral_1 = symbolic::integrate_inf(iintegrand);
        symbolic::poly integral = integral_1 + integral_2;
        
        //          ∞
        //         ⌠
        //  Idir = |  j(x) j'(x) φ(x) dx
        //         ⌡
        //        0
        
        auto integrand = [&](double x) -> double
        {
            double j = ric_j(l,k*x);
            double jp = ric_j(lp,kp*x);
            
            return symbolic::eval(integral, x) * j * jp;
        };
        
        GaussKronrod<decltype(integrand)> Q(integrand);
        Q.integrate(0., Inf);
        
        if (not Q.ok())
        {
            std::cerr << format
            (
                "ComputeJdir failed for λ=%d, xi=%s, l=%d, k=%g, xip=%s, lp=%d, kp=%g (\"%s\").\n",
                lambda, symbolic::tostring(xi).c_str(), l, k, symbolic::tostring(xip).c_str(), lp, k, Q.status().c_str()
            );
        }
        
        return Q.result();
    }
    else
    {
        //          ∞
        //         ⌠              -λ-1
        //  φ(x) = |  P(y) P'(y) y     dy
        //         ⌡
        //        x
        
        symbolic::poly iintegrand_1 = xixip;
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
        
        symbolic::poly iintegrand_2 = xixip;
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
        
        symbolic::poly phipsi = phi + psi;
        
        // collect exponential factor
        // TODO
        
        auto integrand = [&](double x) -> double
        {
            double j = ric_j(l,k*x);
            double jp = ric_j(lp,kp*x);
            
            return symbolic::eval(phipsi, x) * j * jp;
        };
        GaussKronrod<decltype(integrand)> Q(integrand);
        Q.integrate(0., Inf);
        
        if (not Q.ok())
        {
            std::cerr << format
            (
                "ComputeJdir failed for λ=%d, xi=%s, l=%d, k=%g, xip=%s, lp=%d, kp=%g (\"%s\").\n",
                lambda, symbolic::tostring(xi).c_str(), l, k, symbolic::tostring(xip).c_str(), lp, k, Q.status().c_str()
            );
        }
        
        return Q.result();
    }
}
