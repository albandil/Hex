//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2014, Jakub Benda, Charles University in Prague                    //
//                                                                                   //
// MIT License:                                                                      //
//                                                                                   //
//  Permission is hereby granted, free of charge, to any person obtaining a          //
// copy of this software and associated documentation files (the "Software"),        //
// to deal in the Software without restriction, including without limitation         //
// the rights to use, copy, modify, merge, publish, distribute, sublicense,          //
// and/or sell copies of the Software, and to permit persons to whom the             //
// Software is furnished to do so, subject to the following conditions:              //
//                                                                                   //
//  The above copyright notice and this permission notice shall be included          //
// in all copies or substantial portions of the Software.                            //
//                                                                                   //
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS          //
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       //
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE       //
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, //
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF         //
// OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.  //
//                                                                                   //
//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "arrays.h"
#include "complex.h"
#include "gausskronrod.h"
#include "hydrogen.h"
#include "misc.h"
#include "nodeintegrate.h"
#include "radial.h"
#include "special.h"

rArray interpolate_bound_bound_potential
(
    rArray const & x,
    int lambda,
    int Na, int La, int Nb, int Lb
)
{
    // output array
    unsigned N = x.size();
    rArray V (N);
    
    // combined exponential factor
    double c = 1./Na + 1./Nb;
    
    // combined normalization factor
    double Norm = std::sqrt
    (
        std::pow(2./Na,3) * gsl_sf_fact(Na-La-1) / (2 * Na * gsl_sf_fact(Na + La)) *
        std::pow(2./Nb,3) * gsl_sf_fact(Nb-Lb-1) / (2 * Nb * gsl_sf_fact(Nb + Lb))
    );
    
    // compute the integrals (symbolically)
    if (lambda == 0)
    {
        auto potential = [&](double y) -> double
        {
            if (y == 0)
                return 0; // actually Inf, but...
            
            // for all terms of product of Laguerre polynomials
            double suma1 = 0, suma2 = 0;
            for (int ia = 0; ia <= Na - La - 1; ia++)
            for (int ib = 0; ib <= Nb - Lb - 1; ib++)
            {
                double afactor = (ia % 2 == 0 ? 1. : -1.) * std::pow(2./Na,La+ia) * gsl_sf_choose(Na+La,Na-La-1-ia) / gsl_sf_fact(ia);
                double bfactor = (ib % 2 == 0 ? 1. : -1.) * std::pow(2./Nb,Lb+ib) * gsl_sf_choose(Nb+Lb,Nb-Lb-1-ib) / gsl_sf_fact(ib);
                int rpower1 = La + 1 + Lb + 1 + ia + ib - 1;
                int rpower2 = La + 1 + Lb + 1 + ia + ib;
                double integral1 = gsl_sf_gamma(rpower1 + 1) * gsl_sf_gamma_inc_Q(rpower1 + 1, c*y);
                double integral2 = gsl_sf_gamma(rpower2 + 1) * gsl_sf_gamma_inc_Q(rpower2 + 1, c*y);
                suma1 += afactor * bfactor * integral1 / std::pow(c, rpower1 + 1);
                suma2 += afactor * bfactor * integral2 / std::pow(c, rpower2 + 1);
            }
            
            return suma1 - suma2 / y;
        };
        
        for (unsigned i = 0; i < N; i++)
            V[i] = Norm * potential(x[i]);
    }
    else
    {
        auto potential = [&](double y) -> double
        {
            if (y == 0)
                return 0; // actually Inf, but...
            
            // for all terms of product of Laguerre polynomials
            double suma1 = 0, suma2 = 0;
            for (int ia = 0; ia <= Na - La - 1; ia++)
            for (int ib = 0; ib <= Nb - Lb - 1; ib++)
            {
                double afactor = (ia % 2 == 0 ? 1. : -1.) * std::pow(2./Na,La+ia) * gsl_sf_choose(Na+La,Na-La-1-ia) / gsl_sf_fact(ia);
                double bfactor = (ib % 2 == 0 ? 1. : -1.) * std::pow(2./Nb,Lb+ib) * gsl_sf_choose(Nb+Lb,Nb-Lb-1-ib) / gsl_sf_fact(ib);
                int rpower1 = La + 1 + Lb + 1 + ia + ib + lambda;
                int rpower2 = La + 1 + Lb + 1 + ia + ib - lambda - 1;
                double integral1 = gsl_sf_gamma(rpower1 + 1) * gsl_sf_gamma_inc_P(rpower1 + 1, c*y);
                double integral2 = gsl_sf_gamma(rpower2 + 1) * gsl_sf_gamma_inc_Q(rpower2 + 1, c*y);
                suma1 += afactor * bfactor * integral1 / std::pow(c, rpower1 + 1);
                suma2 += afactor * bfactor * integral2 / std::pow(c, rpower2 + 1);
            }
            
            return suma1 * std::pow(y, -lambda-1) + suma2 * std::pow(y, lambda);
        };
        
        for (unsigned i = 0; i < N; i++)
            V[i] = Norm * potential(x[i]);
    }
    
    return V;
}

rArray interpolate_bound_free_potential
(
    rArray const & x,
    int lambda,
    int Na, int La, double Kb, int Lb
)
{
    // array of bound-free potential evaluations
    unsigned N = x.size();
    double rmax = x.back();
    rArray V(N);
    
    // number of Coulomb zeros within the range of grid "x" (initial guess)
    unsigned nzeros = Kb * x.back() / (2 * special::constant::pi);
    
    // calculate the Coulomb zeros
    rArray zeros;
    do
    {
        // double the necessary zero count
        nzeros = 2 * (1 + nzeros);
        zeros.resize(nzeros);
        
        // for large frequencies 'Kb' the integral will be almost zero, so lets check already here
        if (nzeros > N / 2)
            return rArray(N,0.);
        
        // compute zeros of the Coulomb function
        special::coulomb_zeros(-1/Kb,Lb,nzeros,zeros.data());
        zeros /= Kb;
    }
    while (zeros.back() < x.back());
    
    // the classical turning point for non-S-waves
    double rt = (std::sqrt(1 + Kb*Kb*Lb*(Lb+1))-1)/(Kb*Kb);
    
    // Coulomb function scaled by hydrogen orbital
    rArray PF(N);
    for (unsigned i = 0; i < N; i++)
    {
        PF[i] = Hydrogen::P(Na,La,x[i]) * Hydrogen::F(Kb,Lb,x[i]);
    }
    gsl_spline * spline = gsl_spline_alloc (gsl_interp_cspline, N);
    gsl_spline_init (spline, x.data(), PF.data(), N);
    
    // compute the integrals
    if (lambda == 0)
    {
        // setup accelerator for value search (for every thread)
        gsl_interp_accel * acc = gsl_interp_accel_alloc ();
        
        // evaluate potential at grid point 'y'
        auto potential = [Na,La,Kb,Lb,N,rt,rmax,zeros,spline,&acc](double y) -> double
        {
            // use zero potential in the origin (should be Inf, but that would spreads NaNs in PC arithmeric)
            if (y == 0.)
                return 0;
            
            // integrand Pa(r₁) V(r₁,r₂) Fb(r₁) for λ = 0
            auto integrand = [Na,La,y,Kb,Lb,spline,&acc](double x) -> double
            {
                return gsl_spline_eval(spline, x, acc) * (y - x) / (y * x);
            };
            
            // integrate per node of the Coulomb function (precomputed before)
            FixedNodeIntegrator<decltype(integrand),GaussKronrod<decltype(integrand)>,double> Q(integrand,zeros,rt);
            Q.setEpsAbs(0);
            Q.integrate(y, rmax);
            if (not Q.ok())
            {
                std::cerr << format
                (
                    "Bound-free potential V[0]{%d,%d->%g,%d} integration failed for r = %g (\"%s\").",
                    Na, La, Kb, Lb, y, Q.status().c_str()
                ) << std::endl;
                std::terminate();
            }
            return Q.result();
        };
        
        // evaluate potential at all grid points
        for (unsigned i = 0; i < N; i++)
            V[i] = potential(x[i]);
        
        // release memory
        gsl_interp_accel_free(acc);
    }
    else
    {
        // setup accelerator for value search (for every thread)
        gsl_interp_accel * acc = gsl_interp_accel_alloc ();
        
        // evaluate potential at grid point 'y'
        auto potential = [Na,La,Kb,Lb,lambda,zeros,rt,rmax,x,spline,&acc](double y) -> double
        {
            // use zero potential in the origin (should be Inf, but that would spread NaNs in PC arithmeric)
            if (y == 0.)
                return 0;
            
            // integrand Pa(r₁) V(r₁,r₂) Fb(r₁) for λ ≠ 0 and r₁ < r₂
            auto integrand1 = [Na,La,Kb,Lb,lambda,y,spline,&acc](double x) -> double
            {
                return gsl_spline_eval(spline, x, acc) * gsl_sf_pow_int(x/y,lambda);
            };
            
            // integrate per node of the Coulomb function (precomputed before)
            FixedNodeIntegrator<decltype(integrand1),GaussKronrod<decltype(integrand1)>,double> Q1(integrand1,zeros,rt);
            Q1.setEpsAbs(0);
            Q1.integrate(0., y);
            if (not Q1.ok())
            {
                std::cerr << format
                (
                    "Bound-free potential V[%d]{%d,%d->%g,%d} [0,%g] integration failed (\"%s\").",
                    lambda, Na, La, Kb, Lb, y, Q1.status().c_str()
                ) << std::endl;
                std::terminate();
            }
            
            // integrand Pa(r₁) V(r₁,r₂) Fb(r₁) for λ ≠ 0 and r₁ > r₂
            auto integrand2 = [&](double x) -> double
            {
                return gsl_spline_eval(spline, x, acc) * gsl_sf_pow_int(y/x,lambda + 1);
            };
            
            // integrate per node of the Coulomb function (precomputed before)
            FixedNodeIntegrator<decltype(integrand2),GaussKronrod<decltype(integrand2)>,double> Q2(integrand2,zeros,rt);
            Q2.setEpsAbs(0);
            Q2.integrate(y, rmax);
            if (not Q2.ok())
            {
                std::cerr << format
                (
                    "Bound-free potential V[%d]{%d,%d->%g,%d} [%g,inf] integration failed (\"%s\").",
                    lambda, Na, La, Kb, Lb, y, Q2.status().c_str()
                ) << std::endl;
                std::terminate();
            }
            
            return (Q1.result() + Q2.result()) / y;
        };
        
        // evaluate potential in all grid points
        for (unsigned i = 0; i < N; i++)
            V[i] = potential(x[i]);
        
        // release memory
        gsl_interp_accel_free(acc);
    }
    
    // free allocated memory
    gsl_spline_free(spline);
    
    // return the array of evaluations
    return V;
}

rArray interpolate_free_bound_potential
(
    rArray const & x,
    int lambda,
    double Ka, int La, int Nb, int Lb
)
{
    // the potential is real, so just return the transpose
    return interpolate_bound_free_potential (x, lambda, Nb, Lb, Ka, La);
}

rArray interpolate_riccati_bessel_j
(
    rArray const & x,
    int l, double k
)
{
    return x.transform
    (
        [&](double r) -> double
        {
            return special::ric_j(l,k*r);
        }
    );
}

rArray interpolate_riccati_bessel_y
(
    rArray const & x,
    int l, double k
)
{
    return x.transform
    (
        [&](double r) -> double
        {
            if (k == 0. or r == 0.)
                return l == 0 ? 1. : 0.;
            else
                return -special::ric_n(l,k*r);
        }
    );
}

rArray interpolate_riccati_bessel_iscaled
(
    rArray const & x,
    int l, double k
)
{
    return x.transform
    (
        [&](double r) -> double
        {
            return special::ric_i_scaled(l,k*r);
        }
    );
}

rArray interpolate_riccati_bessel_kscaled
(
    rArray const & x,
    int l, double k
)
{
    return x.transform
    (
        [&](double r) -> double
        {
            if (k == 0. or r == 0.)
                return l == 0 ? 1 : 0;
            else
                return special::ric_k_scaled(l,k*r);
        }
    );
}

Complex Idir_allowed
(
    rArray const & grid,
    rArray const & jf, rArray const & Vfn,
    rArray const & jn, rArray const & yn_ji, rArray const & yn_jf,
    rArray const & ji, rArray const & Vni
)
{
    std::size_t N = grid.size();
    
    rArray inner_lower(N), inner_lower_int(N), inner_higher_re(N), inner_higher_re_int(N);
    double sum_outer_re = 0, sum_outer_im = 0;
    
    for (std::size_t i = 1; i < N; i++)
    {
        inner_lower[i] = Vni[i] * jn[i] * ji[i];
        inner_higher_re[i] = Vni[i] * yn_ji[i];
    }
    
    // low integral
    gsl_spline * spline = gsl_spline_alloc(gsl_interp_cspline, N);
    gsl_spline_init(spline, grid.data(), inner_lower.data(), N);
    gsl_interp_accel * acc = gsl_interp_accel_alloc();
    
    // compute forward partial sums
    for (std::size_t i = 0; i < N; i++)
        inner_lower_int[i] = gsl_spline_eval_integ(spline, grid.front(), grid[i], acc);
    
    gsl_interp_accel_free(acc);
    gsl_spline_free(spline);
    
    // high integral
    {
        gsl_spline * spline = gsl_spline_alloc(gsl_interp_cspline, N);
        gsl_spline_init(spline, grid.data(), inner_higher_re.data(), N);
        gsl_interp_accel * acc = gsl_interp_accel_alloc();
        
        // compute backward partial sums
        for (std::size_t i = 0; i < N; i++)
            inner_higher_re_int[i] = gsl_spline_eval_integ(spline, grid[i], grid.back(), acc);
        
        gsl_interp_accel_free(acc);
        gsl_spline_free(spline);
    }
    
    // outer integral (real and imaginary)
    {
        rArray outer_re = Vfn * (inner_lower_int * yn_jf + inner_higher_re_int * jf * jn);
    
        gsl_interp_accel * acc = gsl_interp_accel_alloc();
        gsl_spline * spline = gsl_spline_alloc(gsl_interp_cspline, N);
        gsl_spline_init(spline, grid.data(), outer_re.data(), N);
        
        sum_outer_re = gsl_spline_eval_integ(spline, grid.front(), grid.back(), acc);
        
        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);
    }
    {
        rArray outer_im = Vfn * jf * jn;
        
        gsl_interp_accel * acc = gsl_interp_accel_alloc();
        gsl_spline * spline = gsl_spline_alloc(gsl_interp_cspline, N);
        gsl_spline_init(spline, grid.data(), outer_im.data(), N);
        
        sum_outer_im = inner_lower_int.back() * gsl_spline_eval_integ(spline, grid.front(), grid.back(), acc);
        
        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);
    }
    
    return Complex (sum_outer_re, sum_outer_im);
}

double Idir_forbidden
(
    rArray const & grid,
    rArray const & jf, rArray const & Vfn,
    rArray const & iscaled_n, rArray const & kscaled_n, double kn,
    rArray const & ji, rArray const & Vni
)
{
    unsigned N = grid.size();
    double h = grid.back() / N;
    
    rArray fpart(N),  ipart(N);
    
    for (unsigned i = 0; i < N; i++)
    {
        fpart[i] = jf[i] * Vfn[i];
        ipart[i] = ji[i] * Vni[i];
    }
    
    double suma = 0;
    unsigned last_diagonal = N/2;
    
    // index "d" runs across the contours x - y = konst
    // TODO : use spline integration along and across the diagonals
    for (unsigned d = 0; d < last_diagonal; d++)
    {
        // contribution of these contours
        double contrib = 0;
        
        if (d == 0)
        {
            // index "i" runs along the diagonal [x,y] = [i,i]
            for (unsigned i = 0; i < N; i++)
                contrib += fpart[i] * ipart[i] * iscaled_n[i] * kscaled_n[i];
        }
        else
        {
            // index "i" runs along the contour [x,y] = [i,i-d]
            for (unsigned i = d; i < N; i++)
                contrib += fpart[i-d] * ipart[i] * iscaled_n[i-d] * kscaled_n[i];
            // index "i" runs along the contour [x,y] = [i,i+d]
            for (unsigned i = 0; i < N - d; i++)
                contrib += fpart[i+d] * ipart[i] * iscaled_n[i] * kscaled_n[i+d];
        }
        
        // add properly scaled contribution
        contrib *= std::exp(-kn*h*d);
        suma += contrib;
        
        // check convergence (do not blindly integrate everything)
        if (std::abs(contrib) < 1e-10 * std::abs(suma))
            last_diagonal = d;
    }
    
    return suma * h * h;
}

Complex Idir_nBound_allowed
(
    rArray const & grid, int L,
    int Nf, int Lf, double kf, int lf,
    int Nn, int Ln, double kn, int ln,
    int Ni, int Li, double ki, int li,
    std::ostream & log
)
{
    Complex result = 0;
    
    double nrm = std::sqrt(1. / kn);
    
    int lambdaf_min = std::max(std::abs(lf-ln), std::abs(Lf-Ln));
    int lambdaf_max = std::min(lf+ln, Lf+Ln);
    
    int lambdai_min = std::max(std::abs(li-ln), std::abs(Li-Ln));
    int lambdai_max = std::min(li+ln, Li+Ln);
    
    for (int lambdaf = lambdaf_min; lambdaf <= lambdaf_max; lambdaf++)
    for (int lambdai = lambdai_min; lambdai <= lambdai_max; lambdai++)
    {
        double ff = special::computef(lambdaf, Lf, lf, Ln, ln, L);
        double fi = special::computef(lambdai, Ln, ln, Li, li, L);
        
        // skip non-contributing transfers
        if (ff == 0. or fi == 0.)
            continue;
        
        // evaluate initial / final radial part of the projectile partial wave
        rArray ji = std::move(interpolate_riccati_bessel_j(grid, li, ki));
        rArray jf = std::move(interpolate_riccati_bessel_j(grid, lf, kf));
        
        // evaluate initial / final multipole potential
        rArray Vfn = std::move(interpolate_bound_bound_potential(grid, lambdaf, Nf, Lf, Nn, Ln));
        rArray Vni = std::move(interpolate_bound_bound_potential(grid, lambdai, Nn, Ln, Ni, Li));
        
        // evaluate the Green function terms (Bessel functions "jn" and "yn")
        // NOTE : Assuming kn > 0.
        rArray jn = interpolate_riccati_bessel_j(grid, ln, kn) * nrm;
        rArray yn = interpolate_riccati_bessel_y(grid, ln, kn) * nrm;
        
        // check finiteness (use asymptotic form if out)
        rArray yn_ji = yn * ji;
        for (unsigned i = 1; i < yn_ji.size(); i++)
        {
            if (not std::isfinite(yn_ji[i]))
                yn_ji[i] = gsl_sf_doublefact(2*ln-1) / gsl_sf_doublefact(2*li+1) * std::pow(kn,-ln) * std::pow(ki,li+1) * std::pow(grid[i],li+1-ln);
        }
        rArray yn_jf = yn * jf;
        for (unsigned i = 1; i < yn_jf.size(); i++)
        {
            if (not std::isfinite(yn_jf[i]))
                yn_jf[i] = gsl_sf_doublefact(2*ln-1) / gsl_sf_doublefact(2*lf+1) * std::pow(kn,-ln) * std::pow(kf,lf+1) * std::pow(grid[i],lf+1-ln);
        }
        
        // finally integrate
        Complex inte = -Idir_allowed (grid, jf, Vfn, jn, yn_ji, yn_jf, ji, Vni);
        
        // comment this result
        log << format
        (
            "\t\ttransfer [%d %d] initial (%d %d, %g %d), intermediate (%d %d, %g %d) final (%d %d, %g %d) : (%g,%g)",
            lambdai, lambdaf, Ni, Li, ki, li, Nn, Ln, kn, ln, Nf, Lf, kf, lf, inte.real(), inte.imag()
        ) << std::endl;
        
        // update result
        result += ff * fi * inte;
    }
    
    return result;
}

double Idir_nBound_forbidden
(
    rArray const & grid, int L,
    int Nf, int Lf, double kf, int lf,
    int Nn, int Ln, double kappan, int ln,
    int Ni, int Li, double ki, int li,
    std::ostream & log
)
{
    double result = 0;
    
    double nrm = std::sqrt(1. / kappan);
    
    int lambdaf_min = std::max(std::abs(lf-ln), std::abs(Lf-Ln));
    int lambdaf_max = std::min(lf+ln, Lf+Ln);
    
    int lambdai_min = std::max(std::abs(li-ln), std::abs(Li-Ln));
    int lambdai_max = std::min(li+ln, Li+Ln);
    
    for (int lambdaf = lambdaf_min; lambdaf <= lambdaf_max; lambdaf++)
    for (int lambdai = lambdai_min; lambdai <= lambdai_max; lambdai++)
    {
        double ff = special::computef(lambdaf, Lf, lf, Ln, ln, L);
        double fi = special::computef(lambdai, Ln, ln, Li, li, L);
        
        // skip non-contributing transfers
        if (ff == 0. or fi == 0.)
            continue;
        
        // evaluate initial / final radial part of the projectile partial wave
        rArray ji = std::move(interpolate_riccati_bessel_j(grid, li, ki));
        rArray jf = std::move(interpolate_riccati_bessel_j(grid, lf, kf));
        
        // evaluate initial / final multipole potential
        rArray Vfn = std::move(interpolate_bound_bound_potential(grid, lambdaf, Nf, Lf, Nn, Ln));
        rArray Vni = std::move(interpolate_bound_bound_potential(grid, lambdai, Nn, Ln, Ni, Li));
        
        // evaluate the Green function terms (Bessel functions "in" and "kn")
        // NOTE : Assuming kappan > 0.
        rArray iscaled_n = interpolate_riccati_bessel_iscaled(grid, ln, kappan) * nrm;
        rArray kscaled_n = interpolate_riccati_bessel_kscaled(grid, ln, kappan) * nrm;
        
        // finally integrate
        double inte = -Idir_forbidden(grid, jf, Vfn, iscaled_n, kscaled_n, kappan, ji, Vni);
        
        // comment this result
        log << format
        (
            "\t\ttransfer [%d %d] initial (%d %d, %g %d) intermediate (%d %d, %g %d) final (%d %d, %g %d) : (%g,%g)",
            lambdai, lambdaf, Ni, Li, ki, li, Nn, Ln, kappan, ln, Nf, Lf, kf, lf, inte, 0.
        ) << std::endl;
        
        // update result
        result += ff * fi * inte;
    }
    
    return result;
}

Complex Idir_nFree_allowed
(
    rArray const & grid, int L,
    int Nf, int Lf, double kf, int lf,
    double Kn, int Ln, double kn, int ln,
    int Ni, int Li, double ki, int li,
    std::ostream & log
)
{
    Complex result = 0;
    
    double nrm = std::sqrt(1. / kn);
    
    // check if kn is small enough for asymptotics
    bool kn_is_small = (kn * grid.back() < 1e-6);
    
    int lambdaf_min = std::max(std::abs(lf-ln), std::abs(Lf-Ln));
    int lambdaf_max = std::min(lf+ln, Lf+Ln);
    
    int lambdai_min = std::max(std::abs(li-ln), std::abs(Li-Ln));
    int lambdai_max = std::min(li+ln, Li+Ln);
    
    for (int lambdaf = lambdaf_min; lambdaf <= lambdaf_max; lambdaf++)
    for (int lambdai = lambdai_min; lambdai <= lambdai_max; lambdai++)
    {
        double ff = special::computef(lambdaf, Lf, lf, Ln, ln, L);
        double fi = special::computef(lambdai, Ln, ln, Li, li, L);
        
        // skip non-contributing transfers
        if (ff == 0. or fi == 0.)
            continue;
        
        // evaluate initial / final radial part of the projectile partial wave
        rArray ji = std::move(interpolate_riccati_bessel_j(grid, li, ki));
        rArray jf = std::move(interpolate_riccati_bessel_j(grid, lf, kf));
        
        // evaluate initial / final multipole potential
        rArray Vfn = std::move(interpolate_bound_free_potential(grid, lambdaf, Nf, Lf, Kn, Ln));
        rArray Vni = std::move(interpolate_free_bound_potential(grid, lambdai, Kn, Ln, Ni, Li));
        
        // evaluate the Green function terms (Bessel functions "jn" and "yn")
        // and cancel dependence on "kn" if too small (using Bessel function asymptotics)
        rArray jn = (
            kn_is_small ? // can we use the "jn" Bessel function asymptotic form?
            pow(grid,ln+1) / gsl_sf_doublefact(2*ln+1) :                // yes
            interpolate_riccati_bessel_j(grid, ln, kn) * nrm            // no
        );
        rArray yn = (
            kn_is_small ?  // can we use the "yn" Bessel function asymptotic form?
            pow(grid,-ln) * gsl_sf_doublefact(std::max(0,2*ln-1)) :     // yes
            interpolate_riccati_bessel_y(grid, ln, kn) * nrm            // no
        );
        
        // check finiteness of product (use asymptotic form of the product if not finite)
        rArray yn_ji = yn * ji;
        for (unsigned i = 1; i < yn_ji.size(); i++)
        {
            if (not std::isfinite(yn_ji[i]))
                yn_ji[i] = gsl_sf_doublefact(2*ln-1) / gsl_sf_doublefact(2*li+1) * std::pow(kn,-ln) * std::pow(ki,li+1) * std::pow(grid[i],li+1-ln);
        }
        rArray yn_jf = yn * jf;
        for (unsigned i = 1; i < yn_jf.size(); i++)
        {
            if (not std::isfinite(yn_jf[i]))
                yn_jf[i] = gsl_sf_doublefact(2*ln-1) / gsl_sf_doublefact(2*lf+1) * std::pow(kn,-ln) * std::pow(kf,lf+1) * std::pow(grid[i],lf+1-ln);
        }
        
        // finally integrate
        Complex inte = -Idir_allowed (grid, jf, Vfn, jn, yn_ji, yn_jf, ji, Vni);
        
        // if we cancelled small kn, we need to throw away nonsensial imaginary part, where this can't be done
        if (kn_is_small)
            inte.imag(0.);
        
        // comment this result
        log << format
        (
            "\t\ttransfer [%d %d] initial (%d %d, %g %d) intermediate (%g %d, %g %d) final (%d %d, %g %d) : (%g,%g)",
            lambdai, lambdaf, Ni, Li, ki, li, Kn, Ln, kn, ln, Nf, Lf, kf, lf, inte.real(), inte.imag()
        ) << std::endl;
        
        // update result
        result += ff * fi * inte;
    }
    
    return result;
}

double Idir_nFree_forbidden
(
    rArray const & grid, int L,
    int Nf, int Lf, double kf, int lf,
    double Kn, int Ln, double kappan, int ln,
    int Ni, int Li, double ki, int li,
    std::ostream & log
)
{
    double result = 0;
    
    double nrm = std::sqrt(1. / kappan);
    
    // check if kn is small enough for asymptotics
    bool kn_is_small = (kappan * grid.back() < 1e-6);
    
    int lambdaf_min = std::max(std::abs(lf-ln), std::abs(Lf-Ln));
    int lambdaf_max = std::min(lf+ln, Lf+Ln);
    
    int lambdai_min = std::max(std::abs(li-ln), std::abs(Li-Ln));
    int lambdai_max = std::min(li+ln, Li+Ln);
    
    for (int lambdaf = lambdaf_min; lambdaf <= lambdaf_max; lambdaf++)
    for (int lambdai = lambdai_min; lambdai <= lambdai_max; lambdai++)
    {
        double ff = special::computef(lambdaf, Lf, lf, Ln, ln, L);
        double fi = special::computef(lambdai, Ln, ln, Li, li, L);
        
        // skip non-contributing transfers
        if (ff == 0. or fi == 0.)
            continue;
        
        // evaluate initial / final radial part of the projectile partial wave
        rArray ji = std::move(interpolate_riccati_bessel_j(grid, li, ki));
        rArray jf = std::move(interpolate_riccati_bessel_j(grid, lf, kf));
        
        // evaluate initial / final multipole potential
        rArray Vfn = std::move(interpolate_bound_free_potential(grid, lambdaf, Nf, Lf, Kn, Ln));
        rArray Vni = std::move(interpolate_free_bound_potential(grid, lambdai, Kn, Ln, Ni, Li));
        
        // evaluate the Green function terms (Bessel functions "in" and "kn")
        // and cancel dependence on "kn" if too small (using Bessel function asymptotics)
        rArray iscaled_n = (
            kn_is_small ? // can we use the "in" Bessel function asymptotic form?
            pow(grid,ln+1) :                                            // yes
            interpolate_riccati_bessel_iscaled(grid, ln, kappan) * nrm  // no
        );
        rArray kscaled_n = (
            kn_is_small ? // can we use the "kn" Bessel function asymptotic form?
            pow(grid,-ln) :                                             // yes
            interpolate_riccati_bessel_kscaled(grid, ln, kappan) * nrm  // no
        );
        
        // finally integrate
        double inte = -Idir_forbidden(grid, jf, Vfn, iscaled_n, kscaled_n, kappan, ji, Vni);
        
        // comment this result
        log << format
        (
            "\t\ttransfer [%d %d] initial (%d %d, %g %d) intermediate (%g %d, %g %d) final (%d %d, %g %d) : (%g,%g)",
            lambdai, lambdaf, Ni, Li, ki, li, Kn, Ln, kappan, ln, Nf, Lf, kf, lf, inte, 0.
        ) << std::endl;
        
        // update result
        result += ff * fi * inte;
    }
    
    return result;
}
