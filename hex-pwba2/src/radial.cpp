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
            if (y == 0.)
                return 0;
            
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
        
        # pragma omp parallel for
        for (unsigned i = 0; i < N; i++)
            V[i] = Norm * potential(x[i]);
    }
    else
    {
        auto potential = [&](double y) -> double
        {
            if (y == 0.)
                return 0;
            
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
        
        # pragma omp parallel for
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
    rArray V (N);
    
    // precompute Coulomb zeros within the range of "x"
    int nzeros = Kb * x.back() / (2 * special::constant::pi); // initial guess
    rArray zeros;
    do
    {
        // double the necessary zero count
        nzeros = 2 * (1 + nzeros);
        zeros.resize(nzeros);
        
        // compute zeros of the Coulomb function
        special::coulomb_zeros(-1/Kb,Lb,nzeros,zeros.data());
        zeros /= Kb;
    }
    while (zeros.back() < x.back());
    
    // the classical turning point for non-S-waves
    double rt = (std::sqrt(1 + Kb*Kb*Lb*(Lb+1))-1)/(Kb*Kb);
    
    // compute the integrals
    if (lambda == 0)
    {
        auto potential = [&](double y) -> double
        {
            // use zero potential in the origin (should be Inf, but that would spreads NaNs in PC arithmeric)
            if (y == 0.)
                return 0;
            
            // integrand Pa(r₁) V(r₁,r₂) Fb(r₁) for λ = 0
            auto integrand = [&](double x) -> double
            {
                return Hydrogen::P(Na,La,x) * (1./x - 1./y) * Hydrogen::F(Kb,Lb,x);
            };
            
            // integrate per node of the Coulomb function (precomputed before)
            FixedNodeIntegrator<decltype(integrand),GaussKronrod<decltype(integrand)>> Q(integrand,zeros,rt);
            Q.setEpsAbs(0);
            Q.integrate(y, x.back());
            if (not Q.ok())
            {
                throw exception
                (
                    "Bound-free potential V[0]{%d,%d->%g,%d} integration failed for r = %g (\"%s\").",
                    Na, La, Kb, Lb, y, Q.status().c_str()
                );
            }
            return Q.result();
        };
        
        // evaluate potential in all grid points
        # pragma omp parallel for schedule (dynamic)
        for (unsigned i = 0; i < N; i++)
            V[i] = potential(x[i]);
    }
    else
    {
        auto potential = [&](double y) -> double
        {
            // use zero potential in the origin (should be Inf, but that would spreads NaNs in PC arithmeric)
            if (y == 0.)
                return 0;
            
            // integrand Pa(r₁) V(r₁,r₂) Fb(r₁) for λ ≠ 0 and r₁ < r₂
            auto integrand1 = [&](double x) -> double
            {
                return Hydrogen::P(Na,La,x) * std::pow(x/y,lambda) * Hydrogen::F(Kb,Lb,x);
            };
            
            // integrate per node of the Coulomb function (precomputed before)
            FixedNodeIntegrator<decltype(integrand1),GaussKronrod<decltype(integrand1)>> Q1(integrand1,zeros,rt);
            Q1.setEpsAbs(0);
            Q1.integrate(0., y);
            if (not Q1.ok())
            {
                throw exception
                (
                    "Bound-free potential V[%d]{%d,%d->%g,%d} [0,%g] integration failed (\"%s\").",
                    lambda, Na, La, Kb, Lb, y, Q1.status().c_str()
                );
            }
            
            // integrand Pa(r₁) V(r₁,r₂) Fb(r₁) for λ ≠ 0 and r₁ > r₂
            auto integrand2 = [&](double x) -> double
            {
                return Hydrogen::P(Na,La,x) * std::pow(y/x,lambda+1) * Hydrogen::F(Kb,Lb,x);
            };
            
            // integrate per node of the Coulomb function (precomputed before)
            FixedNodeIntegrator<decltype(integrand2),GaussKronrod<decltype(integrand2)>> Q2(integrand2,zeros,rt);
            Q2.setEpsAbs(0);
            Q2.integrate(y, x.back());
            if (not Q2.ok())
            {
                throw exception
                (
                    "Bound-free potential V[%d]{%d,%d->%g,%d} [%g,inf] integration failed (\"%s\").",
                    lambda, Na, La, Kb, Lb, y, Q2.status().c_str()
                );
            }
            
            return (Q1.result() + Q2.result()) / y;
        };
        
        // evaluate potential in all grid points
        # pragma omp parallel for schedule (dynamic)
        for (unsigned i = 0; i < N; i++)
            V[i] = potential(x[i]);
    }
    
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
        [&](double r) -> double { return ric_j(l,k*r); }
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
                return l == 0 ? 1 : 0;
            else
                return -ric_n(l,k*r);
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
        [&](double r) -> double { return ric_i_scaled(l,k*r); }
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
                return ric_k_scaled(l,k*r);
        }
    );
}

Complex Idir_allowed
(
    rArray const & grid,
    rArray const & jf, rArray const & Vfn,
    rArray const & jn, rArray const & yn,
    rArray const & ji, rArray const & Vni
)
{
    size_t N = grid.size();
    double h = grid.back() / N;
    
    rArray inner_lower(N), inner_higher_re(N), inner_higher_im(N);
    
    # pragma omp parallel for
    for (size_t i = 1; i < N; i++)
    {
        inner_lower[i] = Vni[i] * jn[i] * ji[i];
        inner_higher_im[i] = inner_lower[i];
        inner_higher_re[i] = Vni[i] * yn[i] * ji[i];
    }
    
    // forward partial trapezoidal sum (for low integral)
    /*
    double prev = 0., curr = inner_lower[1];
    inner_lower[0] = 0.;
    inner_lower[1] = 0.5 * (curr + prev);
    for (size_t i = 2; i < N; i++)
    {
        prev = curr; curr = inner_lower[i];
        inner_lower[i] = inner_lower[i-1] + 0.5 * (curr + prev);
    }
    */
    {
        gsl_interp_accel * acc = gsl_interp_accel_alloc ();
        gsl_spline * spline = gsl_spline_alloc (gsl_interp_cspline, N);
        gsl_spline_init (spline, grid.data(), inner_lower.data(), N);
        
        # pragma omp parallel for
        for (size_t i = 0; i < N; i++)
            inner_lower[i] = gsl_spline_eval_integ (spline, grid.front(), grid[i], acc);
    }
    
    // reverse partial trapezoidal sums (for high integral)
    /*
    double prevRe = inner_higher_re[N-1], currRe = inner_higher_re[N-2];
    double prevIm = inner_higher_im[N-1], currIm = inner_higher_im[N-2];
    inner_higher_re[N-1] = 0.; inner_higher_re[N-2] = 0.5 * (currRe + prevRe);
    inner_higher_im[N-1] = 0.; inner_higher_im[N-2] = 0.5 * (currIm + prevIm);
    for (size_t i = 2; i < N; i++)
    {
        prevRe = currRe; currRe = inner_higher_re[N-i-1];
        prevIm = currIm; currIm = inner_higher_im[N-i-1];
        inner_higher_re[N-i-1] = inner_higher_re[N-i] + 0.5 * (currRe + prevRe);
        inner_higher_im[N-i-1] = inner_higher_im[N-i] + 0.5 * (currIm + prevIm);
    }
    */
    {
        gsl_interp_accel * acc = gsl_interp_accel_alloc ();
        gsl_spline * spline = gsl_spline_alloc (gsl_interp_cspline, N);
        gsl_spline_init (spline, grid.data(), inner_higher_re.data(), N);
        
        # pragma omp parallel for
        for (size_t i = 0; i < N; i++)
            inner_higher_re[i] = gsl_spline_eval_integ (spline, grid[i], grid.back(), acc);
    }
    {
        gsl_interp_accel * acc = gsl_interp_accel_alloc ();
        gsl_spline * spline = gsl_spline_alloc (gsl_interp_cspline, N);
        gsl_spline_init (spline, grid.data(), inner_higher_im.data(), N);
        
        # pragma omp parallel for
        for (size_t i = 0; i < N; i++)
            inner_higher_im[i] = gsl_spline_eval_integ (spline, grid[i], grid.back(), acc);
    }
    
    rArray outer_re(N), outer_im(N);
    double sum_outer_re = 0, sum_outer_im = 0;
    
    /*
    # pragma omp parallel for reduction (+:sum_outer_re,sum_outer_im)
    for (size_t i = 1; i < N; i++)
    {
        sum_outer_re += Vfn[i] * jf[i] * (inner_lower[i] * yn[i] + inner_higher_re[i] * jn[i]);
        sum_outer_im += Vfn[i] * jf[i] * (inner_lower[i] * jn[i] + inner_higher_im[i] * jn[i]);
    }
    */
    outer_re = Vfn * jf * (inner_lower * yn + inner_higher_re * jn);
    outer_im = Vfn * jf * (inner_lower * jn + inner_higher_im * jn);
    {
        gsl_interp_accel * acc = gsl_interp_accel_alloc ();
        gsl_spline * spline = gsl_spline_alloc (gsl_interp_cspline, N);
        gsl_spline_init (spline, grid.data(), outer_re.data(), N);
        
        sum_outer_re = gsl_spline_eval_integ (spline, grid.front(), grid.back(), acc);
    }
    {
        gsl_interp_accel * acc = gsl_interp_accel_alloc ();
        gsl_spline * spline = gsl_spline_alloc (gsl_interp_cspline, N);
        gsl_spline_init (spline, grid.data(), outer_im.data(), N);
        
        sum_outer_im = gsl_spline_eval_integ (spline, grid.front(), grid.back(), acc);
    }
    
    return Complex (sum_outer_re, sum_outer_im) * h * h;
}

double Idir_forbidden
(
    rArray const & grid,
    rArray const & jf, rArray const & Vfn,
    rArray const & iscaled_n, rArray const & kscaled_n,
    rArray const & ji, rArray const & Vni
)
{
    unsigned N = grid.size();
    double h = grid.back() / N;
    
    rArray fpart = jf * Vfn;
    rArray ipart = ji * Vni;
    
    double suma = 0;
    
    // index "d" runs across the contours x - y = konst
    for (unsigned d = 0; d < N/2; d++)
    {
        // contribution of these contours
        double contrib = 0;
        
        if (d == 0)
        {
            // index "i" runs along the diagonal i = x = y
            for (unsigned i = 0; i < N; i++)
                contrib += fpart[i] * ipart[i] * iscaled_n[i] * kscaled_n[i];
        }
        else
        {
            // index "i" runs along the contour i = x = y + d
            for (unsigned i = d; i < N; i++)
                contrib += fpart[i-d] * ipart[i] * iscaled_n[i-d] * kscaled_n[i];
            // index "i" runs along the contour i = x = y - d
            for (unsigned i = 0; i < N - d; i++)
                contrib += fpart[i+d] * ipart[i] * iscaled_n[i] * kscaled_n[i+d];
        }
        
        // add properly scaled contribution
        contrib *= std::exp(-h*d);
        suma += contrib;
        
        // check convergence (do not blindly integrate everything)
        if (std::abs(contrib) < 1e-10 * std::abs(suma))
            break;
    }
    
    return suma * h * h;
}

Complex Idir_nBound_allowed
(
    rArray const & grid, int L,
    int Nf, int Lf, double kf, int lf,
    int Nn, int Ln, double kn, int ln,
    int Ni, int Li, double ki, int li
)
{
    Complex result = 0;
    
    int lambdaf_min = std::max(std::abs(lf-ln), std::abs(Lf-Ln));
    int lambdaf_max = std::min(lf+ln, Lf+Ln);
    
    int lambdai_min = std::max(std::abs(li-ln), std::abs(Li-Ln));
    int lambdai_max = std::min(li+ln, Li+Ln);
    
    for (int lambdaf = lambdaf_min; lambdaf <= lambdaf_max; lambdaf++)
    for (int lambdai = lambdai_min; lambdai <= lambdai_max; lambdai++)
    {
        double ff = computef(lambdaf, Li, li, Ln, ln, L);
        double fi = computef(lambdai, Ln, ln, Li, li, L);
        
        rArray ji = std::move(interpolate_riccati_bessel_j(grid, li, ki));
        rArray jf = std::move(interpolate_riccati_bessel_j(grid, lf, kf));
        rArray Vfn = std::move(interpolate_bound_bound_potential(grid, lambdaf, Nf, Lf, Nn, Ln));
        rArray Vni = std::move(interpolate_bound_bound_potential(grid, lambdai, Nn, Ln, Ni, Li));
        rArray jn = std::move(interpolate_riccati_bessel_j(grid, ln, kn));
        rArray yn = std::move(interpolate_riccati_bessel_y(grid, ln, kn));
        
        Complex inte = Idir_allowed (grid, jf, Vfn, jn, yn, ji, Vni);
        
        std::cout << format
        (
            "\t\ttransfer [%d %d] initial (%d %d, %g %d), intermediate (%d %d, %g %d) final (%d %d, %g %d) : (%g,%g)",
            lambdai, lambdaf, Ni, Li, ki, li, Nn, Ln, kn, ln, Nf, Lf, kf, lf, inte.real(), inte.imag()
        ) << std::endl;
        
        result += ff * fi * inte;
    }
    
    return result;
}

Complex Idir_nBound_forbidden
(
    rArray const & grid, int L,
    int Nf, int Lf, double kf, int lf,
    int Nn, int Ln, double kappan, int ln,
    int Ni, int Li, double ki, int li
)
{
    Complex result = 0;
    
    int lambdaf_min = std::max(std::abs(lf-ln), std::abs(Lf-Ln));
    int lambdaf_max = std::min(lf+ln, Lf+Ln);
    
    int lambdai_min = std::max(std::abs(li-ln), std::abs(Li-Ln));
    int lambdai_max = std::min(li+ln, Li+Ln);
    
    for (int lambdaf = lambdaf_min; lambdaf <= lambdaf_max; lambdaf++)
    for (int lambdai = lambdai_min; lambdai <= lambdai_max; lambdai++)
    {
        double ff = computef(lambdaf, Li, li, Ln, ln, L);
        double fi = computef(lambdai, Ln, ln, Li, li, L);
        
        rArray ji = std::move(interpolate_riccati_bessel_j(grid, li, ki));
        rArray jf = std::move(interpolate_riccati_bessel_j(grid, lf, kf));
        rArray Vfn = std::move(interpolate_bound_bound_potential(grid, lambdaf, Nf, Lf, Nn, Ln));
        rArray Vni = std::move(interpolate_bound_bound_potential(grid, lambdai, Nn, Ln, Ni, Li));
        rArray iscaled_n = std::move(interpolate_riccati_bessel_iscaled(grid, ln, kappan));
        rArray kscaled_n = std::move(interpolate_riccati_bessel_kscaled(grid, ln, kappan));
        
        Complex inte = 2. / special::constant::pi * Complex(0., -Idir_forbidden(grid, jf, Vfn, iscaled_n, kscaled_n, ji, Vni));
        
        std::cout << format
        (
            "\t\ttransfer [%d %d] initial (%d %d, %g %d) intermediate (%d %d, %g %d) final (%d %d, %g %d) : (%g,%g)",
            lambdai, lambdaf, Ni, Li, ki, li, Nn, Ln, kappan, ln, Nf, Lf, kf, lf, inte.real(), inte.imag()
        ) << std::endl;
        
        result += ff * fi * inte;
    }
    
    return result;
}

Complex Idir_nFree_allowed
(
    rArray const & grid, int L,
    int Nf, int Lf, double kf, int lf,
    double Kn, int Ln, double kn, int ln,
    int Ni, int Li, double ki, int li
)
{
    Complex result = 0;
    
    int lambdaf_min = std::max(std::abs(lf-ln), std::abs(Lf-Ln));
    int lambdaf_max = std::min(lf+ln, Lf+Ln);
    
    int lambdai_min = std::max(std::abs(li-ln), std::abs(Li-Ln));
    int lambdai_max = std::min(li+ln, Li+Ln);
    
    for (int lambdaf = lambdaf_min; lambdaf <= lambdaf_max; lambdaf++)
    for (int lambdai = lambdai_min; lambdai <= lambdai_max; lambdai++)
    {
        double ff = computef(lambdaf, Li, li, Ln, ln, L);
        double fi = computef(lambdai, Ln, ln, Li, li, L);
        
        rArray ji = std::move(interpolate_riccati_bessel_j(grid, li, ki));
        rArray jf = std::move(interpolate_riccati_bessel_j(grid, lf, kf));
        rArray Vfn = std::move(interpolate_bound_free_potential(grid, lambdaf, Nf, Lf, Kn, Ln));
        rArray Vni = std::move(interpolate_free_bound_potential(grid, lambdai, Kn, Ln, Ni, Li));
        rArray jn = std::move(interpolate_riccati_bessel_j(grid, ln, kn));
        rArray yn = std::move(interpolate_riccati_bessel_y(grid, ln, kn));
        
        Complex inte = Idir_allowed (grid, jf, Vfn, jn, yn, ji, Vni);
        
        std::cout << format
        (
            "\t\ttransfer [%d %d] initial (%d %d, %g %d) intermediate (%g %d, %g %d) final (%d %d, %g %d) : (%g,%g)",
            lambdai, lambdaf, Ni, Li, ki, li, Kn, Ln, kn, ln, Nf, Lf, kf, lf, inte.real(), inte.imag()
        ) << std::endl;
        
        result += ff * fi * inte;
    }
    
    return result;
}

Complex Idir_nFree_forbidden
(
    rArray const & grid, int L,
    int Nf, int Lf, double kf, int lf,
    double Kn, int Ln, double kappan, int ln,
    int Ni, int Li, double ki, int li
)
{
    Complex result = 0;
    
    int lambdaf_min = std::max(std::abs(lf-ln), std::abs(Lf-Ln));
    int lambdaf_max = std::min(lf+ln, Lf+Ln);
    
    int lambdai_min = std::max(std::abs(li-ln), std::abs(Li-Ln));
    int lambdai_max = std::min(li+ln, Li+Ln);
    
    for (int lambdaf = lambdaf_min; lambdaf <= lambdaf_max; lambdaf++)
    for (int lambdai = lambdai_min; lambdai <= lambdai_max; lambdai++)
    {
        double ff = computef(lambdaf, Li, li, Ln, ln, L);
        double fi = computef(lambdai, Ln, ln, Li, li, L);
        
        rArray ji = std::move(interpolate_riccati_bessel_j(grid, li, ki));
        rArray jf = std::move(interpolate_riccati_bessel_j(grid, lf, kf));
        rArray Vfn = std::move(interpolate_bound_free_potential(grid, lambdaf, Nf, Lf, Kn, Ln));
        rArray Vni = std::move(interpolate_free_bound_potential(grid, lambdai, Kn, Ln, Ni, Li));
        rArray iscaled_n = std::move(interpolate_riccati_bessel_iscaled(grid, ln, kappan));
        rArray kscaled_n = std::move(interpolate_riccati_bessel_kscaled(grid, ln, kappan));
        
        Complex inte = 2. / special::constant::pi * Complex(0., -Idir_forbidden(grid, jf, Vfn, iscaled_n, kscaled_n, ji, Vni));
        
        std::cout << format
        (
            "\t\ttransfer [%d %d] initial (%d %d, %g %d) intermediate (%g %d, %g %d) final (%d %d, %g %d) : (%g,%g)",
            lambdai, lambdaf, Ni, Li, ki, li, Kn, Ln, kappan, ln, Nf, Lf, kf, lf, inte.real(), inte.imag()
        ) << std::endl;
        
        result += ff * fi * inte;
    }
    
    return result;
}
