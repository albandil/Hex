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
#include <iostream>

#include <gsl/gsl_errno.h>

#include "cmdline.h"
#include "complex.h"
#include "radial.h"
#include "version.h"

#include "diophantine.h"
#include "gausskronrod.h"
#include "hydrogen.h"
#include "specf.h"

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
        nzeros *= 2;
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
        # pragma omp parallel for
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
        # pragma omp parallel for
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
    
    rArray inner_lower = Vni * jn * ji;
    rArray inner_higher_re = Vni * yn * ji;
    rArray inner_higher_im = inner_lower;
    
    for (size_t i = 1; i < N; i++)
    {
        // forward partial sum for low integral
        inner_lower[i] += inner_lower[i-1];
        
        // reverse partial sum for high integral
        inner_higher_re[N-i-1] += inner_higher_re[N-i];
        inner_higher_im[N-i-1] += inner_higher_im[N-i];
    }
    
    rArray inner_lower_re = inner_lower * yn * h;
    rArray inner_lower_im = inner_lower * jn * h;
    inner_higher_re *= jn * h;
    inner_higher_im *= jn * h;
    
    rArray outer_re = inner_lower_re + inner_higher_re;
    rArray outer_im = inner_lower_im + inner_higher_im;
    
    outer_re *= Vfn * jf;
    outer_im *= Vfn * jf;
    
    return Complex (sum(outer_re) * h, sum(outer_im) * h);
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

//
// Convert string to a specific number type.
//

// general template
template <class T> T string_to (std::string str)
{
    throw exception
    (
        "Conversion to \"%s\" not implemented.", typeid(T).name()
    );
}

// specialization for int
template <> int string_to (std::string str)
{
    // convert to int
    char* tail; long val = strtol (str.c_str(), &tail, 10);
    
    // throw or return
    if (*tail != 0x0)
        throw exception ("The string \"%s\" cannot be converted to integer number.", str.c_str());
    else
        return val;
}

// specialization for double
template <> double string_to (std::string str)
{
    // convert to float
    char* tail; double val = strtod (str.c_str(), &tail);
    
    // throw or return
    if (*tail != 0x0)
        throw exception ("The string \"%s\" cannot be converted to real number.", str.c_str());
    else
        return val;
}

// read next entry from input stream
template <class T> T read_next (std::ifstream & f)
{
    // text buffer
    std::string s;
    
    // while there is something in the file
    while (not f.eof())
    {
        // read string (it won't start with white character)
        f >> s;
        
        // check length (skip empty reads)
        if (s.size() == 0)
            continue;
        
        // check if it is a beginning of a comment
        if (s.front() == '#')
        {
            // get the rest of the line
            std::getline(f, s);
            continue;
        }
        
        // otherwise exit the loop (a valid entry was found)
        break;
    }
    
    // check for special character; exit if found
    if (s == "*")
        throw true;
    
    // convert entry to type T
    T val = string_to<T>(s.c_str());
    
    // return
    return val;
}

void parse_input_file
(
    const char * filename,
    int & L, int & Pi,
    int & Ni, int & Li, int & Nf, int & Lf, double & Ei,
    double & Rmax, int & N,
    int & maxNn, int & maxLn, double & Enmax,
    int & maxlevel_allowed, int & maxlevel_forbidden
)
{
    std::ifstream inf (filename);
    
    // check if file exists
    if (not inf.good())
        throw exception ("Can't open input file \"%s\".", filename);
    std::cout << "Reading input file \"" << filename << "\"." << std::endl << std::endl;
    
    // read all data
    L = read_next<int>(inf);  Pi = read_next<int>(inf);
    Ni = read_next<int>(inf); Li = read_next<int>(inf);
    Nf = read_next<int>(inf); Lf = read_next<int>(inf);
    Ei = read_next<double>(inf);
    Rmax = read_next<double>(inf);
    N = read_next<int>(inf);
    maxNn = read_next<int>(inf);
    maxLn = read_next<int>(inf);
    Enmax = read_next<double>(inf);
    maxlevel_allowed = read_next<int>(inf);
    maxlevel_forbidden = read_next<int>(inf);
}

namespace PWBA2
{
cArrays PartialWave_direct
(
    rArray grid,
    int L, int Pi,
    int Ni, int Li, double ki,
    int Nf, int Lf, double kf,
    int nL, int maxNn, double Enmax,
    int maxlevel_allowed, int maxlevel_forbidden
)
{
    cArrays Tdir;
    double Etot = ki*ki - 1./(Ni*Ni);
    
    // compute all angular contributions to the T-matrix
    for (int lf = std::abs(L - Lf); lf <= L + Lf; lf++)
    {
        // conserve parity
        if ((L + Lf + lf) % 2 != Pi)
        {
            std::cout << "Skipping lf = " << lf << " due to parity conservation." << std::endl;
            continue;
        }
        
        std::cout << "---------- lf = " << lf << " ----------" << std::endl << std::endl;
        
        cArray Tdir_lf ((2*Li+1)*(2*Lf+1));
        
        for (int li = std::abs(L - Li); li <= L + Li; li++)
        {
            Complex Tdir_lf_li = 0;
            
            // conserve parity
            if ((L + Li + li) % 2 != Pi)
            {
                std::cout << "Skipping li = " << li << " due to parity conservation." << std::endl;
                continue;
            }
            
//             for (int ell = 0; ell <= nL; ell++)
//             for (int Ln = ell; Ln <= ell + L + Pi; Ln++)
            {
                int Ln = 26, ln = 29;
//                 int ln = 2 * ell + L + Pi - Ln;
                std::cout << "\nli = " << li << ", Ln = " << Ln << ", ln = " << ln << std::endl << std::endl;
                
                // sum over bound states
                std::cout << "\tBound intermediate states" << std::endl;
                for (int Nn = Ln + 1; Nn <= maxNn; Nn++)
                {
                    // compute energy of the intermediate projectile state
                    double en = ki*ki - 1./(Ni*Ni) + 1./(Nn*Nn);
                    
                    // check energy
                    if (en > 0)
                    {
                        double kn = std::sqrt(en);
                        
                        // integrate
                        Tdir_lf_li += Idir_nBound_allowed
                        (
                            grid, L,
                            Nf, Lf, kf, lf,
                            Nn, Ln, kn, ln,
                            Ni, Li, ki, li
                        ) * (-2. / kn);
                    }
                    else
                    {
                        double kappan = std::sqrt(-en);
                        
                        // integrate
                        Tdir_lf_li += Idir_nBound_forbidden
                        (
                            grid, L,
                            Nf, Lf, kf, lf,
                            Nn, Ln, kappan, ln,
                            Ni, Li, ki, li
                        ) * (-2. / kappan);
                    }
                }
                
                // integrate over allowed free states (Kn^2 < ki^2 - 1/Ni^2)
                cArray allowed_integrals = { 0. };
                cArray allowed_contributions = { 0. };
                for (int Nlevel = 2; Nlevel <= maxlevel_allowed; Nlevel *= 2)
                {
                    std::cout << "\n\tAllowed intermediate states (n = " << Nlevel << ")" << std::endl;
                    
                    // get energy samples for this level
                    double H = std::min(Enmax,Etot) / Nlevel;
                    rArray allowed_energies = linspace (0., std::min(Enmax,Etot), Nlevel + 1);
                    
                    // for all samples
                    for (int ie = 1; ie <= Nlevel; ie++)
                    {
                        // skip those already computed
                        if (ie % 2 == 0)
                            continue;
                        
                        // get energy and momentum of the intermediate hydrogen continuum state
                        double En = allowed_energies[ie];
                        double Kn = std::sqrt(En);
                        
                        // get momentum of the projectile
                        double kn = std::sqrt(ki*ki - 1./(Ni*Ni) - En);
                        
                        // compute the radial integral
                        Complex I = Idir_nFree_allowed
                        (
                            grid, L,
                            Nf, Lf, kf, lf,
                            Kn, Ln, kn, ln,
                            Ni, Li, ki, li
                        ) * (-2. / kn);
                        
                        // insert the contribution
                        allowed_contributions.insert (allowed_contributions.begin() + ie, En * I);
                    }
                    
                    // integrate
                    allowed_integrals.push_back(special::integral::trapz(H, allowed_contributions));
                    
                    std::cout << "\t\tsum: " << allowed_integrals.back() << ", romberg: " << special::integral::romberg(allowed_integrals).back() << std::endl;
                    
                    // check convergence
                    const double eps = 1e-6;
                    if (allowed_integrals.size() > 1)
                    {
                        if (std::abs(allowed_integrals.back(0)) < eps * std::abs(allowed_integrals.back(1)))
                            break;
                    }
                }
                
                Tdir_lf_li += allowed_integrals.back();
                
                // integrate over forbidden free states (Kn^2 > ki^2 - 1/Ni^2)
                cArray forbidden_integrals = { 0. };
                cArray forbidden_contributions = { 0. };
                if (Etot < Enmax)
                for (int Nlevel = 2; Nlevel <= maxlevel_forbidden; Nlevel *= 2)
                {
                    std::cout << "\n\tForbidden intermediate states (n = " << Nlevel << ")" << std::endl;
                    
                    // get energy samples for this level
                    double H = (Enmax - Etot) / Nlevel;
                    rArray forbidden_energies = linspace (Etot, Enmax, Nlevel + 1);
                    
                    // for all samples
                    for (int ie = 1; ie <= Nlevel; ie++)
                    {
                        // skip those already computed
                        if (ie % 2 == 0)
                            continue;
                        
                        // get energy and momentum of the intermediate hydrogen continuum state
                        double En = forbidden_energies[ie];
                        double Kn = std::sqrt(En);
                        
                        // get momentum of the projectile
                        double kappan = std::sqrt(En - ki*ki + 1./(Ni*Ni));
                        
                        // compute the radial integral
                        Complex I = Idir_nFree_forbidden
                        (
                            grid, L,
                            Nf, Lf, kf, lf,
                            Kn, Ln, kappan, ln,
                            Ni, Li, ki, li
                        ) * (-2. / kappan);
                        
                        // insert the contribution
                        forbidden_contributions.insert (forbidden_contributions.begin() + ie, En * I);
                    }
                    
                    // integrate
                    forbidden_integrals.push_back(special::integral::trapz(H, forbidden_contributions));
                    
                    std::cout << "\t\tsum: " << forbidden_integrals.back() << ", romberg: " << special::integral::romberg(forbidden_integrals).back() << std::endl;
                    
                    // check convergence
                    const double eps = 1e-6;
                    if (allowed_integrals.size() > 1)
                    {
                        if (std::abs(forbidden_integrals.back(0)) < eps * std::abs(forbidden_integrals.back(1)))
                            break;
                    }
                }
                
                Tdir_lf_li += forbidden_integrals.back();
            }
            
            Complex factor = std::pow(Complex(0.,1.),li-lf) * std::pow(4*special::constant::pi, 1.5) * std::sqrt(2*li + 1.);
            
            for (int Mi = -Li; Mi <= Li; Mi++)
            for (int Mf = -Lf; Mf <= Lf; Mf++)
            {
                double Cf = ClebschGordan(Lf, Mf, lf, Mi-Mf, L, Mi);
                double Ci = ClebschGordan(Li, Mi, li, 0, L, Mi);
                
                std::cout << "Contribution of li = " << li << ": " << factor * Cf * Ci * Tdir_lf_li << std::endl;
                
                Tdir_lf[(Mi + Li) * (2 * Lf + 1) + Mf + Lf] += factor * Cf * Ci * Tdir_lf_li;
            }
        }
        
        Tdir.push_back(Tdir_lf / (ki * kf));
    }
    
    return Tdir;
}

cArrays FullTMatrix_direct
(
    rArray grid,
    int Ni, int Li, double ki,
    int Nf, int Lf, double kf,
    int maxNn, int maxLn, double maxEn,
    int maxlevel_allowed, int maxlevel_forbidden
)
{
    cArray Tdir;
    
    // for all bound intermediate states
    std::cout << "\tBound intermediate states" << std::endl;
    for (int Nn = 1; Nn <= maxNn; Nn++)
    for (int Ln = 0; Ln <= std::min(maxLn, Nn-1); Ln++)
    {
        // compute energy of the intermediate projectile state
        Complex en = ki*ki - 1./(Ni*Ni) + 1./(Nn*Nn);
        
        // get momentum of the projectile in the intermediate state
        Complex kn = std::sqrt(en);
    }
    
    // for all free intermediate states
    // TODO
    
    return cArrays ({ Tdir });
}

}; // end of namespace "PWBA2"

const std::string sample_input =
    "# ------------------------\n"
    "# Quantum state numbers\n"
    "# ------------------------\n"
    "\n"
    "# L  Pi\n"
    "  0  0\n"
    "\n"
    "# initial state\n"
    "# ni li\n"
    "   1  0\n"
    "\n"
    "# final state\n"
    "# nf lf\n"
    "   1  0\n"
    "\n"
    "# impact energy\n"
    "# Ei\n"
    "   4\n"
    "\n"
    "# ------------------------\n"
    "# Grid parameters\n"
    "# ------------------------\n"
    "\n"
    "# maximal radius\n"
    "# Rmax\n"
    "   100\n"
    "\n"
    "# linear samples\n"
    "# N\n"
    "  1000\n"
    "\n"
    "# ------------------------\n"
    "# Intermediate states\n"
    "# ------------------------\n"
    "\n"
    "# maximal quantum numbers\n"
    "# maxNn  nL     maxEn\n"
    "  8       3     20\n"
    "\n"
    "# continuum integration samples\n"
    "# allowed forbidden\n"
    "  64      32\n";


int main (int argc, char* argv[])
{
    // print program logo
    std::cout << logo_raw() << std::endl;
    std::cout << "=== Plane wave second Born approximation ===" << std::endl << std::endl;
    
    // disable fatal GSL errors
    gsl_set_error_handler_off();
    
    // parse command line
    bool partial_wave = false;
    ParseCommandLine
    (
        argc, argv,
        
        "example", "e", 0, [&](std::string optarg) -> bool
            {
                std::cout << "Writing sample input file to \"example.inp\".\n\n";
                
                // produce sample input file
                std::ofstream out("example.inp");
                if (out.bad())
                    throw exception ("Error: Cannot write to \"example.inp\"\n");
                
                out << sample_input;
                    
                out.close();
                exit(0);
            },
        "help", "h", 0, [](std::string optarg) -> bool
            {
                // print usage information
                std::cout << "\n"
                    "Available switches (short forms in parentheses):                                                                  \n"
                    "                                                                                                                  \n"
                    "\t--example                 (-e)  create sample input file                                                        \n"
                    "\t--help                    (-h)  display this help                                                               \n"
                    "\t--partial-wave            (-w)  compute only contribution of single partial wave                                \n"
                    "                                                                                                                  \n"
                ;
                exit(0);
            },
        "partial-wave", "w", 0, [&](std::string optarg) -> bool
            {
                // compute only contribution of a single partial wave
                partial_wave = true;
                return true;
            },
        
        [](std::string opt, std::string optarg) -> bool
            {
                throw exception
                (
                    "Unknown option \"%s\" with argument \"%s\".",
                    opt.c_str(), optarg.c_str()
                );
            }
    );
    
    // grid parameters
    int N;
    double Rmax;
    
    // quantum numbers
    int Pi, L, maxNn, nL;
    int Ni, Li, Nf, Lf;
    double Ei, Enmax;
    int maxlevel_allowed;
    int maxlevel_forbidden;
    
    parse_input_file
    (
        "pwba2.inp",
        L, Pi,
        Ni, Li, Nf, Lf, Ei,
        Rmax, N,
        maxNn, nL, Enmax,
        maxlevel_allowed, maxlevel_forbidden
    );
    
    // compute other variables from input
    double ki = std::sqrt(Ei);
    double Etot = ki*ki - 1./(Ni*Ni);
    rArray grid = linspace(0., Rmax, N);
    
    // echo input data
    std::cout << "Quantum state parameters:" << std::endl;
    std::cout << "\t- total angular momentum: L = " << L << (partial_wave ? " (not used)" : "") << std::endl;
    std::cout << "\t- total parity: Π = " << Pi << (partial_wave ? " (not used)" : "") << std::endl;
    std::cout << "\t- initial atomic state: Ni = " << Ni << ", Li = " << Li << std::endl;
    std::cout << "\t- final atomic state: Nf = " << Nf << ", Lf = " << Lf << std::endl;
    std::cout << "\t- impact energy: Ei = " << ki * ki << std::endl;
    std::cout << "\t- total energy: Etot = " << Etot << std::endl;
    std::cout << std::endl;
    std::cout << "Grid parameters:" << std::endl;
    std::cout << "\t- grid length: Rmax = " << Rmax << std::endl;
    std::cout << "\t- grid total samples: N = " << N << std::endl;
    std::cout << "\t- grid spacing: h = " << Rmax / N << std::endl;
    std::cout << std::endl;
    std::cout << "Intermediate atomic states:" << std::endl;
    std::cout << "\t- maximal bound state principal quantum number: maxNn = " << maxNn << std::endl;
    std::cout << "\t- maximal intermediate angular momentum sum (- L): nL = " << nL << std::endl;
    std::cout << "\t- maximal energy: Enmax = " << Enmax << std::endl;
    std::cout << "\t- how many allowed to integrate: " << maxlevel_allowed << std::endl;
    std::cout << "\t- how many forbidden to integrate: " << maxlevel_forbidden << std::endl;
    
    if (partial_wave)
    {
        // write all intermediate states' angular momenta pairs
        for (int ell = 0; ell <= nL; ell++)
        {
            std::cout << "\t- angular momenta [" << ell << "]: ";
            for (int Ln = ell; Ln <= ell + L + Pi; Ln++)
            {
                int ln = 2 * ell + L + Pi - Ln;
                std::cout << "(" << Ln << "," << ln << ") ";
            }
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;
    
    // final energy
    double ef = ki*ki - 1./(Ni*Ni) + 1./(Nf*Nf);
    
    // check energy
    if (ef < 0)
    {
        std::cout << "Excitation from Ni = " << Ni << " to Nf = " << Nf << " is not possible at given energy.";
        return 1;
    }
    
    // final momentum
    double kf = std::sqrt(ef);
    
    // check grid spacing
    double min_wavelength = 2 * special::constant::pi / std::sqrt(Enmax);
    std::cout << "There are " << min_wavelength / (Rmax / N) << " grid samples per shortest wavelength." << std::endl;
    if (Rmax / N > min_wavelength / 10)
        std::cout << "Warning: Grid is not sufficiently fine!" << std::endl;
    std::cout << std::endl;
    
    // outgoing electron partial T-matrices
    cArrays Tdir =
        partial_wave ?
        PWBA2::PartialWave_direct(grid, L, Pi, Ni, Li, ki, Nf, Lf, kf, nL, maxNn, Enmax, maxlevel_allowed, maxlevel_forbidden) :
        PWBA2::FullTMatrix_direct(grid, Ni, Li, ki, Nf, Lf, kf, maxNn, nL, Enmax, maxlevel_allowed, maxlevel_forbidden);
    
    std::cout << "Tdir = " << Tdir << std::endl;
    std::cout << "Tdir sums = " << sums(Tdir) << std::endl;
    std::cout << std::endl;
    
    return 0;
}
