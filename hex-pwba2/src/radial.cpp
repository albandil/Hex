//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2015, Jakub Benda, Charles University in Prague                    //
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

#include <mpc.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "arrays.h"
#include "gausskronrod.h"
#include "hydrogen.h"
#include "misc.h"
#include "nodeintegrate.h"
#include "radial.h"
#include "special.h"

// include Michel's hypergeometric function implementation
#include "hyp_2F1.H"

double integ_pow_exp_coulomb
(
    unsigned long a, double b, double k, unsigned long L, double r,
    int bits = 256, int max_iter = 10000, double tolerance = 1e-5, unsigned * iter = nullptr
)
{
    // update coordinate power
    a += L + 1;
    
    // initialize some auxiliary variables
    mpc_t lambda; mpc_init2(lambda, bits); mpc_set_d_d(lambda, b,   k,       MPC_RNDNN); // lambda = b + ik
    mpc_t rho;    mpc_init2(rho,    bits); mpc_set_d_d(rho,    k*r, 0.0,     MPC_RNDNN); // rho = k*r
    mpc_t R;      mpc_init2(R,      bits); mpc_set_d_d(R,      r,   0.0,     MPC_RNDNN); // R = r
    mpc_t xi;     mpc_init2(xi,     bits); mpc_set_d_d(xi,     b*r, k*r,     MPC_RNDNN); // xi = lambda*r = (b + ik)*r
    mpc_t T;      mpc_init2(T,      bits); mpc_set_d_d(T,      0.0, 2.0 * k, MPC_RNDNN); // T = 2ik
    mpc_t U;      mpc_init2(U,      bits); // U = L + i/k + n
    
    // initialize inner sum
    mpc_t inner_term; // = 1
        mpc_init2(inner_term, bits);
        mpc_set_ui_ui(inner_term, 1, 0, MPC_RNDNN);
    mpc_t inner_sum; // exp[(b+ik)r] - 1
        mpc_init2(inner_sum,  bits);
        mpc_exp(inner_sum, xi, MPC_RNDNN);
        mpc_sub(inner_sum, inner_sum, inner_term, MPC_RNDNN);
    for (unsigned long j = 1; j <= a; j++)
    {
        // inner_term *= (b+ik)*r/j;
        mpc_mul(inner_term, inner_term, xi, MPC_RNDNN);
        mpc_div_ui(inner_term, inner_term, j, MPC_RNDNN);
        
        // inner_sum -= inner_term;
        mpc_sub(inner_sum, inner_sum, inner_term, MPC_RNDNN);
    }
    
    // outer sum
    mpc_t factor; // = 1
        mpc_init2(factor, bits);
        mpc_set_ui_ui(factor, 1, 0, MPC_RNDNN);
    mpc_t term; // = inner_sum
        mpc_init2(term, bits);
        mpc_set(term, inner_sum, MPC_RNDNN);
    mpc_t sum; // = term
        mpc_init2(sum, bits);
        mpc_set(sum, term, MPC_RNDNN);
    mpfr_t sqrmag_term; // = |term|^2
        mpfr_init2(sqrmag_term, bits);
        mpc_norm(sqrmag_term, term, MPFR_RNDN);
    mpfr_t sqrmag_sum; // = |sum|^2
        mpfr_init2(sqrmag_sum, bits);
        mpc_norm(sqrmag_sum, sum, MPFR_RNDN);
    mpfr_t contrib; // = sqrmag_term / sqrmag_sum
        mpfr_init2(contrib, bits);
        mpfr_div(contrib, sqrmag_term, sqrmag_sum, MPFR_RNDN);
        mpfr_sqrt(contrib, contrib, MPFR_RNDN);
    mpfr_t tmp;
        mpfr_init2(tmp, bits);
    int n = 0;
    for (n = 1; n < max_iter and mpfr_get_d(contrib, MPFR_RNDN) > tolerance; n++)
    {
        /// DEBUG
//         std::cout << "n = " << n-1 << std::endl;
//         std::cout << "\tcontrib = " << mpfr_get_d(contrib, MPFR_RNDN) << std::endl;
        mpc_abs(tmp, inner_sum, MPFR_RNDN);
//         std::cout << "\t|inner_sum| = " << mpfr_get_d(tmp, MPFR_RNDN) << std::endl;
        
        // inner_term *= (b + ik) * r / (a + n)
        mpc_mul(inner_term, inner_term, xi, MPC_RNDNN);
        mpc_div_ui(inner_term, inner_term, a + n, MPC_RNDNN);
        
        // inner_sum -= inner_term;
        mpc_sub(inner_sum, inner_sum, inner_term, MPC_RNDNN);
        
        // factor *= (L + n + i/k) / (2*L + 1 + n)
        mpc_set_d_d(U, L + n, 1.0 / k, MPC_RNDNN);
        mpc_mul(factor, factor, U, MPC_RNDNN);
        mpc_div_ui(factor, factor, 2*L + 1 + n, MPC_RNDNN);
        
        // factor *= 2ik / (b+ik)
        mpc_mul(factor, factor, T, MPC_RNDNN);
        mpc_div(factor, factor, lambda, MPC_RNDNN);
        
        // factor *= (a + n) / n
        mpc_mul_ui(factor, factor, a + n, MPC_RNDNN);
        mpc_div_ui(factor, factor, n, MPC_RNDNN);
        
        // term = factor * inner_sum;
        mpc_mul(term, factor, inner_sum, MPC_RNDNN);
        
        // sum += term;
        mpc_add(sum, sum, term, MPC_RNDNN);
        
        // get magnitude of the last added term and of the sum
        mpc_norm(sqrmag_term, term, MPFR_RNDN);
        mpc_norm(sqrmag_sum, sum, MPFR_RNDN);
        mpfr_div(contrib, sqrmag_term, sqrmag_sum, MPFR_RNDN);
        mpfr_sqrt(contrib, contrib, MPFR_RNDN);
    }
    
    // calculate Coulomb wave normalization factor
    gsl_sf_result Cl;
    gsl_sf_coulomb_CL_e(L, -1/k, &Cl);
    
    // store iteration count
    if (iter != nullptr)
        *iter = n;
    
    // convert result to common 16-byte complex type
    Complex csum
    (
        mpfr_get_d(mpc_realref(sum),MPFR_RNDN),
        mpfr_get_d(mpc_imagref(sum),MPFR_RNDN)
    );
    
    // return result
    return /* Cl.val * std::pow(k,L+1) * */ gsl_sf_fact(a) * ( csum / std::pow(Complex(b,k),a+1) * std::exp(-Complex(b*r,k*r)) ).real();
}

rArray interpolate_bound_bound_potential (rArray const & x, int lambda, int Na, int La, int Nb, int Lb)
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

rArray interpolate_bound_free_potential_0 (rArray const & x, int lambda, int Na, int La, double Kb, int Lb)
{
    // array of bound-free potential evaluations
    rArray V(x.size());
    
    // bound state normalization constant
    double Nnl = std::sqrt(std::pow(2./Na,3) * gsl_sf_fact(Na-La-1) / (2 * Na * std::pow(gsl_sf_fact(Na+La),3)));
    
    // free state normalization constant
    gsl_sf_result Cl;
    gsl_sf_coulomb_CL_e(Lb, -1/Kb, &Cl);
    Cl.val *= special::constant::sqrt_two / (Kb * special::constant::sqrt_pi);
    
    double a = 0, inte_0_inf_L = 0, inte_0_inf_mLm1 = 0, inte_0_r_mLm1 = 0, sign = 0, laguerre = 0, inte_0_r_L = 0, eps = 1e-5;
    unsigned bits = 1024, max_iter = 1000, iter;
    
    // compute the integrals
    if (lambda == 0)
    {
        // sum all terms
        for (int j = 2*La+1; j <= Na+La; j++)
        {
            // Laguerre polynomial factor
            sign = -(j % 2 == 0 ? 1. : -1.); // FIXME Which sign ?
            laguerre = (gsl_sf_pow_int(gsl_sf_fact(Na+La),2) * sign * gsl_sf_pow_int(2./Na,j-La-1)) / (gsl_sf_fact(Na+La-j) * gsl_sf_fact(j) * gsl_sf_fact(j-2*La-1));
            
            // integral of P(r) F(r) / r, r = 0 .. infinity
            a = (j - La - 1) + 1 + (Lb + 1) - 1;
            inte_0_inf_mLm1 = gsl_sf_fact(a) * (hyp_2F1(Complex(Lb+1.,1./Kb), a+1, 2.*Lb+2., Complex(0.,2.*Kb)/Complex(1./Na,Kb)) / std::pow(Complex(1./Na,Kb),a+1)).real();
            
            // integral of P(r) F(r), r = 0 .. infinity
            a = (j - La - 1) + 1 + (Lb + 1);
            inte_0_inf_L = (La == Lb ? 0. : gsl_sf_fact(a) * (hyp_2F1(Complex(Lb+1.,1./Kb), a+1, 2.*Lb+2., Complex(0.,2.*Kb)/Complex(1./Na,Kb)) / std::pow(Complex(1./Na,Kb),a+1)).real());
            
            // calculate from Taylor series (small distances)
            unsigned i;
            for (i = 1; i < x.size(); i++)
            {
                // lower integral of P(r) F(r) / r
                inte_0_r_mLm1 = integ_pow_exp_coulomb((j - La - 1) + 1 - 1, 1./Na, Kb, Lb, x[i], bits, max_iter, eps, &iter);
                if (iter == max_iter)
                    break;
                
                // lower integral of P(r) F(r)
                inte_0_r_L = integ_pow_exp_coulomb((j - La - 1) + 1, 1./Na, Kb, Lb, x[i], bits, max_iter, eps, &iter);
                if (iter == max_iter)
                    break;
                
                std::cout << "Taylor: " << x[i] << " " << inte_0_inf_mLm1 << " " << inte_0_r_mLm1 << " " << inte_0_inf_L << " " << inte_0_r_L << std::endl;
                
                // calculate potential for all grid points
                V[i] += laguerre * ((inte_0_inf_mLm1 - inte_0_r_mLm1) - (inte_0_inf_L - inte_0_r_L) / x[i]);
            }
            
            // calculate from asymptotic series (large distances)
            for (; i < x.size(); i++)
            {
                // TODO
            }
            
            std::exit(0);
        }
    }
    else
    {
        // sum all terms
        for (int j = 2*La+1; j <= Na+La; j++)
        {
            // Laguerre polynomial factor  // FIXME Which sign ?
            sign = -(j % 2 == 0 ? 1. : -1.);
            laguerre = (gsl_sf_pow_int(gsl_sf_fact(Na+La),2) * sign * gsl_sf_pow_int(2./Na,j-La-1)) / (gsl_sf_fact(Na+La-j) * gsl_sf_fact(j) * gsl_sf_fact(j-2*La-1));
            
            // full integral of P(r) F(r) / r^(lambda + 1)
            a = (j - La - 1) + 1 + (Lb + 1) - lambda - 1;
            inte_0_inf_mLm1 = gsl_sf_fact(a) * (hyp_2F1(Complex(Lb+1.,1./Kb), a+1, 2.*Lb+2., Complex(0.,2.*Kb)/Complex(1./Na,Kb)) / std::pow(Complex(1./Na,Kb),a+1)).real();
            
            // for all grid points
            for (unsigned i = 1; i < x.size(); i++)
            {
                // lower integral of P(r) F(r) / r^(lambda + 1)
                a = (j - La - 1) + 1 + (Lb + 1) - lambda - 1;
                inte_0_r_mLm1 = integ_pow_exp_coulomb(a, 1./Na, Kb, Lb, x[i], bits, max_iter, eps, &iter);
                
                // lower integral of P(r) F(r) * r^lambda
                a = (j - La - 1) + 1 + (Lb + 1) + lambda;
                inte_0_r_L = integ_pow_exp_coulomb(a, 1./Na, Kb, Lb, x[i], bits, max_iter, eps, &iter);
                
                // calculate potential for all grid points
                V[i] += laguerre * (gsl_sf_pow_int(x[i],lambda) * (inte_0_inf_mLm1 - inte_0_r_mLm1) + inte_0_r_L / gsl_sf_pow_int(x[i],lambda+1));
            }
        }
    }
    
    // correct norm and value in the origin
    V *= Nnl * Cl.val * gsl_sf_pow_int(Kb,Lb+1);
    V[0] = 0;
    
    write_array(V, "V-new.dat");
    std::exit(0);
    
    // return the array of evaluations
    return V;
}

rArray interpolate_bound_free_potential (rArray const & grid, int lambda, int Na, int La, double Kb, int Lb)
{
    // array of bound-free potential evaluations
    unsigned N = grid.size();
    double rmax = grid.back();
    rArray V(N);
    
    // number of Coulomb zeros within the range of grid "x" (initial guess)
    unsigned nzeros = Kb * rmax / (2 * special::constant::pi);
    
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
    while (zeros.back() < rmax);
    
    // the classical turning point for non-S-waves
    double rt = (std::sqrt(1 + Kb*Kb*Lb*(Lb+1))-1)/(Kb*Kb);
    
    // Coulomb function scaled by hydrogen orbital
    rArray PF(N);
    for (unsigned i = 0; i < N; i++)
        PF[i] = Hydrogen::P(Na,La,grid[i]) * Hydrogen::F(Kb,Lb,grid[i]);
    
    // cubic spline interpolation of the above
    gsl_spline * spline = gsl_spline_alloc (gsl_interp_cspline, N);
    gsl_spline_init (spline, grid.data(), PF.data(), N);
    gsl_interp_accel * acc = gsl_interp_accel_alloc ();
    
    // evaluate potential at all grid points 'y'
    for (unsigned i = 1; i < N; i++)
    {
        // grid point
        double y = grid[i];
        
        // integrand Pa(r₁) V(r₁,r₂) Fb(r₁) for r₁ < r₂ and r₁ > r₂, resp.
        auto integrand1 = [&](double x) -> double { return gsl_spline_eval(spline, x, acc) * gsl_sf_pow_int(x/y,lambda); };
        auto integrand2 = [&](double x) -> double { return gsl_spline_eval(spline, x, acc) * gsl_sf_pow_int(y/x,lambda + 1); };
        
        // Coulomb function node integrators
        FixedNodeIntegrator<decltype(integrand1),GaussKronrod<decltype(integrand1)>,double> Q1(integrand1,zeros,rt); Q1.setEpsAbs(0);
        FixedNodeIntegrator<decltype(integrand2),GaussKronrod<decltype(integrand2)>,double> Q2(integrand2,zeros,rt); Q2.setEpsAbs(0);
        
        // integrate
        if (not Q1.integrate(0., y))
            HexException("Bound-free potential V[%d]{%d,%d->%g,%d} [0,%g] integration failed (\"%s\").", lambda, Na, La, Kb, Lb, y, Q1.status().c_str());
        if (not Q2.integrate(y, rmax))
            HexException("Bound-free potential V[%d]{%d,%d->%g,%d} [%g,inf] integration failed (\"%s\").", lambda, Na, La, Kb, Lb, y, Q2.status().c_str());
        
        // return sum
        V[i] = (Q1.result() + Q2.result()) / y;
    }
    
    // free allocated memory
    gsl_interp_accel_free(acc);
    gsl_spline_free(spline);
    
    // return the array of evaluations
    return V;
}

rArray interpolate_free_bound_potential (rArray const & x, int lambda, double Ka, int La, int Nb, int Lb)
{
    // the potential is real, so just return the transpose (= complex conjugate)
    return interpolate_bound_free_potential (x, lambda, Nb, Lb, Ka, La);
}

rArray interpolate_riccati_bessel_j (rArray const & x, int l, double k)
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
                return l == 0 ? 1. : 0.; // NOTE: zero should be inf, actually
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
                return l == 0 ? 1 : 0; // NOTE: zero should be inf, actually
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
        
        // also, the real part may have exploded
        if (not std::isfinite(inte.real()))
            inte.real(0.);
        
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
