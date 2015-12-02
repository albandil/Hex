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

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "hex-arrays.h"
#include "hex-gausskronrod.h"
#include "hex-hydrogen.h"
#include "hex-misc.h"
#include "hex-nodeintegrate.h"
#include "hex-special.h"

#include "radial.h"

// include Michel's hypergeometric function implementation
#include "hyp_2F1.H"

struct ScaledFunction
{
    // behaviour near zero
    int power;
    double factor;
    
    // last point scaled by (factor*x)^power
    int lScaled;
    
    // evaluated points
    rArray eval;
};

struct MultipolePotential : ScaledFunction
{
    // has the potential monopole component?
    bool hasMonopole;
    
    // asymptotic multipole potential
    double asy_factor;
};

struct RiccatiBesselJ : public ScaledFunction
{
    /// Constructor.
    RiccatiBesselJ (int l_, double k_, rArray const & grid)
    {
        power = l_ + 1;
        factor = k_;
        lScaled = -1;
        eval.resize(grid.size());
        l = l_;
        k = k_;
        rt = std::sqrt((l_+0.5)*(l_+1.5))/k_;
        
        // number of grid points
        unsigned Npt = grid.size();
        
//         std::ofstream out (format("j-%d-%g.out",l,k).c_str());
        
        // for all grid points
        for (unsigned ipt = 0; ipt < Npt; ipt++)
        {
            // get current radial coordinate
            double r = grid[ipt];
            
            // how far is the grid point from the origin?
            if (k*r < 1.)
            {
                lScaled = ipt;
                // Evaluate the regular Riccati-Bessel function using power series.
                // - The scaling factor (kr)^(l+1) is explicitly put out.
                double term = 1. / gsl_sf_doublefact(2*l+1), sum = term;
                for (unsigned j = 1; std::abs(term) > 1e-15 * std::abs(sum); j++)
                {
                    term *= -0.5*(k*r)*(k*r) / (j * (2*l+2*j+1));
                    sum += term;
                }
                eval[ipt] = sum;
                
//                 out << r << " " << eval[ipt] << " " << eval[ipt] * gsl_sf_pow_int(k*r,power) << std::endl;
            }
            else if (r < rt)
            {
                lScaled = ipt;
                // Evaluate the function within the rest of the classically forbidden region.
                // - Use library function and just divide the result by the scale factor.
                eval[ipt] = special::ric_j(l,k*r) / gsl_sf_pow_int(k*r,l+1);
                
//                 out << r << " " << eval[ipt] << " " << eval[ipt] << std::endl;
            }
            else
            {
                // Calculate the function using a library routine.
                // - No scaling is done (function is oscillatory).
                eval[ipt] = special::ric_j(l,k*r);
                
//                 out << r << " " << eval[ipt] << " " << eval[ipt] << std::endl;
            }
        }
    }
    
    /// Angular momentum.
    int l;
    
    /// Linear momentum.
    double k;
    
    /// Turning point.
    double rt;
};

struct RiccatiBesselI : public ScaledFunction
{
    /// Constructor.
    RiccatiBesselI (int l_, double k_, rArray const & grid)
    {
        power = l_ + 1;
        factor = k_;
        lScaled = -1;
        eval.resize(grid.size());
        l = l_;
        k = k_;
        rt = std::sqrt((l_+0.5)*(l_+1.5))/k_;
        
        // number of grid points
        unsigned Npt = grid.size();
        
//         std::ofstream out (format("i-%d-%g.out",l,k).c_str());
        
        // for all grid points
        for (unsigned ipt = 0; ipt < Npt; ipt++)
        {
            // get current radial coordinate
            double r = grid[ipt];
            
            // how far is the grid point from the origin?
            if (k*r < 1.)
            {
                lScaled = ipt;
                // Evaluate the regular Riccati-Bessel function using power series.
                // - The scaling factor (kr)^(l+1) is explicitly put out.
                double term = 1. / gsl_sf_doublefact(2*l+1), sum = term;
                for (unsigned j = 1; term > 1e-15 * sum; j++)
                {
                    term *= 0.5*(k*r)*(k*r) / (j * (2*l+2*j+1));
                    sum += term;
                }
                eval[ipt] = sum * std::exp(-k*r);
                
//                 out << r << " " << eval[ipt] << " " << eval[ipt] * gsl_sf_pow_int(k*r,power) << std::endl;
            }
            /*else if (r < rt)
            {
                lScaled = ipt;
                // Evaluate the function within the rest of the classically forbidden region.
                // - Use library function and just divide the result by the scale factor.
                eval[ipt] = special::ric_i_scaled(l,k*r) / gsl_sf_pow_int(k*r,l+1);
                
                out << r << " " << eval[ipt] << " " << eval[ipt] << std::endl;
            }*/
            else
            {
                // Calculate the function using a library routine.
                // - No scaling is done (function is oscillatory).
                eval[ipt] = special::ric_i_scaled(l,k*r);
                
//                 out << r << " " << eval[ipt] << " " << eval[ipt] << std::endl;
            }
        }
    }
    
    /// Angular momentum.
    int l;
    
    /// Linear momentum.
    double k;
    
    /// Turning point.
    double rt;
};

struct RiccatiBesselY : public ScaledFunction
{
    /// Constructor.
    RiccatiBesselY (int l_, double k_, rArray const & grid)
    {
        power = -l_;
        factor = k_;
        lScaled = -1;
        eval.resize(grid.size());
        l = l_;
        k = k_;
        rt = std::sqrt((l_+0.5)*(l_+1.5))/k_;
        
        // number of grid points
        unsigned Npt = grid.size();
        
//         std::ofstream out (format("y-%d-%g.out",l,k).c_str());
        
        // for all grid points
        for (unsigned ipt = 0; ipt < Npt; ipt++)
        {
            // get current radial coordinate
            double r = grid[ipt];
            
            // how far is the grid point from the origin?
            if (k*r == 0.)
            {
                // Use the limit of x^l y_l(k*r) at k*r = 0.
                lScaled = ipt;
                eval[ipt] = (l == 0 ? -1. : -gsl_sf_doublefact(2*l - 1));
                
//                 out << r << " " << eval[ipt] << " " << eval[ipt] * gsl_sf_pow_int(k*r,power) << std::endl;
            }
            else if (k*r < 1.)
            {
                // Evaluate the regular Riccati-Bessel function using power series.
                // - The scaling factor is explicitly put out.
                Complex term = std::pow(Complex(0,1),-l-1) * gsl_sf_pow_int(k*r,l), sum = term;
                for (int m = 1; m <= l; m++)
                {
                    term *= Complex(0,1) * double(l + m) * double(l + 1 - m) / (2 * m * k * r);
                    sum += term;
                }
                lScaled = ipt;
                eval[ipt] = (sum * Complex(std::cos(k*r),std::sin(k*r))).imag();
                
//                 out << r << " " << eval[ipt] << " " << eval[ipt] * gsl_sf_pow_int(k*r,power) << std::endl;
            }
            /*else if (r < rt)
            {
                // Evaluate the function within the rest of the classically forbidden region.
                // - Use library function and just divide the result by the scale factor.
                lScaled = ipt;
                eval[ipt] = special::ric_n(l,k*r) * gsl_sf_pow_int(k*r,l);
                
                out << r << " " << eval[ipt] << " " << eval[ipt] * gsl_sf_pow_int(k*r,power) << std::endl;
            }*/
            else
            {
                // Calculate the function using a library routine.
                // - No scaling is done (function is oscillatory).
                eval[ipt] = special::ric_n(l,k*r);
                
//                 out << r << " " << eval[ipt] << " " << eval[ipt] << std::endl;
            }
        }
        
        // change sign
        eval = -eval;
    }
    
    /// Angular momentum.
    int l;
    
    /// Linear momentum.
    double k;
    
    /// Turning point.
    double rt;
};

struct RiccatiBesselK : public ScaledFunction
{
    /// Constructor.
    RiccatiBesselK (int l_, double k_, rArray const & grid)
    {
        power = -l_;
        factor = k_;
        lScaled = -1;
        eval.resize(grid.size());
        l = l_;
        k = k_;
        rt = std::sqrt((l_+0.5)*(l_+1.5))/k_;
        
        // number of grid points
        unsigned Npt = grid.size();
        
//         std::ofstream out (format("k-%d-%g.out",l,k).c_str());
        
        // for all grid points
        for (unsigned ipt = 0; ipt < Npt; ipt++)
        {
            // get current radial coordinate
            double r = grid[ipt];
            
            // how far is the grid point from the origin?
            if (k*r == 0.)
            {
                // Use the limit of x^l K_l(k*r) at k*r = 0.
                lScaled = ipt;
                eval[ipt] = (l == 0 ? 1. : gsl_sf_doublefact(2*l - 1));
                
//                 out << r << " " << eval[ipt] << " " << eval[ipt] * gsl_sf_pow_int(k*r,power) << std::endl;
            }
            else if (k*r < 1.)
            {
                // Evaluate the regular Riccati-Bessel function using power series.
                // - The scaling factor is explicitly put out.
                double term = gsl_sf_pow_int(k*r,l), sum = term;
                for (int m = 1; m <= l; m++)
                {
                    term *= (l + m) * (l + 1. - m) / (2 * m * k * r);
                    sum += term;
                }
                lScaled = ipt;
                eval[ipt] = sum;
                
//                 out << r << " " << eval[ipt] << " " << eval[ipt] * gsl_sf_pow_int(k*r,power) << std::endl;
            }
            /*else if (r < rt)
            {
                // Evaluate the function within the rest of the classically forbidden region.
                // - Use library function and just divide the result by the scale factor.
                lScaled = ipt;
                eval[ipt] = special::ric_n(l,k*r) * gsl_sf_pow_int(k*r,l);
                
                out << r << " " << eval[ipt] << " " << eval[ipt] * gsl_sf_pow_int(k*r,power) << std::endl;
            }*/
            else
            {
                // Calculate the function using a library routine.
                // - No scaling is done (function is oscillatory).
                eval[ipt] = special::ric_k_scaled(l,k*r);
                
//                 out << r << " " << eval[ipt] << " " << eval[ipt] << std::endl;
            }
        }
    }
    
    /// Angular momentum.
    int l;
    
    /// Linear momentum.
    double k;
    
    /// Turning point.
    double rt;
};

struct BoundOrbital : ScaledFunction
{
    BoundOrbital (int n, int l, rArray const & grid)
    {
        power = l + 1;
        factor = 2. / n;
        lScaled = -1;
        eval.resize(grid.size());
        
        // normalization factor
        double Norm = std::sqrt(std::pow(2./n,3) * gsl_sf_fact(n-l-1) / (2 * n * gsl_sf_fact(n + l)));
        
        // grid size
        unsigned N = grid.size();
        
        // evaluate bound orbital on grid
        for (unsigned ipt = 0; ipt < N; ipt++)
        {
            // get current grid point
            double r = grid[ipt];
            
            // evaluate (add one inverse factor)
            lScaled = ipt;
            eval[ipt] = Norm * gsl_sf_laguerre_n(n - l - 1, 2 * l + 1, 2 * r / n) * std::exp(-r/n) * n / 2.;
        }
    }
};

struct FreeOrbital : ScaledFunction
{
    FreeOrbital (double k, int l, rArray const & grid)
    {
        power = l + 1;
        factor = k;
        lScaled = -1;
        eval.resize(grid.size());
        
        // normalization factor
        gsl_sf_result Norm;
        gsl_sf_coulomb_CL_e(l, -1./k, &Norm);
        Norm.val *= special::constant::sqrt_two / (k * special::constant::sqrt_pi);
        
        // grid size
        unsigned N = grid.size();
        
        // get Coulomb wave phase shift
        double sigma = special::coul_F_sigma(l, k);
        
        // evaluate bound orbital on grid
        for (unsigned ipt = 0; ipt < N; ipt++)
        {
            // get current grid point
            double r = grid[ipt];
            
            // how far is the grid point from the origin?
            if (k*r < 1)
            {
                // Evaluate from power series.
                lScaled = ipt;
                Complex term = 1., sum = term;
                for (unsigned n = 1; std::abs(term) > 1e-15 * std::abs(sum); n++)
                {
                    term *= Complex(l + n, 1./k) * Complex(0., 2*k*r) / (Complex(2*l + 1 + n, 0.) * double(n));
                    sum += term;
                }
                eval[ipt] = Norm.val * (Complex(std::cos(-k*r),std::sin(-k*r)) * sum).real();
            }
            else
            {
                // Use the library function.
                eval[ipt] = Hydrogen::F(k, l, r, sigma);
            }
        }
    }
};

struct MultipolePotentialPP : MultipolePotential
{
    MultipolePotentialPP (int lambda, int Na, int La, int Nb, int Lb, rArray const & grid)
    {
        power = lambda;
        factor = 1.;
        lScaled = 0;
        eval.resize(grid.size());
        hasMonopole = (Na == Nb and La == Lb and lambda == 0);
        
        // grid size
        unsigned N = grid.size();
        
        // normalization constants
        double Norma = std::sqrt(gsl_sf_pow_int(2./Na,3) * gsl_sf_fact(Na-La-1) / (2 * Na * gsl_sf_fact(Na + La)));
        double Normb = std::sqrt(gsl_sf_pow_int(2./Nb,3) * gsl_sf_fact(Nb-Lb-1) / (2 * Nb * gsl_sf_fact(Nb + Lb)));
        
        // exponential scale
        double c = 1./Na + 1./Nb;
        
//         std::ofstream out (format("V[%d]-%d-%d-%d-%d.out", lambda, Na, La, Nb, Lb).c_str());
        
        // evaluate potential in all grid points
        for (unsigned ipt = 0; ipt < N; ipt++)
        {
            // get current grid point
            double r = grid[ipt];
            
            // how far is the grid point from the origin?
            if (r < 1.)
            {
                // low integral will be computed numerically, with the scaled integrand
                auto integrand1 = [&](double x) -> double
                {
                    return gsl_sf_pow_int(x, La + Lb + 2 + lambda)
                         * gsl_sf_laguerre_n(Na - La - 1, 2*La + 1, 2*r*x/Na)
                         * gsl_sf_laguerre_n(Nb - Lb - 1, 2*Lb + 1, 2*r*x/Nb)
                         * gsl_sf_pow_int(2./Na,La) * gsl_sf_pow_int(2./Nb,Lb)
                         * std::exp(-r*x*(1./Na + 1./Nb));
                };
                GaussKronrod<decltype(integrand1)> Q(integrand1);
                if (not Q.integrate(0,1)) // ~ x = 0 .. r
                    HexException("Failed to evaluate bound-bound potential (low contrib).");
                double low = Q.result() * gsl_sf_pow_int(r, La + Lb + 2 - lambda);
                
                // high integral will be computed analytically
                double high = 0;
                for (int ia = 0; ia <= Na - La - 1; ia++)
                for (int ib = 0; ib <= Nb - Lb - 1; ib++)
                {
                    double afactor = (ia % 2 == 0 ? 1. : -1.) * gsl_sf_pow_int(2./Na,La+ia) * gsl_sf_choose(Na+La,Na-La-1-ia) / gsl_sf_fact(ia);
                    double bfactor = (ib % 2 == 0 ? 1. : -1.) * gsl_sf_pow_int(2./Nb,Lb+ib) * gsl_sf_choose(Nb+Lb,Nb-Lb-1-ib) / gsl_sf_fact(ib);
                    int rpower = La + 1 + Lb + 1 + ia + ib - lambda - 1;
                    double integral = gsl_sf_gamma(rpower + 1) * gsl_sf_gamma_inc_Q(rpower + 1, c*r);
                    high += afactor * bfactor * integral / gsl_sf_pow_int(c, rpower + 1);
                }
                
                // sum the contributions
                lScaled = ipt;
                eval[ipt] = Norma * Normb * (low + high);
                
//                 out << r << " " << eval[ipt] << " " << eval[ipt] * gsl_sf_pow_int(factor*r,power) << std::endl;
            }
            else
            {
                // calculate the potential without scaling; both integrals analytically
                double suma1 = 0, suma2 = 0;
                for (int ia = 0; ia <= Na - La - 1; ia++)
                for (int ib = 0; ib <= Nb - Lb - 1; ib++)
                {
                    double afactor = (ia % 2 == 0 ? 1. : -1.) * gsl_sf_pow_int(2./Na,La+ia) * gsl_sf_choose(Na+La,Na-La-1-ia) / gsl_sf_fact(ia);
                    double bfactor = (ib % 2 == 0 ? 1. : -1.) * gsl_sf_pow_int(2./Nb,Lb+ib) * gsl_sf_choose(Nb+Lb,Nb-Lb-1-ib) / gsl_sf_fact(ib);
                    int rpower1 = La + 1 + Lb + 1 + ia + ib + lambda;
                    int rpower2 = La + 1 + Lb + 1 + ia + ib - lambda - 1;
                    double integral1 = gsl_sf_gamma(rpower1 + 1) * gsl_sf_gamma_inc_P(rpower1 + 1, c*r);
                    double integral2 = gsl_sf_gamma(rpower2 + 1) * gsl_sf_gamma_inc_Q(rpower2 + 1, c*r);
                    suma1 += afactor * bfactor * integral1 / gsl_sf_pow_int(c, rpower1 + 1);
                    suma2 += afactor * bfactor * integral2 / gsl_sf_pow_int(c, rpower2 + 1);
                }
                
                eval[ipt] = Norma * Normb * (suma1 * gsl_sf_pow_int(r, -lambda-1) + suma2 * gsl_sf_pow_int(r, lambda));
                
//                 out << r << " " << eval[ipt] << " " << eval[ipt] << std::endl;
            }
        }
        
        // evaluate asymptotic potential
        asy_factor = 0;
        for (int ia = 0; ia <= Na - La - 1; ia++)
        for (int ib = 0; ib <= Nb - Lb - 1; ib++)
        {
            double afactor = (ia % 2 == 0 ? 1. : -1.) * gsl_sf_pow_int(2./Na,La+ia) * gsl_sf_choose(Na+La,Na-La-1-ia) / gsl_sf_fact(ia);
            double bfactor = (ib % 2 == 0 ? 1. : -1.) * gsl_sf_pow_int(2./Nb,Lb+ib) * gsl_sf_choose(Nb+Lb,Nb-Lb-1-ib) / gsl_sf_fact(ib);
            int rpower = La + 1 + Lb + 1 + ia + ib + lambda;
            asy_factor += afactor * bfactor * gsl_sf_gamma(rpower + 1) / gsl_sf_pow_int(c, rpower + 1);
        }
        asy_factor *= Norma * Normb;
    }
};

struct MultipolePotentialPF : MultipolePotential
{
    MultipolePotentialPF (int lambda, int Na, int La, double Kb, int Lb, rArray const & grid)
    {
        power = lambda;
        factor = 1.;
        lScaled = 0;
        eval.resize(grid.size());
        hasMonopole = false;
        
        // FIXME : Safely handle origin.
        eval = interpolate_bound_free_potential(grid, lambda, Na, La, Kb, Lb);
        
        // FIXME : evaluate asymptotic potential
        asy_factor = 0;
    }
    
    rArray interpolate_bound_free_potential (rArray const & grid, int lambda, int Na, int La, double Kb, int Lb)
    {
        // array of bound-free potential evaluations
        unsigned N = grid.size();
        double rmax = grid.back();
        rArray V(N);
        
        // the classical turning point for non-S-waves
        double rt = (std::sqrt(1 + Kb*Kb*Lb*(Lb+1))-1)/(Kb*Kb);
        
        // Coulomb function scaled by hydrogen orbital and its nodes
        rArray PF (N), zeros;
        for (unsigned i = 0; i < N; i++)
        {
            // evaluate the product of P and F
            PF[i]  = Hydrogen::P(Na,La,grid[i]) * Hydrogen::F(Kb,Lb,grid[i]);
            
            // update nodes
            if (PF[i] == 0.)
                zeros.push_back(grid[i]);
            else if (PF[i-1] != 0. and signum(PF[i-1]) != signum(PF[i]))
                zeros.push_back((grid[i-1] * PF[i] - grid[i] * PF[i-1]) / (PF[i] - PF[i-1]));
        }
        
        // cubic spline interpolation of P*F
        gsl_spline * spline = gsl_spline_alloc(gsl_interp_cspline, N);
        gsl_spline_init(spline, grid.data(), PF.data(), N);
        gsl_interp_accel * acc = gsl_interp_accel_alloc();
        
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
};

double Idir_allowed_re
(
    rArray const & grid,
    RiccatiBesselJ const & jf, MultipolePotential const & Vfn, RiccatiBesselJ const & jn,
    RiccatiBesselY const & yn, MultipolePotential const & Vni, RiccatiBesselJ const & ji
)
{
    // grid size
    int N = grid.size();
    
    static int w = 0; w++;
    
    // evaluate inner low integrand
    rArray jn_Vni_ji(N);
    for (int i = 0; i < N; i++)
    {
        // get true (unscaled) values, except for the possible monopole term in the potential
        double Vni_val = (i <= Vni.lScaled ? gsl_sf_pow_int(Vni.factor*grid[i], Vni.power) : 1.) * Vni.eval[i];
        double jn_val = (i <= jn.lScaled ? gsl_sf_pow_int(jn.factor*grid[i], jn.power) : 1.) * jn.eval[i];
        double ji_val = (i <= ji.lScaled ? gsl_sf_pow_int(ji.factor*grid[i], ji.power) : 1.) * ji.eval[i];
        
        // start with non-monopole contribution
        jn_Vni_ji[i] = jn_val * Vni_val * ji_val;
        
        // do we miss the monopole term?
        if (Vni.hasMonopole)
        {
            if (i <= std::min(jn.lScaled, ji.lScaled))
            {
                // near the origin we need to carefully combine the asymptotics
                jn_Vni_ji[i] -= gsl_sf_pow_int(jn.factor,jn.power) * gsl_sf_pow_int(ji.factor,ji.power) * gsl_sf_pow_int(grid[i],jn.power+ji.power-1) * jn.eval[i] * ji.eval[i];
            }
            else
            {
                // far from the origin we just add the monopole term
                jn_Vni_ji[i] -= jn_val * ji_val / grid[i];
            }
        }
    }
    
    // evalute inner high integrand (the real part)
    rArray yn_Vni_ji(N);
    for (int i = 0; i < N; i++)
    {
        // get true (unscaled) values, except for the possible monopole term in the potential
        double Vni_val = (i <= Vni.lScaled ? gsl_sf_pow_int(Vni.factor*grid[i], Vni.power) : 1.) * Vni.eval[i];
        double yn_val = (i <= yn.lScaled ? gsl_sf_pow_int(yn.factor*grid[i], yn.power) : 1.) * yn.eval[i];
        double ji_val = (i <= ji.lScaled ? gsl_sf_pow_int(ji.factor*grid[i], ji.power) : 1.) * ji.eval[i];
        
        // near the origin carefully combine the asymptotics
        if (i <= mmin(yn.lScaled,Vni.lScaled,ji.lScaled))
        {
            yn_Vni_ji[i] = gsl_sf_pow_int(yn.factor,yn.power) * yn.eval[i]
                         * gsl_sf_pow_int(ji.factor,ji.power) * ji.eval[i]
                         * gsl_sf_pow_int(Vni.factor,Vni.power) * Vni.eval[i]
                         * gsl_sf_pow_int(grid[i],yn.power+Vni.power+ji.power);
            
            // add monopole term (implying lambdai = 0)
            if (Vni.hasMonopole)
            {
                yn_Vni_ji[i] -= gsl_sf_pow_int(yn.factor,yn.power) * yn.eval[i]
                              * gsl_sf_pow_int(ji.factor,ji.power) * ji.eval[i]
                              * gsl_sf_pow_int(grid[i],yn.power+ji.power-1);
            }
        }
        
        // far from the origin just simply multiply the values
        else
        {
            yn_Vni_ji[i] = yn_val * Vni_val * ji_val;
            
            // add monopole term  (implying lambdai = 0)
            if (Vni.hasMonopole)
                yn_Vni_ji[i] -= yn_val * ji_val / grid[i];
        }
    }
    
    // evaluate outer low integrand
    rArray jf_Vfn_jn(N);
    for (int i = 0; i < N; i++)
    {
        // get true (unscaled) values, except for the possible monopole term in the potential
        double Vfn_val = (i <= Vfn.lScaled ? gsl_sf_pow_int(Vfn.factor*grid[i], Vfn.power) : 1.) * Vfn.eval[i];
        double jn_val = (i <= jn.lScaled ? gsl_sf_pow_int(jn.factor*grid[i], jn.power) : 1.) * jn.eval[i];
        double jf_val = (i <= jf.lScaled ? gsl_sf_pow_int(jf.factor*grid[i], jf.power) : 1.) * jf.eval[i];
        
        // start with non-monopole contribution
        jf_Vfn_jn[i] = jf_val * Vfn_val * jn_val;
        
        // do we miss the monopole term?
        if (Vfn.hasMonopole)
        {
            if (i <= std::min(jn.lScaled, jf.lScaled))
            {
                // near the origin we need to carefully combine the asymptotics
                jf_Vfn_jn[i] -= gsl_sf_pow_int(jf.factor,jf.power) * gsl_sf_pow_int(jn.factor,jn.power) * gsl_sf_pow_int(grid[i],jf.power+jn.power-1) * jf.eval[i] * jn.eval[i];
            }
            else
            {
                // far from the origin we just add the monopole term
                jf_Vfn_jn[i] -= jf_val * jn_val / grid[i];
            }
        }
    }
    
    // evalute outer high integrand (the real part)
    rArray jf_Vfn_yn(N);
    for (int i = 0; i < N; i++)
    {
        // get true (unscaled) values, except for the possible monopole term in the potential
        double Vfn_val = (i <= Vfn.lScaled ? gsl_sf_pow_int(Vfn.factor*grid[i], Vfn.power) : 1.) * Vfn.eval[i];
        double yn_val = (i <= yn.lScaled ? gsl_sf_pow_int(yn.factor*grid[i], yn.power) : 1.) * yn.eval[i];
        double jf_val = (i <= jf.lScaled ? gsl_sf_pow_int(jf.factor*grid[i], jf.power) : 1.) * jf.eval[i];
        
        // near the origin carefully combine the asymptotics
        if (i <= mmin(yn.lScaled,Vfn.lScaled,jf.lScaled))
        {
            jf_Vfn_yn[i] = gsl_sf_pow_int(yn.factor,yn.power) * yn.eval[i]
                               * gsl_sf_pow_int(jf.factor,jf.power) * jf.eval[i]
                               * gsl_sf_pow_int(Vfn.factor,Vfn.power) * Vfn.eval[i]
                               * gsl_sf_pow_int(grid[i],yn.power+Vfn.power+jf.power);
            
            // add monopole term (implying lambdaf = 0)
            if (Vfn.hasMonopole)
            {
                jf_Vfn_yn[i] -= gsl_sf_pow_int(yn.factor,yn.power) * yn.eval[i]
                              * gsl_sf_pow_int(jf.factor,jf.power) * jf.eval[i]
                              * gsl_sf_pow_int(grid[i],yn.power+jf.power-1);
            }
        }
        
        // far from the origin just simply multiply the values
        else
        {
            jf_Vfn_yn[i] = yn_val * Vfn_val * jf_val;
            
            // add monopole term  (implying lambdaf = 0)
            if (Vfn.hasMonopole)
                jf_Vfn_yn[i] -= yn_val * jf_val / grid[i];
        }
    }
    
    // working arrays and other auxiliary variables
    rArray yn_Vfn_jf_r_inf(N), yn_Vni_ji_r_inf(N);
    double sum_outer_re = 0;
    
    // calculate all high integrals over intervals (grid[i],grid.back())
    {
        gsl_spline * spline = gsl_spline_alloc(gsl_interp_cspline, N);
        gsl_spline_init(spline, grid.data(), jf_Vfn_yn.data(), N);
        gsl_interp_accel * acc = gsl_interp_accel_alloc();
        for (int i = N - 1; i >= 0; i--)
        {
            if (i == N - 1)
                yn_Vfn_jf_r_inf[i] = 0;
            else
                yn_Vfn_jf_r_inf[i] = yn_Vfn_jf_r_inf[i+1] + gsl_spline_eval_integ(spline, grid[i], grid[i+1], acc);
        }
        gsl_interp_accel_free(acc);
        gsl_spline_free(spline);
        
        // Add asymptotic multipole contribution = integral of
        //     yn jf A / r^(lambda' + 1) ~ cos(kn*r-pi*ln/2) sin(kf*r-pi*lf/2) A / r^(lambda'+1)
        // for r = grid.back() .. inf, where A = <psi_f|r^lambda'|psi_n>.
        if (Vfn.power != 0)
        {
            double U_phi = -special::constant::pi_half * (jf.l + yn.l);
            Complex U_k (0., -(jf.k + yn.k));
            Complex U = special::cis(U_phi) * std::pow(U_k, Vfn.power) * special::cfgamma(-Vfn.power, U_k * grid.back());
            
            double V_phi = -special::constant::pi_half * (jf.l - yn.l);
            Complex V_k (0., -(jf.k - yn.k));
            Complex V = (V_k.imag() == 0. ? 0. : special::cis(V_phi) * std::pow(V_k, Vfn.power) * special::cfgamma(-Vfn.power, V_k * grid.back()));
            
            yn_Vfn_jf_r_inf += 0.5 * Vfn.asy_factor * (U.imag() + V.imag());
        }
    }
    {
        gsl_spline * spline = gsl_spline_alloc(gsl_interp_cspline, N);
        gsl_spline_init(spline, grid.data(), yn_Vni_ji.data(), N);
        gsl_interp_accel * acc = gsl_interp_accel_alloc();
        for (int i = N - 1; i >= 0; i--)
        {
            if (i == N - 1)
                yn_Vni_ji_r_inf[i] = 0;
            else
                yn_Vni_ji_r_inf[i] = yn_Vni_ji_r_inf[i+1] + gsl_spline_eval_integ(spline, grid[i], grid[i+1], acc);
        }
        gsl_interp_accel_free(acc);
        gsl_spline_free(spline);
        
        // Add asymptotic multipole contribution = integral of
        //     yn ji A / r^(lambda + 1) ~ cos(kn*r-pi*ln/2) sin(ki*r-pi*li/2) A / r^(lambda+1)
        // for r = grid.back() .. inf, where A = <psi_i|r^lambda'|psi_n>.
        if (Vni.power != 0)
        {
            double U_phi = -special::constant::pi_half * (ji.l + yn.l);
            Complex U_k (0., -(ji.k + yn.k));
            Complex U = special::cis(U_phi) * std::pow(U_k, Vni.power) * special::cfgamma(-Vni.power, U_k * grid.back());
            
            double V_phi = -special::constant::pi_half * (ji.l - yn.l);
            Complex V_k (0., -(ji.k - yn.k));
            Complex V = (V_k.imag() == 0 ? 0. : special::cis(V_phi) * std::pow(V_k, Vni.power) * special::cfgamma(-Vni.power, V_k * grid.back()));
            
            yn_Vni_ji_r_inf += 0.5 * Vni.asy_factor * (U.imag() + V.imag());
        }
    }
    
    // real outer integral
    {
        rArray outer_re = jn_Vni_ji * yn_Vfn_jf_r_inf + jf_Vfn_jn * yn_Vni_ji_r_inf;
        
        gsl_interp_accel * acc = gsl_interp_accel_alloc();
        gsl_spline * spline = gsl_spline_alloc(gsl_interp_cspline, N);
        gsl_spline_init(spline, grid.data(), outer_re.data(), N);
        
        sum_outer_re = gsl_spline_eval_integ(spline, grid.front(), grid.back(), acc);
        
        // TODO : Add asymptotic contribution from region III.
        
        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);
    }
    
    return sum_outer_re;
}

double Idir_allowed_im
(
    rArray const & grid,
    RiccatiBesselJ const & jf, MultipolePotential const & Vfn, RiccatiBesselJ const & jn,
    RiccatiBesselY const & yn, MultipolePotential const & Vni, RiccatiBesselJ const & ji
)
{
    // grid size
    int N = grid.size();
    
    static int w = 0; w++;
    
    // evaluate inner low integrand
    rArray jn_Vni_ji(N);
    for (int i = 0; i < N; i++)
    {
        // get true (unscaled) values, except for the possible monopole term in the potential
        double Vni_val = (i <= Vni.lScaled ? gsl_sf_pow_int(Vni.factor*grid[i], Vni.power) : 1.) * Vni.eval[i];
        double jn_val = (i <= jn.lScaled ? gsl_sf_pow_int(jn.factor*grid[i], jn.power) : 1.) * jn.eval[i];
        double ji_val = (i <= ji.lScaled ? gsl_sf_pow_int(ji.factor*grid[i], ji.power) : 1.) * ji.eval[i];
        
        // start with non-monopole contribution
        jn_Vni_ji[i] = jn_val * Vni_val * ji_val;
        
        // do we miss the monopole term?
        if (Vni.hasMonopole)
        {
            if (i <= std::min(jn.lScaled, ji.lScaled))
            {
                // near the origin we need to carefully combine the asymptotics
                jn_Vni_ji[i] -= gsl_sf_pow_int(jn.factor,jn.power) * gsl_sf_pow_int(ji.factor,ji.power) * gsl_sf_pow_int(grid[i],jn.power+ji.power-1) * jn.eval[i] * ji.eval[i];
            }
            else
            {
                // far from the origin we just add the monopole term
                jn_Vni_ji[i] -= jn_val * ji_val / grid[i];
            }
        }
    }
    
    // evaluate outer low integrand
    rArray jf_Vfn_jn(N);
    for (int i = 0; i < N; i++)
    {
        // get true (unscaled) values, except for the possible monopole term in the potential
        double Vfn_val = (i <= Vfn.lScaled ? gsl_sf_pow_int(Vfn.factor*grid[i], Vfn.power) : 1.) * Vfn.eval[i];
        double jn_val = (i <= jn.lScaled ? gsl_sf_pow_int(jn.factor*grid[i], jn.power) : 1.) * jn.eval[i];
        double jf_val = (i <= jf.lScaled ? gsl_sf_pow_int(jf.factor*grid[i], jf.power) : 1.) * jf.eval[i];
        
        // start with non-monopole contribution
        jf_Vfn_jn[i] = jf_val * Vfn_val * jn_val;
        
        // do we miss the monopole term?
        if (Vfn.hasMonopole)
        {
            if (i <= std::min(jn.lScaled, jf.lScaled))
            {
                // near the origin we need to carefully combine the asymptotics
                jf_Vfn_jn[i] -= gsl_sf_pow_int(jf.factor,jf.power) * gsl_sf_pow_int(jn.factor,jn.power) * gsl_sf_pow_int(grid[i],jf.power+jn.power-1) * jf.eval[i] * jn.eval[i];
            }
            else
            {
                // far from the origin we just add the monopole term
                jf_Vfn_jn[i] -= jf_val * jn_val / grid[i];
            }
        }
    }
    
    // working arrays and other auxiliary variables
    double sum_outer_im = 0;
    
    // imaginary outer integral
    {
        gsl_interp_accel * acc = gsl_interp_accel_alloc();
        gsl_spline * spline = gsl_spline_alloc(gsl_interp_cspline, N);
        
        gsl_spline_init(spline, grid.data(), jf_Vfn_jn.data(), N);
        double integ_fn = gsl_spline_eval_integ(spline, grid.front(), grid.back(), acc);
        
        // Add asymptotic multipole contribution = integral of
        //     jf jn A / r^(lambda' + 1)
        // for r = grid.back() .. inf.
        if (Vfn.power != 0)
        {
            double U_phi = -special::constant::pi_half * (jn.l - jf.l);
            Complex U_k (0., -(jn.k - jf.k));
            Complex U = (U_k.imag() == 0. ? 0. : special::cis(U_phi) * std::pow(U_k, Vfn.power) * special::cfgamma(-Vfn.power, U_k * grid.back()));
            
            double V_phi = -special::constant::pi_half * (jn.l + jf.l);
            Complex V_k (0., -(jn.k + jf.k));
            Complex V = special::cis(V_phi) * std::pow(V_k, Vfn.power) * special::cfgamma(-Vfn.power, V_k * grid.back());
            
            integ_fn += 0.5 * Vfn.asy_factor * (U.real() - V.real());
        }
        
        gsl_spline_init(spline, grid.data(), jn_Vni_ji.data(), N);
        double integ_ni = gsl_spline_eval_integ(spline, grid.front(), grid.back(), acc);
        
        // Add asymptotic multipole contribution = integral of
        //     jn ji A / r^(lambda + 1)
        // for r = grid.back() .. inf.
        if (Vni.power != 0)
        {
            double U_phi = -special::constant::pi_half * (jn.l - ji.l);
            Complex U_k (0., -(jn.k - ji.k));
            Complex U = (U_k.imag() == 0. ? 0. : special::cis(U_phi) * std::pow(U_k, Vni.power) * special::cfgamma(-Vni.power, U_k * grid.back()));
            
            double V_phi = -special::constant::pi_half * (jn.l + ji.l);
            Complex V_k (0., -(jn.k + ji.k));
            Complex V = special::cis(V_phi) * std::pow(V_k, Vni.power) * special::cfgamma(-Vni.power, V_k * grid.back());
            integ_ni += 0.5 * Vni.asy_factor * (U.real() - V.real());
        }
        
        sum_outer_im = integ_fn * integ_ni;
        
        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);
    }
    
    return sum_outer_im;
}

double Idir_forbidden
(
    rArray const & grid,
    RiccatiBesselJ const & jf, MultipolePotential const & Vfn,
    RiccatiBesselI const & in, RiccatiBesselK const & kn,
    RiccatiBesselJ const & ji, MultipolePotential const & Vni
)
{
    int N = grid.size();
    double h = grid.back() / N;
    
    double suma = 0;
    int last_diagonal = N/2;
    double scale, ji_ev, Vfn_ev, jf_ev, Vni_ev, in_ev, kn_ev;
    
    // index "d" runs across the contours x - y = konst
    // TODO : use spline integration along and across the diagonals
    for (int d = 0; d < last_diagonal; d++)
    {
        // contribution of these contours
        double contrib = 0;
        
        if (d == 0)
        {
            // index "i" runs along the diagonal [x,y] = [i,i]
            for (int i = 1; i < N; i++)
            {
                if (mmin(jf.lScaled,Vfn.lScaled,ji.lScaled,Vni.lScaled,in.lScaled,kn.lScaled) >= i)
                {
                    // combine scale factors
                    scale = gsl_sf_pow_int(jf.factor, jf.power) * gsl_sf_pow_int(Vfn.factor, Vfn.power) *
                            gsl_sf_pow_int(ji.factor, ji.power) * gsl_sf_pow_int(Vni.factor, Vni.power) *
                            gsl_sf_pow_int(in.factor, in.power) * gsl_sf_pow_int(kn.factor, kn.power) *
                            gsl_sf_pow_int(grid[i], jf.power + Vfn.power + ji.power + Vni.power + in.power + kn.power);
                    ji_ev = ji.eval[i]; jf_ev = jf.eval[i]; Vfn_ev = Vfn.eval[i]; Vni_ev = Vni.eval[i]; in_ev = in.eval[i]; kn_ev = kn.eval[i];
                }
                else
                {
                    // use scale factors
                    scale = 1;
                    jf_ev  = (i <= jf.lScaled  ? gsl_sf_pow_int(jf.factor  * grid[i], jf.power)  * jf.eval[i]  : jf.eval[i]);
                    Vfn_ev = (i <= Vfn.lScaled ? gsl_sf_pow_int(Vfn.factor * grid[i], Vfn.power) * Vfn.eval[i] : Vfn.eval[i]);
                    in_ev  = (i <= in.lScaled  ? gsl_sf_pow_int(in.factor  * grid[i], in.power)  * in.eval[i]  : in.eval[i]);
                    kn_ev  = (i <= kn.lScaled  ? gsl_sf_pow_int(kn.factor  * grid[i], kn.power)  * kn.eval[i]  : kn.eval[i]);
                    Vni_ev = (i <= Vni.lScaled ? gsl_sf_pow_int(Vni.factor * grid[i], Vni.power) * Vni.eval[i] : Vni.eval[i]);
                    ji_ev  = (i <= ji.lScaled  ? gsl_sf_pow_int(ji.factor  * grid[i], ji.power)  * ji.eval[i]  : ji.eval[i]);
                }
                contrib += scale * jf_ev * Vfn_ev * ji_ev * Vni_ev * in_ev * kn_ev;
            }
        }
        else
        {
            // index "i" runs along the contour [x,y] = [i,i-d]
            for (int i = d + 1; i < N; i++)
            {
                if (mmin(jf.lScaled+d,Vfn.lScaled+d,ji.lScaled,Vni.lScaled,in.lScaled+d,kn.lScaled) >= i)
                {
                    // combine scale factors
                    scale = gsl_sf_pow_int(jf.factor, jf.power) * gsl_sf_pow_int(Vfn.factor, Vfn.power) *
                            gsl_sf_pow_int(ji.factor, ji.power) * gsl_sf_pow_int(Vni.factor, Vni.power) *
                            gsl_sf_pow_int(in.factor, in.power) * gsl_sf_pow_int(kn.factor, kn.power) *
                            gsl_sf_pow_int(grid[i]/grid[i-d], ji.power + Vni.power + kn.power) *
                            gsl_sf_pow_int(grid[i-d], jf.power + Vfn.power + in.power + kn.power + Vni.power + ji.power);
                    ji_ev = ji.eval[i]; jf_ev = jf.eval[i-d]; Vfn_ev = Vfn.eval[i-d]; Vni_ev = Vni.eval[i]; in_ev = in.eval[i-d]; kn_ev = kn.eval[i];
                }
                else
                {
                    // use scale factors
                    scale = 1;
                    jf_ev  = (i-d <= jf.lScaled  ? gsl_sf_pow_int(jf.factor  * grid[i-d], jf.power)  * jf.eval[i-d]  : jf.eval[i-d]);
                    Vfn_ev = (i-d <= Vfn.lScaled ? gsl_sf_pow_int(Vfn.factor * grid[i-d], Vfn.power) * Vfn.eval[i-d] : Vfn.eval[i-d]);
                    in_ev  = (i-d <= in.lScaled  ? gsl_sf_pow_int(in.factor  * grid[i-d], in.power)  * in.eval[i-d]  : in.eval[i-d]);
                    kn_ev  = (i   <= kn.lScaled  ? gsl_sf_pow_int(kn.factor  * grid[i],   kn.power)  * kn.eval[i]    : kn.eval[i]);
                    Vni_ev = (i   <= Vni.lScaled ? gsl_sf_pow_int(Vni.factor * grid[i],   Vni.power) * Vni.eval[i]   : Vni.eval[i]);
                    ji_ev  = (i   <= ji.lScaled  ? gsl_sf_pow_int(ji.factor  * grid[i],   ji.power)  * ji.eval[i]    : ji.eval[i]);
                }
                contrib += scale * jf_ev * Vfn_ev * ji_ev * Vni_ev * in_ev * kn_ev;
            }
            // index "i" runs along the contour [x,y] = [i,i+d]
            for (int i = 1; i < N - d - 1; i++)
            {
                if (mmin(jf.lScaled-d,Vfn.lScaled-d,ji.lScaled,Vni.lScaled,in.lScaled,kn.lScaled-d) >= i)
                {
                    // combine scale factors
                    scale = gsl_sf_pow_int(jf.factor, jf.power) * gsl_sf_pow_int(Vfn.factor, Vfn.power) *
                            gsl_sf_pow_int(ji.factor, ji.power) * gsl_sf_pow_int(Vni.factor, Vni.power) *
                            gsl_sf_pow_int(in.factor, in.power) * gsl_sf_pow_int(kn.factor, kn.power) *
                            gsl_sf_pow_int(grid[i+d]/grid[i], jf.power + Vfn.power + kn.power) *
                            gsl_sf_pow_int(grid[i], jf.power + Vfn.power + ji.power + Vni.power + in.power + kn.power);
                    ji_ev = ji.eval[i]; jf_ev = jf.eval[i+d]; Vfn_ev = Vfn.eval[i+d]; Vni_ev = Vni.eval[i]; in_ev = in.eval[i]; kn_ev = kn.eval[i+d];
                }
                else
                {
                    // discard scale factor
                    scale = 1;
                    jf_ev  = (i+d <= jf.lScaled  ? gsl_sf_pow_int(jf.factor  * grid[i+d], jf.power)  * jf.eval[i+d]  : jf.eval[i+d]);
                    Vfn_ev = (i+d <= Vfn.lScaled ? gsl_sf_pow_int(Vfn.factor * grid[i+d], Vfn.power) * Vfn.eval[i+d] : Vfn.eval[i+d]);
                    kn_ev  = (i+d <= kn.lScaled  ? gsl_sf_pow_int(kn.factor  * grid[i+d], kn.power)  * kn.eval[i+d]  : kn.eval[i+d]);
                    in_ev  = (i   <= in.lScaled  ? gsl_sf_pow_int(in.factor  * grid[i],   in.power)  * in.eval[i]    : in.eval[i]);
                    Vni_ev = (i   <= Vni.lScaled ? gsl_sf_pow_int(Vni.factor * grid[i],   Vni.power) * Vni.eval[i]   : Vni.eval[i]);
                    ji_ev  = (i   <= ji.lScaled  ? gsl_sf_pow_int(ji.factor  * grid[i],   ji.power)  * ji.eval[i]    : ji.eval[i]);
                }
                contrib += scale * jf_ev * Vfn_ev * ji_ev * Vni_ev * in_ev * kn_ev;
            }
        }
        
        // add properly scaled contribution
        contrib *= std::exp(-in.k*h*d);
        suma += contrib;
        
        // check convergence (do not blindly integrate everything)
        if (std::abs(contrib) < 1e-15 * std::abs(suma))
            last_diagonal = d;
    }
    
    return suma * h * h;
}

double Idir_forbidden_old
(
    rArray const & grid,
    RiccatiBesselJ const & jf, MultipolePotential const & Vfn,
    RiccatiBesselI const & in, RiccatiBesselK const & kn,
    RiccatiBesselJ const & ji, MultipolePotential const & Vni
)
{
    int N = grid.size();
    double h = grid.back() / N;
    
    rArray fpart(N),  ipart(N);
    
    for (int i = 0; i < N; i++)
    {
        double jf_ev  = (i <= jf.lScaled  ? gsl_sf_pow_int(jf.factor  * grid[i], jf.power)  * jf.eval[i]  : jf.eval[i]);
        double Vfn_ev = (i <= Vfn.lScaled ? gsl_sf_pow_int(Vfn.factor * grid[i], Vfn.power) * Vfn.eval[i] : Vfn.eval[i]);
        double Vni_ev = (i <= Vni.lScaled ? gsl_sf_pow_int(Vni.factor * grid[i], Vni.power) * Vni.eval[i] : Vni.eval[i]);
        double ji_ev  = (i <= ji.lScaled  ? gsl_sf_pow_int(ji.factor  * grid[i], ji.power)  * ji.eval[i]  : ji.eval[i]);
        
        fpart[i] = jf_ev * Vfn_ev;
        ipart[i] = ji_ev * Vni_ev;
    }
    
    double suma = 0;
    int last_diagonal = N/2;
    
    // index "d" runs across the contours x - y = konst
    // TODO : use spline integration along and across the diagonals
    for (int d = 0; d < last_diagonal; d++)
    {
        // contribution of these contours
        double contrib = 0;
        
        if (d == 0)
        {
            // index "i" runs along the diagonal [x,y] = [i,i]
            for (int i = 0; i < N; i++)
            {
                double in_ev  = (i <= in.lScaled  ? gsl_sf_pow_int(in.factor  * grid[i], in.power)  * in.eval[i]  : in.eval[i]);
                double kn_ev  = (i <= kn.lScaled  ? gsl_sf_pow_int(kn.factor  * grid[i],   kn.power)  * kn.eval[i]    : kn.eval[i]);
                contrib += fpart[i] * ipart[i] * in_ev * kn_ev;
            }
        }
        else
        {
            // index "i" runs along the contour [x,y] = [i,i-d]
            for (int i = d; i < N; i++)
            {
                double in_ev  = (i-d <= in.lScaled  ? gsl_sf_pow_int(in.factor  * grid[i-d], in.power)  * in.eval[i-d]  : in.eval[i-d]);
                double kn_ev  = (i <= kn.lScaled  ? gsl_sf_pow_int(kn.factor  * grid[i],   kn.power)  * kn.eval[i]    : kn.eval[i]);
                contrib += fpart[i-d] * ipart[i] * in_ev * kn_ev;
            }
            // index "i" runs along the contour [x,y] = [i,i+d]
            for (int i = 0; i < N - d; i++)
            {
                double in_ev  = (i <= in.lScaled  ? gsl_sf_pow_int(in.factor  * grid[i], in.power)  * in.eval[i]  : in.eval[i]);
                double kn_ev  = (i+d <= kn.lScaled  ? gsl_sf_pow_int(kn.factor  * grid[i+d],   kn.power)  * kn.eval[i+d]    : kn.eval[i+d]);
                contrib += fpart[i+d] * ipart[i] * in_ev * kn_ev;
            }
        }
        
        // add properly scaled contribution
        contrib *= std::exp(-kn.k*h*d);
        suma += contrib;
        
        // check convergence (do not blindly integrate everything)
        if (std::abs(contrib) < 1e-15 * std::abs(suma))
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
    // resulting value of the double integral
    Complex result = 0;
    
    // final angular momentum transfer bounds
    int lambdaf_min = std::max(std::abs(lf-ln), std::abs(Lf-Ln));
    int lambdaf_max = std::min(lf+ln, Lf+Ln);
    
    // initial angular momentum transfer bounds
    int lambdai_min = std::max(std::abs(li-ln), std::abs(Li-Ln));
    int lambdai_max = std::min(li+ln, Li+Ln);
    
    // for all angular momentum transfers
    for (int lambdaf = lambdaf_min; lambdaf <= lambdaf_max; lambdaf++)
    for (int lambdai = lambdai_min; lambdai <= lambdai_max; lambdai++)
    {
        // calculate angular integrals
        double ff = special::computef(lambdaf, Lf, lf, Ln, ln, L);
        double fi = special::computef(lambdai, Ln, ln, Li, li, L);
        
        // skip non-contributing transfers
        if (ff == 0. or fi == 0.)
            continue;
        
        // evaluate initial and final radial part of the projectile partial wave
        RiccatiBesselJ ji (li, ki, grid);
        RiccatiBesselJ jf (lf, kf, grid);
        
        // evaluate bound states
        BoundOrbital Pi (Ni, Li, grid);
        BoundOrbital Pn (Nn, Ln, grid);
        BoundOrbital Pf (Nf, Lf, grid);
        
        // evaluate initial and final multipole potential
        MultipolePotentialPP Vfn (lambdaf, Nf, Lf, Nn, Ln, grid);
        MultipolePotentialPP Vni (lambdai, Nn, Ln, Ni, Li, grid);
        
        // evaluate the Green function terms - the Bessel functions "jn" and "yn"
        RiccatiBesselJ jn (ln, kn, grid);
        RiccatiBesselY yn (ln, kn, grid);
        
        // integrate
        Complex inte
        (
            -Idir_allowed_re(grid, jf, Vfn, jn, yn, Vni, ji),
            -Idir_allowed_im(grid, jf, Vfn, jn, yn, Vni, ji)
        );
        
        // update result (NOTE : assuming kn > 0)
        result += ff * fi * inte / kn;
        
        // comment this result
        log << format
        (
            "\t\ttransfer [%d %d] initial (%d %d, %g %d), intermediate (%d %d, %g %d) final (%d %d, %g %d) : (%g,%g)",
            lambdai, lambdaf, Ni, Li, ki, li, Nn, Ln, kn, ln, Nf, Lf, kf, lf, inte.real(), inte.imag()
        ) << std::endl;
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
    // resulting value of the double integral
    double result = 0;
    
    // final angular momentum transfer bounds
    int lambdaf_min = std::max(std::abs(lf-ln), std::abs(Lf-Ln));
    int lambdaf_max = std::min(lf+ln, Lf+Ln);
    
    // initial angular momentum transfer bounds
    int lambdai_min = std::max(std::abs(li-ln), std::abs(Li-Ln));
    int lambdai_max = std::min(li+ln, Li+Ln);
    
    // for all angular momentum transfers
    for (int lambdaf = lambdaf_min; lambdaf <= lambdaf_max; lambdaf++)
    for (int lambdai = lambdai_min; lambdai <= lambdai_max; lambdai++)
    {
        // calculate angular integrals
        double ff = special::computef(lambdaf, Lf, lf, Ln, ln, L);
        double fi = special::computef(lambdai, Ln, ln, Li, li, L);
        
        // skip non-contributing transfers
        if (ff == 0. or fi == 0.)
            continue;
        
        // evaluate initial / final radial part of the projectile partial wave
        RiccatiBesselJ ji (li, ki, grid);
        RiccatiBesselJ jf (lf, kf, grid);
        
        // evaluate initial / final multipole potential
        MultipolePotentialPP Vfn (lambdaf, Nf, Lf, Nn, Ln, grid);
        MultipolePotentialPP Vni (lambdai, Nn, Ln, Ni, Li, grid);
        
        // evaluate the Green function terms (Bessel functions "in" and "kn")
        RiccatiBesselI iscaled_n (ln, kappan, grid);
        RiccatiBesselK kscaled_n (ln, kappan, grid);
        
        // finally integrate
        double inte = -Idir_forbidden(grid, jf, Vfn, iscaled_n, kscaled_n, ji, Vni);
        
        // update result (assume kappan > 0)
        result += ff * fi * inte / kappan;
        
        // comment this result
        log << format
        (
            "\t\ttransfer [%d %d] initial (%d %d, %g %d) intermediate (%d %d, %g %d) final (%d %d, %g %d) : (%g,%g)",
            lambdai, lambdaf, Ni, Li, ki, li, Nn, Ln, kappan, ln, Nf, Lf, kf, lf, inte, 0.
        ) << std::endl;
    }
    
    return result;
}

double Idir_nFree_allowed_re
(
    rArray const & grid, int L,
    int Nf, int Lf, double kf, int lf,
    double Kn, int Ln, double kn, int ln,
    int Ni, int Li, double ki, int li,
    std::ostream & log
)
{
    double result = 0;
    
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
        RiccatiBesselJ ji (li, ki, grid);
        RiccatiBesselJ jf (lf, kf, grid);
        
        // evaluate initial / final multipole potential
        MultipolePotentialPF Vfn (lambdaf, Nf, Lf, Kn, Ln, grid);
        MultipolePotentialPF Vni (lambdai, Ni, Li, Kn, Ln, grid);
        
        // evaluate the Green function terms (Bessel functions "jn" and "yn")
        RiccatiBesselJ jn (ln, kn, grid);
        RiccatiBesselY yn (ln, kn, grid);
        
        // integrate
        double inte = -Idir_allowed_re (grid, jf, Vfn, jn, yn, Vni, ji) / kn;
        
        // if we cancelled small kn, we need to throw away nonsensial imaginary part, where this can't be done
//         if (kn_is_small)
//             inte.imag(0.);
        
        // also, the real part may have exploded
//         if (not std::isfinite(inte.real()))
//             inte.real(0.);
        
        // comment this result
        log << format
        (
            "\t\ttransfer [%d %d] initial (%d %d, %g %d) intermediate (%g %d, %g %d) final (%d %d, %g %d) : real %g%s",
            lambdai, lambdaf, Ni, Li, ki, li, Kn, Ln, kn, ln, Nf, Lf, kf, lf,
            std::isfinite(inte) ? inte : 0.,
            std::isfinite(inte) ? "" : " *"
        ) << std::endl;
        
        // update result
        result += ff * fi * inte;
    }
    
    return result;
}

double Idir_nFree_allowed_im
(
    rArray const & grid, int L,
    int Nf, int Lf, double kf, int lf,
    double Kn, int Ln, double kn, int ln,
    int Ni, int Li, double ki, int li,
    std::ostream & log
)
{
    double result = 0;
    
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
        RiccatiBesselJ ji (li, ki, grid);
        RiccatiBesselJ jf (lf, kf, grid);
        
        // evaluate initial / final multipole potential
        MultipolePotentialPF Vfn (lambdaf, Nf, Lf, Kn, Ln, grid);
        MultipolePotentialPF Vni (lambdai, Ni, Li, Kn, Ln, grid);
        
        // evaluate the Green function terms (Bessel functions "jn" and "yn")
        RiccatiBesselJ jn (ln, kn, grid);
        RiccatiBesselY yn (ln, kn, grid);
        
        // integrate
        double inte = -Idir_allowed_im (grid, jf, Vfn, jn, yn, Vni, ji) / kn;
        
        // if we cancelled small kn, we need to throw away nonsensial imaginary part, where this can't be done
//         if (kn_is_small)
//             inte.imag(0.);
        
        // also, the real part may have exploded
//         if (not std::isfinite(inte.real()))
//             inte.real(0.);
        
        // comment this result
        log << format
        (
            "\t\ttransfer [%d %d] initial (%d %d, %g %d) intermediate (%g %d, %g %d) final (%d %d, %g %d) : imag %g%s",
            lambdai, lambdaf, Ni, Li, ki, li, Kn, Ln, kn, ln, Nf, Lf, kf, lf,
            std::isfinite(inte) ? inte : 0.,
            std::isfinite(inte) ? "" : " *"
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
        RiccatiBesselJ ji (li, ki, grid);
        RiccatiBesselJ jf (lf, kf, grid);
        
        // evaluate initial / final multipole potential
        MultipolePotentialPF Vfn (lambdaf, Nf, Lf, Kn, Ln, grid);
        MultipolePotentialPF Vni (lambdai, Ni, Li, Kn, Ln, grid);
        
        // evaluate the Green function terms (Bessel functions "in" and "kn")
        RiccatiBesselI iscaled_n (ln, kappan, grid);
        RiccatiBesselK kscaled_n (ln, kappan, grid);
        
        // finally integrate
        double inte = -Idir_forbidden(grid, jf, Vfn, iscaled_n, kscaled_n, ji, Vni) / kappan;
        
        // comment this result
        log << format
        (
            "\t\ttransfer [%d %d] initial (%d %d, %g %d) intermediate (%g %d, %g %d) final (%d %d, %g %d) : %g",
            lambdai, lambdaf, Ni, Li, ki, li, Kn, Ln, kappan, ln, Nf, Lf, kf, lf, inte
        ) << std::endl;
        
        // update result
        result += ff * fi * inte;
    }
    
    return result;
}
