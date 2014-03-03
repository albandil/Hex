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

#include <iostream>

#include "arrays.h"
#include "radial.h"
#include "romberg.h"
#include "specf.h"

double compute_Idir (int li, int lf, int lambda, int Ni, int Li, double ki, int Nf, int Lf, double kf)
{
    if (lambda == 0)
    {
        //
        // r1 > r2
        //
        
        double skew = 2;
        auto kernel = [&](rArray const & u) -> double
        {
            // get subgrid size
            unsigned n = u.size();
            
            // unscale coordinates
            rArray x = u.transform([&](double t) -> double { return skew * t / (1.-t); });
            
            // precompute radial functions
            rArray xi  = x.transform ([&](double r) -> double { return hydro_P(Ni, Li, r); });
            rArray xf  = x.transform ([&](double r) -> double { return hydro_P(Nf, Lf, r); });
            rArray ji  = x.transform ([&](double r) -> double { return ric_j(li, ki*r);    });
            rArray jf  = x.transform ([&](double r) -> double { return ric_j(lf, kf*r);    });
            
            // sum of evaluations
            double suma = 0.;
            
            // evaluate integrand on the carthesian product x × x
            for (unsigned ir1 = 0; ir1 < n; ir1++)
            for (unsigned ir2 = 0; ir2 < ir1; ir2++)
                suma += xi[ir1]*xf[ir1]*ji[ir2]*jf[ir2]*(1./x[ir1] - 1./x[ir2]) / ((1.-u[ir1])*(1.-u[ir1])*(1.-u[ir2])*(1.-u[ir2]));
            
            // returh the aggregated result
            return skew * skew * suma;
        };
        
        UnitSquareRomberg<double,decltype(kernel)> R(kernel);
        R.setEpsAbs(1e-10);
        R.setEpsRel(1e-5);
        R.setMinLevel(3);
        R.setMaxLevel(10);
        R.setMaxRombLevel(1);
        R.integrate_extern();
        
        if (not R.ok())
        {
            std::cerr << format
            (
                "compute_Idir failed for λ=%d, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d (\"%s\").\n\tresult = %g\n",
                lambda, Ni, Li, ki, li, Nf, Lf, kf, lf, R.status().c_str(), R.result()
            );
        }
        
        return R.result();
    }
    else
    {
        double skew = 2;
        auto kernel = [&](rArray const & u) -> double
        {
            // get subgrid size
            unsigned n = u.size();
            
            // unscale coordinates
            rArray x = u.transform([&](double t) -> double { return skew * t / (1.-t); });
            
            // precompute radial functions
            rArray xi  = x.transform ([&](double r) -> double { return hydro_P(Ni, Li, r); });
            rArray xf  = x.transform ([&](double r) -> double { return hydro_P(Nf, Lf, r); });
            rArray ji  = x.transform ([&](double r) -> double { return ric_j(li, ki*r);    });
            rArray jf  = x.transform ([&](double r) -> double { return ric_j(lf, kf*r);    });
            
            // sum of evaluations
            double suma = 0.;
            
            // evaluate integrand on the carthesian product x × x
            for (unsigned ir1 = 0; ir1 < n; ir1++)
            {
                for (unsigned ir2 = 0; ir2 < ir1; ir2++)
                    suma += xi[ir1]*xf[ir1]*ji[ir2]*jf[ir2]*pow(x[ir2]/x[ir1],lambda)/x[ir1] / ((1.-u[ir1])*(1.-u[ir1])*(1.-u[ir2])*(1.-u[ir2]));
                for (unsigned ir2 = ir1; ir2 < n; ir2++)
                    suma += xi[ir1]*xf[ir1]*ji[ir2]*jf[ir2]*pow(x[ir1]/x[ir2],lambda)/x[ir2] / ((1.-u[ir1])*(1.-u[ir1])*(1.-u[ir2])*(1.-u[ir2]));
            }
            
            // returh the aggregated result
            return skew * skew * suma;
        };
        
        UnitSquareRomberg<double,decltype(kernel)> R(kernel);
        R.setEpsAbs(1e-10);
        R.setEpsRel(1e-5);
        R.setMinLevel(3);
        R.setMaxLevel(10);
        R.setMaxRombLevel(1);
        R.integrate_extern();
        
        if (not R.ok())
        {
            std::cerr << format
            (
                "compute_Idir failed for λ=%d, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d (\"%s\").\n\tresult = %g\n",
                lambda, Ni, Li, ki, li, Nf, Lf, kf, lf, R.status().c_str(), R.result()
            );
        }
        
        return R.result();
    }
}

double compute_Iexc (int li, int lf, int lambda, int Ni, int Li, double ki, int Nf, int Lf, double kf)
{
    if (lambda == 0)
    {
        //
        // r1 > r2
        //
        
        double skew = 2;
        auto kernel = [&](rArray const & u) -> double
        {
            // get subgrid size
            unsigned n = u.size();
            
            // unscale coordinates
            rArray x = u.transform([&](double t) -> double { return skew * t / (1.-t); });
            
            // precompute radial functions
            rArray xi  = x.transform ([&](double r) -> double { return hydro_P(Ni, Li, r); });
            rArray xf  = x.transform ([&](double r) -> double { return hydro_P(Nf, Lf, r); });
            rArray ji  = x.transform ([&](double r) -> double { return ric_j(li, ki*r);    });
            rArray jf  = x.transform ([&](double r) -> double { return ric_j(lf, kf*r);    });
            
            // sum of evaluations
            double suma = 0.;
            
            // evaluate integrand on the carthesian product x × x
            for (unsigned ir1 = 0; ir1 < n; ir1++)
            for (unsigned ir2 = 0; ir2 < ir1; ir2++)
                suma += xi[ir1]*xf[ir2]*ji[ir2]*jf[ir1]*(1./x[ir1] - 1./x[ir2]) / ((1.-u[ir1])*(1.-u[ir1])*(1.-u[ir2])*(1.-u[ir2]));
            
            // returh the aggregated result
            return skew * skew * suma;
        };
        
        UnitSquareRomberg<double,decltype(kernel)> R(kernel);
        R.setEpsAbs(1e-10);
        R.setEpsRel(1e-5);
        R.setMinLevel(3);
        R.setMaxLevel(10);
        R.setMaxRombLevel(1);
        R.integrate_extern();
        
        if (not R.ok())
        {
            std::cerr << format
            (
                "compute_Iexc failed for λ=%d, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d (\"%s\").\nresult = %g\n",
                lambda, Ni, Li, ki, li, Nf, Lf, kf, lf, R.status().c_str(), R.result()
            );
        }
        
        return R.result();
    }
    else
    {
        double skew = 2;
        auto kernel = [&](rArray const & u) -> double
        {
            // get subgrid size
            unsigned n = u.size();
            
            // unscale coordinates
            rArray x = u.transform([&](double t) -> double { return skew * t / (1.-t); });
            
            // precompute radial functions
            rArray xi  = x.transform ([&](double r) -> double { return hydro_P(Ni, Li, r); });
            rArray xf  = x.transform ([&](double r) -> double { return hydro_P(Nf, Lf, r); });
            rArray ji  = x.transform ([&](double r) -> double { return ric_j(li, ki*r);    });
            rArray jf  = x.transform ([&](double r) -> double { return ric_j(lf, kf*r);    });
            
            // sum of evaluations
            double suma = 0.;
            
            // evaluate integrand on the carthesian product x × x
            for (unsigned ir1 = 0; ir1 < n; ir1++)
            {
                for (unsigned ir2 = 0; ir2 < ir1; ir2++)
                    suma += xi[ir1]*xf[ir2]*ji[ir2]*jf[ir1]*pow(x[ir2]/x[ir1],lambda)/x[ir1] / ((1.-u[ir1])*(1.-u[ir1])*(1.-u[ir2])*(1.-u[ir2]));
                for (unsigned ir2 = ir1; ir2 < n; ir2++)
                    suma += xi[ir1]*xf[ir2]*ji[ir2]*jf[ir1]*pow(x[ir1]/x[ir2],lambda)/x[ir2] / ((1.-u[ir1])*(1.-u[ir1])*(1.-u[ir2])*(1.-u[ir2]));
            }
            
            // returh the aggregated result
            return skew * skew * suma;
        };
        
        UnitSquareRomberg<double,decltype(kernel)> R(kernel);
        R.setEpsAbs(1e-10);
        R.setEpsRel(1e-5);
        R.setMinLevel(3);
        R.setMaxLevel(10);
        R.setMaxRombLevel(1);
        R.integrate_extern();
        
        if (not R.ok())
        {
            std::cerr << format
            (
                "compute_Iexc failed for λ=%d, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d (\"%s\").\n\tresult = %g\n",
                lambda, Ni, Li, ki, li, Nf, Lf, kf, lf, R.status().c_str(), R.result()
            );
        }
        
        return R.result();
    }
}
