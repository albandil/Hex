/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
*                                                                           *
*                       / /   / /    __    \ \  / /                         *
*                      / /__ / /   / _ \    \ \/ /                          *
*                     /  ___  /   | |/_/    / /\ \                          *
*                    / /   / /    \_\      / /  \ \                         *
*                                                                           *
*                         Jakub Benda (c) 2013                              *
*                     Charles University in Prague                          *
*                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <o2scl/ode_funct.h>
#include <o2scl/ode_iv_solve.h>

#include <gsl/gsl_sf.h>

#include "arrays.h"
#include "complex.h"
#include "wave_irreg.h"

// -------------------------------------------------------------------------- //
//                                                                            //
//  IrregularWave class members                                               //
//                                                                            //
// -------------------------------------------------------------------------- //

IrregularWave::IrregularWave(double _kn, int _ln, DistortingPotential const & _U)
    : Evaluations(0), U(_U), kn(_kn), ln(_ln), r0(sqrt(ln*(ln+1)) / kn), rf(U.getFarRadius())
{
    // determine discretization in the classically allowed region
    int N = 1000;               // N samples per wave length
    this->h = std::min (        // grid step
        2*M_PI/(N*kn),          //  -> N samples per wave length and
        (rf-r0)/1000            //  -> at least 1000 samples total
    );
    samples = (rf - r0)/h + 1;  // with both boundaries
    
    // determine discretization in the classically forbidden region
    samples0 = r0 * 100;         // FIXME some heuristic?
    
    cArray array, array0; // auxiliary arrays for HDF I/O
    
    // create grids
    grid = rArray(samples);
    grid0 = rArray(samples0);
    for (int i = 0; i < samples; i++)
        grid[i] = r0 + h * i;
    for (int i = 0; i < samples0; i++)
        grid0[i] = h0 * i;
    
    char etaname[255], etaname0[255];
    sprintf(etaname,  "dwi_a-N%d-K%g-l%d-k%g.dwf", U.n, U.k, ln, kn);
    sprintf(etaname0, "dwi_f-N%d-K%g-l%d-k%g.dwf", U.n, U.k, ln, kn);
    
    if (not array.hdfload(etaname) or not array0.hdfload(etaname0))
    {
        //
        // solve in the classically allowed region
        //
        
        // data arrays
        o2scl::ovector xg(samples);
        o2scl::omatrix y1g(samples,2), y1pg(samples,2), y1errg(samples,2);	// Re η
        o2scl::omatrix y2g(samples,2), y2pg(samples,2), y2errg(samples,2);	// Im η
        
        // derivatives callback
        o2scl::ode_funct_mfptr<IrregularWave> derivs(this, & IrregularWave::derivs);
        
        // faster adaptive Cash-Karp Runge-Kutta stepper
        o2scl::gsl_astep<decltype(derivs)> adapt_stepper;
        o2scl::gsl_rkck_fast<2,decltype(derivs)> stepper;
        adapt_stepper.set_step(stepper);
        
        // create grid
        for (int i = 0; i < samples; i++)
            xg[i] = r0 + h * (samples - i - 1);
        
        // set (right) boundary conditions from asymptotics of Riccati-Bessel function
        y1g[0][0] = -ric_n(ln,kn*xg[0]);        // Re η(R)
        y1g[0][1] = -kn*dric_n(ln,kn*xg[0]);    // Re η'(R)
        y2g[0][0] = ric_j(ln,kn*xg[0]);         // Im η(R)
        y2g[0][1] = kn*dric_j(ln,kn*xg[0]);     // Im η'(R)
        
        int n1 = solve2(xg, samples + 1, -h, y1g, y1pg, y1errg, adapt_stepper, derivs, RETURN_ON_OVERFLOW);
        int n2 = solve2(xg, samples + 1, -h, y2g, y2pg, y2errg, adapt_stepper, derivs, RETURN_ON_OVERFLOW);
        
        if (std::min(n1,n2) < samples)
            throw exception ("[IrregularWave] Unable to solve the wave function even in the classically allowed domain.");
        
        // combine real and imaginary part
        array.resize(samples);
        array_re.resize(samples);
        array_im.resize(samples);
        for (int i = 0; i < samples; i++)
        {
            array_re[i] = y1g[samples-i-1][0];
            if (not finite(array_re[i]))
                array_re[i] = 0.;
            
            array_im[i] = y2g[samples-i-1][0];
            if (not finite(array_im[i]))
                array_im[i] = 0.;
            
            array[i] = Complex(array_re[i], array_im[i]);
        }
        
        array.hdfsave(etaname);
        
        //
        // solve in the classically forbidden region (if any)
        //
        
        if (r0 > 0.)
        {
            // data arrays
            o2scl::ovector x0g(samples0);
            o2scl::omatrix y01g(samples0,2), y01pg(samples0,2), y01errg(samples0,2);	// Re η
            o2scl::omatrix y02g(samples0,2), y02pg(samples0,2), y02errg(samples0,2);	// Im η
            
            // derivatives callback
            o2scl::ode_funct_mfptr<IrregularWave> derivs0(this, & IrregularWave::derivs0);
            
            // faster adaptive Cash-Karp Runge-Kutta stepper
            o2scl::gsl_astep<decltype(derivs0)> adapt_stepper0;
            o2scl::gsl_rkck_fast<2,decltype(derivs0)> stepper0;
            adapt_stepper0.set_step(stepper0);
            
            // create grid
            for (int i = 0; i < samples0; i++)
                x0g[i] = h0 * (samples0 - i - 1);
            
            // set (right) boundary conditions from left value of the allowed computation
            y01g[0][0] = y1g[samples-1][0];    // Re η(r₀)
            y01g[0][1] = y1g[samples-1][1];    // Re η'(r₀)
            y02g[0][0] = y2g[samples-1][0];    // Im η(r₀)
            y02g[0][1] = y2g[samples-1][1];    // Im η'(r₀)
            
            int n01 = solve2(x0g, samples0, -h0, y01g, y01pg, y01errg, adapt_stepper0, derivs0, RETURN_ON_OVERFLOW);
            int n02 = solve2(x0g, samples0, -h0, y02g, y02pg, y02errg, adapt_stepper0, derivs0, RETURN_ON_OVERFLOW);
            
            if (std::min(n01,n02) < samples0)
                throw exception ("[IrregularWave] Unable to solve the wave function in the classically forbidden domain.");
            
            // combine real and imaginary part
            array0.resize(samples0);
            array0_re.resize(samples0);
            array0_im.resize(samples0);
            for (int i = 0; i < samples0; i++)
            {
                array0_re[i] = y01g[samples0-i-1][0];
                if (not finite(array0_re[i]))
                    array0_re[i] = 0.;
                
                array0_im[i] = y02g[samples0-i-1][0];
                if (not finite(array0_im[i]))
                    array0_im[i] = 0.;
                
                array0[i] = Complex(array0_re[i], array0_im[i]);
            }
            array0.hdfsave(etaname0);
        }
    }
    else
    {
        array_re = rArray(samples);
        array_im = rArray(samples);
        for (int i = 0; i < samples; i++)
        {
            array_re[i] = array[i].real();
            array_im[i] = array[i].imag();
        }
        
        array0_im = rArray(samples0);
        array0_re = rArray(samples0);
        for (int i = 0; i < samples0; i++)
        {
            array0_re[i] = array0[i].real();
            array0_im[i] = array0[i].imag();
        }
    }
    
    // setup the interpolators
    this->interpolator_re.set(grid, array_re, array_re.size());
    this->interpolator_im.set(grid, array_im, array_im.size());
    this->interpolator0_re.set(grid0, array0_re, array0_re.size());
    this->interpolator0_im.set(grid0, array0_im, array0_im.size());
}

IrregularWave IrregularWave::operator=(const IrregularWave& W)
{
    U = W.U;
    kn = W.kn;
    ln = W.ln;
    samples = W.samples;
    h = W.h;
    array_im = W.array_im;
    array_re = W.array_re;
    array0_im = W.array0_im;
    array0_re = W.array0_re;
    
    interpolator_im.set(W.grid, W.array_im, W.array_im.size());
    interpolator_re.set(W.grid, W.array_re, W.array_re.size());
    
    interpolator0_im.set(W.grid0, W.array0_im, W.array0_im.size());
    interpolator0_re.set(W.grid0, W.array0_re, W.array0_re.size());
    
    return *this;
}

int IrregularWave::derivs(double x, size_t nv, const o2scl::ovector_base& y, o2scl::ovector_base& dydx)
{
    dydx[0] = y[1];
    dydx[1] = (x == 0.) ? 0. : (-kn*kn + 2.*U(x) + ln*(ln+1)/(x*x)) * y[0];
    return 0;
}

int IrregularWave::derivs0(double x, size_t nv, const o2scl::ovector_base& y, o2scl::ovector_base& dydx)
{
    dydx[0] = y[1];
    dydx[1] = (x == 0.) ? 0. : (-kn*kn + 2.*U(x)) * y[0] + 2*ln/x * y[1];
    return 0;
}

Complex IrregularWave::operator()(double x) const
{
    std::cout << "Evaluationm, x = " << x << "\n";
    
    Evaluations++;
    
    if (x <= 0.)
        return 0.;
    
    // extrapolate using asymptotic form if far away
    else if (x > rf)
        return Complex(-ric_n(ln,kn*x),ric_j(ln,kn*x));
    
    // interpolate using cspline interpolator in the allowed region
    else if (x >= r0)
    {
        std::cout << "x = " << x << "\n";
        return Complex(interpolator_re.interp(x), interpolator_im.interp(x));
    }
    
    // interpolate using cspline interpolator in the forbidden region
    else
    {
        std::cout << "x0 = " << x << "\n";
        return Complex(interpolator0_re.interp(x), interpolator0_im.interp(x)) / pow(x, ln);
    }
}

void IrregularWave::toFile(const char* filename) const
{
    int N = array_re.size();
    cArray array(N);
    for (int i = 0; i < N; i++)
        array[i] = Complex(array_re[i], array_im[i]);
    
    write_array(grid, array, filename);
}
