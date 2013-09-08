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
 * \* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <o2scl/ode_funct.h>
#include <o2scl/ode_iv_solve.h>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_interp.h>

#include "arrays.h"
#include "complex.h"
#include "wave_forbid.h"

// -------------------------------------------------------------------------- //
//                                                                            //
//  ForbiddenWave class members                                               //
//                                                                            //
// -------------------------------------------------------------------------- //

ForbiddenWave::ForbiddenWave(double _kn, int _ln, DistortingPotential const & _U)
: Evaluations(0), U(_U), kn(_kn), ln(_ln), Scaled(false), r0(sqrt(ln*(ln+1))/kn), rf(U.getFarRadius()/kn)
{
    // determine discretization
    samples = 10001;		// with both boundaries
    h = (rf - r0) / (samples - 1);	// grid step
    
    samples0 = 100;
    h0 = r0 / (samples0 - 1);
    
    if (ln == 0)
        samples0 = 0;
    
    // create grids
    grid.resize(samples);
    for (int i = 0; i < samples; i++)
        grid[i] = r0 + i * h;
    grid0.resize(samples0);
    for (int i = 0; i < samples0; i++)
        grid0[i] = i * h0;
    
    // look for relevant data
    char filename[255], filename0[255];
    sprintf(filename, "dwf_a-N%d-K%g-l%d-k%g.dwf", U.n, U.k, ln, kn);
    sprintf(filename0, "dwf_f-N%d-K%g-l%d-k%g.dwf", U.n, U.k, ln, kn);
    if (not array.hdfload(filename) or (ln > 0 and not array0.hdfload(filename0)))
    {
        // data arrays
        o2scl::ovector x0g(samples0);
        o2scl::omatrix y0g(samples0,2), y0pg(samples0,2), y0errg(samples0,2);
        
        //
        // solve in the classically allowed region
        //
        
        // prepare derivative callback MFPTR-wrapper
        o2scl::ode_funct_mfptr<ForbiddenWave> derivs(this, &ForbiddenWave::derivs);
        
        // prepare faster adaptive Cash-Karp Runge-Kutta stepper
        o2scl::gsl_astep<decltype(derivs)> adapt_stepper;
        o2scl::gsl_rkck_fast<2,decltype(derivs)> stepper;
        adapt_stepper.set_step(stepper);
        
        // data arrays
        o2scl::ovector xg(samples);
        o2scl::omatrix yg(samples,2), ypg(samples,2), yerrg(samples,2);
        
        // create grid
        for (int i = 0; i < samples; i++)
            xg[i] = r0 + h * (samples - i - 1);
        
        // set FAR RIGHT boundary condition
        yg[0][0] =  ric_k_scaled(ln,kn*rf);
        yg[0][1] = dric_k_scaled(ln,kn*rf) * kn;
        
        // solve
        solve2(xg, samples, -h, yg, ypg, yerrg, adapt_stepper, derivs, ABORT_ON_OVERFLOW);
        
        //
        // solve in the classically forbidden region (if any)
        //
        
        if (ln > 0)
        {
            // prepare derivative callback MFPTR-wrapper
            o2scl::ode_funct_mfptr<ForbiddenWave> derivs0(this, &ForbiddenWave::derivs0);
            
            // prepare faster adaptive Cash-Karp Runge-Kutta stepper
            o2scl::gsl_astep<decltype(derivs0)> adapt_stepper0;
            o2scl::gsl_rkck_fast<2,decltype(derivs0)> stepper0;
            adapt_stepper0.set_step(stepper0);
            
            // create grid
            for (int i = 0; i < samples0; i++)
                x0g[i] = h0 * (samples0 - i - 1);
            
            // set RIGHT boundary condition
            y0g[0][0] = pow(r0,ln) * yg[samples-1][0];
            y0g[0][1] = ln*pow(r0,ln-1)*yg[samples-1][0] + pow(r0,ln)*yg[samples-1][1];
            
            // solve
            solve2(x0g, samples0, -h0, y0g, y0pg, y0errg, adapt_stepper0, derivs0, ABORT_ON_OVERFLOW);
        }
        
        // copy array to class storage
        array = rArray(samples);
        for (int i = 0; i < samples; i++)
            array[i] = yg[samples-i-1][0];
        
        // copy also the forbidden section
        array0 = rArray(samples0);
        for (int i = 0; i < samples0; i++)
            array0[i] = y0g[samples0-i-1][0];
        
        // save for recyclation
        array.hdfsave(filename);
        if (ln > 0)
            array0.hdfsave(filename0);
    }
    
    // setup the interpolators
    interpolator = gsl_interp_alloc(gsl_interp_cspline, array.size());
    gsl_interp_init(interpolator, grid.data(), array.data(), array.size());
    if (ln > 0)
    {
        interpolator0 = gsl_interp_alloc(gsl_interp_cspline, array0.size());
        gsl_interp_init(interpolator0, grid0.data(), array0.data(), array0.size());
    }
}

ForbiddenWave ForbiddenWave::operator=(ForbiddenWave const & W)
{
    gsl_interp_free(interpolator);
    if (ln > 0)
        gsl_interp_free(interpolator0);
    
    kn = W.kn;
    ln = W.ln;
    U = W.U;
    r0 = W.r0;
    rf = W.rf;
    
    grid = W.grid;
    array = W.array;
    samples = W.samples;
    h = W.h;
    
    grid0 = W.grid0;
    array0 = W.array0;
    samples0 = W.samples0;
    h0 = W.h0;
    
    interpolator = gsl_interp_alloc(gsl_interp_cspline, array.size());
    gsl_interp_init(interpolator, grid.data(), array.data(), array.size());
    if (ln > 0)
    {
        interpolator0 = gsl_interp_alloc(gsl_interp_cspline, array0.size());
        gsl_interp_init(interpolator0, grid0.data(), array0.data(), array0.size());
    }
    
    return *this;
}

int ForbiddenWave::derivs(double x, size_t nv, const o2scl::ovector_base& y, o2scl::ovector_base& dydx)
{
    dydx[0] = y[1];
    dydx[1] = (x == 0.) ? 0. : (2.*U(x) + ln*(ln+1)/(x*x)) * y[0] + 2*kn*y[1];
    
    return 0;
}

int ForbiddenWave::derivs0(double x, size_t nv, const o2scl::ovector_base& y, o2scl::ovector_base& dydx)
{
    dydx[0] = y[1];
    dydx[1] = (x == 0.) ? 0. : (2.*U(x) - 2*kn*ln/x) * y[0] + (2*kn + 2*ln/x)*y[1];
    
    return 0;
}

double ForbiddenWave::operator()(double x) const
{
    if (x <= 0.)
        return 0.;
    
    // extrapolate
    else if (x > grid.back())
        return Scaled ? ric_k_scaled(ln,kn*x) : ric_k(ln,kn*x);
    
    // interpolate allowed region
    else if (x >= grid.front())
    {
        double y = gsl_interp_eval(interpolator, grid.data(), array.data(), x, nullptr);
        return Scaled ? y : y * exp(-kn*x);
    }
    
    // interpolate forbidden region
    else
    {
        double y = gsl_interp_eval(interpolator0, grid0.data(), array0.data(), x, nullptr) * pow(x,-ln);
        return Scaled ? y : y * exp(-kn*x);
    }
}

void ForbiddenWave::toFile(const char* filename) const
{
    write_array(grid, array, filename);
}

void ForbiddenWave::scale(bool s)
{
    Scaled = s;
}

ForbiddenWave::~ForbiddenWave()
{
    gsl_interp_free(interpolator);
    if (ln > 0)
        gsl_interp_free(interpolator0);
}

std::pair<double,int> ForbiddenWave::getZeroAsymptotic(double x) const
{
    double y;
    
    if (x <= 0.)
        y = 0.;
    
    // extrapolate
    else if (x > grid.back())
        y = Scaled ? ric_k_scaled(ln,kn*x) * pow(x,ln) : ric_k(ln,kn*x) * pow(x,ln);
    
    // interpolate allowed region
    else if (x >= grid.front())
    {
        y = Scaled ?
            gsl_interp_eval(interpolator, grid.data(), array.data(), x, nullptr) * pow(x,ln) :
            gsl_interp_eval(interpolator, grid.data(), array.data(), x, nullptr) * pow(x,ln) * exp(-kn*x);
    }
    
    // interpolate forbidden region
    else
    {
        y = Scaled ?
            gsl_interp_eval(interpolator0, grid0.data(), array0.data(), x, nullptr) :
            gsl_interp_eval(interpolator0, grid0.data(), array0.data(), x, nullptr) * exp(-kn*x);
    }
    
    return std::make_pair(y, -ln);
}
