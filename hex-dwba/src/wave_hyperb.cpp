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
#include "wave_hyperb.h"

// -------------------------------------------------------------------------- //
//                                                                            //
//  HyperbolicWave class members                                              //
//                                                                            //
// -------------------------------------------------------------------------- //

HyperbolicWave::HyperbolicWave(double _kn, int _ln, DistortingPotential const & _U)
  : Evaluations(0), U(_U), kn(_kn), ln(_ln), Scaled(false),
    r0(sqrt(ln*(ln+1))/kn), rf(U.getFarRadius()/kn)
{
    // determine discretization
    samples = 10001;        // with both boundaries; FIXME some heuristic ?
    h = (rf - r0) / (samples - 1); // average grid step
    
    samples0 = 100;
    h0 = r0 / (samples0 - 1);;
    
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
    char etaname[255], etaname0[255];
    sprintf(etaname, "dwh_a-N%d-K%g-l%d-k%g.dwf", U.n, U.k, ln, kn);
    sprintf(etaname0, "dwh_f-N%d-K%g-l%d-k%g.dwf", U.n, U.k, ln, kn);
    if (not array.hdfload(etaname) or (ln > 0 and not array0.hdfload(etaname0)))
    {
        //
        // solve in the classically forbidden region (if any)
        //
        
        // data arrays
        o2scl::ovector x0g(samples0);
        o2scl::omatrix y0g(samples0,2), y0pg(samples0,2), y0errg(samples0,2);
        
        if (ln > 0)
        {
            // derivatives callback
            o2scl::ode_funct_mfptr<HyperbolicWave> derivs0(this, & HyperbolicWave::derivs0);
            
            // faster adaptive Cash-Karp Runge-Kutta stepper
            o2scl::gsl_astep<decltype(derivs0)> adapt_stepper0;
            o2scl::gsl_rkck_fast<2,decltype(derivs0)> stepper0;
            adapt_stepper0.set_step(stepper0);
            
            // create grid
            for (int i = 0; i < samples0; i++)
                x0g[i] = h0 * i;
            
            // set boundary conditions
            y0g[0][0] = 0;
            y0g[0][1] = 1;
            
            // solve the equation
            solve2(x0g, samples0 + 1, h0, y0g, y0pg, y0errg, adapt_stepper0, derivs0, NORMALIZE_ON_OVERFLOW);
        }
        
        //
        // solve in the classically allowed region
        //
        
        // data arrays
        o2scl::ovector xg(samples);
        o2scl::omatrix yg(samples,2), ypg(samples,2), yerrg(samples,2);
        
        // derivatives callback
        o2scl::ode_funct_mfptr<HyperbolicWave> derivs(this, & HyperbolicWave::derivs);
        
        // faster adaptive Cash-Karp Runge-Kutta stepper
        o2scl::gsl_astep<decltype(derivs)> adapt_stepper;
        o2scl::gsl_rkck_fast<2,decltype(derivs)> stepper;
        adapt_stepper.set_step(stepper);
        
        // create grid
        for (int i = 0; i < samples; i++)
            xg[i] = r0 + h * i;
        
        // set left boundary conditions
        if (ln > 0)
        {
            // match to the forbidden region
            yg[0][0] = pow(r0,ln+1) * y0g[samples0-1][0]; // Θ (h)
            yg[0][1] = (ln+1) * pow(r0,ln) * y0g[samples0-1][0] + pow(r0,ln+1) * y0g[samples0-1][1]; // Θ'(h)
        }
        else
        {
            // start regular soltion
            yg[0][0] = 0;
            yg[0][1] = 1;
        }
        
        // solve the equation (do not re-normalize!)
        solve2(xg, samples + 1, h, yg, ypg, yerrg, adapt_stepper, derivs, ABORT_ON_OVERFLOW);
        
        // save with correct normalization
        array = rArray(samples);
        double norm = ric_i_scaled(ln,kn*rf) / yg[samples-1][0];
        for (int i = 0; i < samples; i++)
            array[i] = yg[i][0] * norm;
        
        // save also the forbidden part
        array0 = rArray(samples0);
        for (int i = 0; i < samples0; i++)
            array0[i] = y0g[i][0] * norm;
        
        array.hdfsave(etaname);
        if (ln > 0)
            array0.hdfsave(etaname0);
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

HyperbolicWave HyperbolicWave::operator=(HyperbolicWave const & W)
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
    h = W.h;
    samples = W.samples;
    
    grid0 = W.grid0;
    array0 = W.array0;
    h0 = W.h0;
    samples0 = W.samples0;
    
    interpolator = gsl_interp_alloc(gsl_interp_cspline, array.size());
    gsl_interp_init(interpolator, grid.data(), array.data(), array.size());
    
    if (ln > 0)
    {
        interpolator0 = gsl_interp_alloc(gsl_interp_cspline, array0.size());
        gsl_interp_init(interpolator0, grid0.data(), array0.data(), array0.size());
    }
    
    return *this;
}

int HyperbolicWave::derivs(double x, size_t nv, const o2scl::ovector_base& y, o2scl::ovector_base& dydx)
{
    dydx[0] = y[1];
    dydx[1] = (x == 0.) ? 0. : (2.*U(x) + ln*(ln+1)/(x*x)) * y[0] - 2*kn*y[1];
    return 0;
}

int HyperbolicWave::derivs0(double x, size_t nv, const o2scl::ovector_base& y, o2scl::ovector_base& dydx)
{
    dydx[0] = y[1];
    dydx[1] = (x == 0.) ? 0. : (2.*U(x) - 2*kn*(ln+1)/x) * y[0] - (2*kn + 2*(ln+1)/x)*y[1];
    return 0;
}

double HyperbolicWave::operator()(double x) const
{
    if (x <= 0)
        return 0;
    
    // extrapolate
    else if (x > grid.back())
        return Scaled ? ric_i_scaled(ln,kn*x) : ric_i(ln,kn*x);
    
    // interpolate allowed region
    else if (x >= grid.front())
    {
        double y = gsl_interp_eval(interpolator, grid.data(), array.data(), x, nullptr);
        return Scaled ? y : y * exp(kn*x);
    }
    
    // interpolate forbidden region
    else
    {
        double y = gsl_interp_eval(interpolator0, grid0.data(), array0.data(), x, nullptr) * pow(x,ln+1);
        return Scaled ? y : y * exp(kn*x);
    }
}

void HyperbolicWave::toFile(const char* filename) const
{
    write_array(grid, array, filename);
}

void HyperbolicWave::scale(bool s)
{
    Scaled = s;
}

HyperbolicWave::~HyperbolicWave()
{
    gsl_interp_free(interpolator);
    if (ln > 0)
        gsl_interp_free(interpolator0);
}

std::pair<double,int> HyperbolicWave::getZeroAsymptotic(double x) const
{
    double y;
    
    if (x <= 0.)
        y = 0.;
    
    // extrapolate
    else if (x > grid.back())
        y = Scaled ? ric_i_scaled(ln,kn*x) / pow(x,ln+1) : ric_i(ln,kn*x) / pow(x,ln+1);
    
    // interpolate allowed region
    else if (x >= grid.front())
    {
        y = Scaled ?
            gsl_interp_eval(interpolator, grid.data(), array.data(), array.size(), nullptr) / pow(x,ln+1) :
            gsl_interp_eval(interpolator, grid.data(), array.data(), array.size(), nullptr) / pow(x,ln+1) * exp(kn*x);
    }
    
    // interpolate forbidden region
    else
    {
        y = Scaled ?
            gsl_interp_eval(interpolator0, grid0.data(), array0.data(), array0.size(), nullptr) :
            gsl_interp_eval(interpolator0, grid0.data(), array0.data(), array0.size(), nullptr) * exp(kn*x);
    }
    
    return std::make_pair(y, ln + 1);
}
