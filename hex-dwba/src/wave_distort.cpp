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

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <gsl/gsl_sf.h>

#include "arrays.h"
#include "ode.h"
#include "potential.h"
#include "wave_distort.h"

// -------------------------------------------------------------------------- //
//                                                                            //
//  DistortedWave class members                                               //
//                                                                            //
// -------------------------------------------------------------------------- //

DistortedWave DistortedWave::operator= (DistortedWave const& W)
{
    gsl_interp_free(interpolator);
    if (ln > 0)
        gsl_interp_free(interpolator0);
    
    U = W.U;
    kn = W.kn;
    ln = W.ln;
    r0 = W.r0;
    rf = W.rf;
    
    phase = W.phase;
    
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

int derivs (double x, const double y[2], double dydx[2], void * data)
{
    // retype parameter
    std::tuple<double,int,DistortingPotential const *> *pdata =
        (std::tuple<double,int,DistortingPotential const *> *)data;
    
    // get individual sub-parameters
    double kn = std::get<0>(*pdata);
    int ln = std::get<1>(*pdata);
    DistortingPotential const & U = *std::get<2>(*pdata);
    
    // compute derivatives
    dydx[0] = y[1];
    dydx[1] = (x == 0.) ? 0. : (-kn*kn + 2.*U(x) + ln*(ln+1)/(x*x)) * y[0];
    
    return GSL_SUCCESS;
}

int derivs0 (double x, const double y[2], double dydx[2], void * data)
{
    // retype parameter
    std::tuple<double,int,DistortingPotential const *> *pdata =
        (std::tuple<double,int,DistortingPotential const *> *)data;
    
    // get individual sub-parameters
    double kn = std::get<0>(*pdata);
    int ln = std::get<1>(*pdata);
    DistortingPotential const & U = *std::get<2>(*pdata);
    
    // compute derivatives
    dydx[0] = y[1];
    dydx[1] = (x == 0.) ? 0. : (-kn*kn + 2.*U(x)) * y[0] - 2*(ln+1)/x * y[1];
    
    return GSL_SUCCESS;
}

DistortedWave::DistortedWave(double _kn, int _ln, DistortingPotential const & _U)
    : Evaluations(0), U(_U), kn(_kn), ln(_ln), r0(sqrt(ln*(ln+1)) / kn), rf(std::max(2*r0,U.getFarRadius()))
{
    // determine discretization, use at least 1000 samples
    int N = 100;                // N samples per wave length
    h = std::min (              // grid step
        2*M_PI/(N*kn),          //  -> N samples per wave length and
        (rf-r0)/1000            //  -> at least 1000 samples totally
    );
    samples = (rf-r0)/h + 1;    // with both boundaries
    
    // discretization of the classically forbidden region
    samples0 = 1000;
    h0 = r0 / (samples0 - 1);
    
    if (ln == 0)
        samples0 = 0;
    
    // create grids
    grid.resize(samples);
    grid0.resize(samples0);
    for (int i = 0; i < samples; i++)
        grid[i] = r0 + h * i;
    for (int i = 0; i < samples0; i++)
        grid0[i] = h0 * i;
    
    // look for relevant data
    char filename[255], filename0[255];
    sprintf(filename, "dwr_a-N%d-K%g-l%d-k%g.dwf", U.n(), U.k(), ln, kn);
    sprintf(filename0, "dwr_f-N%d-K%g-l%d-k%g.dwf", U.n(), U.k(), ln, kn);
    if (array.hdfload(filename) and (ln == 0 or array0.hdfload(filename0)))
    {
        // the last element is the phase shift
        phase = array.pop_back();
    }
    else
    {
        // assemble data needed to compute the derivatives
        std::tuple<double, int, DistortingPotential const *> data = std::make_tuple (kn, ln, &U);
        
        //
        // solve in the classically forbidden region (if any)
        //
        
        // data arrays
        double x0g[samples0], y0g[samples0][2];
        
        if (ln > 0)
        {
            // create grid
            for (int i = 0; i < samples0; i++)
                x0g[i] = h0 * i;
            
            // set left boundary conditions (normalized regular wave divided by r^(l+1))
            y0g[0][0] = 1;
            y0g[0][1] = -1./(ln+1);
            
            // solve the equation
            solve2 (x0g, samples0 + 1, h0, y0g, derivs0, &data, NORMALIZE_ON_OVERFLOW);
        }
        
        //
        // solve in the classically allowed region
        //
        
        // data arrays
        double xg[samples], yg[samples][2];
        
        // create grid
        for (int i = 0; i < samples; i++)
            xg[i] = r0 + h * i;
        
        // set left boundary conditions
        if (ln > 0)
        {
            // match to the forbidden region
            yg[0][0] = y0g[samples0-1][0];
            yg[0][1] = (ln+1)/r0 * y0g[samples0-1][0] + y0g[samples0-1][1];
        }
        else
        {
            // start general regular function
            yg[0][0] = 0;
            yg[0][1] = 1;
        }
        
        // solve
        solve2 (xg, samples + 1, h, yg, derivs, &data, NORMALIZE_ON_OVERFLOW);
        
        // values y0g[last] and yg[first] shoud be equal -- normalize the classically forbidden part (if any)
        if (ln > 0)
        {
            double scale_factor = yg[0][0] / y0g[samples0-1][0];
            
            for (int i = 0; i < samples0; i++)
                y0g[i][0] *= scale_factor;
        }
        
        // compute phase shift; use last point
        double R   = xg[samples-2];
        double val =  yg[samples-2][0];
        double der = yg[samples-2][1];
        double D = der / val;
        double pilhalf = M_PI * ln / 2.;
        double tg_delta = (kn * cos(kn*R - pilhalf) - D * sin(kn*R - pilhalf)) /
                          (kn * sin(kn*R - pilhalf) + D * cos(kn*R - pilhalf));
        phase = atan(tg_delta);
        
        // normalize ( far away should be η ≈ exp(iδ)sin(kr-πl/2+δ) )
        array = rArray(samples);
        double inverse_norm = 1./sqrt(val*val + der*der/(kn*kn));
        for (int ir = 0; ir < samples; ir++)
            array[ir] = yg[ir][0] * inverse_norm;
        
        // normalize "forbidden" segment as well
        array0 = rArray(samples0);
        for (int ir = 0; ir < samples0; ir++)
            array0[ir] = y0g[ir][0] * inverse_norm;
        
        // save for recyclation, append last element on tail
        if (ln > 0)
            array0.hdfsave(filename0);
        array.push_back(phase);
        array.hdfsave(filename);
        array.pop_back();
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

DistortedWave::~DistortedWave()
{
    gsl_interp_free(interpolator);
    if (ln > 0)
        gsl_interp_free(interpolator0);
}

double DistortedWave::operator() (double x) const
{
    Evaluations++;
    
    if (x <= 0.)
        return 0.;
    
    // extrapolate using asymptotic form if far away
    else if (x > grid.back())
        return sin(kn*x - ln*M_PI*0.5 + phase);
    
    // or interpolate using cspline interpolator in the classically allowed region
    else if (x >= grid.front())
        return gsl_interp_eval(interpolator, grid.data(), array.data(), x, nullptr);
    
    // or interpolate using cspline interpolator in the classically forbidden region
    else
        return gsl_interp_eval(interpolator0, grid0.data(), array0.data(), x, nullptr) * pow(x/r0,ln+1);
}

std::pair<double,int> DistortedWave::getZeroAsymptotic (double x) const
{
    double y;
    
    if (x <= 0.)
        y = 0.;
    
    // extrapolate using asymptotic form if far away
    else if (x > grid.back())
        y = sin(kn*x - ln*M_PI*0.5 + phase) / pow(x, ln + 1);
    
    // or interpolate using cspline interpolator in the classically allowed region
    else if (x >= grid.front())
        y = gsl_interp_eval(interpolator, grid.data(), array.data(), x, nullptr) / pow(x, ln + 1);
    
    // or interpolate using cspline interpolator in the classically forbidden region
    else
        y = gsl_interp_eval(interpolator0, grid0.data(), array0.data(), x, nullptr) / pow(r0,ln+1);
    
    return std::make_pair(y, ln + 1);
}
