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
#include <gsl/gsl_interp.h>

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
    h = std::min (              // grid step
        2*M_PI/(N*kn),          //  -> N samples per wave length and
        (rf-r0)/1000            //  -> at least 1000 samples total
    );
    samples = (rf - r0)/h + 1;  // with both boundaries
    
    // determine discretization in the classically forbidden region
    samples0 = 100;
    h0 = r0 / (samples0 - 1);
    
    if (ln == 0)
        samples0 = 0;
    
    // auxiliary arrays for HDF I/O
    cArray array, array0;
    
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
    
    if (not array.hdfload(etaname) or (ln > 0 and not array0.hdfload(etaname0)))
    {
        // data arrays
        o2scl::ovector x0g(samples0);
        o2scl::omatrix y01g(samples0,2), y01pg(samples0,2), y01errg(samples0,2);    // Re η
        o2scl::omatrix y02g(samples0,2), y02pg(samples0,2), y02errg(samples0,2);    // Im η
        
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
        
        // integrate to the left bound only if is is not zero
        int n = (r0 == 0.) ? 0 : 1;
        
        // integrate the differential equations (real and imaginary)
        solve2(xg, samples + n, -h, y1g, y1pg, y1errg, adapt_stepper, derivs, ABORT_ON_OVERFLOW);
        solve2(xg, samples + n, -h, y2g, y2pg, y2errg, adapt_stepper, derivs, ABORT_ON_OVERFLOW);
        
        // fix left boundary (extrapolate to zero, if necessary)
        if (r0 == 0.)
        {
            y1g[samples-1][0] = y1g[samples-2][0];
            y2g[samples-1][0] = y2g[samples-2][0];
        }
        
        //
        // solve in the classically forbidden region (if any)
        //
        
        if (ln > 0)
        {
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
            y01g[0][0] = pow(r0,ln)*y1g[samples-1][0];
            y01g[0][1] = ln*pow(r0,ln-1)*y1g[samples-1][0] + pow(r0,ln)*y1g[samples-1][1];
            y02g[0][0] = pow(r0,ln)*y2g[samples-1][0];
            y02g[0][1] = ln*pow(r0,ln-1)*y2g[samples-1][0] + pow(r0,ln)*y2g[samples-1][1];
            
            // integrate the differential equations (real and imaginary)
            solve2(x0g, samples0 + 1, -h0, y01g, y01pg, y01errg, adapt_stepper0, derivs0, ABORT_ON_OVERFLOW);
            solve2(x0g, samples0 + 1, -h0, y02g, y02pg, y02errg, adapt_stepper0, derivs0, ABORT_ON_OVERFLOW);
        }
        
        // combine real and imaginary part
        array.resize(samples);
        array_re.resize(samples);
        array_im.resize(samples);
        for (int i = 0; i < samples; i++)
        {
            array_re[i] = y1g[samples-i-1][0];
            array_im[i] = y2g[samples-i-1][0];
            array[i] = Complex(array_re[i], array_im[i]);
        }
        
        // combine real and imaginary part (forbidden segment)
        array0.resize(samples0);
        array0_re.resize(samples0);
        array0_im.resize(samples0);
        for (int i = 0; i < samples0; i++)
        {
            array0_re[i] = y01g[samples0-i-1][0];
            array0_im[i] = y02g[samples0-i-1][0];
            array0[i] = Complex(array0_re[i], array0_im[i]);
        }
        
        // check origin
        if (not finite(array[0].real()) or not finite(array[0].imag()))
            array[0] = array[1];
        
        // check origin
        if (ln > 0 and (not finite(array[0].real()) or not finite(array[0].imag())))
            array0[0] = array0[1];
        
        // write results to file
        if (ln > 0)
            array0.hdfsave(etaname0);
        array.hdfsave(etaname);
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
    interpolator_re = gsl_interp_alloc(gsl_interp_cspline, array_re.size());
    interpolator_im = gsl_interp_alloc(gsl_interp_cspline, array_im.size());
    gsl_interp_init(interpolator_re, grid.data(), array_re.data(), array_re.size());
    gsl_interp_init(interpolator_im, grid.data(), array_im.data(), array_im.size());
    if (ln > 0)
    {
        interpolator0_re = gsl_interp_alloc(gsl_interp_cspline, array0_re.size());
        interpolator0_im = gsl_interp_alloc(gsl_interp_cspline, array0_im.size());
        gsl_interp_init(interpolator0_re, grid0.data(), array0_re.data(), array0_re.size());
        gsl_interp_init(interpolator0_im, grid0.data(), array0_im.data(), array0_im.size());
    }
}

IrregularWave IrregularWave::operator=(const IrregularWave& W)
{
    gsl_interp_free(interpolator_re);
    gsl_interp_free(interpolator_im);
    if (ln > 0)
    {
        gsl_interp_free(interpolator0_re);
        gsl_interp_free(interpolator0_im);
    }
    
    U = W.U;
    kn = W.kn;
    ln = W.ln;
    r0 = W.r0;
    rf = W.rf;
    
    grid = W.grid;
    array_re = W.array_re;
    array_im = W.array_im;
    h = W.h;
    samples = W.samples;
    
    grid0 = W.grid0;
    array0_re = W.array0_re;
    array0_im = W.array0_im;
    h0 = W.h0;
    samples0 = W.samples0;
    
    interpolator_re = gsl_interp_alloc(gsl_interp_cspline, array_re.size());
    gsl_interp_init(interpolator_re, grid.data(), array_re.data(), array_re.size());
    
    interpolator_im = gsl_interp_alloc(gsl_interp_cspline, array_im.size());
    gsl_interp_init(interpolator_im, grid.data(), array_im.data(), array_im.size());
    
    if (ln > 0)
    {
        interpolator0_re = gsl_interp_alloc(gsl_interp_cspline, array0_re.size());
        gsl_interp_init(interpolator0_re, grid0.data(), array0_re.data(), array0_re.size());
        
        interpolator0_im = gsl_interp_alloc(gsl_interp_cspline, array0_im.size());
        gsl_interp_init(interpolator0_im, grid0.data(), array0_im.data(), array0_im.size());
    }
    
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
    Evaluations++;
    
    if (x < 0.)
        return 0.;
    
    // extrapolate using asymptotic form if far away
    else if (x > grid.back())
        return Complex(-ric_n(ln,kn*x),ric_j(ln,kn*x));
    
    // interpolate using cspline interpolator in the allowed region
    else if (x >= grid.front())
    {
        return Complex (
            gsl_interp_eval(interpolator_re, grid.data(), array_re.data(), x, nullptr),
            gsl_interp_eval(interpolator_im, grid.data(), array_im.data(), x, nullptr)
        );
    }
    
    // interpolate using cspline interpolator in the forbidden region
    else
    {
        return Complex (
            gsl_interp_eval(interpolator0_re, grid0.data(), array0_re.data(), x, nullptr),
            gsl_interp_eval(interpolator0_im, grid0.data(), array0_im.data(), x, nullptr)
        ) / pow(x, ln);
    }
}

std::pair<Complex,int> IrregularWave::getZeroAsymptotic(double x) const
{
    Complex y;
    
    if (x < 0.)
        y = 0.;
    
    // extrapolate using asymptotic form if far away
    else if (x > grid.back())
        y = Complex(-ric_n(ln,kn*x),ric_j(ln,kn*x)) * pow(x, ln);
    
    // interpolate using cspline interpolator in the allowed region
    else if (x >= grid.front())
    {
        y = Complex (
            gsl_interp_eval(interpolator_re, grid.data(), array_re.data(), x, nullptr),
            gsl_interp_eval(interpolator_im, grid.data(), array_im.data(), x, nullptr)
        ) * pow(x, ln);
    }
    
    // interpolate using cspline interpolator in the forbidden region
    else
    {
        y = Complex (
            gsl_interp_eval(interpolator0_re, grid0.data(), array0_re.data(), x, nullptr),
            gsl_interp_eval(interpolator0_im, grid0.data(), array0_im.data(), x, nullptr)
        );
    }
    
    return std::make_pair (y, -ln);
}

void IrregularWave::toFile(const char* filename) const
{
    int N = array_re.size();
    cArray array(N);
    for (int i = 0; i < N; i++)
        array[i] = Complex(array_re[i], array_im[i]);
    
    write_array(grid, array, filename);
}

IrregularWave::~IrregularWave()
{
    gsl_interp_free(interpolator_re);
    gsl_interp_free(interpolator_im);
    if (ln > 0)
    {
        gsl_interp_free(interpolator0_re);
        gsl_interp_free(interpolator0_im);
    }
}
