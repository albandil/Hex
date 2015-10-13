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

#ifndef HEX_ODE
#define HEX_ODE

#include <gsl/gsl_odeiv2.h>

#define DIVERGENCE_THRESHOLD    100

#define ABORT_ON_OVERFLOW       0
#define RETURN_ON_OVERFLOW      1
#define NORMALIZE_ON_OVERFLOW   2

#include "hex-misc.h"

/**
 * @brief Second-order differential equation solver.
 * 
 * This routine is adapted from the O₂scl file "ode_iv_solve.h".
 * The modification it the following:
 * - When the solution reaches some dangerously high absolute value, the function
 *   either exits (returning the last solved index that was kept inside the
 *   limits) or renormalizes the up-to-now solution so that the overflow is
 *   avoided.
 * @param xg Independent variable grid.
 * @param N Size of the grid.
 * @param h Grid spacing.
 * @param yg (out) Solution.
 * @param derivs Second derivative callback of the signature
 * @code
 *   int derivs(double x, const double y[2], double dydx[2], void* params)
 * @endcode
 * @param data Custom data pointer to pass to the "derivs" routine.
 * @param flag Behaviour on overflow, one of { @ref ABORT_ON_OVERFLOW | @ref RETURN_ON_OVERFLOW | @ref NORMALIZE_ON_OVERFLOW }.
 * 
 * @return N for successful run of n < N for forced terminantion on overflow,
 *         where n is the index of last valid field in yg and ypg.
 */
int solve2
(
    double * xg, int N, double h, double ** yg,
    int (*derivs) (double, const double*, double*, void*),
    void * data, int flag
)
{
    // initialize auxiliary variables
    size_t ntrial = 100000000;	// 10⁷
    double x0 = xg[0], x1 = xg[N - 2];
    double x = x0, xnext;
    int status = GSL_SUCCESS, first_status = GSL_SUCCESS;
    size_t nsteps = 0;
    xg[0] = x0;
    
    // setup the stepper
    const gsl_odeiv2_step_type * step_type = gsl_odeiv2_step_rkck;
    gsl_odeiv2_step * step = gsl_odeiv2_step_alloc(step_type, 2);
    gsl_odeiv2_control * control = gsl_odeiv2_control_y_new(1e-6, 0.0);
    gsl_odeiv2_evolve * evolve = gsl_odeiv2_evolve_alloc(2);
    
    // setup the system
    gsl_odeiv2_system sys;
    sys.function = derivs;
    sys.jacobian = nullptr;
    sys.dimension = 2;
    sys.params = data;
    
    // for all steps
    for(int i = 1; i < N - 1 and status == GSL_SUCCESS; i++)
    {
        xnext = x0 + (x1 - x0) * ((double)i)/((double)(N - 2));
        double * y_row = yg[i];
        y_row[0] = yg[i-1][0];
        y_row[1] = yg[i-1][1];
        
        // Step until we reach the next grid point
        bool done = false;
        while (not done and status == GSL_SUCCESS)
        {
            status = gsl_odeiv2_evolve_apply
            (
                evolve,
                control,
                step,
                &sys,
                &x, xnext,
                &h,
                y_row
            );
            
            // check finiteness
            if (not std::isfinite(y_row[0]))
                HexException("[solve2] Infinite result (%g) for i = %d", y_row[0], i);
            
            // advance number of steps
            nsteps++;
            if (status != GSL_SUCCESS and first_status != GSL_SUCCESS)
                first_status = status;
            if (nsteps > ntrial)
                std::cerr << "[solve2] Too many steps required (ntrial=" << ntrial << ", x=" << x << ").";
            
            // avoid overflow (normalize)
            if (fabs(y_row[0]) > DIVERGENCE_THRESHOLD)
            {
                if (flag == NORMALIZE_ON_OVERFLOW)
                {
                    double inverse_norm = 1./fabs(y_row[0]);
                    for (int j = 0; j <= i; j++)
                    {
                        yg[j][0] *= inverse_norm;
                        yg[j][1] *= inverse_norm;
                    }
                }
                if (flag == RETURN_ON_OVERFLOW)
                {
                    return i - 1;
                }
                if (flag == ABORT_ON_OVERFLOW)
                {
                    HexException("[solve2] Overflow error.");
                }
            }
            
            // decide if we have reached the grid point
            if (x1 > x0)
            {
                if (x >= xnext)
                    done = true;
            }
            else
            {
                if (x <= xnext)
                    done = true;
            }
            
        }
    }
    
    gsl_odeiv2_evolve_free(evolve);
    gsl_odeiv2_control_free(control);
    gsl_odeiv2_step_free(step);
    
    return N;
}

#endif
