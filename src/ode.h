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

#ifndef HEX_ODE
#define HEX_ODE

#include <o2scl/ovector_tlate.h>
#include <o2scl/omatrix_tlate.h>
#include <o2scl/ode_funct.h>
#include <o2scl/ode_iv_solve.h>
#include <o2scl/gsl_rkf45.h>
#include <o2scl/gsl_rk8pd.h>

#define DIVERGENCE_THRESHOLD	100

#define RETURN_ON_OVERFLOW 		1
#define NORMALIZE_ON_OVERFLOW	2

template <class AdaptiveStepper, class DerivativeCallback>
int solve2(
	o2scl::ovector xg, int N, double h,								// grid info
	o2scl::omatrix& yg, o2scl::omatrix& ypg, o2scl::omatrix& yerrg, // solution
	AdaptiveStepper& adapt_stepper, 								// adaptive stepper
	DerivativeCallback& derivs,										// derivative callback
	int flag														// overflow behaviour
) {
	// solve (taken from O₂scl header "ode_iv_solve.h" ------------------------
	//
	size_t ntrial = 100000000;	// 10⁷
	double x0 = xg[0], x1 = xg[N - 2];
	double x = x0, xnext;
	int ret = 0, first_ret = 0;
	size_t nsteps = 0;
	xg[0] = x0;
	//
	o2scl::ovector ystart(2), dydx_start(2);
	ystart[0] = yg[0][0];
	ystart[1] = yg[0][1];
	derivs(x0,2,ystart,dydx_start);
	ypg[0][0] = dydx_start[0];	yerrg[0][0]=0.;
	ypg[0][1] = dydx_start[1];	yerrg[0][1]=0.;
	//
	for(int i = 1; i < N - 1 and ret == 0; i++)
	{
		xnext = x0 + (x1 - x0) * ((double)i)/((double)(N - 2));
		o2scl::omatrix_row y_row(yg,i);
		o2scl::omatrix_row dydx_row(ypg,i);
		o2scl::omatrix_row yerr_row(yerrg,i);
		//
		// Step until we reach the next grid point
		bool done = false;
		while (not done and ret == 0)
		{
			ret = adapt_stepper.astep_full(
				x,xnext,xg[i],h,2,ystart,dydx_start,
				y_row,yerr_row,dydx_row,derivs
			);
			nsteps++;
			if (ret != 0 and first_ret != 0)
					first_ret = ret;
			if (nsteps > ntrial)
			{
				std::string str="Too many steps required (ntrial="+o2scl::itos(ntrial)+
					", x="+o2scl::dtos(x)+") in ode_iv_solve::solve_grid().";
				std::cerr << str << std::endl;
			}
			//
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
			}
			//
			// Copy the results of the last step to the starting point for the next step
			x = xg[i];
			ystart[0] = y_row[0];	dydx_start[0] = dydx_row[0];
			ystart[1] = y_row[1];	dydx_start[1] = dydx_row[1];
			//
			// Decide if we have reached the grid point
			if (x1 > x0)
			{
				if (x >= xnext) done = true;
			}
			else
			{
				if (x <= xnext) done = true;
			}
		
	}
	}
	//------------------------------------------------------------------------
	return N;
}


#endif