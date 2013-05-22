/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2012                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <algorithm>
#include <cstdio>
#include <complex>
#include <numeric>

#include "arrays.h"
#include "bspline.h"

extern rArrays GaussLegendreRoots, GaussLegendreWeights;
bool GaussLegendreRoots_init();
bool GaussLegendreWeights_init();

int p_roots(int points, const double* &vx, const double* &vw)
{
	static bool init1 = GaussLegendreRoots_init();
	static bool init2 = GaussLegendreWeights_init();
	
	// available Gauss-Legendre data
	int N = std::min(GaussLegendreRoots.size(), GaussLegendreWeights.size());
	if (points < 3)
	{
		points = 3;
	}
	if (points >= N)
	{
		std::cerr << "[quad] Warning: insufficent Gauss-Legendre points, using " << N - 1 << std::endl;
		points = N - 1;
	}
	
	// get correct pointer
	vx = GaussLegendreRoots[points].data();
	vw = GaussLegendreWeights[points].data();
	
	return points;
}


cArray p_points(int& points, Complex x1, Complex x2)
{
	// get the Gauss-Legendre nodes and weights
	const double *vx, *vw;
	points = p_roots(points, vx, vw);
	
	// prepare centre and half-width of the interval
	Complex hw = 0.5 * (x2 - x1);
	Complex ct = 0.5 * (x2 + x1);
	
	// prepare evaluation nodes
	cArray xs(points);
	for (int dat = 0; dat < points; dat++)
	{
		// compute evaluation points
		xs[dat] = ct + vx[dat] * hw;
	}
	
	return xs;
}


cArray p_weights(int& points, Complex x1, Complex x2)
{
	// get the Gauss-Legendre nodes and weights
	const double *vx, *vw;
	points = p_roots(points, vx, vw);
	
	// prepare halfwidth of the interval
	Complex hw = 0.5 * (x2 - x1);
	
	// prepare weights
	cArray ws(points);
	for (int dat = 0; dat < points; dat++)
	{
		// scale weights
		ws[dat] = vw[dat] * hw;
	}
	
	return ws;
}


Complex quad (
	void (*f)(int, Complex*, Complex*, void*), void *data,
	int points, int iknot, Complex x1, Complex x2
){
	// check boundaries
	if (x1.real() < Bspline::ECS().t(iknot).real() or Bspline::ECS().t(iknot+1).real() < x1.real() or
		x2.real() < Bspline::ECS().t(iknot).real() or Bspline::ECS().t(iknot+1).real() < x2.real())
	{
		std::fprintf(stderr, "[quad] Error: boundaries not for this iknot!\n");
		exit(1);
	}
	
	// get evaluation points and weights
	cArray xs = p_points(points, x1, x2);
	cArray ws = p_weights(points, x1, x2);
	
	// evaluate the function
	Complex values[points];
	f(points, xs.data(), values, data);
	
	// sum the results
	Complex result = std::inner_product (
		values, values + points,
		ws.begin(),
		Complex(0.)
	);
	
	return result;
}
