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

#include <gsl/gsl_sf.h>

#include "arrays.h"
#include "potential.h"

DistortingPotential DistortingPotential::operator= (DistortingPotential const& V)
{
	n = V.n;
	k = V.k;
	
	return *this;
}

bool DistortingPotential::operator== (DistortingPotential const & V) const
{
	return n == V.n and k == V.k;
}

double DistortingPotential::operator() (double x) const
{
	if (n == 0 and k == 0)
		return 0.;
	
	if (n == 1)
		return -(1.+1./x)*exp(-2.*x);
	
	if (n == 2)
		return -(0.125*x*x+0.25*x+0.75+1./x)*exp(-x);
	
	printf("U not implemented for n = %d\n", n);	// FIXME
	abort();
}

double DistortingPotential::plusMonopole(double x) const
{
	if (x == 0.)
		return getConstant();
	
	return 1./x + (*this)(x);
}

double DistortingPotential::getConstant() const
{
	if (n == 1)
		return -1.;
	
	printf("U not implemented for n = %d\n", n);	// FIXME
	abort();
}

double DistortingPotential::getFarRadius() const
{
	// determine from parameters of the potential
	#define U_THRESHOLD    1e-50
	#define U_MAX_ITERS    100
	
	// hunt for low value
	double far = 1., far_value;
	while ((far_value = fabs((*this)(far *= 2))) > U_THRESHOLD)
		/* continue hunting */;
	
	// bisect for exact value
	double near = 1.;
	for (int i = 0; i < U_MAX_ITERS and near != far; i++)
	{
		double middle = (near + far) * 0.5;
		double middle_val = fabs((*this)(middle));
		
		if (U_THRESHOLD > middle_val)
			far = middle;
			
		if (middle_val > U_THRESHOLD)
			near = middle;
	}
	
	return (near + far) * 0.5;
}

void DistortingPotential::toFile(const char* filename) const
{
	rArray grid = linspace(0., getFarRadius(), 1000);
	rArray vals(grid.size());
	for (int i = 0; i < (int)grid.size(); i++)
		vals[i] = (*this)(grid[i]);
	write_array(grid, vals, filename);
}

DistortedWave DistortingPotential::getDistortedWave(double kn, int ln) const
{
	return DistortedWave(kn, ln, *this);
}

IrregularWave DistortingPotential::getIrregularWave(double kn, int ln) const
{
	return IrregularWave(kn, ln, *this);
}

ForbiddenWave DistortingPotential::getForbiddenWave(double kn, int ln) const
{
	return ForbiddenWave(kn, ln, *this);
}

HyperbolicWave DistortingPotential::getHyperbolicWave(double kn, int ln) const
{
	return HyperbolicWave(kn, ln, *this);
}
