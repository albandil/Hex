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
    n_ = V.n_;
    k_ = V.k_;
    
    return *this;
}

bool DistortingPotential::operator== (DistortingPotential const & V) const
{
    return n_ == V.n_ and k_ == V.k_;
}

double DistortingPotential::operator() (double x) const
{
    if (n_ == 0 and k_ == 0)
        return 0.;
    
    if (n_ == 1)
        return -(1.+1./x)*exp(-2.*x);
    
    if (n_ == 2)
        return -((0.125*x+0.25)*x+0.75+1./x)*exp(-x);
    
    if (n_ == 3)
        return -(((((4*x-12)*x+108)*x+324)*x+1215)/2187 + 1/x)*exp(-2*x/3);
    
    printf("U not implemented for n = %d\n", n_);    // FIXME
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
    if (n_ == 1)
        return -1.;
    
    printf("U not implemented for n = %d\n", n_);    // FIXME
    abort();
}

double DistortingPotential::getFarRadius() const
{
    return 400; // FIXME
    
/*
    // determine from parameters of the potential
    #define U_THRESHOLD    1e-50
    #define U_MAX_ITERS    100
    
    // hunt for low value
    double far = 1., far_value;
    while ((far_value = fabs((*this)(far *= 2))) > U_THRESHOLD)
        // continue hunting ;
    
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
*/
}

void DistortingPotential::toFile(const char* filename) const
{
    rArray grid = linspace(0., getFarRadius(), 1000);
    rArray vals(grid.size());
    for (int i = 0; i < (int)grid.size(); i++)
        vals[i] = (*this)(grid[i]);
    write_array(grid, vals, filename);
}
