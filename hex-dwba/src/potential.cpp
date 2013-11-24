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

const rArrays Ucoeffs = {
  // empty
    {},
  // U(1s)
    { -1. },
  // U(2s)
    { -0.125, -0.25, -0.25 },
  // U(3s)
    { -0.001828989483310471, 0.005486968449931412, -0.04938271604938271, -0.1481481481481481, -0.5555555555555556 },
  // U(4s)
    { -6.781684027777778E-6, 1.3563368055555555E-4, -0.00146484375, 0.00390625, -0.0234375, -0.09375, -0.4375 },
  // U(5s)
    { -9.102222222222223E-9, 5.006222222222222E-7, -1.2288E-5, 1.4563555555555555E-4, -0.001024, 0.00256, -0.0128, -0.064, -0.36 }
  // ... TODO the rest ...
};

DistortingPotential DistortingPotential::operator= (DistortingPotential const& V)
{
    // copy data
    n_ = V.n_;
    k_ = V.k_;
    return *this;
}

bool DistortingPotential::operator== (DistortingPotential const & V) const
{
    // check equality
    return n_ == V.n_ and k_ == V.k_;
}

double DistortingPotential::operator () (double x) const
{
    // skip empty potential
    if (n_ == 0 and k_ == 0)
        return 0.;
    
    // disallow distorsion by free states
    if (k_ != 0)
        throw exception ("U not implemented for k != 0.");
    
    // stop if not enough precomputed coefficients
    if (n_ >= (int)Ucoeffs.size())
        throw exception ("U not implemented for n = %d > %d.", n_, Ucoeffs.size() - 1);
    
    // get correct polynomial coefficients
    rArray const & coeffs = Ucoeffs[n_];
    
    // evaluate (P(x) - 1/x) * exp(-2*x/n)
    double eval = 0;
    for (double factor : coeffs)
        eval = eval * x + factor;
    return (eval - 1./x) * exp(-2.*x/n_);
}

double DistortingPotential::plusMonopole (double x) const
{
    // in origin returm true asymptotic
    if (x == 0.)
        return getConstant();
    
    // otherwise compute
    return 1./x + (*this)(x);
}

double DistortingPotential::getConstant () const
{
    // skip empty potential
    if (n_ == 0 and k_ == 0)
        return 0.;
    
    // disallow distorsion by free states
    if (k_ != 0)
        throw exception ("U not implemented for k != 0.");
    
    // stop if not enough precomputed coefficients
    if (n_ >= (int)Ucoeffs.size())
        throw exception ("U not implemented for n = %d > %d.", n_, Ucoeffs.size() - 1);
    
    // return the x -> 0 limiting value withou the monopole term
    return Ucoeffs[n_].back();
}

double DistortingPotential::getFarRadius () const
{
    // if the rmax has been overriden, use the supplied value
    if (rmax_ > 0.)
        return rmax_;
    
    //
    // otherwise compute rmax using hunt & bisect run
    //
    
    // determine from parameters of the potential
    #define U_THRESHOLD    1e-100
    #define U_MAX_ITERS    1000
    
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

void DistortingPotential::toFile (const char* filename) const
{
    // generate evaluation grid
    rArray grid = linspace(0., getFarRadius(), 1000);
    rArray vals(grid.size());
    
    // evaluate
    for (int i = 0; i < (int)grid.size(); i++)
        vals[i] = (*this)(grid[i]);
    
    // write to text file
    write_array(grid, vals, filename);
}
