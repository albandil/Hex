/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2014                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <algorithm>
#include <cstdio>
#include <complex>
#include <numeric>

#include <gsl/gsl_integration.h>

#include "arrays.h"
#include "gauss.h"
#include "misc.h"

std::vector<std::pair<double*,double*>> GaussLegendreData::data_ = {
    std::make_pair(nullptr, nullptr),  // n = 0
    std::make_pair(nullptr, nullptr)   // n = 1
};

int GaussLegendreData::gauss_nodes_and_weights (int points, const double* & vx, const double* & vw) const
{
    // enforce at least second order rule
    if (points < 2)
        throw exception ("[gauss_nodes_and_weights] Nor implemented for orders less than 2. Your input: %d.", points);
    
    // first of all generate any missing data
    # pragma omp critical
    if (points >= (int)data_.size())
    {
        for (int n = data_.size(); n <= points; n++)
        {
            double* nodes = new double [n];
            double* weights = new double [n];
            
            gsl_integration_glfixed_table *t = gsl_integration_glfixed_table_alloc(n);
            
            // add anti/symmetrically the nodes and weights
            for (int i = 0; i < n/2; i++)
            {
                nodes[n/2-i-1] = -t->x[i+n%2];
                nodes[n/2+i+n%2] = t->x[i+n%2];
                
                weights[n/2-i-1] = t->w[i+n%2];
                weights[n/2+i+n%2] = t->w[i+n%2];
            }
            
            // for odd 'n' add also the node and weight in zero
            if (n % 2 != 0)
            {
                nodes[n/2] = 0.;
                weights[n/2] = t->w[0];
            }
            
            data_.push_back(std::make_pair(nodes,weights));
        }
    }
    
    // return the arrays
    vx = data_[points].first;
    vw = data_[points].second;
    return points;
}

cArray GaussLegendre::p_points (int points, Complex x1, Complex x2) const
{
    // get the Gauss-Legendre nodes and weights
    const double *vx, *vw;
    points = gauss_nodes_and_weights (points, vx, vw);
    
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

cArray GaussLegendre::p_weights (int points, Complex x1, Complex x2) const
{
    // get the Gauss-Legendre nodes and weights
    const double *vx, *vw;
    points = gauss_nodes_and_weights (points, vx, vw);
    
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
