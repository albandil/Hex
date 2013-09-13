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

#include <algorithm>
#include <cstdio>
#include <complex>
#include <numeric>

#include <gsl/gsl_integration.h>

#include "arrays.h"
#include "bspline.h"

std::vector<std::pair<double*,double*>> gauss_data = {
    std::make_pair(nullptr, nullptr),  // n = 0
    std::make_pair(nullptr, nullptr)   // n = 1
};

int gauss_nodes_and_weights(int points, const double* & vx, const double* & vw)
{
    // first of all generate any missing data
    # pragma omp critical
    if (points >= gauss_data.size())
    {
        for (int n = gauss_data.size(); n <= points; n++)
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
            
            gauss_data.push_back(std::make_pair(nodes,weights));
        }
    }
    
    // return the arrays
    vx = gauss_data[points].first;
    vw = gauss_data[points].second;
    return points;
}

cArray p_points(int& points, Complex x1, Complex x2)
{
    // get the Gauss-Legendre nodes and weights
    const double *vx, *vw;
    points = gauss_nodes_and_weights(points, vx, vw);
    
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
    points = gauss_nodes_and_weights(points, vx, vw);
    
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
        throw exception ("[quad] Error: boundaries not for this iknot!");
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
