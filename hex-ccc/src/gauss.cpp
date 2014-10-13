//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2014, Jakub Benda, Charles University in Prague                    //
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

#include <algorithm>
#include <cstdio>
#include <complex>
#include <numeric>

#include <gsl/gsl_integration.h>

#include "arrays.h"
#include "gauss.h"
#include "misc.h"

//
// Forward explicit specializatons
//

template rArray GaussLegendre<double>::nodes (int points, double x1, double x2) const;
template rArray GaussLegendre<double>::weights (int points, double x1, double x2) const;

template cArray GaussLegendre<Complex>::nodes (int points, Complex x1, Complex x2) const;
template cArray GaussLegendre<Complex>::weights (int points, Complex x1, Complex x2) const;

//
// Templates
//

template <class T> std::vector<std::pair<double*,double*>> GaussLegendre<T>::data_ =
{
    std::make_pair(nullptr, nullptr),  // n = 0
    std::make_pair(nullptr, nullptr)   // n = 1
};

template <class T> int GaussLegendre<T>::nodes_and_weights (int points, const double* & vx, const double* & vw) const
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

template <class T> NumberArray<T> GaussLegendre<T>::nodes (int points, T x1, T x2) const
{
    // get the Gauss-Legendre nodes and weights
    const double *vx, *vw;
    points = nodes_and_weights (points, vx, vw);
    
    // prepare centre and half-width of the interval
    T hw = 0.5 * (x2 - x1);
    T ct = 0.5 * (x2 + x1);
    
    // prepare evaluation nodes
    NumberArray<T> xs(points);
    for (int dat = 0; dat < points; dat++)
    {
        // compute evaluation points
        xs[dat] = ct + vx[dat] * hw;
    }
    
    return xs;
}

template <class T> NumberArray<T> GaussLegendre<T>::weights (int points, T x1, T x2) const
{
    // get the Gauss-Legendre nodes and weights
    const double *vx, *vw;
    points = nodes_and_weights (points, vx, vw);
    
    // prepare halfwidth of the interval
    T hw = 0.5 * (x2 - x1);
    
    // prepare weights
    NumberArray<T> ws(points);
    for (int dat = 0; dat < points; dat++)
    {
        // scale weights
        ws[dat] = vw[dat] * hw;
    }
    
    return ws;
}
