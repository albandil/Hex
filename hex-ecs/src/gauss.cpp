//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2016, Jakub Benda, Charles University in Prague                    //
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

#include "hex-arrays.h"
#include "hex-misc.h"

#include "gauss.h"

std::vector<std::pair<Real*,Real*>> GaussLegendreData::data_ = {
    std::make_pair(nullptr, nullptr),  // n = 0
    std::make_pair(nullptr, nullptr)   // n = 1
};

GaussLegendreData::~GaussLegendreData ()
{
    for (std::pair<Real*,Real*> & p : data_)
    {
        if (p.first != nullptr)
            delete p.first;
        if (p.second != nullptr)
            delete p.second;
        
        p.first = p.second = nullptr;
    }
}

void GaussLegendreData::precompute_nodes_and_weights (int points) const
{
    # pragma omp critical
    for (int n = data_.size(); n <= points; n++)
    {
        Real* nodes = new Real [n];
        Real* weights = new Real [n];
        
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

void GaussLegendreData::gauss_nodes_and_weights (int points, const Real* & vx, const Real* & vw) const
{
    // enforce at least second order rule
    if (points < 2)
        HexException("[gauss_nodes_and_weights] Nor implemented for orders less than 2. Your input: %d.", points);
    
    // first of all generate any missing data
    if (points >= (int)data_.size())
        HexException("[gauss_nodes_and_weights] Not enough precomputed data; requested order = %d, available = %d", points, data_.size());
    
    // choose arrays
    vx = data_[points].first;
    vw = data_[points].second;
}

void GaussLegendre::scaled_nodes_and_weights (int points, Complex x1, Complex x2, Complex* xs, Complex* ws) const
{
    // get the Gauss-Legendre nodes and weights
    const Real *vx, *vw;
    gauss_nodes_and_weights(points, vx, vw);
    
    // prepare centre and half-width of the interval
    Complex hw = 0.5_r * (x2 - x1);
    Complex ct = 0.5_r * (x2 + x1);
    
    // prepare evaluation nodes and weights
    for (int dat = 0; dat < points; dat++)
    {
        xs[dat] = ct + vx[dat] * hw;
        ws[dat] = vw[dat] * hw;
    }
}
