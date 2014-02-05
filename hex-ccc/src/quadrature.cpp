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

#include <tuple>

#include <gsl/gsl_integration.h>

#include "basis.h"
#include "quadrature.h"

QuadratureRule::QuadratureRule (LaguerreBasis const & basis, double E)
    : basis_(basis)
{
    // quadrature order
    int order = 20; // ( FIXME : we want it more dynamic )
    
    // precompute nodes and weights on the interval (-1,1)
    gsl_integration_glfixed_table *t = gsl_integration_glfixed_table_alloc (order);
    
    // maximal projectile (or propagator) angular momentum
    maxpell_ = basis.size() - 1; // ( FIXME : equal to "maxell" now )
    
    // get iterator
    size_t pos = 0;
    
    // for all angular momenta of the basis
    for (int ell = 0; ell < basis.size(); ell++)
    {
        // for all angular momenta of the projectile (or propagator)
        for (int l = 0; l < maxpell_; l++)
        {
            // for all eigenstates of the basis
            for (int n = 1; n <= basis.size(ell); n++)
            {
                // compute the position of the singularity
                double sgE = E - basis.energy(ell, n);
                double k0 = (sgE > 0.) ? sqrt(sgE) : 0.;
                
                // half-interval around the singularity to handle separately
                double dk = std::min(0.05, k0); // FIXME : more flexible
                
                // how far after the singularity interval to integrate
                double Dk = sqrt(E); // FIXME : change to something more reasonable
                
                // compute evaluation nodes and weights for the interval (0, k0-dk), if any
                if (k0 != 0.)
                {
                    // add index
                    indices_[std::make_tuple(ell, l, n)] = pos;
                    
                    // allocate space for next quadrature coefficients
                    npoints_.push_back(3 * order);
                    nodes_.resize(nodes_.size() + npoints_.back());
                    weights_.resize(weights_.size() + npoints_.back());
                    
                    linearWeightsAndNodes ( E, ell, n, t,
                        rArrayView(order, nodes_.data() + pos),
                        rArrayView(order, weights_.data() + pos),
                        0., k0 - dk
                    );
                    pos += order;
                    
                    // compute evaluation nodes and weights for the interval (k0-d0, k0+dk)
                    linearWeightsAndNodes ( E, ell, n, t,
                        rArrayView(order, nodes_.data() + pos),
                        rArrayView(order, weights_.data() + pos),
                        k0 - dk, k0 + dk
                    );
                    pos += order;
                    
                    // compute evaluation nodes and weights for the interval (sgE+dE, infty)
                    linearWeightsAndNodes ( E, ell, n, t,
                        rArrayView(order, nodes_.data() + pos),
                        rArrayView(order, weights_.data() + pos),
                        k0 + dk, k0 + dk + Dk
                    );
                    pos += order;
                }
                else
                {
                    // add index
                    indices_[std::make_tuple(ell, l, n)] = pos;
                    
                    // allocate space for next quadrature coefficients
                    npoints_.push_back(3 * order);
                    nodes_.resize(nodes_.size() + npoints_.back());
                    weights_.resize(weights_.size() + npoints_.back());
                    
                    // compute evaluation nodes and weights for the interval (sgE+dE, infty)
                    linearWeightsAndNodes ( E, ell, n, t,
                        rArrayView(order, nodes_.data() + pos),
                        rArrayView(order, weights_.data() + pos),
                        0., Dk
                    );
                    pos += order;
                }
            }
        }
    }
    
    // sanity check
#ifndef NDEBUG
    size_t npoints_total = std::accumulate(npoints_.begin(), npoints_.end(), 0);
    assert (npoints_total == nodes_.size());
    assert (npoints_total == weights_.size());
#endif
}

void QuadratureRule::linearWeightsAndNodes (
    double E, int ell, int n,
    gsl_integration_glfixed_table const * t,
    rArrayView nodes, rArrayView weights,
    double a, double b
){
    // available energy to the propagator
    double Eg = 2 * E - basis_.energy(ell, n); // a.u.
    
    // parameters of the interval
    double m = 0.5 * (b + a);   // centre of the interval
    double d = 0.5 * (b - a);   // half-width of the interval
    
    // auxiliary variable (the evaluation linear momentum)
    double k;
    
    // for all available points
    for (int i = 0; i < t->n/2; i++)
    {
        k = m - d * t->x[i];
        nodes[t->n/2 - i - 1]   = k;
        weights[t->n/2 - i - 1] = d * t->w[i] / (Eg - 0.5*k*k);
        
        k = m + d * t->x[i];
        nodes[t->n/2 + i]       = k;
        weights[t->n/2 + i]     = d * t->w[i] / (Eg - 0.5*k*k);
    }
}

const rArrayView QuadratureRule::nodes (int ell, int l, int n) const
{
    assert (ell < basis_.size());
    assert (l < basis_.size()); // FIXME
    assert (ell < 0 or n < basis_.size(ell));
    assert (n != 0);
    
    if (ell < 0)
    {
        return nodes_;
    }
    
    else if (l < 0)
    {
        return nodes_.slice (
            indices_.at(std::make_tuple(ell, 0, 1)),
            (ell + 1 == basis_.size()) ? nodes_.size() : indices_.at(std::make_tuple(ell + 1, 0, 1))
        );
    }
    
    else if (n < 0)
    {
        return nodes_.slice (
            indices_.at(std::make_tuple(ell, l, 1)),
            (l + 1 == basis_.size()) ? nodes_.size() : indices_.at(std::make_tuple(ell, l + 1, 1)) // FIXME
        );
    }
    
    else return nodes_.slice (
        indices_.at(std::make_tuple(ell, l, n)),
        (n + 1 == basis_.size(ell)) ? nodes_.size() : indices_.at(std::make_tuple(ell, l, n + 1))
    );
}

const rArrayView QuadratureRule::weights (int ell, int l, int n) const
{
    assert (ell < basis_.size());
    assert (l < basis_.size()); // FIXME
    assert (ell < 0 or n < basis_.size(ell));
    assert (n != 0);
    
    if (ell < 0)
    {
        return weights_;
    }
    
    else if (l < 0)
    {
        return weights_.slice (
            indices_.at(std::make_tuple(ell, 0, 1)),
            (ell + 1 == basis_.size()) ? nodes_.size() : indices_.at(std::make_tuple(ell + 1, 0, 1))
        );
    }
    
    else if (n < 0)
    {
        return weights_.slice (
            indices_.at(std::make_tuple(ell, l, 1)),
            (l + 1 == basis_.size()) ? nodes_.size() : indices_.at(std::make_tuple(ell, l + 1, 1)) // FIXME
        );
    }
    
    else return weights_.slice (
        indices_.at(std::make_tuple(ell, l, n)),
        (n + 1 == basis_.size(ell)) ? nodes_.size() : indices_.at(std::make_tuple(ell, l, n + 1))
    );
}
