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

#ifndef HEX_CCC_GAUSS
#define HEX_CCC_GAUSS

#include <vector>

#include "arrays.h"

/**
* @brief Gauss-Legendre quadrature.
* 
* This class's purpose is aid in fixed-order Gauss-Legendre quadrature.
* The functions "p_points" and "p_weights" return evaluation nodes
* and weights, respectively, and the functions "quad" do the fixed-order
* quadrature itself. Every call to p_points or p_weights resuts in call
* to gauss_nodes_and_weights, which uses GSL to get the requested data.
* The computed nodes and weights are stored in a cache table to allow
* faster subsequent computations.
*/
template <class T> class GaussLegendre
{
    public:
        
        /**
        * @brief Retrieve Gauss-Legendre data.
        * 
        * Precompute Gauss-Legendre quadrature data (if not already done
        * is some previous call) and return pointers to the chache table.
        * 
        * @param points Gauss-Legendre points half-count. If too low/high, the return value
        *               will contain the (used) nearest implemented value.
        * @param vx     On return, the Gauss-Legendre nodes (nonnegative half of them).
        * @param vw     On return, the corresponding Gauss-Legendre weights.
        */
        int nodes_and_weights (int points, const double* &vx, const double* &vw) const;
        
        /**
        * @brief Get Gauss-Legendre quadrature points in interval.
        * 
        * Map Gauss-Legendre points provided by @ref gauss_nodes_and_weights
        * to a complex interval @f$ (x_1,x_2) @f$.
        * 
        * @param points Number of Gauss-Legendre points.
        * @param x1 Left boundary of the interval.
        * @param x2 Right boundary of the interval.
        */
        NumberArray<T> nodes (int points, T x1, T x2) const;
        
        /**
        * @brief Get Gauss-Legendre quadrature weights in interval.
        * 
        * Map Gauss-Legendre weights provided by @ref gauss_nodes_and_weights
        * to a complex interval @f$ (x_1,x_2) @f$.
        * 
        * @param points Number of Gauss-Legendtre points.
        * @param x1 Left boundary of the interval.
        * @param x2 Right boundary of the interval.
        */
        NumberArray<T> weights (int points, T x1, T x2) const;
        
    private:
        
        // precomputed nodes and weights, common to all instances
        static std::vector<std::pair<double*,double*>> data_;
};

#endif
