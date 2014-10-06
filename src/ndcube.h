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

#ifndef HEX_NDCUBE
#define HEX_NDCUBE

#include <algorithm>
#include <bitset>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "special.h"

namespace geom
{

/**
 * @brief N-dimensional cube.
 * 
 * This class wraps some utility functions for work with n-dimensional cubes.
 * The attributes subdivision history of the original unit hypercube. The core routine
 * used in Hex is the member function @ref subdivide, which decomposes the cube
 * into a set of 2^dim smaller cubes by cutting every edge in the middle. This
 * is used in the adaptive sparse grid integrator @ref spgrid::SparseGrid.
 */
template <int dim> class ndCube
{
    public:
        
        // throw a compile-time error if the class is used with an unsupported template parameter
        static_assert
        (
            dim > 0,
            "Dimension of the n-cube has to be at least 1."
        );
        
        //
        // Constructors
        //
        
        ndCube () : history_({0}) {}
        
        ndCube (std::vector<std::bitset<dim>> history) : history_(history) {}
        
        //
        // Subdivision
        //
        
        /**
         * @brief Subdivision of the cube.
         * 
         * This function will return a set of smaller cubes that arise
         * when splitting every edge of the original cube into two equal parts.
         */
        std::vector<ndCube<dim>> subdivide () const
        {
            std::vector<ndCube<dim>> subcubes;
            
            // copy the subdivision history and append a new element for this subdivision
            std::vector<std::bitset<dim>> new_history(history_.size() + 1);
            std::copy(history_.begin(), history_.end(), new_history.begin());
            std::bitset<dim> & subhistory = new_history.back();
            
            for (int i = 0; i < nVertex(); i++)
            {
                // get subdivision history of the sub-cube
                for (int j = 0; j < dim; j++)
                    subhistory[j] = (i bitand special::pow2(j));
                
                // add sub-cube
                subcubes.push_back(ndCube<dim>(new_history));
            }
            
            return subcubes;
        }
        
        //
        // Getters
        //
        
        /// Coordinates of the origin of the cube.
        std::vector<double> const origin () const
        {
            std::vector<double> coords(dim,0.0);
            double edge_length = 1.0;
            for (int i = 1; i < history_.size(); i++)
            {
                // update edge length for this subdivision level
                edge_length *= 0.5;
                
                // for all coordinates
                for (int j = 0; j < dim; j++)
                {
                    // check if this coordinate is shifted with respect to sub-origin
                    if (history_[i][j] == 1)
                        coords[j] += edge_length;
                }
            }
            return coords;
        }
        
        /// Length of edge of the cube.
        double edge () const { return std::pow(0.5, history_.size() - 1); }
        
        /// Volume of the hypercube.
        double volume () const { return std::pow(0.5, dim * (history_.size() - 1)); }
        
        /// Subdivision history.
        std::vector<std::bitset<dim>> history () const { return history_; }
        
        /// Coordinates of the centre of the cube.
        std::vector<double> centre () const
        {
            std::vector<double> coords = origin();
            for (double & coord : coords)
                coord += 0.5 * edge();
            return coords;
        }
        
        /// Number of vertices for 'dim'-dimensional cube.
        static int nVertex () { return special::pow2(dim); }
        
    private:
        
        /**
         * @brief Subdivision history.
         * 
         * This vector contains one element for every refinement that took place, starting from some
         * original cube. Every element has exactly 'D' fields (equal to zero or one each), where 'D'
         * is the dimension count. For a given dimension, zero indicates equal coordinate of origin and
         * suborigin, whereas one indicates that the suborigin has been shifted to half of the edge.
         * Strings of zeros and ones from stacked history entries compose binary numbers (less than one)
         * which uniquely define position of the suborigin.
         */
        std::vector<std::bitset<dim>> history_;
};

/**
 * @brief Stream output of the ndCube class.
 * 
 * This is an overload of the stream output operator. The representation of the
 * class is a list (enclosed in brackets "[" and "]") of positions (which are
 * coordinates enclosed in parentheses "(" and ")").
 */
template <int dim> std::ostream & operator << (std::ostream & os, ndCube<dim> const & dCube)
{
    os << "[(";
    for (int i = 0; i < dim; i++)
        if (i == 0) os << dCube.origin()[i]; else os << "," << dCube.origin()[i];
    os << "),(";
    for (int i = 0; i < dim; i++)
        if (i == 0) os << dCube.origin()[i] + dCube.edge(); else os << "," << dCube.origin()[i] + dCube.edge();
    os << ")]";
    return os;
}

}; // namespace geom

#endif // HEX_NDCUBE
