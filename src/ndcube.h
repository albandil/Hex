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
 * The attributes keep the origin and edge length of the cube. The core routine
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
        
        ndCube (double edge = 1.)
            : origin_(dim),
              edge_(edge)
        {}
        
        ndCube (std::vector<double> const & origin, double edge)
            : origin_(origin),
              edge_(edge)
        {}
        
        ndCube (std::initializer_list<double> list)
            : origin_(dim),
              edge_(0.)
        {
            // first n-tuple is position of origin, then one number is edge length
            assert(list.size() == dim + 1);
            
            // copy initializer data
            auto iter = list.begin();
            for (int i = 0; i < dim; i++)
                origin_[i] = (*iter++);
            edge_ = (*iter++);
        }
        
        /**
         * @brief Subdivision of the cube.
         * 
         * This function will return a set of smaller cubes that arise
         * when splitting every edge of the original cube into two equal parts.
         */
        std::vector<ndCube<dim>> subdivide () const
        {
            std::vector<ndCube<dim>> subcubes;
            std::vector<double> suborigin(dim);
            
            for (int i = 0; i < nVertex(); i++)
            {
                // get coordinates of the sub-origin
                for (int j = 0; j < dim; j++)
                    suborigin[j] = origin_[j] + ((i & special::pow2(j)) != 0 ? 1 : 0) * 0.5 * edge_;
                
                // add sub-cube
                subcubes.push_back(ndCube<dim>(suborigin, 0.5 * edge_));
            }
            
            return subcubes;
        }
        
        // getters
        
        std::vector<double> const & origin () const { return origin_; }
        double edge () const { return edge_; }
        double volume () const { return std::pow(edge_, dim); }
        static int nVertex () { return special::pow2(dim); }
        
    private:
        
        /// Coordinates of the cube's origin.
        std::vector<double> origin_;
        
        /// Length of an edge.
        double edge_;
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
    os << "[";
    for (int i = 0; i < dCube.nVertex(); i++)
    {
        os << "(";
        for (int j = 0; j < dim; j++)
        {
            os << dCube.vertex(i)[j];
            if (j < dim - 1)
                os << " ";
        }
        os << ")";
        if (i < dCube.nVertex() - 1)
            os << ",";
    }
    os << "]";
    return os;
}

}; // namespace geom

// some abbreviations
#define Unit_2Cube geom::ndCube<2>()
#define Unit_3Cube geom::ndCube<3>()
#define Unit_4Cube geom::ndCube<4>()
#define Unit_5Cube geom::ndCube<5>()
#define Unit_6Cube geom::ndCube<6>()

#endif // HEX_NDCUBE
