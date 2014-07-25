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
              edge_(edge),
              vertices_(dim * nVertex()),
              subvertices_(dim * nSubvertex())
        {
            // pre-compute vertices
            init_();
        }
        
        ndCube (std::vector<double> origin, double edge)
            : origin_(origin),
              edge_(edge),
              vertices_(dim * nVertex()),
              subvertices_(dim * nSubvertex())
        {
            // pre-compute vertices
            init_();
        }
        
        ndCube (std::initializer_list<double> list)
            : origin_(dim),
              edge_(0.),
              vertices_(dim * nVertex()),
              subvertices_(dim * nSubvertex())
        {
            // first n-tuple is position of origin, then one number is edge length
            assert(list.size() == dim + 1);
            
            // copy initializer data
            auto iter = list.begin();
            for (int i = 0; i < dim; i++)
                origin_[i] = (*iter++);
            edge_ = (*iter++);
            
            // pre-compute vertices
            init_();
        }
        
        std::vector<double> const & origin () const
        {
            return origin_;
        }
        
        double edge () const
        {
            return edge_;
        }
        
        double volume () const
        {
            return std::pow(edge_, dim);
        }
        
        static int nVertex ()
        {
            return special::pow2(dim);
        }
        
        static int nSubvertex ()
        {
            return special::pow3(dim);
        }
        
        double const * vertex (int i) const
        {
            assert(i < nVertex());
            return vertices_.data() + dim * i;
        }
        
        double const * subvertex (int i) const
        {
            assert(i < nSubvertex());
            return subvertices_.data() + dim * i;
        }
        
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
        
    private:
        
        void init_ ()
        {
            init_vertices_();
            init_subvertices_();
        }
        
        void init_vertices_ ()
        {
            // start with origin
            for (int i = 0; i < dim; i++)
                vertices_[i] = origin_[i];
            
            // extrude already existing n-cube into next dimensions
            for (int d = 0; d < dim; d++)
            {
                // current position in 'vertices'
                int pos = dim * special::pow2(d);
                
                // extrude existing points
                for (int p = 0; p < special::pow2(d); p++)
                {
                    // copy existing point and shift the 'd'-th coordinate of this point
                    for (int i = 0; i < dim; i++)
                        vertices_[pos + p * dim + i] = vertices_[p * dim + i];
                    vertices_[pos + p * dim + d] += edge_;
                }
            }
        }
        
        void init_subvertices_ ()
        {
            // start with origin
            for (int i = 0; i < dim; i++)
                subvertices_[i] = origin_[i];
            
            // extrude already existing n-cube into next dimensions
            for (int d = 0; d < dim; d++)
            {
                // current position in 'vertices'
                int pos = dim * special::pow3(d);
                
                // extrude existing points
                for (int p = 0; p < special::pow3(d); p++)
                {
                    // copy existing point and shift the 'd'-th coordinate of this point (halfway)
                    for (int i = 0; i < dim; i++)
                        subvertices_[pos + p * dim + i] = subvertices_[p * dim + i];
                    subvertices_[pos + p * dim + d] += 0.5 * edge_;
                    
                    // copy existing point and shift the 'd'-th coordinate of this point (full)
                    for (int i = 0; i < dim; i++)
                        subvertices_[pos + (p + 1) * dim + i] = subvertices_[p * dim + i];
                    subvertices_[pos + (p + 1) * dim + d] += edge_;
                }
            }
        }
        
        std::vector<double> origin_;
        
        double edge_;
        
        std::vector<double> vertices_;
        
        std::vector<double> subvertices_;
};

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

#define Unit_5Cube geom::ndCube<5>()
#define Unit_6Cube geom::ndCube<6>()

#endif // HEX_NDCUBE
