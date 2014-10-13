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
        static_assert
        (
            dim <= 8,
            "Dimension of the n-cube has to be at most 8."
        );
        
        //
        // Constructors & destructor
        //
        
        ndCube () : lvl_(0), history_(nullptr) {}
        
        ndCube (ndCube<dim> const & cube) : ndCube(cube.lvl_, cube.history_) {}
        
        ndCube (ndCube<dim> && cube) : ndCube()
        {
            std::swap(lvl_, cube.lvl_);
            std::swap(history_, cube.history_);
        }
        
        ndCube (int lvl, char const * history) : lvl_(lvl), history_(nullptr)
        {
            if (lvl > 0)
            {
                history_ = new char [lvl_];
                std::memcpy(history_, history, lvl_ * sizeof(char));
            }
        }
        
        ~ndCube ()
        {
            if (history_ != nullptr)
                delete [] history_;
        }
        
        //
        // Assignment operator
        //
        
        ndCube<dim> & operator= (ndCube<dim> const & cube)
        {
            // delete current data
            if (history_ != nullptr)
            {
                assert(history_ != cube.history_);
                delete [] history_;
            }
            
            if (cube.lvl_ == 0)
            {
                // assignment of empty structure
                lvl_ = 0;
                history_ = nullptr;
            }
            else
            {
                // valid assignment
                lvl_ = cube.lvl_;
                history_ = new char [lvl_];
                std::memcpy(history_, cube.history_, lvl_);
            }
            
            return *this;
        }
        
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
            
            // new subdivision history for sub-cubes
            char new_history[lvl_ + 1];
            
            // copy the current subdivision history and append a new element for new subdivision
            if (lvl_ > 0)
                std::memcpy(new_history, history_, lvl_ * sizeof(char));
            
            // the last subdivision info
            char & subhistory = *(new_history + lvl_);
            
            // for all vertices (there is a sub-cube adjacent to each of them)
            for (int i = 0; i < nVertex(); i++)
            {
                // get subdivision history of this sub-cube
                subhistory = i;
                
                // add sub-cube
                subcubes.push_back(ndCube<dim>(lvl_ + 1, new_history));
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
            for (int i = 0; i < lvl_; i++)
            {
                // update edge length for this subdivision level
                edge_length *= 0.5;
                
                // for all coordinates
                for (int j = 0; j < dim; j++)
                {
                    // check if this coordinate is shifted with respect to sub-origin
                    if ((history_[i] & special::pow2(j)) != 0)
                        coords[j] += edge_length;
                }
            }
            return coords;
        }
        
        /// Length of edge of the cube.
        double edge () const { return std::pow(0.5, lvl_); }
        
        /// Volume of the hypercube.
        double volume () const { return std::pow(0.5, dim * lvl_); }
        
        /// Subdivision level.
        int level () const { return lvl_; }
        
        /// Subdivision history.
        char const * history () const { return history_; }
        
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
         * @brief Subdivision level.
         */
        char lvl_;
        
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
        char * history_;
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
    if (dCube.level() > 0)
        os << int(*(dCube.history())) << " ";
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
