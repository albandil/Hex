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

#ifndef HEX_VEC3D
#define HEX_VEC3D

#include <iostream>

namespace geom
{

/// Simple 3D vector class.
class vec3d
{
    public:
        
        double x, y, z;
        
        /// Dot product.
        static double dot (vec3d const & u, vec3d const & v);

        /// Cross product.
        static vec3d cross (vec3d const & u, vec3d const & v);
        
        /// Vector length.
        static double norm (vec3d const & v);
        
        /// Return normal vector.
        static vec3d normalize (vec3d const & v);
};

/// Write to stream.
std::ostream & operator << (std::ostream & os, vec3d const & v);

/// Read from stream.
std::istream & operator >> (std::istream & is, vec3d & v);

/// Vector difference.
vec3d operator - (vec3d const & u, vec3d const & v);

/// Vector scaling.
vec3d operator * (vec3d const & u, double a);
vec3d operator * (double a, vec3d const & u);

}; // namespace geom

#endif
