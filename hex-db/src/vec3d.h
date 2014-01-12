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

/// Simple 3D vector class.
struct vec3d
{
    double x;
    double y;
    double z;
};

/// Write to stream.
std::ostream & operator << (std::ostream & os, vec3d const & v);

/// Read from stream.
std::istream & operator >> (std::istream & is, vec3d & v);

/// Dot product.
double dot (vec3d const & u, vec3d const & v);

/// Cross product.
vec3d cross (vec3d const & u, vec3d const & v);

/// Vector difference.
vec3d operator - (vec3d const & u, vec3d const & v);

/// Vector scaling.
vec3d operator * (vec3d const & u, double a);

/// Vector length.
double norm (vec3d const & v);

/// Return normal vector.
vec3d normalize (vec3d const & v);

#endif
