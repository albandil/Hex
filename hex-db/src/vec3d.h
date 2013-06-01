/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2013                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HEX_VEC3D
#define HEX_VEC3D

#include <iostream>

struct vec3d
{
	double x;
	double y;
	double z;
};

// write to stream
std::ostream & operator << (std::ostream & os, vec3d const & v);

// read from stream
std::istream & operator >> (std::istream & is, vec3d & v);

#endif
