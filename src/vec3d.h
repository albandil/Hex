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
