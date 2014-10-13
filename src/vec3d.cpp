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

#include <cctype>
#include <iostream>
#include <sstream>

#include "misc.h"
#include "vec3d.h"

namespace geom
{

std::ostream & operator << (std::ostream & os, vec3d const & v)
{
    os << "( " << v.x << " " << v.y << " " << v.z << " )";
    return os;
}

std::istream & operator >> (std::istream & is, vec3d & v)
{
    // read buffer
    char c;
    
    // whole string
    std::string vec;
    
    // read characters
    while (true)
    {
        // read white characters as well
        is >> std::noskipws >> c;
        
        // throw away leading spaces and the opening bracket
        if (vec.empty() and isspace(c))
            continue;
        
        // check that we start with the opening parenthesis
        if (vec.empty() and c != '(')
            throw exception("A specification of a vector has to start with '('!");
        
        // add character to the whole string
        vec.push_back(c);
        
        // exit on right parenthesis
        if (c == ')')
            break;
    }
    
    // strip parentheses
    vec.erase(vec.begin());
    vec.erase(vec.end()-1);
    
    // read components
    std::istringstream iss(vec);
    iss >> v.x >> v.y >> v.z;
    
    return is;
}

double vec3d::dot (vec3d const & u, vec3d const & v)
{
    return u.x*v.x + u.y*v.y + u.z*v.z;
}

vec3d vec3d::cross (vec3d const & u, vec3d const & v)
{
    return vec3d
    (
        {
            u.y * v.z - u.z * v.y,
            u.z * v.x - u.x * v.z,
            u.x * v.y - u.y * v.x
        }
    );
}

vec3d operator - (vec3d const & u, vec3d const & v)
{
    return vec3d
    (
        {
            u.x - v.x,
            u.y - v.y,
            u.z - v.z
        }
    );
}

vec3d operator * (vec3d const & u, double a)
{
    return vec3d
    (
        {
            a * u.x,
            a * u.y,
            a * u.z
        }
    );
}

vec3d operator * (double a, vec3d const & u)
{
    return u * a;
}

double vec3d::norm (vec3d const & v)
{
    return std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

vec3d vec3d::normalize (vec3d const & v)
{
    return v * (1/norm(v));
}

}; // namespace geom
