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

#include <cctype>
#include <iostream>
#include <sstream>

#include "misc.h"
#include "vec3d.h"

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

double dot (vec3d const & u, vec3d const & v)
{
    return u.x*v.x + u.y*v.y + u.z*v.z;
}

vec3d cross (vec3d const & u, vec3d const & v)
{
    return vec3d (
        {
            u.y * v.z - u.z * v.y,
            u.z * v.x - u.x * v.z,
            u.x * v.y - u.y * v.x
        }
    );
}

vec3d operator - (vec3d const & u, vec3d const & v)
{
    return vec3d (
        {
            u.x - v.x,
            u.y - v.y,
            u.z - v.z
        }
    );
}

vec3d operator * (vec3d const & u, double a)
{
    return vec3d (
        {
            a * u.x,
            a * u.y,
            a * u.z
        }
    );
}

double norm (vec3d const & v)
{
    return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

vec3d normalize(vec3d const & v)
{
    return v * (1/norm(v));
}
