//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2015, Jakub Benda, Charles University in Prague                    //
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

#include "hex-misc.h"
#include "hex-vec3d.h"

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
        
        // throw away leading spaces
        if (vec.empty() and std::isspace(c))
            continue;
        
        // check that we start with the opening parenthesis
        if (vec.empty() and c != '(')
            HexException("A specification of a vector has to start with '('!");
        
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

} // namespace geom
