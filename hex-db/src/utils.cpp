//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2016, Jakub Benda, Charles University in Prague                    //
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

#include <algorithm>

#include "hex-special.h"
#include "hex-vec3d.h"
#include "hex-version.h"

#include "utils.h"

double change_units (eUnit A, eUnit B)
{
    // no change
    if (A == B)
        return 1.;
    
    double ufactor = 1.;
    
    // transform to Rydbergs
    if (A == eUnit_au)
        ufactor *= 2.;
    if (A == eUnit_eV)
        ufactor *= 1./13.605692;
    
    // tranform from Rydbergs
    if (B == eUnit_au)
        ufactor *= 0.5;
    if (B == eUnit_eV)
        ufactor *= 13.605692;
    
    return ufactor;
}

double change_units (lUnit A, lUnit B)
{
    // no change
    if (A == B)
        return 1.;
    
    double ufactor = 1.;
    
    // transform to a.u.
    if (A == lUnit_cgs)
        ufactor *= 1./5.29177211e-9;
    
    // tranform from a.u.
    if (B == lUnit_cgs)
        ufactor *= 5.29177211e-9;
    
    return ufactor;
}

double change_units (aUnit A, aUnit B)
{
    // no change
    if (A == B)
        return 1.;
    
    double ufactor = 1.;
    
    // transform to radians
    if (A == aUnit_deg)
        ufactor *= special::constant::pi / 180;
    
    // transform from radians
    if (B == aUnit_deg)
        ufactor *= 180 / special::constant::pi;
    
    return ufactor;
}

std::string unit_name (eUnit u)
{
    switch (u)
    {
        case eUnit_au:
            return std::string("a.u.");
        case eUnit_eV:
            return std::string("eV");
        case eUnit_Ry:
            return std::string("Ry");
        default:
            return std::string("");
    }
}

std::string unit_name (lUnit u)
{
    switch (u)
    {
        case lUnit_au:
            return std::string("a.u.");
        case lUnit_cgs:
            return std::string("CGS");
        default:
            return std::string("");
    }
}

std::string unit_name (aUnit u)
{
    switch (u)
    {
        case aUnit_deg:
            return std::string("deg");
        case aUnit_rad:
            return std::string("rad");
        default:
            return std::string("");
    }
}
