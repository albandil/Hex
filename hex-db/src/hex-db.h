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

#ifndef HEX_HEX_DB
#define HEX_HEX_DB

#include <map>
#include <string>
#include <vector>

#include "arrays.h"
#include "complex.h"
#include "interfaces.h"

/// Energy units
enum eUnit
{
    eUnit_Ry,    // Rydberg (13.605692 eV, default)
    eUnit_au,    // Hartree (2 Ry)
    eUnit_eV    // electron-Volt
};

/// Output (length) units
enum lUnit
{
    lUnit_au,    // atomic units (Bohr radius a₀=5.29x10⁻⁹ cm)
    lUnit_cgs    // centimeters (1 cm = (1cm/a₀) a₀)
};

// Angular units
enum aUnit
{
    aUnit_deg,    // degrees
    aUnit_rad    // radians
};

// global unit settings
extern eUnit Eunits;
extern lUnit Lunits;
extern aUnit Aunits;

/**
 * Run the computations.
 */
int hex_run
(
    std::vector<std::string> const & vars,
    std::map<std::string,std::string> const & sdata
);

#endif
