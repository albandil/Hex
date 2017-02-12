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

#include <map>
#include <string>
#include <vector>

// --------------------------------------------------------------------------------- //

#include "hex-interpolate.h"
#include "hex-special.h"
#include "hex-version.h"

// --------------------------------------------------------------------------------- //

#include "../interfaces.h"
#include "../quantities.h"
#include "../utils.h"

// --------------------------------------------------------------------------------- //

createNewScatteringQuantity(SpinAsymmetry, "asy")

// --------------------------------------------------------------------------------- //

std::string SpinAsymmetry::description ()
{
    return "Spin asymetry.";
}

std::vector<std::string> SpinAsymmetry::dependencies ()
{
    return std::vector<std::string>
    {
        "dcs"
    };
}

std::vector<std::pair<std::string,std::string>> SpinAsymmetry::params ()
{
    return std::vector<std::pair<std::string,std::string>>
    {
        {"ni", "Initial atomic principal quantum number."},
        {"li", "Initial atomic orbital quantum number."},
        {"mi", "Initial atomic magnetic quantum number."},
        {"nf", "Final atomic principal quantum number."},
        {"lf", "Final atomic orbital quantum number."},
        {"mf", "Final atomic magnetic quantum number."},
        {"Ei", "Projectile impact energy (Rydberg)."},
        {"theta", "Scattering angles for which to compute the spin asymetry."}
    };
}

std::vector<std::string> SpinAsymmetry::vparams ()
{
    return std::vector<std::string>
    {
        "theta"
    };
}

// --------------------------------------------------------------------------------- //

bool SpinAsymmetry::initialize (sqlitepp::session & db)
{
    return ScatteringQuantity::initialize(db);
}

bool SpinAsymmetry::createTable ()
{
    return ScatteringQuantity::createTable();
}

bool SpinAsymmetry::updateTable ()
{
    return ScatteringQuantity::updateTable();
}

// --------------------------------------------------------------------------------- //

bool SpinAsymmetry::run (std::map<std::string,std::string> const & sdata)
{
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
    double afactor = change_units(Aunits, aUnit_rad);
    
    // scattering event parameters
    int ni = Conv<int>(sdata, "ni", name());
    int li = Conv<int>(sdata, "li", name());
    int mi0= Conv<int>(sdata, "mi", name());
    int nf = Conv<int>(sdata, "nf", name());
    int lf = Conv<int>(sdata, "lf", name());
    int mf0= Conv<int>(sdata, "mf", name());
    double E = Conv<double>(sdata, "Ei", name()) * efactor;
    
    // use mi >= 0; if mi < 0, flip both signs
    int mi = (mi0 < 0 ? -mi0 : mi0);
    int mf = (mi0 < 0 ? -mf0 : mf0);
    
    // angles
    rArray angles;
    
    // get angle / angles
    try {
        
        // is there a single angle specified using command line ?
        angles.push_back(Conv<double>(sdata, "theta", name()));
        
    } catch (std::exception e) {
        
        // are there more angles specified using the STDIN ?
        angles = readStandardInput<double>();
    }
    
    // compute cross sections
    rArray scaled_angles = angles * afactor, dcs0(angles.size()), dcs1(angles.size());
    hex_differential_cross_section (ni,li,mi, nf,lf,mf, 0, 1,&E, angles.size(),scaled_angles.data(), dcs0.data(),nullptr);
    hex_differential_cross_section (ni,li,mi, nf,lf,mf, 1, 1,&E, angles.size(),scaled_angles.data(), dcs1.data(),nullptr);
    
    // compute spin asymetry
    rArray asy = (dcs0 - dcs1 / 3.) / (dcs0 + dcs1);
    
    // substitute possible "nan"-s and "inf"-s by zeros
    for (double & x : asy)
    {
        if (not std::isfinite(x))
            x = 0.;
    }
    
    // write out
    std::cout << logo("#") <<
        "# Spin asymetry for \n" <<
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi0<< ",\n" <<
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf0<< ",\n" <<
        "#     E = " << E/efactor << " " << unit_name(Eunits)
                     << " ordered by angle in " << unit_name(Aunits) << "\n";
    OutputTable table;
    table.setWidth(15, 15);
    table.setAlignment(OutputTable::left);
    table.write("# angle    ", "asymetry ");
    table.write("# ---------", "---------");
    
    for (std::size_t i = 0; i < angles.size(); i++)
        table.write(angles[i], asy[i]);
    
    return true;
}
