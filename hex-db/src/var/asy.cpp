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

#include <map>
#include <string>
#include <vector>

#include "../interfaces.h"
#include "../interpolate.h"
#include "../special.h"
#include "../variables.h"
#include "../version.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


const std::string SpinAsymmetry::Id = "asy";
const std::string SpinAsymmetry::Description = "Spin asymetry.";
const std::vector<std::pair<std::string,std::string>> SpinAsymmetry::Dependencies = {
    {"ni", "Initial atomic principal quantum number."},
    {"li", "Initial atomic orbital quantum number."},
    {"mi", "Initial atomic magnetic quantum number."},
    {"nf", "Final atomic principal quantum number."},
    {"lf", "Final atomic orbital quantum number."},
    {"mf", "Final atomic magnetic quantum number."},
    {"Ei", "Projectile impact energy (Rydberg)."},
    {"theta", "Scattering angles for which to compute the spin asymetry."}
};
const std::vector<std::string> SpinAsymmetry::VecDependencies = { "theta" };

bool SpinAsymmetry::initialize (sqlitepp::session & db) const
{
    return true;
}

std::vector<std::string> const & SpinAsymmetry::SQL_CreateTable () const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

std::vector<std::string> const & SpinAsymmetry::SQL_Update () const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

bool SpinAsymmetry::run (std::map<std::string,std::string> const & sdata) const
{
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
    double afactor = change_units(Aunits, aUnit_rad);
    
    // scattering event parameters
    int ni = Conv<int>(sdata, "ni", Id);
    int li = Conv<int>(sdata, "li", Id);
    int mi = Conv<int>(sdata, "mi", Id);
    int nf = Conv<int>(sdata, "nf", Id);
    int lf = Conv<int>(sdata, "lf", Id);
    int mf = Conv<int>(sdata, "mf", Id);
    double E = Conv<double>(sdata, "Ei", Id) * efactor;
    
    // angles
    rArray angles;
    
    // get angle / angles
    try {
        
        // is there a single angle specified using command line ?
        angles.push_back(Conv<double>(sdata, "theta", Id));
        
    } catch (std::exception e) {
        
        // are there more angles specified using the STDIN ?
        angles = readStandardInput<double>();
    }
    
    // compute cross sections
    rArray scaled_angles = angles * afactor, dcs0(angles.size()), dcs1(angles.size());
    hex_differential_cross_section (ni,li,mi, nf,lf,mf, 0, E, angles.size(), scaled_angles.data(), dcs0.data());
    hex_differential_cross_section (ni,li,mi, nf,lf,mf, 1, E, angles.size(), scaled_angles.data(), dcs1.data());
    
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
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n" <<
        "#     E = " << E/efactor << " " << unit_name(Eunits)
                     << " ordered by angle in " << unit_name(Aunits) << "\n" <<
        "# θ\t dσ\n";
    for (std::size_t i = 0; i < angles.size(); i++)
        std::cout << angles[i] << "\t" << asy[i] << "\t" << dcs0[i] << "\t" << dcs1[i] << "\n";
    
    return true;
}
