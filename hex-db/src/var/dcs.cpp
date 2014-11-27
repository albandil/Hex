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

#include <cmath>
#include <map>
#include <string>
#include <vector>

#include <sqlitepp/sqlitepp.hpp>

#include "../chebyshev.h"
#include "../interpolate.h"
#include "../special.h"
#include "../variables.h"
#include "../version.h"

const std::string DifferentialCrossSection::Id = "dcs";
const std::string DifferentialCrossSection::Description = "Differential cross section.";
const std::vector<std::pair<std::string,std::string>> DifferentialCrossSection::Dependencies = {
    {"ni", "Initial atomic principal quantum number."},
    {"li", "Initial atomic orbital quantum number."},
    {"mi", "Initial atomic magnetic quantum number."},
    {"nf", "Final atomic principal quantum number."},
    {"lf", "Final atomic orbital quantum number."},
    {"mf", "Final atomic magnetic quantum number."},
    {"S", "Total spin of atomic + projectile electron."},
    {"Ei", "Projectile impact energy (Rydberg)."},
    {"theta", "Scattering angles for which to compute the cross section."}
};
const std::vector<std::string> DifferentialCrossSection::VecDependencies = { "theta" };

bool DifferentialCrossSection::initialize (sqlitepp::session & db) const
{
    return true;
}

std::vector<std::string> const & DifferentialCrossSection::SQL_CreateTable () const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

std::vector<std::string> const & DifferentialCrossSection::SQL_Update () const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

void hex_differential_cross_section_
(
    int * ni, int * li, int * mi,
    int * nf, int * lf, int * mf,
    int * S, double * E, int * N,
    double * angles, double * dcs
)
{
    double ki = std::sqrt((*E));
    double kf = std::sqrt((*E) - 1./((*ni)*(*ni)) + 1./((*nf)*(*nf)));
    
    // check if this is an allowed transition
    if (not std::isfinite(ki) or not std::isfinite(kf))
        return;
    
    // get scattering amplitudes
    cArray amplitudes(*N);
    hex_scattering_amplitude_
    (
        ni, li, mi, nf, lf, mf, S, E, N, angles,
        reinterpret_cast<double*>(amplitudes.data())
    );
    
    // calculate differential cross section
    rArrayView(*N, dcs) = sqrabs(amplitudes) * (kf/ki * 0.25 * (2 * (*S) + 1));
}

bool DifferentialCrossSection::run (std::map<std::string,std::string> const & sdata) const
{
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
    double lfactor = change_units(lUnit_au, Lunits);
    double afactor = change_units(Aunits, aUnit_rad);
    
    // scattering event parameters
    int ni = Conv<int>(sdata, "ni", Id);
    int li = Conv<int>(sdata, "li", Id);
    int mi = Conv<int>(sdata, "mi", Id);
    int nf = Conv<int>(sdata, "nf", Id);
    int lf = Conv<int>(sdata, "lf", Id);
    int mf = Conv<int>(sdata, "mf", Id);
    int  S = Conv<int>(sdata, "S", Id);
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
    
    // compute cross section
    rArray scaled_angles = angles * afactor, dcs(angles.size());
    hex_differential_cross_section(ni,li,mi, nf,lf,mf, S, E, angles.size(), scaled_angles.data(), dcs.data());
    
    // write out
    std::cout << logo("#") <<
        "# Differential cross section in " << unit_name(Lunits) << " for \n" <<
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n" <<
        "#     S = " << S << ", E = " << E/efactor << " " << unit_name(Eunits)
                     << " ordered by angle in " << unit_name(Aunits) << "\n";
    OutputTable table;
    table.setWidth(15);
    table.setAlignment(OutputTable::left);
    table.write("# angle    ", "dcs      ");
    table.write("# ---------", "---------");
    for (std::size_t i = 0; i < angles.size(); i++)
        table.write(angles[i], dcs[i]*lfactor*lfactor);
    
    return true;
}
