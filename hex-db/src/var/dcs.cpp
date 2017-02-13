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

#include <cmath>
#include <map>
#include <string>
#include <vector>

// --------------------------------------------------------------------------------- //

#include "hex-chebyshev.h"
#include "hex-interpolate.h"
#include "hex-special.h"
#include "hex-version.h"

// --------------------------------------------------------------------------------- //

#include "../quantities.h"
#include "../utils.h"

// --------------------------------------------------------------------------------- //

createNewScatteringQuantity(DifferentialCrossSection, "dcs")

// --------------------------------------------------------------------------------- //

std::string DifferentialCrossSection::description ()
{
    return "Differential cross section.";
}

std::vector<std::string> DifferentialCrossSection::dependencies ()
{
    return std::vector<std::string>
    {
        "scatamp"
    };
}

std::vector<std::pair<std::string,std::string>> DifferentialCrossSection::params ()
{
    return std::vector<std::pair<std::string,std::string>>
    {
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
}

std::vector<std::string> DifferentialCrossSection::vparams ()
{
    return std::vector<std::string>
    {
        "theta"
    };
}

// --------------------------------------------------------------------------------- //

bool DifferentialCrossSection::initialize (sqlitepp::session & db)
{
    return ScatteringQuantity::initialize(db);
}

bool DifferentialCrossSection::createTable ()
{
    return ScatteringQuantity::createTable();
}

bool DifferentialCrossSection::updateTable()
{
    return ScatteringQuantity::updateTable();
}

// --------------------------------------------------------------------------------- //

void hex_differential_cross_section
(
    int ni, int li, int mi,
    int nf, int lf, int mf,
    int S,
    int nEnergies, double * energies,
    int nAngles, double * angles,
    double * dcs, double * extra
)
{
    hex_differential_cross_section_
    (
        &ni, &li, &mi,
        &nf, &lf, &mf,
        &S,
        &nEnergies, energies,
        &nAngles, angles,
        dcs, extra
    );
}

void hex_differential_cross_section_
(
    int * ni, int * li, int * mi,
    int * nf, int * lf, int * mf,
    int * S,
    int * nEnergies, double * energies,
    int * nAngles, double * angles,
    double * dcs, double * extra
)
{
    rArrayView E (*nEnergies, energies);
    
    rArray ki = sqrt(E);
    rArray kf = sqrt(E - 1./((*ni)*(*ni)) + 1./((*nf)*(*nf)));
    
    // get scattering amplitudes
    cArray amplitudes ((*nEnergies) * (*nAngles));
    cArray xamplitudes ((*nEnergies) * (*nAngles));
    hex_scattering_amplitude_
    (
        ni, li, mi,
        nf, lf, mf,
        S,
        nEnergies, energies,
        nAngles, angles,
        reinterpret_cast<double*>(amplitudes.data()),
        extra ? reinterpret_cast<double*>(xamplitudes.data()) : nullptr
    );
    
    for (int ie = 0; ie < (*nEnergies); ie++)
    {
        for (int ia = 0; ia < (*nAngles); ia++)
        {
            // calculate differential cross section
            dcs[ie * (*nAngles) + ia] = sqrabs(amplitudes[ie * (*nAngles) + ia]) * (kf[ie]/ki[ie] * 0.25 * (2 * (*S) + 1));
            
            // calculate extrapolated cross section
            if (extra)
                extra[ie * (*nAngles) + ia] = sqrabs(xamplitudes[ie * (*nAngles) + ia]) * (kf[ie]/ki[ie] * 0.25 * (2 * (*S) + 1));
        }
    }
    
}

bool DifferentialCrossSection::run (std::map<std::string,std::string> const & sdata)
{
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
    double lfactor = change_units(lUnit_au, Lunits);
    double afactor = change_units(Aunits, aUnit_rad);
    
    // scattering event parameters
    int ni = Conv<int>(sdata, "ni", name());
    int li = Conv<int>(sdata, "li", name());
    int mi0= Conv<int>(sdata, "mi", name());
    int nf = Conv<int>(sdata, "nf", name());
    int lf = Conv<int>(sdata, "lf", name());
    int mf0= Conv<int>(sdata, "mf", name());
    int  S = Conv<int>(sdata, "S", name());
    double E = Conv<double>(sdata, "Ei", name()) * efactor;
    
    // use mi >= 0; if mi < 0, flip both signs
    int mi = (mi0 < 0 ? -mi0 : mi0);
    int mf = (mi0 < 0 ? -mf0 : mf0);
    
    // angles
    rArray angles;
    
    // get angle / angles
    try
    {
        // is there a single angle specified using command line ?
        angles.push_back(Conv<double>(sdata, "theta", name()));
    }
    catch (std::exception e)
    {
        // are there more angles specified using the STDIN ?
        angles = readStandardInput<double>();
    }
    
    // compute cross section
    rArray scaled_angles = angles * afactor, dcs (angles.size()), dcs_ex (angles.size());
    hex_differential_cross_section(ni,li,mi, nf,lf,mf, S, 1, &E, angles.size(), scaled_angles.data(), dcs.data(), dcs_ex.data());
    
    // write out
    std::cout << logo("#") <<
        "# Differential cross section in " << unit_name(Lunits) << " for \n" <<
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi0 << ",\n" <<
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf0 << ",\n" <<
        "#     S = " << S << ", E = " << E/efactor << " " << unit_name(Eunits)
                     << " ordered by angle in " << unit_name(Aunits) << "\n";
    OutputTable table;
    table.setWidth(15);
    table.setAlignment(OutputTable::left);
    table.write("# angle    ", "dcs      ", "dcs [ex] ");
    table.write("# ---------", "---------", "---------");
    for (std::size_t i = 0; i < angles.size(); i++)
    {
        table.write
        (
            angles[i],
            dcs[i] * lfactor * lfactor,
            dcs_ex[i] * lfactor * lfactor
        );
    }
    
    return true;
}
