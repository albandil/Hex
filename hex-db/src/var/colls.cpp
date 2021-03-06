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
#include <map>
#include <string>
#include <vector>

// --------------------------------------------------------------------------------- //

#include "hex-interpolate.h"
#include "hex-version.h"

// --------------------------------------------------------------------------------- //

#include "../quantities.h"
#include "../utils.h"

// --------------------------------------------------------------------------------- //

createNewScatteringQuantity(CollisionStrength, "colls")

// --------------------------------------------------------------------------------- //

std::string CollisionStrength::description ()
{
    return "Collision strength (energy scaled integral cross section).";
}

std::vector<std::string> CollisionStrength::dependencies ()
{
    return std::vector<std::string>
    {
        "ics"
    };
}

std::vector<std::pair<std::string,std::string>> CollisionStrength::params ()
{
    return std::vector<std::pair<std::string,std::string>>
    {
        {"ni", "Initial atomic principal quantum number."},
        {"li", "Initial atomic orbital quantum number."},
//         {"mi", "Initial atomic magnetic quantum number."},
        {"nf", "Final atomic principal quantum number."},
        {"lf", "Final atomic orbital quantum number."},
//         {"mf", "Final atomic magnetic quantum number."},
        {"S", "Total spin of atomic + projectile electron."},
        {"Ei", "Projectile impact energy (Rydberg)."}
    };
}

std::vector<std::string> CollisionStrength::vparams ()
{
    return std::vector<std::string>
    {
        "Ei"
    };
}

// --------------------------------------------------------------------------------- //

bool CollisionStrength::initialize (sqlitepp::session & db)
{
    return ScatteringQuantity::initialize(db);
}

bool CollisionStrength::createTable ()
{
    return ScatteringQuantity::createTable();
}

bool CollisionStrength::updateTable ()
{
    return ScatteringQuantity::updateTable();
}

// --------------------------------------------------------------------------------- //

bool CollisionStrength::run (std::map<std::string,std::string> const & sdata)
{
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);

    // scattering event parameters
    int ni = Conv<int>(sdata, "ni", name());
    int li = Conv<int>(sdata, "li", name());
    int nf = Conv<int>(sdata, "nf", name());
    int lf = Conv<int>(sdata, "lf", name());

    // initial magnetic quantum number
    int mi0= 0; bool all_mi = false;
    if (sdata.find("mi") == sdata.end() or sdata.at("mi") == "*")
        all_mi = true;
    else
        mi0= Conv<int>(sdata, "mi", name());

    // final magnetic quantum number
    int mf0= 0; bool all_mf = false;
    if (sdata.find("mf") == sdata.end() or sdata.at("mf") == "*")
        all_mf = true;
    else
        mf0= Conv<int>(sdata, "mf", name());

    // use mi >= 0; if mi < 0, flip both signs
    int mi = (mi0 < 0 ? -mi0 : mi0);
    int mf = (mi0 < 0 ? -mf0 : mf0);

    // check if we have the optional pw list
    iArray pws;
    if (sdata.find("pws") != sdata.end())
        pws = Conv<iArray>(sdata, "pws", name());

    // energies and cross sections
    rArray energies;

    // get energy / energies
    try
    {
        // is there a single energy specified using command line ?
        energies.push_back(Conv<double>(sdata, "Ei", name()));
    }
    catch (std::exception e)
    {
        // are there more energies specified using the STDIN ?
        energies = readStandardInput<double>();
    }

    // energies in Rydbergs
    rArray scaled_energies = energies * efactor;

    // complete cross section array
    rArray ccs (energies.size()), sum_ccs (energies.size());

    // resize arrays if all energies are requested
    if (energies[0] < 0)
    {
        int N;

        // get number of energies
        hex_complete_cross_section(ni, li, mi, nf, lf, mf, energies.size(), energies.data(), nullptr, &N, pws.size(), pws.empty() ? nullptr : pws.data());

        // resize
        scaled_energies.resize(N);
        ccs.resize(N);
        if (N > 0)
            scaled_energies[0] = -1;

        // get list of energies
        hex_complete_cross_section(ni, li, mi, nf, lf, mf, scaled_energies.size(), scaled_energies.data(), nullptr, &N, pws.size(), pws.empty() ? nullptr : pws.data());
    }

    // retrieve requested energies
    if (scaled_energies.size() > 0)
    {
        // sum all cross sections from the initial state to all final magnetic sub-levels
        sum_ccs.resize(ccs.size());
        for (int mi0 = -li; mi0 <= li; mi0++)
        for (int mf0 = -lf; mf0 <= lf; mf0++)
        {
            if ((all_mi or mi0 == mi) and (all_mf or mf0 == mf))
                hex_complete_cross_section(ni, li, mi0, nf, lf, mf0, scaled_energies.size(), scaled_energies.data(), ccs.data(), nullptr, pws.size(), pws.empty() ? nullptr : pws.data());
            sum_ccs += ccs;
            ccs.fill(0);
        }
        ccs = sum_ccs;

        // average initial state
        if (all_mi)
            ccs /= (2. * li + 1.);
    }

    // write header
    std::cout << logo("#") <<
        "# Collision strengths in " << unit_name(Lunits) << " for\n" <<
        "#     ni = " << ni << ", li = " << li << ", mi = " << (all_mi ? "*" : std::to_string(mi0)) << ",\n" <<
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << (all_mf ? "*" : std::to_string(mf0)) << ",\n";
    if (not pws.empty())
    {
        std::cout <<
        "#     limited to partial waves ell = " << pws << "\n";
    }
        std::cout <<
        "# ordered by energy in " << unit_name(Eunits) << "\n" <<
        "# \n";
    OutputTable table;
    table.setWidth(15, 15);
    table.setAlignment(OutputTable::left);
    table.write("# E        ", "Omega    ");
    table.write("# ---------", "---------");

    if (scaled_energies.size() == 0)
        std::cout << "#\n# The database contains no relevant data." << std::endl;

    // write data
    for (std::size_t i = 0; i < scaled_energies.size(); i++)
    {
        table.write
        (
            scaled_energies[i] / efactor,
            std::isfinite(ccs[i]) ? scaled_energies[i] * ccs[i] : 0.
        );
    }

    return true;
}
