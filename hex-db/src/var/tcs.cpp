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
#include "hex-special.h"
#include "hex-version.h"

// --------------------------------------------------------------------------------- //

#include "../interfaces.h"
#include "../quantities.h"
#include "../utils.h"

// --------------------------------------------------------------------------------- //

createNewScatteringQuantity(TotalCrossSection, "tcs")

// --------------------------------------------------------------------------------- //

std::string TotalCrossSection::description ()
{
    return "Total cross section.";
}

std::vector<std::string> TotalCrossSection::dependencies ()
{
    return std::vector<std::string>
    {
        "ics"
    };
}

std::vector<std::pair<std::string,std::string>> TotalCrossSection::params ()
{
    return std::vector<std::pair<std::string,std::string>>
    {
        {"ni", "Initial atomic principal quantum number."},
        {"li", "Initial atomic orbital quantum number."},
        {"mi", "Initial atomic magnetic quantum number."},
        {"nf", "Final atomic principal quantum number."},
        {"lf", "Final atomic orbital quantum number."},
        {"mf", "Final atomic magnetic quantum number."},
        {"Ei", "Projectile impact energy (Rydberg)."}
    };
}

std::vector<std::string> TotalCrossSection::vparams ()
{
    return std::vector<std::string>
    {
        "Ei"
    };
}

// --------------------------------------------------------------------------------- //

bool TotalCrossSection::initialize (sqlitepp::session & db)
{
    return ScatteringQuantity::initialize(db);
}

bool TotalCrossSection::createTable ()
{
    return ScatteringQuantity::createTable();
}

bool TotalCrossSection::updateTable ()
{
    return ScatteringQuantity::updateTable();
}

// --------------------------------------------------------------------------------- //

bool TotalCrossSection::run (std::map<std::string,std::string> const & sdata)
{
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
    double lfactor = change_units(lUnit_au, Lunits);

    // scattering event parameters
    int ni = Conv<int>(sdata, "ni", name());
    int li = Conv<int>(sdata, "li", name());
    int mi = Conv<int>(sdata, "mi", name());

    if (mi < 0)
        mi = -mi;

    // energies and cross sections
    double E;
    rArray energies, sigma_arr;

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

    // get all impact energies
    if (not energies.empty() and energies.front() == -1)
    {
        energies.clear();
        sqlitepp::statement st (session());
        st << "SELECT DISTINCT Ei FROM ics WHERE ni = :ni AND li = :li AND mi = :mi",
            sqlitepp::into(E), sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi);
        while (st.exec())
        {
            energies.push_back(E);
        }
    }

    // get all final states
    iArray nfs, lfs, mfs;
    int Nf = 0, nf, lf, mf;
    sqlitepp::statement st(session());
    st << "SELECT DISTINCT nf, lf, mf FROM ics WHERE ni = :ni AND li = :li AND mi = :mi",
        sqlitepp::into(nf), sqlitepp::into(lf), sqlitepp::into(mf),
        sqlitepp::use(ni),  sqlitepp::use(li),  sqlitepp::use(mi);
    while (st.exec())
    {
        nfs.push_back(nf);
        lfs.push_back(lf);
        mfs.push_back(mf);
        Nf++;
    }

    // sum integral cross sections for all transitions
    if (not energies.empty())
    {
        rArray cs (energies.size());

        sigma_arr.resize(energies.size());

        for (int i = 0; i < Nf; i++)
        {
            hex_complete_cross_section
            (
                ni, li, mi,
                nfs[i], lfs[i], mfs[i],
                energies.size(), &energies[0], &cs[0],
                nullptr, -1, nullptr
            );

            sigma_arr += cs;
        }
    }

    // calculate total cross section from the optical theorem
    cArray singlet (energies.size()), triplet (energies.size());
    double zero = 0;
    hex_scattering_amplitude(ni, li, mi, ni, li, mi, 0, energies.size(), &energies[0], 1, &zero, reinterpret_cast<double*>(&singlet[0]), nullptr);
    hex_scattering_amplitude(ni, li, mi, ni, li, mi, 1, energies.size(), &energies[0], 1, &zero, reinterpret_cast<double*>(&triplet[0]), nullptr);
    rArray csopt = -4 * special::constant::pi * (0.25 * imagpart(singlet) + 0.75 * imagpart(triplet)) / sqrt(energies);

    // write header
    std::cout << logo("#") <<
        "# Total cross section in " << unit_name(Lunits) << " for\n"
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << "\n" <<
        "# ordered by energy in " << unit_name(Eunits) << "\n" <<
        "# \n";
    OutputTable table;
    table.setWidth(15);
    table.setAlignment(OutputTable::left);
    table.write("# E        ", "sigma (sum)", "sigma (opt)");
    table.write("# ---------", "-----------", "-----------");

    // terminate if no data
    if (energies.empty())
    {
        std::cout << "No data for this transition." << std::endl;
        return true;
    }

    // output corrected cross section
    for (std::size_t i = 0; i < energies.size(); i++)
    {
        table.write
        (
            energies[i] / efactor,
            sigma_arr[i] * lfactor * lfactor,
            csopt[i] * lfactor * lfactor
        );
    }

    return true;
}
