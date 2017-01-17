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
    
    // energies and cross sections
    double E, sigma;
    rArray energies, E_arr, sigma_arr;
    
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
    
    //
    // load partial wave constributions
    //
    
    // compose query
    sqlitepp::statement st (session());
    st << "SELECT Ei, sum(sigma) FROM 'ics' "
            "WHERE ni = :ni "
            "  AND li = :li "
            "  AND mi = :mi "
            "GROUP BY Ei "
            "ORDER BY Ei ASC",
        sqlitepp::into(E), sqlitepp::into(sigma),
        sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi);
    
    // retrieve data
    while (st.exec())
    {
        E_arr.push_back(E);
        sigma_arr.push_back(sigma);
    }
    
    //
    // compute and write output
    //
    
    // write header
    std::cout << logo("#") <<
        "# Total cross section in " << unit_name(Lunits) << " for\n"
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << "\n" <<
        "# ordered by energy in " << unit_name(Eunits) << "\n" <<
        "# \n";
    OutputTable table;
    table.setWidth(15);
    table.setAlignment(OutputTable::left);
    table.write("# E        ", "sigma    ");
    table.write("# ---------", "---------");
    
    // terminate if no data
    if (E_arr.empty())
    {
        std::cout << "No data for this transition." << std::endl;
        return true;
    }
    
    // threshold for ionization
    double Eion = 1./(ni*ni);
    
    // negative energy indicates output of all available cross sections
    if (energies[0] < 0.)
    {
        // output corrected cross section
        for (std::size_t i = 0; i < E_arr.size(); i++)
            table.write(E_arr[i] / efactor, sigma_arr[i] * lfactor * lfactor);
    }
    else
    {
        // interpolate for given 'energies'
        rArray tcs = (efactor * energies.front() < Eion) ? 
            interpolate_real(E_arr, sigma_arr, energies * efactor, gsl_interp_linear) :
            interpolate_real(E_arr, sigma_arr, energies * efactor, gsl_interp_cspline);
        
        // output
        for (std::size_t i = 0; i < energies.size(); i++)
            table.write(energies[i], (std::isfinite(tcs[i]) ? tcs[i] * lfactor * lfactor : 0.));
    }
    
    return true;
}
