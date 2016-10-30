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
        {"mi", "Initial atomic magnetic quantum number."},
        {"nf", "Final atomic principal quantum number."},
        {"lf", "Final atomic orbital quantum number."},
        {"mf", "Final atomic magnetic quantum number."},
        {"S", "Total spin of atomic + projectile electron."},
        {"Ei", "Projectile impact energy (Rydberg)."},
        {"ell", "Partial wave."}
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
//     double lfactor = change_units(lUnit_au, Lunits);
    
    // scattering event parameters
    int ni = Conv<int>(sdata, "ni", name());
    int li = Conv<int>(sdata, "li", name());
    int mi0= Conv<int>(sdata, "mi", name());
    int nf = Conv<int>(sdata, "nf", name());
    int lf = Conv<int>(sdata, "lf", name());
    int mf0= Conv<int>(sdata, "mf", name());
    int  S = Conv<int>(sdata, "S", name());
    int ell = Conv<int>(sdata, "ell", name());
    
    // use mi >= 0; if mi < 0, flip both signs
    int mi = (mi0 < 0 ? -mi0 : mi0);
    int mf = (mi0 < 0 ? -mf0 : mf0);
    
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
    
    // compose query
    sqlitepp::statement st (session());
    st << "SELECT Ei, sigma FROM 'ics' "
            "WHERE ni = :ni "
            "  AND li = :li "
            "  AND mi = :mi "
            "  AND nf = :nf "
            "  AND lf = :lf "
            "  AND mf = :mf "
            "  AND ell = :ell  "
            "  AND  S = :S  "
            "ORDER BY Ei ASC",
        sqlitepp::into(E), sqlitepp::into(sigma),
        sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
        sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
        sqlitepp::use(ell), sqlitepp::use(S);
    
    // retrieve data
    while (st.exec())
    {
        E_arr.push_back(E);
        sigma_arr.push_back(sigma);
    }
    
    // write header
    std::cout << logo("#") <<
        "# Collision strength (dimensionless) for\n"
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n" <<
        "#     ell = " << ell << ", S = " << S << "\n" <<
        "# ordered by energy in " << unit_name(Eunits) << "\n" <<
        "# \n";
    OutputTable table;
    table.setWidth(15, 15, 15);
    table.setAlignment(OutputTable::left);
    table.write("# E        ", "omega    ");
    table.write("# ---------", "---------");
    
    // terminate if no data
    if (E_arr.size() == 0)
        return true;
    
    if (energies[0] < 0.)
    {
        // negative energy indicates full output
        for (std::size_t i = 0; i < E_arr.size(); i++)
        {
            // kinematic variables to unscale the cross sections
            double ki = std::sqrt(E_arr[i]);
            double kf = std::sqrt(E_arr[i] - 1./(ni*ni) + 1./(nf*nf));
            double cs_prefactor = kf / ki * (2*S+1) / 4.;
            
            // compute collisions strengths
            double Omega = E_arr[i] * (2*li+1) * sigma_arr[i] / cs_prefactor;
            
            // print collisions strengths
            table.write(E_arr[i] / efactor, Omega);
        }
    }
    else
    {
        // threshold for ionization
        double Eion = 1./(ni*ni);
        
        // kinematic variables to unscale the cross sections
        rArray ki = sqrt(energies);
        rArray kf = sqrt(energies - 1./(ni*ni) + 1./(nf*nf));
        rArray cs_prefactor = kf / ki * (2*S+1) / 4.;
        
        // interpolate
        rArray interp = (efactor * energies.front() < Eion) ? 
            interpolate_real(E_arr, sigma_arr, energies * efactor, gsl_interp_linear) :
            interpolate_real(E_arr, sigma_arr, energies * efactor, gsl_interp_cspline);
        
        // compute collision strength
        rArray omegas = energies * efactor * (2*li+1) * interp / cs_prefactor;
        
        // output
        for (std::size_t i = 0; i < energies.size(); i++)
            table.write(energies[i], omegas[i]);
    }
    
    return true;
}
