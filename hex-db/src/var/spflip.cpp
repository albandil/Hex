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

#include <sqlitepp/sqlitepp.hpp>

#include "hex-chebyshev.h"
#include "hex-interpolate.h"
#include "hex-special.h"
#include "hex-version.h"

#include "variables.h"

const std::string SpinFlipCrossSection::Id = "spflip";
const std::string SpinFlipCrossSection::Description = "Spin flip cross section.";
const std::vector<std::pair<std::string,std::string>> SpinFlipCrossSection::Dependencies = {
    {"ni", "Initial atomic principal quantum number."},
    {"li", "Initial atomic orbital quantum number."},
    {"mi", "Initial atomic magnetic quantum number."},
    {"nf", "Final atomic principal quantum number."},
    {"lf", "Final atomic orbital quantum number."},
    {"mf", "Final atomic magnetic quantum number."},
    {"Ei", "Projectile impact energy (Rydberg)."}
};
const std::vector<std::string> SpinFlipCrossSection::VecDependencies = { "Ei" };

bool SpinFlipCrossSection::initialize (sqlitepp::session & db) const
{
    return true;
}

std::vector<std::string> const & SpinFlipCrossSection::SQL_CreateTable () const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

std::vector<std::string> const & SpinFlipCrossSection::SQL_Update () const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

bool SpinFlipCrossSection::run (std::map<std::string,std::string> const & sdata) const
{
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
    double lfactor = change_units(lUnit_au, Lunits);
    
    // scattering event parameters
    int ni = Conv<int>(sdata, "ni", Id);
    int li = Conv<int>(sdata, "li", Id);
    int mi = Conv<int>(sdata, "mi", Id);
    int nf = Conv<int>(sdata, "nf", Id);
    int lf = Conv<int>(sdata, "lf", Id);
    int mf = Conv<int>(sdata, "mf", Id);
    
    // energies and cross sections
    rArray energies;
    
    // get energy / energies
    try {
        
        // is there a single energy specified using command line ?
        energies.push_back(Conv<double>(sdata, "Ei", Id));
        
    } catch (std::exception e) {
        
        // are there more energies specified using the STDIN ?
        energies = readStandardInput<double>();
    }
    
    // energies in Rydbergs
    rArray scaled_energies = energies * efactor;
    
    // compute cross section
    double E, sigma;
    sqlitepp::statement st(db);
    st << "SELECT singlet.Ei AS Ei, SUM(singlet.ReT*singlet.ReT+singlet.ImT*singlet.ImT+triplet.ReT*triplet.ReT+triplet.ImT*triplet.ImT-2*singlet.ReT*triplet.ReT-2*singlet.ImT*triplet.ImT)/157.91367 "
              "FROM "
              "( "
                "SELECT Ei, ell, SUM(Re_T_ell) AS ReT, SUM(Im_T_ell) As ImT "
                "FROM tmat "
                "WHERE ni = :ni AND li = :li AND mi = :mi AND nf = :nf AND lf = :lf AND mf = :mf AND S = 0 "
                "GROUP BY Ei, ell, L "
                "ORDER BY Ei, ell ASC "
              ") AS singlet "
              "INNER JOIN "
              "( "
                "SELECT Ei, ell, SUM(Re_T_ell) AS ReT, SUM(Im_T_ell) AS ImT "
                "FROM tmat "
                "WHERE ni = :ni AND li = :li AND mi = :mi AND nf = :nf AND lf = :lf AND mf = :mf AND S = 1 "
                "GROUP BY Ei, ell, L "
                "ORDER BY Ei, ell ASC "
              ") AS triplet "
              "ON singlet.Ei = triplet.Ei AND singlet.ell = triplet.ell "
              "GROUP BY singlet.Ei "
              "ORDER BY singlet.Ei",
            sqlitepp::into(E), sqlitepp::into(sigma),
            sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi), sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf);
    
    // retrieve data
    rArray E_arr, sigma_arr;
    while (st.exec())
    {
        E_arr.push_back(E);
        sigma_arr.push_back(sigma);
    }
    
    // terminate if no data
    if (E_arr.empty() or energies.empty())
        return true;
    
    rArray spflip;
    double Eion = 1./(ni*ni);
    if (energies[0] < 0)
    {
        // all available
        energies = E_arr;
        spflip = sigma_arr;
    }
    else
    {
        // interpolate
        spflip = (energies[0] < Eion) ? 
            interpolate_real(E_arr, sigma_arr, energies, gsl_interp_linear) :
            interpolate_real(E_arr, sigma_arr, energies, gsl_interp_cspline);
    }
    
    // write out
    std::cout << logo("#") <<
        "# Spin flip cross section in " << unit_name(Lunits) << " for \n" <<
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n" <<
        "#     ordered by energy in " << unit_name(Eunits) << "\n";
    OutputTable table;
    table.setWidth(15);
    table.setAlignment(OutputTable::left);
    table.write("# Ei       ", "spin flip");
    table.write("# ---------", "---------");
    for (std::size_t i = 0; i < energies.size(); i++)
        table.write(energies[i], spflip[i]*lfactor*lfactor);
    
    return true;
}
