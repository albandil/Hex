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

#include "variables.h"

// --------------------------------------------------------------------------------- //

createNewScatteringQuantity(ExtrapolatedCrossSection);

// --------------------------------------------------------------------------------- //

const std::string ExtrapolatedCrossSection::Id = "xcs";
const std::string ExtrapolatedCrossSection::Description = "Extrapolated cross section (using Aitken Δ²-process).";
const std::vector<std::pair<std::string,std::string>> ExtrapolatedCrossSection::Dependencies = {
    {"ni", "Initial atomic principal quantum number."},
    {"li", "Initial atomic orbital quantum number."},
    {"mi", "Initial atomic magnetic quantum number."},
    {"nf", "Final atomic principal quantum number."},
    {"lf", "Final atomic orbital quantum number."},
    {"mf", "Final atomic magnetic quantum number."},
    {"Ei", "Projectile impact energy (Rydberg)."}
};
const std::vector<std::string> ExtrapolatedCrossSection::VecDependencies = { "Ei" };

bool ExtrapolatedCrossSection::initialize (sqlitepp::session & db) const
{
    return true;
}

std::vector<std::string> const & ExtrapolatedCrossSection::SQL_CreateTable () const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

std::vector<std::string> const & ExtrapolatedCrossSection::SQL_Update () const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

bool ExtrapolatedCrossSection::run (std::map<std::string,std::string> const & sdata) const
{
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
    double lfactor = change_units(lUnit_au, Lunits);
    
    // scattering event parameters
    int ni = Conv<int>(sdata, "ni", Id);
    int li = Conv<int>(sdata, "li", Id);
    int mi0= Conv<int>(sdata, "mi", Id);
    int nf = Conv<int>(sdata, "nf", Id);
    int lf = Conv<int>(sdata, "lf", Id);
    int mf0= Conv<int>(sdata, "mf", Id);
    
    // use mi >= 0; if mi < 0, flip both signs
    int mi = (mi0 < 0 ? -mi0 : mi0);
    int mf = (mi0 < 0 ? -mf0 : mf0);
    
    // energies
    rArray energies;
    
    // get energy / energies
    try {
        
        // is there a single energy specified using command line ?
        energies.push_back(Conv<double>(sdata, "Ei", Id));
        
    } catch (std::exception e) {
        
        // are there more energies specified using the STDIN ?
        energies = readStandardInput<double>();
    }
    
    // SQL interface variables
    double E, sigma; int delta;
    
    // last contributions to the complete integral cross section
    // for every energy from E0
    rArray E0, dsigma0;
    
    // last but one contributions to the complete integral cross section
    // for every energy from E0
    rArray E1, dsigma1;
    
    // compose query
    sqlitepp::statement st(db);
    st << "SELECT Dat.Ei AS Ei, Lim.ell - Dat.ell AS Delta, Dat.sigma                "
          "FROM (                                                                    "
          "  SELECT Ei, L, SUM(sigma) AS sigma FROM " + IntegralCrossSection::Id + " "
          "  WHERE ni = :ni AND li = :li AND mi = :mi                                "
          "    AND nf = :nf AND lf = :lf AND mf = :mf                                "
          "  GROUP BY Ei, L                                                          "
          "  ORDER BY Ei ASC                                                         "
          ") AS Dat CROSS JOIN (                                                     "
          "  SELECT Ei, MAX(L) AS L FROM " + IntegralCrossSection::Id + "            "
          "  WHERE ni = :Ni AND li = :Li AND mi = :Mi                                "
          "    AND nf = :Nf AND lf = :Lf AND mf = :Mf                                "
          "  GROUP BY Ei                                                             "
          "  ORDER BY Ei ASC                                                         "
          ") AS Lim ON Dat.Ei = Lim.Ei                                               "
          "WHERE Dat.ell = Lim.ell OR Dat.ell = Lim.ell - 1                          "
          "ORDER BY Ei ASC, Delta DESC",
        sqlitepp::into(E), sqlitepp::into(delta), sqlitepp::into(sigma),
        sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
        sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
        sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
        sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf);
    
    // retrieve last and last but one contributions to the cross sections
    while (st.exec())
    {
        if (delta == 0)
        {
            E0.push_back(E);
            dsigma0.push_back(sigma);
        }
        else /* delta == 1 */
        {
            E1.push_back(E);
            dsigma1.push_back(sigma);
        }
    }
    
    // full complete cross section
    rArray Efull, sigmafull;
    
    // compose query
    sqlitepp::statement st1(db);
    st1 << "SELECT Ei, sigma FROM " + CompleteCrossSection::Id + " "
            "WHERE ni = :ni "
            "  AND li = :li "
            "  AND mi = :mi "
            "  AND nf = :nf "
            "  AND lf = :lf "
            "  AND mf = :mf "
            "ORDER BY Ei ASC",
        sqlitepp::into(E), sqlitepp::into(sigma),
        sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
        sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf);
    
    // retrieve data
    while (st1.exec())
    {
        Efull.push_back(E);
        sigmafull.push_back(sigma);
    }
    
    // terminate if no data
    if (Efull.empty())
        return true;
    
    // write header
    std::cout << logo("#") <<
        "# Extrapolated and complete cross section in " << unit_name(Lunits) << " for\n" <<
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n" <<
        "# ordered by energy in " << unit_name(Eunits) << "\n" <<
        "#\n" <<
        "# E\t σx\t σc\n";
    
    // threshold for ionization
    double Eion = 1./(ni*ni);
    
    // compute complete cross section
    rArray ccs = (energies[0] < 0) ? sigmafull : ((efactor * energies.front() < Eion) ? 
        interpolate_real(Efull, sigmafull, energies * efactor, gsl_interp_linear) :
        interpolate_real(Efull, sigmafull, energies * efactor, gsl_interp_cspline));
    
    // reshape arrays for compatible indexing by energies
    rArray dsigma0empty(dsigma0.size());    // zero-filled ghost
    rArray dsigma1empty(dsigma1.size());    // zero-filled ghost
    merge (Efull, sigmafull, E0, dsigma0empty);    // inflate Efull and sigmafull to accomodate E0 and dsigma0
    merge (Efull, sigmafull, E1, dsigma1empty);    // inflate Efull and sigmafull to accomodate E1 and dsigma1
    rArray sigmafullempty(sigmafull.size());    // // zero-filled ghost
    merge (E0, dsigma0, Efull, sigmafullempty);    // inflate E0 and dsigma0
    merge (E1, dsigma1, Efull, sigmafullempty);    // inflate E1 and dsigma1
    
    // extrapolate using Aitken Δ²-method
    sigmafull -= dsigma0 * dsigma0 / (dsigma0 - dsigma1);
    
    // interpolate
    rArray xcs = (energies[0] < 0) ? sigmafull : ((efactor * energies.front() < 1.) ? 
        interpolate_real(Efull, sigmafull, energies * efactor, gsl_interp_linear) :
        interpolate_real(Efull, sigmafull, energies * efactor, gsl_interp_cspline));
    
    if (energies[0] < 0.)
    {
        // negative energy indicates full output
        for (std::size_t i = 0; i < Efull.size(); i++)
            std::cout << Efull[i] / efactor << "\t" << xcs[i]*lfactor*lfactor << "\t" << ccs[i]*lfactor*lfactor << "\n";
    }
    else
    {
        // output
        for (std::size_t i = 0; i < energies.size(); i++)
            std::cout << energies[i] << "\t" << xcs[i]*lfactor*lfactor << "\t" << ccs[i]*lfactor*lfactor << "\n";
    }
    
    return true;
}
