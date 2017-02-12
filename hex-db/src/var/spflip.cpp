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

createNewScatteringQuantity(SpinFlipCrossSection, "spflip")

// --------------------------------------------------------------------------------- //

extern void hex_tmat_pw_transform
(
    int ni, int li, int mi,
    int nf, int lf, int mf,
    int nEnergies,
    double * energies,
    bool extra,
    std::function
    <
        void
        (
            int ell,
            iArrays const & converged,
            iArrays const & complete,
            cArrays const & tmatrices
        )
    > f
);

// --------------------------------------------------------------------------------- //

std::string SpinFlipCrossSection::description ()
{
    return "Spin flip cross section.";
}

std::vector<std::string> SpinFlipCrossSection::dependencies ()
{
    return std::vector<std::string>
    {
        "tmat"
    };
}

std::vector<std::pair<std::string,std::string>> SpinFlipCrossSection::params ()
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

std::vector<std::string> SpinFlipCrossSection::vparams ()
{
    return std::vector<std::string>
    {
        "Ei"
    };
}

// --------------------------------------------------------------------------------- //

bool SpinFlipCrossSection::initialize (sqlitepp::session & db)
{
    return ScatteringQuantity::initialize(db);
}

bool SpinFlipCrossSection::createTable ()
{
    return ScatteringQuantity::createTable();
}

bool SpinFlipCrossSection::updateTable ()
{
    return ScatteringQuantity::updateTable();
}

// --------------------------------------------------------------------------------- //

bool SpinFlipCrossSection::run (std::map<std::string,std::string> const & sdata)
{
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
    double lfactor = change_units(lUnit_au, Lunits);
    
    // scattering event parameters
    int ni = Conv<int>(sdata, "ni", name());
    int li = Conv<int>(sdata, "li", name());
    int mi0= Conv<int>(sdata, "mi", name());
    int nf = Conv<int>(sdata, "nf", name());
    int lf = Conv<int>(sdata, "lf", name());
    int mf0= Conv<int>(sdata, "mf", name());
    
    // use mi >= 0; if mi < 0, flip both signs
    int mi = (mi0 < 0 ? -mi0 : mi0);
    int mf = (mi0 < 0 ? -mf0 : mf0);
    
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
    
    // get available energies if requested
    if (not scaled_energies.empty() and scaled_energies.front() == -1)
    {
        scaled_energies.clear();
        double E;
        sqlitepp::statement st (session());
        st << "SELECT DISTINCT Ei FROM 'tmat' "
              "WHERE ni = :ni "
              "  AND li = :li "
              "  AND mi = :mi "
              "  AND nf = :nf "
              "  AND lf = :lf "
              "  AND mf = :mf "
              "ORDER BY Ei ASC",
              sqlitepp::into(E),
              sqlitepp::use(ni), sqlitepp::use(li),  sqlitepp::use(mi),
              sqlitepp::use(nf), sqlitepp::use(lf),  sqlitepp::use(mf);
        while (st.exec())
        {
            scaled_energies.push_back(E);
        }
    }
    
    // results
    rArray spflip (scaled_energies.size()), spflip_ex (scaled_energies.size());
    
    // This function will be called by "hex_tmat_pw_transform" for every partial wave
    //    -  "ell" is the angular momentum of the partial wave
    //    -  "converged" indicates whether the partial wave expansion is converged (0 = not yet, 1 = yes, 2 = convergence failed)
    //    -  "complete" indicates whether all T-matrix components were found in the database (0 = no, 1 = yes)
    //    -  "tmatrices" are singlet and triplet T-matrix arrays at given energies
    auto fun = [&]
    (
        int ell,
        iArrays const & converged,
        iArrays const & complete,
        cArrays const & tmatrices
    )
    {
        // all energies
        for (unsigned j = 0; j < scaled_energies.size(); j++)
        if (converged[0][j] == 0 and converged[1][j] == 0)
        {
            double contrib = sqrabs(tmatrices[0][j])
                           + sqrabs(tmatrices[1][j])
                           + 2.0 * (tmatrices[0][j] * tmatrices[1][j]).real();
            
            // update scattering amplitude
            if (complete[0][j] and complete[1][j])
                spflip[j] += contrib;
            
            // update the extrapolated scattering amplitude
            spflip_ex[j] += contrib;
            
        }
    };
    
    // call the T-matrix retrieval driver
    hex_tmat_pw_transform
    (
        ni, li, mi, nf, lf, mf,
        scaled_energies.size(), scaled_energies.data(),
        true,
        fun
    );
    
    // momenta
    rArray ki = sqrt(scaled_energies);
    rArray kf = sqrt(scaled_energies - 1.0 / (ni * ni) + 1.0 / (nf * nf));
    
    // add prefactor
    spflip *= kf / ki / 157.9136704174297; // 16pi^2
    spflip_ex *= kf / ki / 157.9136704174297;
    
    // write out
    std::cout << logo("#") <<
        "# Spin flip cross section in " << unit_name(Lunits) << " for \n" <<
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi0 << ",\n" <<
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf0 << ",\n" <<
        "#     ordered by energy in " << unit_name(Eunits) << "\n";
    OutputTable table;
    table.setWidth(15);
    table.setAlignment(OutputTable::left);
    table.write("# Ei       ", "spin flip     ", "spin flip [ex]");
    table.write("# ---------", "--------------", "--------------");
    for (std::size_t i = 0; i < scaled_energies.size(); i++)
    {
        table.write
        (
            scaled_energies[i]*efactor,
            spflip[i]*lfactor*lfactor,
            spflip_ex[i]*lfactor*lfactor
        );
    }
    
    return true;
}
