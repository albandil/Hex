//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2017, Jakub Benda, Charles University in Prague                    //
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

#include "hex-chebyshev.h"
#include "hex-version.h"

// --------------------------------------------------------------------------------- //

#include "../quantities.h"
#include "../utils.h"

// --------------------------------------------------------------------------------- //

createNewScatteringQuantity(SingleDifferentialCrossSectionWrtEnergyShare, "sdcsE")

// --------------------------------------------------------------------------------- //

std::string SingleDifferentialCrossSectionWrtEnergyShare::description ()
{
    return "Single differential ionization cross section with respect to the energy share.";
}

std::vector<std::string> SingleDifferentialCrossSectionWrtEnergyShare::dependencies ()
{
    return std::vector<std::string>
    {
        "ionf"
    };
}

std::vector<std::pair<std::string,std::string>> SingleDifferentialCrossSectionWrtEnergyShare::params ()
{
    return std::vector<std::pair<std::string,std::string>>
    {
        {"ni", "Initial atomic principal quantum number."},
        {"li", "Initial atomic orbital quantum number."},
        {"mi", "Initial atomic magnetic quantum number."},
        {"S", "Total spin of atomic + projectile electron."},
        {"Ei", "Impact energy in Ry."},
        {"Eshare", "Energy share (= min(E1,E2)/(E1 + E2))."}
    };
}

std::vector<std::string> SingleDifferentialCrossSectionWrtEnergyShare::vparams ()
{
    return std::vector<std::string>
    {
        "Eshare"
    };
}

// --------------------------------------------------------------------------------- //

bool SingleDifferentialCrossSectionWrtEnergyShare::initialize (sqlitepp::session & db)
{
    return ScatteringQuantity::initialize(db);
}

bool SingleDifferentialCrossSectionWrtEnergyShare::createTable ()
{
    return true;
}

bool SingleDifferentialCrossSectionWrtEnergyShare::updateTable ()
{
    return ScatteringQuantity::updateTable();
}

// --------------------------------------------------------------------------------- //

bool SingleDifferentialCrossSectionWrtEnergyShare::run (std::map<std::string,std::string> const & sdata)
{
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
    double lfactor = change_units(lUnit_au, Lunits);

    // scattering event parameters
    int ni = Conv<int>(sdata, "ni", name());
    int li = Conv<int>(sdata, "li", name());
    int mi = Conv<int>(sdata, "mi", name());
    int  S = Conv<int>(sdata, "S",  name());
    double Ei = Conv<double>(sdata, "Ei", name()) * efactor;

    // energy shares and cross sections
    int l1, l2, L;
    std::string blob;
    std::vector<Chebyshev<double,Complex>> CB;
    rArray energy_shares;

    // get energy / energies
    try
    {
        // is there a single energy specified using command line ?
        energy_shares.push_back(Conv<double>(sdata, "Eshare", name()));
    }
    catch (std::exception e)
    {
        // are there more energies specified using the STDIN ?
        energy_shares = readStandardInput<double>();
    }

    // write header
    std::cout << logo("#") <<
        "# Single differential cross section (wrt energy share) in " << unit_name(Lunits) << " for\n" <<
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
        "#     S = " << S << ", Ei = " << Ei << " " << unit_name(Eunits) << "\n" <<
        "# ordered by energy share\n" <<
        "#\n";
    OutputTable table;
    table.setWidth(15);
    table.setAlignment(OutputTable::left);
    table.write("# Eshare   ", "dsigma   ");
    table.write("# ---------", "---------");

    // compose query
    sqlitepp::statement st (session());
    st << "SELECT L,l1,l2,QUOTE(cheb) FROM 'ionf' "
            "WHERE ni = :ni "
            "  AND li = :li "
            "  AND mi = :mi "
            "  AND  S = :S  "
            "  AND Ei = :Ei "
            "ORDER BY Ei ASC",
        sqlitepp::into(L), sqlitepp::into(l1), sqlitepp::into(l2), sqlitepp::into(blob),
        sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
        sqlitepp::use(S), sqlitepp::use(Ei);

    // retrieve data (terminate if no data)
    while (st.exec())
    {
        if (not blob.empty())
        {
            // translate blob to Chebyshev coefficients
            cArray coeffs;
            coeffs.fromBlob(blob);
            CB.emplace_back(coeffs, 0, 1);
        }
    }

    // process all energy shares
    for (double x : energy_shares)
    {
        // Given the definition of ionf in hex-ecs, the sdcs is
        //
        //     1
        //    ---- |f(k₁,k₂)|²
        //    k₁k₂
        //
        // summed over all partial waves.
        double dsigma = 0;

        // translate energy share to hyper-angle and momenta
        double Etot = Ei - 1./(ni*ni);
        double kmax = std::sqrt(Etot);
        double k1 = std::sqrt(Etot * x);
        double k2 = std::sqrt(Etot * (1 - x));
        double alpha = k1/kmax;

        // process all partial wave contributions
        for (Chebyshev<double,Complex> const & cb : CB)
            dsigma += sqrabs(cb.clenshaw(alpha, cb.tail(1e-10)));

        // write line to table
        table.write(x, dsigma * lfactor * lfactor / (k1 * k2));
    }

    return true;
}
