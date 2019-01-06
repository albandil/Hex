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

#include "hex-chebyshev.h"
#include "hex-clenshawcurtis.h"
#include "hex-interpolate.h"
#include "hex-version.h"

// --------------------------------------------------------------------------------- //

#include "../quantities.h"
#include "../utils.h"

// --------------------------------------------------------------------------------- //

createNewScatteringQuantity(CompleteCrossSection, "ccs")

// --------------------------------------------------------------------------------- //

std::string CompleteCrossSection::description ()
{
    return "Complete cross section (L- and S-summed integral cross section).";
}

std::vector<std::string> CompleteCrossSection::dependencies ()
{
    return std::vector<std::string>
    {
        "ics"
    };
}

std::vector<std::pair<std::string,std::string>> CompleteCrossSection::params ()
{
    return std::vector<std::pair<std::string,std::string>>
    {
        {"ni", "Initial atomic principal quantum number."},
        {"li", "Initial atomic orbital quantum number."},
//         {"mi", "Initial atomic magnetic quantum number."},
        {"nf", "Final atomic principal quantum number."},
        {"lf", "Final atomic orbital quantum number."},
//         {"mf", "Final atomic magnetic quantum number."},
        {"Ei", "Projectile impact energy (Rydberg)."}
    };
}

std::vector<std::string> CompleteCrossSection::vparams ()
{
    return std::vector<std::string>
    {
        "Ei"
    };
}

// --------------------------------------------------------------------------------- //

bool CompleteCrossSection::initialize (sqlitepp::session & db)
{
    return ScatteringQuantity::initialize(db);
}

bool CompleteCrossSection::createTable ()
{
    return ScatteringQuantity::createTable();
}

bool CompleteCrossSection::updateTable ()
{
    return ScatteringQuantity::updateTable();
}

// --------------------------------------------------------------------------------- //

void hex_complete_cross_section
(
    int ni, int li, int mi,
    int nf, int lf, int mf,
    int N, double * energies,
    double * ccs, int * Nall,
    int npws, int * pws
)
{
    hex_complete_cross_section_
    (
        &ni, &li, &mi,
        &nf, &lf, &mf,
        &N, energies, ccs, Nall,
        &npws, pws
    );
}

void hex_complete_cross_section_
(
    int * ni, int * li, int * pmi,
    int * nf, int * lf, int * pmf,
    int * N, double * energies,
    double * ccs, int * n,
    int * npws, int * pws
)
{
    CompleteCrossSection * CCS = dynamic_cast<CompleteCrossSection*>(get_quantity("ccs"));

    double E, sigma;
    rArray E_arr, sigma_arr;

    // use mi >= 0; if mi < 0, flip both signs
    int mi = ((*pmi) < 0 ? -(*pmi) : (*pmi));
    int mf = ((*pmi) < 0 ? -(*pmf) : (*pmf));

    // if there is nothing to compute, return
    if (*N == 0)
        return;

    // compose the list of wanted partial waves
    std::string pwlimit;
    if (npws != nullptr and pws != nullptr)
    {
        std::ostringstream oss;
        oss << iArrayView(*npws, pws);
        pwlimit = oss.str();
        pwlimit.front() = '(';
        pwlimit.back() = ')';
        pwlimit = " AND ell IN " + pwlimit;
    }

    // get number of all energies
    if (*energies < 0 and n != nullptr)
    {
        sqlitepp::statement st (CCS->session());
        st << "SELECT COUNT(DISTINCT Ei) FROM 'ics' "
                "WHERE ni = :ni "
                "  AND li = :li "
                "  AND mi = :mi "
                "  AND nf = :nf "
                "  AND lf = :lf "
                "  AND mf = :mf " + pwlimit,
            sqlitepp::into(*n),
            sqlitepp::use(*ni), sqlitepp::use(*li), sqlitepp::use(mi),
            sqlitepp::use(*nf), sqlitepp::use(*lf), sqlitepp::use(mf);
        st.exec();

        // if a reallocation is needed, return
        if (*N != *n)
            return;

        // fill the available energies
        if (*N >= *n)
        {
            double E;
            sqlitepp::statement st (CCS->session());
            st << "SELECT DISTINCT Ei FROM 'ics' "
                    "WHERE ni = :ni "
                    "  AND li = :li "
                    "  AND mi = :mi "
                    "  AND nf = :nf "
                    "  AND lf = :lf "
                    "  AND mf = :mf " + pwlimit +
                    "ORDER BY Ei ASC",
                sqlitepp::into(E),
                sqlitepp::use(*ni), sqlitepp::use(*li), sqlitepp::use(mi),
                sqlitepp::use(*nf), sqlitepp::use(*lf), sqlitepp::use(mf);
            for (int i = 0; i < *N and st.exec(); i++)
                *(energies + i) = E;
        }

        // if no space for cross section has been given, return
        if (ccs == nullptr)
            return;
    }

    // get number of partial waves stored in the database
    int max_available_ell;
    sqlitepp::statement ell_st (CCS->session());
    ell_st << "SELECT MAX(ell) FROM 'ics' "
              "WHERE ni = :ni AND li = :li AND mi = :mi"
              "  AND nf = :nf AND lf = :lf AND mf = :mf",
        sqlitepp::into(max_available_ell),
        sqlitepp::use(*ni), sqlitepp::use(*li), sqlitepp::use(mi),
        sqlitepp::use(*nf), sqlitepp::use(*lf), sqlitepp::use(mf);
    ell_st.exec();

    // retrieve cross sections for all requested (and available) partial waves
    rArrays E_data, sigma_data;
    for (int ell = 0; ell <= max_available_ell; ell++)
    {
        // compose query
        sqlitepp::statement st (CCS->session());
        st << "SELECT Ei, SUM(sigma) FROM 'ics' "
              "WHERE ni = :ni AND li = :li AND mi = :mi AND nf = :nf AND lf = :lf AND mf = :mf AND ell = :ell "
              "GROUP BY Ei "
              "ORDER BY Ei ASC",
            sqlitepp::into(E), sqlitepp::into(sigma),
            sqlitepp::use(*ni), sqlitepp::use(*li), sqlitepp::use(mi),
            sqlitepp::use(*nf), sqlitepp::use(*lf), sqlitepp::use(mf), sqlitepp::use(ell);

        // prepare storage
        E_data.push_back(rArray());
        sigma_data.push_back(rArray());

        // retrieve data (if requested)
        if (pws == nullptr or std::find(pws, pws + *npws, ell) != pws + *npws)
        {
            while (st.exec())
            {
                E_data.back().push_back(E);
                sigma_data.back().push_back(sigma);
            }
        }
    }

    // merge energies of all datasets
    E_arr = join(E_data);
    std::sort(E_arr.begin(), E_arr.end());
    double* ptr = std::unique(E_arr.begin(), E_arr.end());
    E_arr.resize(ptr - E_arr.begin());

    // terminate if no data
    if (E_arr.empty())
        return;

    //
    // interpolate
    //

    // negative energy indicates output of all available cross sections
    if (*energies < 0.)
    {
        rArrayView(*N,energies) = E_arr;
        rArrayView(*N,ccs).fill(0);
        for (int ell = 0; ell <= max_available_ell; ell++)
        if (pws == nullptr or std::find(pws, pws + *npws, ell) != pws + *npws)
        {
            rArrayView(*N,ccs) += interpolate_real(E_data[ell], E_data[ell] * sigma_data[ell], E_arr, gsl_interp_akima) / E_arr;
        }
    }
    else
    {
        rArrayView(*N,ccs).fill(0);
        for (int ell = 0; ell <= max_available_ell; ell++)
        if (pws == nullptr or std::find(pws, pws + *npws, ell) != pws + *npws)
        {
            rArrayView(*N,ccs) += interpolate_real(E_data[ell], E_data[ell] * sigma_data[ell], rArrayView(*N,energies), gsl_interp_akima) / rArrayView(*N,energies);
        }
    }

    // erase data for energies below threshold
    for (int i = 0; i < *N; i++)
    {
        if (energies[i] < 1./(*ni * *ni) - 1./(*nf * *nf))
            ccs[i] = 0;
    }
}

// --------------------------------------------------------------------------------- //

bool CompleteCrossSection::run (std::map<std::string,std::string> const & sdata)
{
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
    double lfactor = change_units(lUnit_au, Lunits);

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
        "# Complete cross section in " << unit_name(Lunits) << " for\n" <<
        "#     ni = " << ni << ", li = " << li << ", mi = " << (all_mi ? "*" : std::to_string(mi0)) << ",\n" <<
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << (all_mi ? "*" : std::to_string(mi0)) << ",\n";
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
    table.write("# E        ", "sigma    ");
    table.write("# ---------", "---------");

    if (scaled_energies.size() == 0)
        std::cout << "#\n# The database contains no relevant data." << std::endl;

    // write data
    for (std::size_t i = 0; i < scaled_energies.size(); i++)
        table.write(scaled_energies[i]/efactor, (std::isfinite(ccs[i]) ? ccs[i] * lfactor * lfactor : 0.));

    return true;
}
