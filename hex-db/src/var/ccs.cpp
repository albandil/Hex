//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2015, Jakub Benda, Charles University in Prague                    //
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

#include <sqlite3.h>

#include "../chebyshev.h"
#include "../clenshawcurtis.h"
#include "../interpolate.h"
#include "../variables.h"
#include "../version.h"

// -------------------------------------------------------------------------- //

//
// custom function for integration of BLOB-represented Chebyshev
// expansion of the Born T-matrix
//

void db_bornICS (sqlite3_context* pdb, int n, sqlite3_value** val)
{
    // get blob data as text; reinterpret_cast is save as we are using
    // the low ASCII only
    std::string blob = reinterpret_cast<const char*>(sqlite3_value_text(*val));
    
    // convert text data to binary array
    cArray coeffs;
    coeffs.fromBlob(blob);
    
    // construct Chebyshev approximation object from the data
    Chebyshev<double,Complex> CB (coeffs, -1., 1.);
    
    // integrate
    //
    //   1
    //   ⌠
    //   ⎮
    //   ⎮ |T(cos θ)|² dcos θ
    //   ⎮
    //   ⌡
    //  -1
    //
    int tail = CB.tail(1e-10);
    auto fsqr = [&](double cosTheta) -> double { return sqrabs(CB.clenshaw(cosTheta, tail)); };
    ClenshawCurtis<decltype(fsqr),double> integrator(fsqr);
    double result = 2 * special::constant::pi * integrator.integrate(-1, 1);
    
    // use result of the integration
    sqlite3_result_double(pdb, result);
}

// -------------------------------------------------------------------------- //

const std::string CompleteCrossSection::Id = "ccs";
const std::string CompleteCrossSection::Description = "Complete cross section (L- and S-summed integral cross section).";
const std::vector<std::pair<std::string,std::string>> CompleteCrossSection::Dependencies = {
    {"ni", "Initial atomic principal quantum number."},
    {"li", "Initial atomic orbital quantum number."},
    {"mi", "Initial atomic magnetic quantum number."},
    {"nf", "Final atomic principal quantum number."},
    {"lf", "Final atomic orbital quantum number."},
    {"mf", "Final atomic magnetic quantum number."},
    {"Ei", "Projectile impact energy (Rydberg)."}
};
const std::vector<std::string> CompleteCrossSection::VecDependencies = { "Ei" };

bool CompleteCrossSection::initialize (sqlitepp::session & db) const
{
    //
    // define Gauss-Chebyshev integration of Chebyshev expansion
    //
    
    sqlite3_create_function
    (
        db.impl(),
        "borncs",
        1,              // pass single argument
        SQLITE_UTF8,
        nullptr,
        &db_bornICS,
        nullptr,
        nullptr
    );
    
    return true;
}

std::vector<std::string> const & CompleteCrossSection::SQL_CreateTable () const
{
    static const std::vector<std::string> cmd;
    
    return cmd;
}

std::vector<std::string> const & CompleteCrossSection::SQL_Update () const
{
    static std::vector<std::string> cmd;
    
    return cmd;
}

void hex_complete_cross_section_
(
    int * ni, int * li, int * mi,
    int * nf, int * lf, int * mf,
    int * N, double * energies,
    double * ccs, int * n
)
{
    double E, sigma, sigmab, sigmaB;
    rArray E_arr, sigma_arr, sigmab_arr, EB_arr, sigmaB_arr;
    
    // if there is nothing to compute, return
    if (*N == 0)
        return;
    
    // get number of all energies
    if (*energies < 0 and n != nullptr)
    {
        sqlitepp::statement st(db);
        st << "SELECT COUNT(DISTINCT Ei) FROM " + IntegralCrossSection::Id + " "
                "WHERE ni = :ni "
                "  AND li = :li "
                "  AND mi = :mi "
                "  AND nf = :nf "
                "  AND lf = :lf "
                "  AND mf = :mf ",
            sqlitepp::into(*n),
            sqlitepp::use(*ni), sqlitepp::use(*li), sqlitepp::use(*mi),
            sqlitepp::use(*nf), sqlitepp::use(*lf), sqlitepp::use(*mf);
        st.exec();
        
        // if a reallocation is needed, return
        if (*N != *n or ccs == nullptr)
            return;
    }
    
    //
    // get contributions from the partial wave expansion
    //
    
    // compose query
    sqlitepp::statement st(db);
    st << "SELECT Ei, SUM(sigma), SUM(sigmaB) FROM " + IntegralCrossSection::Id + " "
            "WHERE ni = :ni "
            "  AND li = :li "
            "  AND mi = :mi "
            "  AND nf = :nf "
            "  AND lf = :lf "
            "  AND mf = :mf "
            "GROUP BY Ei "
            "ORDER BY Ei ASC",
        sqlitepp::into(E), sqlitepp::into(sigma), sqlitepp::into(sigmab),
        sqlitepp::use(*ni), sqlitepp::use(*li), sqlitepp::use(*mi),
        sqlitepp::use(*nf), sqlitepp::use(*lf), sqlitepp::use(*mf);
    
    // retrieve data
    while (st.exec())
    {
        E_arr.push_back(E);
        sigma_arr.push_back(sigma);
        sigmab_arr.push_back(sigmab);
    }
    
    // terminate if no data
    if (E_arr.empty())
        return;
    
    //
    // get the whole Born cross section
    //
    
    // compose query
    sqlitepp::statement stb(db);
    stb << "SELECT Ei, sqrt(Ei-1./(ni*ni)+1./(nf*nf))/sqrt(Ei)*borncs(cheb)/39.478418 FROM " + BornFullTMatrix::Id + " " // 4π²
            "WHERE ni = :ni "
            "  AND li = :li "
            "  AND mi = :mi "
            "  AND nf = :nf "
            "  AND lf = :lf "
            "  AND mf = :mf "
            "ORDER BY Ei ASC",
        sqlitepp::into(E), sqlitepp::into(sigmaB),
        sqlitepp::use(*ni), sqlitepp::use(*li), sqlitepp::use(*mi),
        sqlitepp::use(*nf), sqlitepp::use(*lf), sqlitepp::use(*mf);
    
    // retrieve data
    while (stb.exec())
    {
        EB_arr.push_back(E);
        sigmaB_arr.push_back(sigmaB);
    }
    
    //
    // interpolate
    //
    
    // threshold for ionization
    double Eion = 1./((*ni)*(*ni));
    
    // negative energy indicates output of all available cross sections
    if (*energies < 0.)
    {
        // interpolate Born cross section to partial waves' energies
        rArray sigmaBorn = (*energies < Eion) ? 
            interpolate_real(EB_arr, sigmaB_arr, E_arr, gsl_interp_linear) :
            interpolate_real(EB_arr, sigmaB_arr, E_arr, gsl_interp_cspline);
        
        // correct energies
        rArrayView(*N,energies) = E_arr;
        
        // correct cross sections by Born
        rArrayView(*N,ccs) = sigmaBorn + sigma_arr - sigmab_arr;
    }
    else
    {
        // interpolate for given 'energies'
        rArray ccs0 = (*energies < Eion) ? 
            interpolate_real(E_arr, sigma_arr, rArray(*N,energies), gsl_interp_linear) :
            interpolate_real(E_arr, sigma_arr, rArray(*N,energies), gsl_interp_cspline);
        rArray ccsb = (*energies < Eion) ? 
            interpolate_real(E_arr, sigmab_arr, rArray(*N,energies), gsl_interp_linear) :
            interpolate_real(E_arr, sigmab_arr, rArray(*N,energies), gsl_interp_cspline);
        rArray ccsB = (*energies < Eion) ? 
            interpolate_real(EB_arr, sigmaB_arr, rArray(*N,energies), gsl_interp_linear) :
            interpolate_real(EB_arr, sigmaB_arr, rArray(*N,energies), gsl_interp_cspline);
        
        // corrected cross section
        rArrayView(*N,ccs) = ccsB + (ccs0 - ccsb);
    }
}

bool CompleteCrossSection::run (std::map<std::string,std::string> const & sdata) const
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
    
    // complete cross section array
    rArray ccs(energies.size());
    
    // resize arrays if all energies are requested
    if (energies[0] < 0)
    {
        int N;
        
        // get number of energies
        hex_complete_cross_section (ni, li, mi, nf, lf, mf, energies.size(), energies.data(), nullptr, &N);
        
        // resize
        scaled_energies.resize(N);
        ccs.resize(N);
        if (N > 0)
            scaled_energies[0] = -1;
    }
    
    // retrieve requested energies
    if (scaled_energies.size() > 0)
        hex_complete_cross_section (ni, li, mi, nf, lf, mf, scaled_energies.size(), scaled_energies.data(), ccs.data(), nullptr);
    
    // write header
    std::cout << logo("#") <<
        "# Complete cross section in " << unit_name(Lunits) << " for\n" <<
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n" <<
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
