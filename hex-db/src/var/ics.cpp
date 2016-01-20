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

#include <map>
#include <string>
#include <vector>

#include <sqlite3.h>

#include "hex-chebyshev.h"
#include "hex-clenshawcurtis.h"
#include "hex-interpolate.h"
#include "hex-version.h"

#include "variables.h"

// -------------------------------------------------------------------------- //

//
// custom function for evaluation of square root within a SQL statement
//

void db_sqrt (sqlite3_context* pdb, int n, sqlite3_value** val)
{
    assert(n == 1);
    sqlite3_result_double(pdb, std::sqrt(sqlite3_value_double(*val)));
}

//
// custom function for integration of BLOB-represented Chebyshev
// expansion of the ionization amplitude
//

void db_ioncs (sqlite3_context* pdb, int n, sqlite3_value** val)
{
    assert(n == 1);
    
    // get blob data as text; reinterpret_cast is save as we are using
    // the low ASCII only
    std::string blob = reinterpret_cast<const char*>(sqlite3_value_text(*val));
    
    // convert text data to binary array
    cArray coeffs;
    coeffs.fromBlob(blob);
    
    // construct Chebyshev approximation object from the data
    Chebyshev<double,Complex> CB(coeffs, 0, 1);
    
    // integrate
    //
    // 1/√2                π/4
    //  ⌠                   ⌠
    //  ⎮            dκ     ⎮
    //  ⎮ |f(κ)|² ------- = ⎮ |f(sin β)|² dβ
    //  ⎮         √(1-κ²)   ⎮
    //  ⌡                   ⌡
    //  0                   0
    //
    int tail = CB.tail(1e-10);
    auto fsqr = [&](double beta) -> double { return sqrabs(CB.clenshaw(std::sin(beta), tail)); };
    ClenshawCurtis<decltype(fsqr),double> integrator(fsqr);
    double result = integrator.integrate(0, special::constant::pi_quart);
    
    // use result of the integration
    sqlite3_result_double(pdb, result);
}

//
// custom function for linear interpolation of two points
//

void db_interpolate (sqlite3_context* pdb, int n, sqlite3_value** val)
{
    assert(n == 5);
    
    // extract parameters
    double x1 = sqlite3_value_double(*(val + 0));
    double y1 = sqlite3_value_double(*(val + 1));
    double x2 = sqlite3_value_double(*(val + 2));
    double y2 = sqlite3_value_double(*(val + 3));
    double x  = sqlite3_value_double(*(val + 4));
    
    // handle degenerate case
    if (x1 == x2)
    {
        if (x == x1)
            sqlite3_result_double(pdb, y1);
        else
            sqlite3_result_double(pdb, special::constant::Nan);
    }
    
    // interpolate
    else
    {
        sqlite3_result_double(pdb, ((x - x1) * y2 + (x2 - x) * y1) / (x2 - x1));
    }
}

// -------------------------------------------------------------------------- //

const std::string IntegralCrossSection::Id = "ics";
const std::string IntegralCrossSection::Description = "Integral cross section.";
const std::vector<std::pair<std::string,std::string>> IntegralCrossSection::Dependencies = {
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
const std::vector<std::string> IntegralCrossSection::VecDependencies = { "Ei" };

bool IntegralCrossSection::initialize (sqlitepp::session & db) const
{
    // define square root function
    sqlite3_create_function(db.impl(), "SQRT", 1, SQLITE_UTF8, nullptr, &db_sqrt, nullptr, nullptr);
    
    // define Gauss-Chebyshev integration of squared Chebyshev expansion
    sqlite3_create_function(db.impl(), "IONCS", 1, SQLITE_UTF8, nullptr, &db_ioncs, nullptr, nullptr);
    
    // define linear interpolation of two data points (5 arguments: x1, y1, x2, y2, x)
    sqlite3_create_function(db.impl(), "INTERPOLATE", 5, SQLITE_UTF8, nullptr, &db_interpolate, nullptr, nullptr);
    
    return true;
}

std::vector<std::string> const & IntegralCrossSection::SQL_CreateTable () const
{
    static const std::vector<std::string> cmd = {
        "CREATE TABLE IF NOT EXISTS '" + IntegralCrossSection::Id + "' "
        "("
            "ni  INTEGER, "
            "li  INTEGER, "
            "mi  INTEGER, "
            "nf  INTEGER, "
            "lf  INTEGER, "
            "mf  INTEGER, "
            "S   INTEGER, "
            "Ei  DOUBLE PRECISION, "
            "ell INTEGER, "
            "sigma DOUBLE PRECISION, "
            "PRIMARY KEY (ni,li,mi,nf,lf,mf,S,Ei,ell)"
        ")"
    };
    
    return cmd;
}

std::vector<std::string> const & IntegralCrossSection::SQL_Update () const
{
    static const std::vector<std::string> cmd = {
        
        // Avoid disk operations.
        //   The SQLite program does normally all operations using a backou journal disk file, which
        // slows down the calculation of the cross sections.
        
        "PRAGMA temp_store = MEMORY",
        
        // Temporary table with merged available energies from all total angular momenta for every discrete transition and partial wave.
        
        "CREATE TEMP TABLE T AS SELECT DISTINCT ni,li,mi,nf,lf,mf,S,Ei,ell FROM tmat",
        
        // Temporary table where every transition, spin, total L, partial wave and required energy from 'T' contains also the highest lower-or-equal energy.
        
        "CREATE TEMP TABLE Low AS SELECT P.ni, P.li, P.mi, P.nf, P.lf, P.mf, P.L, P.S, P.ell, Q.EiReq, Q.EiLow, P.Re_T_ell AS ReTLow, P.Im_T_ell AS ImTLow "
        "FROM tmat AS P NATURAL JOIN "
        "("
            "SELECT T.ni AS ni, T.li AS li, T.mi AS mi, T.nf AS nf, T.lf AS lf, T.mf AS mf, tmat.L AS L, T.S AS S, T.ell AS ell, MAX(tmat.Ei) AS EiLow, T.Ei AS EiReq "
            "FROM T CROSS JOIN tmat USING (ni,li,mi,nf,lf,mf,S,ell) "
            "WHERE tmat.Ei <= T.Ei "
            "GROUP BY T.ni, T.li, T.mi, T.nf, T.lf, T.mf, tmat.L,  T.S,  T.ell, T.Ei "
        ") AS Q WHERE P.Ei = Q.EiLow",
        
        "DELETE FROM Low WHERE EiLow > EiReq ",
        
        // temporary table where every transition, spin, total L, partial wave and required energy from 'T' contains also the lowest higher-or-equal energy.
        
        "CREATE TABLE Upp AS SELECT P.ni, P.li, P.mi, P.nf, P.lf, P.mf, P.L, P.S, P.ell, Q.EiReq, Q.EiUpp, P.Re_T_ell AS ReTUpp, P.Im_T_ell AS ImTUpp "
        "FROM tmat AS P NATURAL JOIN "
        "( "
            "SELECT T.ni AS ni, T.li AS li, T.mi AS mi, T.nf AS nf, T.lf AS lf, T.mf AS mf, tmat.L AS L, T.S AS S, T.ell AS ell, MIN(tmat.Ei) AS EiUpp, T.Ei AS EiReq "
            "FROM T CROSS JOIN tmat USING (ni,li,mi,nf,lf,mf,S,ell) "
            "WHERE T.Ei <= tmat.Ei "
            "GROUP BY T.ni, T.li, T.mi, T.nf, T.lf, T.mf, tmat.L,  T.S,  T.ell, T.Ei "
        ") AS Q WHERE P.Ei = Q.EiUpp",
        
        "DELETE FROM Upp WHERE EiUpp < EiReq ",
        
        // Interpolate partial cross sections for discrete transition.
        
        "INSERT OR REPLACE INTO " + IntegralCrossSection::Id + " "
        "SELECT ni, li, mi, nf, lf, mf, S, Ei, ell, SQRT(Ei - 1./(ni*ni) + 1./(nf*nf)) * (2 * S + 1) * (ReT * ReT + ImT * ImT) / 157.91367 / SQRT(Ei) AS sigma " // 16π²
        "FROM "
        "( "
            "SELECT Low.ni  AS ni, Low.li  AS li,  Low.mi  AS mi, "
                   "Low.nf  AS nf, Low.lf  AS lf,  Low.mf  AS mf, "
                   "Low.S   AS S,  Low.ell AS ell, Low.EiReq AS Ei, "
                   "SUM(INTERPOLATE(Low.EiLow, Low.ReTLow, Upp.EiUpp, Upp.ReTUpp, Low.EiReq)) AS ReT, "
                   "SUM(INTERPOLATE(Low.EiLow, Low.ImTLow, Upp.EiUpp, Upp.ImTUpp, Low.EiReq)) AS ImT "
            "FROM Low NATURAL JOIN Upp "
            "WHERE Ei - 1./(ni*ni) + 1./(nf*nf) >= 0 "
            "GROUP BY ni, li, mi, nf, lf, mf, S, Ei, ell "
        ")",
        
        // Insert ionization (no interpolation needed: L is used as the partial wave number here).
        
        "INSERT OR REPLACE INTO " + IntegralCrossSection::Id + " "
            "SELECT ni, li, mi, "
                   "0,  0,  0,  "
                   "S,  Ei, L,  "
                   "SUM(0.25*(2*S+1)*IONCS(QUOTE(cheb))/SQRT(Ei)) "
                "FROM " + IonizationF::Id + " "
                "GROUP BY ni, li, mi, S, Ei, L"
    };
    return cmd;
}

bool IntegralCrossSection::run (std::map<std::string,std::string> const & sdata) const
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
    int ell= Conv<int>(sdata, "ell",Id);
    int  S = Conv<int>(sdata, "S",  Id);
    
    // use mi >= 0; if mi < 0, flip both signs
    int mi = (mi0 < 0 ? -mi0 : mi0);
    int mf = (mi0 < 0 ? -mf0 : mf0);
    
    // energies and cross sections
    double E, sigma;
    rArray energies, E_arr, sigma_arr;
    
    // get energy / energies
    try {
        
        // is there a single energy specified using command line ?
        energies.push_back(Conv<double>(sdata, "Ei", Id));
        
    } catch (std::exception e) {
        
        // are there more energies specified using the STDIN ?
        energies = readStandardInput<double>();
    }
    
    // compose query
    sqlitepp::statement st(db);
    st << "SELECT Ei, sigma FROM " + IntegralCrossSection::Id + " "
            "WHERE ni = :ni "
            "  AND li = :li "
            "  AND mi = :mi "
            "  AND nf = :nf "
            "  AND lf = :lf "
            "  AND mf = :mf "
            "  AND ell= :ell"
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
        "# Integral cross section in " << unit_name(Lunits) << " for\n" <<
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n" <<
        "#     ell = " << ell << ", S = " << S << "\n" <<
        "# ordered by energy in " << unit_name(Eunits) << "\n" <<
        "#\n";
    OutputTable table;
    table.setWidth(15);
    table.setAlignment(OutputTable::left);
    table.write("# E        ", "sigma    ");
    table.write("# ---------", "---------");
    
    // terminate if no data
    if (E_arr.empty())
        return true;
    
    if (energies[0] < 0.)
    {
        // negative energy indicates full output
        for (std::size_t i = 0; i < E_arr.size(); i++)
            table.write(E_arr[i] / efactor, sigma_arr[i] * lfactor * lfactor);
    }
    else
    {
        // threshold for ionization
        double Eion = 1./(ni*ni);
        
        // interpolate (linear below ionization threshold, cspline above)
        rArray ics = (efactor * energies.front() < Eion) ? 
            interpolate_real(E_arr, sigma_arr, energies * efactor, gsl_interp_linear) :
            interpolate_real(E_arr, sigma_arr, energies * efactor, gsl_interp_cspline);
        
        // output
        for (std::size_t i = 0; i < energies.size(); i++)
            table.write(energies[i], ics[i] * lfactor * lfactor);
    }
    
    return true;
}
