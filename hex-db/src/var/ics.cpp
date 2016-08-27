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

#include <sqlite3.h>

// --------------------------------------------------------------------------------- //

#include <gsl/gsl_interp.h>

// --------------------------------------------------------------------------------- //

#include "hex-chebyshev.h"
#include "hex-clenshawcurtis.h"
#include "hex-interpolate.h"
#include "hex-version.h"

// --------------------------------------------------------------------------------- //

#include "../quantities.h"
#include "../utils.h"

// --------------------------------------------------------------------------------- //

createNewScatteringQuantity(IntegralCrossSection);

// --------------------------------------------------------------------------------- //

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

// --------------------------------------------------------------------------------- //

std::string IntegralCrossSection::name ()
{
    return "ics";
}

std::string IntegralCrossSection::description ()
{
    return "Integral cross section.";
}

std::vector<std::string> IntegralCrossSection::dependencies ()
{
    return std::vector<std::string>
    {
        "tmat",
        "ionf"
    };
}

std::vector<std::pair<std::string,std::string>> IntegralCrossSection::params ()
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

std::vector<std::string> IntegralCrossSection::vparams ()
{
    return std::vector<std::string>
    {
        "Ei"
    };
}

// --------------------------------------------------------------------------------- //

bool IntegralCrossSection::initialize (sqlitepp::session & db)
{
    // define square root function
    sqlite3_create_function(db.impl(), "SQRT", 1, SQLITE_UTF8, nullptr, &db_sqrt, nullptr, nullptr);
    
    // define Gauss-Chebyshev integration of squared Chebyshev expansion
    sqlite3_create_function(db.impl(), "IONCS", 1, SQLITE_UTF8, nullptr, &db_ioncs, nullptr, nullptr);
    
    // define linear interpolation of two data points (5 arguments: x1, y1, x2, y2, x)
    sqlite3_create_function(db.impl(), "INTERPOLATE", 5, SQLITE_UTF8, nullptr, &db_interpolate, nullptr, nullptr);
    
    return ScatteringQuantity::initialize(db);
}

bool IntegralCrossSection::createTable ()
{
    sqlitepp::statement st (session());
    st <<
        "CREATE TABLE IF NOT EXISTS 'ics' "
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
        ")";
    
    try
    {
        st.exec();
    }
    catch (sqlitepp::exception & e)
    {
        std::cerr << "ERROR: Creation of table 'ics' failed!" << std::endl;
        std::cerr << "       code = " << e.code() << " (\"" << e.what() << "\")" << std::endl;
        return false;
    }
    
    return true;
}

bool IntegralCrossSection::updateTable ()
{
    // Get merged available energies from all total angular momenta for every discrete transition and partial wave.
    int ni, li, mi, nf, lf, mf, S, ell;
    sqlitepp::statement st1 (session());
    st1 << "SELECT DISTINCT ni,li,mi,nf,lf,mf,S,ell FROM 'tmat'",
        sqlitepp::into(ni), sqlitepp::into(li), sqlitepp::into(mi),
        sqlitepp::into(nf), sqlitepp::into(lf), sqlitepp::into(mf),
        sqlitepp::into(S),  sqlitepp::into(ell);
    
    // for all transitions and partial waves
    while (st1.exec())
    {
        // get all merged energies for this transition and partial wave
        double E = 0;
        sqlitepp::statement st2 (session());
        st2 << "SELECT DISTINCT Ei FROM 'tmat' "
               "WHERE ni = :ni AND li = :li AND mi = :mi "
               "  AND nf = :nf AND lf = :lf AND mf = :mf "
               "  AND S  = :S  AND ell = :ell ORDER BY Ei ASC",
               sqlitepp::into(E),
               sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
               sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
               sqlitepp::use(S),  sqlitepp::use(ell);
        rArray merged_energies;
        while (st2.exec())
        {
            merged_energies.push_back(E);
        }
        
        // get all total angular momenta and parities for this transition and partial wave
        int L;
        sqlitepp::statement st3 (session());
        st3 << "SELECT DISTINCT L FROM 'tmat' "
            "WHERE ni = :ni AND li = :li AND mi = :mi "
            "  AND nf = :nf AND lf = :lf AND mf = :mf "
            "  AND S  = :S  AND ell = :ell ORDER BY L ASC",
            sqlitepp::into(L),
            sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
            sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
            sqlitepp::use(S),  sqlitepp::use(ell);
        iArray angular_momenta;
        while (st3.exec())
        {
            angular_momenta.push_back(L);
        }
        
        // merge-interpolate all T-matrix data for this transition and partial wave
        cArray merged_T (merged_energies.size());
        for (int L : angular_momenta)
        {
            rArray energies, Re_T, Im_T;
            // retrieve data
            double ret, imt;
            sqlitepp::statement st4 (session());
            st4 << "SELECT Ei,Re_T_ell,Im_T_ell FROM 'tmat' "
                "WHERE ni = :ni AND li = :li AND mi = :mi "
                "  AND nf = :nf AND lf = :lf AND mf = :mf "
                "  AND S  = :S  AND ell = :ell AND L = :L",
                sqlitepp::into(E), sqlitepp::into(ret), sqlitepp::into(imt),
                sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
                sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
                sqlitepp::use(S),  sqlitepp::use(ell), sqlitepp::use(L);
            while (st4.exec())
            {
                energies.push_back(E);
                Re_T.push_back(ret);
                Im_T.push_back(imt);
            }
            
            // skip datasets with too few samples
            if (energies.size() < gsl_interp_type_min_size(gsl_interp_akima))
                continue;
            
            // interpolate data
            gsl_interp* spline_Re = gsl_interp_alloc(gsl_interp_akima, energies.size());
            gsl_interp* spline_Im = gsl_interp_alloc(gsl_interp_akima, energies.size());
            gsl_interp_init(spline_Re, energies.data(), Re_T.data(), energies.size());
            gsl_interp_init(spline_Im, energies.data(), Im_T.data(), energies.size());
            # pragma omp parallel
            {
                gsl_interp_accel* accel_Re = gsl_interp_accel_alloc();
                gsl_interp_accel* accel_Im = gsl_interp_accel_alloc();
                
                double ret, imt;
                
                # pragma omp for
                for (std::size_t i = 0; i < merged_energies.size(); i++)
                {
                    if (energies.front() <= merged_energies[i] and merged_energies[i] <= energies.back())
                    {
                        int err_re = gsl_interp_eval_e(spline_Re, energies.data(), Re_T.data(), merged_energies[i], accel_Re, &ret);
                        int err_im = gsl_interp_eval_e(spline_Im, energies.data(), Im_T.data(), merged_energies[i], accel_Im, &imt);
                        
                        if (err_re != GSL_SUCCESS)
                        {
                            std::cout << "Warning: Failed to interpolate real data for transition (" << ni << "," << li << "," << mi << ") -> (" << nf << "," << lf << "," << mf << "), "
                                "L = " << L << ", ell = " << ell << ", Ei = " << merged_energies[i] << " (" << gsl_strerror(err_re) << ")" << std::endl;
                        }
                        
                        if (err_im != GSL_SUCCESS)
                        {
                            std::cout << "Warning: Failed to interpolate imag data for transition (" << ni << "," << li << "," << mi << ") -> (" << nf << "," << lf << "," << mf << "), "
                                "L = " << L << ", ell = " << ell << ", Ei = " << merged_energies[i] << " (" << gsl_strerror(err_re) << ")" << std::endl;
                        }
                        
                        merged_T[i] += Complex (ret, imt);
                    }
                }
                
                gsl_interp_accel_free(accel_Re);
                gsl_interp_accel_free(accel_Im);
            }
            gsl_interp_free(spline_Re);
            gsl_interp_free(spline_Im);
        }
        
        // write interpolated data to the database
        sqlitepp::statement st5 (session());
        double ics, Ei, Ef;
        std::size_t rows = 0;
        st5 << "INSERT OR REPLACE INTO 'ics' VALUES ";
        for (std::size_t i = 0; i < merged_energies.size(); i++)
        {
            // initial and final energies
            Ei = merged_energies[i];
            Ef = Ei - 1./(ni*ni) + 1./(nf*nf);
            
            // skip forbidden channels
            if (Ef < 0)
                continue;
            
            // calculate partial integral cross section
            ics = std::sqrt(Ef/Ei) * (2 * S + 1) * sqrabs(merged_T[i]) / std::pow(4 * special::constant::pi, 2);
            
            // separate data n-tuples with comma
            if (rows++ > 0)
                st5.q() << ",";
            
            // add another data n-tuple
            st5.q() << " (" << ni << "," << li << "," << mi << "," << nf << "," << lf << "," << mf << "," << S << "," << Ei << "," << ell << "," << ics << ")";
        }
        if (rows > 0)
        {
            // add all data in one transaction
            st5.exec();
        }
    }
    
    // Insert ionization (no interpolation needed: L is used as the partial wave number here.
    
    sqlitepp::statement st6 (session());
    st6 << "INSERT OR REPLACE INTO 'ics' "
            "SELECT ni, li, mi, "
                   "0,  0,  0,  "
                   "S,  Ei, L,  "
                   "SUM(0.25*(2*S+1)*IONCS(QUOTE(cheb))/SQRT(Ei)) "
                "FROM 'ionf' "
                "GROUP BY ni, li, mi, S, Ei, L";
    st6.exec();
    
    return ScatteringQuantity::updateTable();
}

// --------------------------------------------------------------------------------- //

bool IntegralCrossSection::run (std::map<std::string,std::string> const & sdata)
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
    int ell= Conv<int>(sdata, "ell",name());
    int  S = Conv<int>(sdata, "S",  name());
    
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
            interpolate_real(E_arr, E_arr * sigma_arr, energies * efactor, gsl_interp_linear) / (energies * efactor) :
            interpolate_real(E_arr, E_arr * sigma_arr, energies * efactor, gsl_interp_akima ) / (energies * efactor);
        
        // output
        for (std::size_t i = 0; i < energies.size(); i++)
            table.write(energies[i], ics[i] * lfactor * lfactor);
    }
    
    return true;
}
