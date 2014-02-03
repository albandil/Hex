/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2014                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

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
// custom function for evaluation of square root within a SQL statement
//

void db_sqrt(sqlite3_context* pdb, int n, sqlite3_value** val)
{
    sqlite3_result_double(pdb, sqrt(sqlite3_value_double(*val)));
}

//
// custom function for integration of BLOB-represented Chebyshev
// expansion of the ionization amplitude
//

void db_ioncs(sqlite3_context* pdb, int n, sqlite3_value** val)
{
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
    auto fsqr = [&](double beta) -> double { return sqrabs(CB.clenshaw(sin(beta), tail)); };
    ClenshawCurtis<decltype(fsqr),double> integrator(fsqr);
    double result = integrator.integrate(0, 0.25 * M_PI);
    
    // use result of the integration
    sqlite3_result_double(pdb, result);
}

// -------------------------------------------------------------------------- //

const std::string IntegralCrossSection::Id = "ics";
const std::string IntegralCrossSection::Description = "Integral cross section.";
const std::vector<std::string> IntegralCrossSection::Dependencies = {
    "ni", "li", "mi",
    "nf", "lf", "mf",
    "L", "S",
    "Ei"
};
const std::vector<std::string> IntegralCrossSection::VecDependencies = { "Ei" };

bool IntegralCrossSection::initialize(sqlitepp::session & db) const
{
    //
    // define SQRT function
    //
    
    sqlite3_create_function (
        db.impl(),
        "sqrt",
        1,              // pass single argument
        SQLITE_UTF8,
        nullptr,
        &db_sqrt,
        nullptr,
        nullptr
    );
    
    //
    // define Gauss-Chebyshev integration of squared Chebyshev expansion
    //
    
    sqlite3_create_function (
        db.impl(),
        "ioncs",
        1,              // pass single argument
        SQLITE_UTF8,
        nullptr,
        &db_ioncs,
        nullptr,
        nullptr
    );
    
    return true;
}

std::vector<std::string> const & IntegralCrossSection::SQL_CreateTable() const
{
    static const std::vector<std::string> cmd = {
        "CREATE TABLE '" + IntegralCrossSection::Id + "' ("
            "ni INTEGER, "
            "li INTEGER, "
            "mi INTEGER, "
            "nf INTEGER, "
            "lf INTEGER, "
            "mf INTEGER, "
            "L  INTEGER, "
            "S  INTEGER, "
            "Ei DOUBLE PRECISION, "
            "sigma DOUBLE PRECISION, "
            "PRIMARY KEY (ni,li,mi,nf,lf,mf,L,S,Ei)"
        ")"
    };
    
    return cmd;
}

std::vector<std::string> const & IntegralCrossSection::SQL_Update() const
{
    static const std::vector<std::string> cmd = {
        
        // insert discrete transitions
        
        "INSERT OR REPLACE INTO " + IntegralCrossSection::Id + " "
            "SELECT ni, li, mi, nf, lf, mf, L, S, Ei, "
                "sqrt(Ei-1./(ni*ni)+1./(nf*nf))/sqrt(Ei)*(2*S+1)*SUM(Re_T_ell*Re_T_ell+Im_T_ell*Im_T_ell)/157.91367 " // 16π²
            "FROM " + TMatrix::Id + " "
            "WHERE Ei > 0 AND Ei - 1./(ni*ni) + 1./(nf*nf) > 0 "
            "GROUP BY ni, li, mi, nf, lf, mf, L, S, Ei",
        
        // insert ionization
        
        "INSERT OR REPLACE INTO " + IntegralCrossSection::Id + " "
            "SELECT ni, li, mi, 0, 0, 0, L, S, Ei, "
                "SUM(0.25*(2*S+1)*ioncs(QUOTE(cheb))/sqrt(Ei)) "
            "FROM " + IonizationF::Id + " "
            "GROUP BY ni, li, mi, L, S, Ei"
    };
    return cmd;
}

bool IntegralCrossSection::run (
    sqlitepp::session & db,
    std::map<std::string,std::string> const & sdata
) const {
    
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
    double lfactor = change_units(lUnit_au, Lunits);
    
    // scattering event parameters
    int ni = As<int>(sdata, "ni", Id);
    int li = As<int>(sdata, "li", Id);
    int mi = As<int>(sdata, "mi", Id);
    int nf = As<int>(sdata, "nf", Id);
    int lf = As<int>(sdata, "lf", Id);
    int mf = As<int>(sdata, "mf", Id);
    int  L = As<int>(sdata, "L", Id);
    int  S = As<int>(sdata, "S", Id);
    
    // energies and cross sections
    double E, sigma;
    rArray energies, E_arr, sigma_arr;
    
    // get energy / energies
    try {
        
        // is there a single energy specified using command line ?
        energies.push_back(As<double>(sdata, "Ei", Id));
        
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
            "  AND  L = :L  "
            "  AND  S = :S  "
            "ORDER BY Ei ASC",
        sqlitepp::into(E), sqlitepp::into(sigma),
        sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
        sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
        sqlitepp::use(L), sqlitepp::use(S);
    
    // retrieve data
    while (st.exec())
    {
        E_arr.push_back(E);
        sigma_arr.push_back(sigma);
    }
    
    // write header
    std::cout << logo() <<
        "# Integral cross section in " << unit_name(Lunits) << " for\n" <<
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n" <<
        "#     L = " << L << ", S = " << S << "\n" <<
        "# ordered by energy in " << unit_name(Eunits) << "\n" <<
        "#\n" <<
        "# E\t σ\n";
    
    // terminate if no data
    if (E_arr.empty())
        return true;
    
    if (energies[0] < 0.)
    {
        // negative energy indicates full output
        for (size_t i = 0; i < E_arr.size(); i++)
            std::cout << E_arr[i] / efactor << "\t" << sigma_arr[i] * lfactor * lfactor << "\n";
    }
    else
    {
        // threshold for ionization
        double Eion = 1./(ni*ni);
        
        // interpolate
        rArray ics = (efactor * energies.front() < Eion) ? 
            interpolate_real(E_arr, sigma_arr, energies * efactor, gsl_interp_linear) :
            interpolate_real(E_arr, sigma_arr, energies * efactor, gsl_interp_cspline);
        
        // output
        for (size_t i = 0; i < energies.size(); i++)
            std::cout << energies[i] << "\t" << ics[i] * lfactor * lfactor << "\n";
    }
    
    return true;
}
