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

void db_sqrt (sqlite3_context* pdb, int n, sqlite3_value** val)
{
    sqlite3_result_double(pdb, std::sqrt(sqlite3_value_double(*val)));
}

//
// custom function for integration of BLOB-represented Chebyshev
// expansion of the ionization amplitude
//

void db_ioncs (sqlite3_context* pdb, int n, sqlite3_value** val)
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
    auto fsqr = [&](double beta) -> double { return sqrabs(CB.clenshaw(std::sin(beta), tail)); };
    ClenshawCurtis<decltype(fsqr),double> integrator(fsqr);
    double result = integrator.integrate(0, special::constant::pi_quart);
    
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

bool IntegralCrossSection::initialize (sqlitepp::session & db) const
{
    //
    // define SQRT function
    //
    
    sqlite3_create_function
    (
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
    
    sqlite3_create_function
    (
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

std::vector<std::string> const & IntegralCrossSection::SQL_CreateTable () const
{
    static const std::vector<std::string> cmd = {
        "CREATE TABLE IF NOT EXISTS '" + IntegralCrossSection::Id + "' "
        "("
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
            "sigmaB DOUBLE PRECISION DEFAULT 0, "
            "PRIMARY KEY (ni,li,mi,nf,lf,mf,L,S,Ei)"
        ")"
    };
    
    return cmd;
}

std::vector<std::string> const & IntegralCrossSection::SQL_Update () const
{
    static const std::vector<std::string> cmd = {
        
        // insert discrete transitions
        
        "INSERT OR REPLACE INTO " + IntegralCrossSection::Id + " "
            "SELECT a.ni, a.li, a.mi, a.nf, a.lf, a.mf, a.L, a.S, a.Ei, "
                "sqrt(a.Ei-1./(a.ni*a.ni)+1./(a.nf*a.nf))/sqrt(a.Ei)*(2*a.S+1) * SUM(a.Re_T_ell    * b.Re_T_ell     + a.Im_T_ell     * b.Im_T_ell    )/157.91367, " // 16π²
                "sqrt(a.Ei-1./(a.ni*a.ni)+1./(a.nf*a.nf))/sqrt(a.Ei)*(2*a.S+1) * SUM(a.Re_TBorn_ell* b.Re_TBorn_ell + a.Im_TBorn_ell * b.Im_TBorn_ell)/157.91367 " // 16π²
                "FROM "
                "("
                    "SELECT ni,li,mi,nf,lf,mf,L,S,Ei,ell,Re_T_ell,Im_T_ell,Re_TBorn_ell,Im_TBorn_ell FROM " + TMatrix::Id + " "
                    "ORDER BY ell"
                ") AS a "
                "INNER JOIN "
                "("
                    "SELECT ni,li,mi,nf,lf,mf,L,S,Ei,ell,Re_T_ell,Im_T_ell,Re_TBorn_ell,Im_TBorn_ell FROM " + TMatrix::Id + " "
                    "ORDER BY ell"
                ") AS b "
                "ON a.ni = b.ni AND a.li = b.li AND a.mi = b.mi AND "
                "   a.nf = b.nf AND a.lf = b.lf AND a.mf = b.mf AND "
                "   a.S = b.S  AND a.Ei = b.Ei AND a.ell = b.ell "
                "WHERE a.Ei > 0 AND a.Ei - 1./(a.ni*a.ni) + 1./(a.nf*a.nf) > 0 "
                "GROUP BY a.ni, a.li, a.mi, a.nf, a.lf, a.mf, a.L, a.S, a.Ei",
        
        // insert ionization
        
        "INSERT OR REPLACE INTO " + IntegralCrossSection::Id + " "
            "SELECT ni, li, mi, 0, 0, 0, L, S, Ei, "
                "SUM(0.25*(2*S+1)*ioncs(QUOTE(cheb))/sqrt(Ei)), "
                "0 " // TODO Born T-matrix (to be implemented in future)
            "FROM " + IonizationF::Id + " "
            "GROUP BY ni, li, mi, L, S, Ei"
    };
    return cmd;
}

bool IntegralCrossSection::run
(
    sqlitepp::session & db,
    std::map<std::string,std::string> const & sdata
) const
{
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
    double E, sigma, sigmaB;
    rArray energies, E_arr, sigma_arr, sigmaB_arr;
    
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
    st << "SELECT Ei, sigma, sigmaB FROM " + IntegralCrossSection::Id + " "
            "WHERE ni = :ni "
            "  AND li = :li "
            "  AND mi = :mi "
            "  AND nf = :nf "
            "  AND lf = :lf "
            "  AND mf = :mf "
            "  AND  L = :L  "
            "  AND  S = :S  "
            "ORDER BY Ei ASC",
        sqlitepp::into(E), sqlitepp::into(sigma), sqlitepp::into(sigmaB),
        sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
        sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
        sqlitepp::use(L), sqlitepp::use(S);
    
    // retrieve data
    while (st.exec())
    {
        E_arr.push_back(E);
        sigma_arr.push_back(sigma);
        sigmaB_arr.push_back(sigmaB);
    }
    
    // write header
    std::cout << logo("#") <<
        "# Integral cross section in " << unit_name(Lunits) << " for\n" <<
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n" <<
        "#     L = " << L << ", S = " << S << "\n" <<
        "# ordered by energy in " << unit_name(Eunits) << "\n" <<
        "#\n" <<
        "# E\tσ\tσBorn\n";
    
    // terminate if no data
    if (E_arr.empty())
        return true;
    
    if (energies[0] < 0.)
    {
        // negative energy indicates full output
        for (std::size_t i = 0; i < E_arr.size(); i++)
            std::cout << E_arr[i] / efactor << "\t" << sigma_arr[i] * lfactor * lfactor << "\n";
    }
    else
    {
        // threshold for ionization
        double Eion = 1./(ni*ni);
        
        // interpolate (linear below ionization threshold, cspline above)
        rArray ics = (efactor * energies.front() < Eion) ? 
            interpolate_real(E_arr, sigma_arr, energies * efactor, gsl_interp_linear) :
            interpolate_real(E_arr, sigma_arr, energies * efactor, gsl_interp_cspline);
        rArray icsB = (efactor * energies.front() < Eion) ? 
            interpolate_real(E_arr, sigmaB_arr, energies * efactor, gsl_interp_linear) :
            interpolate_real(E_arr, sigmaB_arr, energies * efactor, gsl_interp_cspline);
        
        // output
        for (std::size_t i = 0; i < energies.size(); i++)
            std::cout << energies[i] << "\t" << ics[i] * lfactor * lfactor << "\t" << icsB[i] * lfactor * lfactor << std::endl;
    }
    
    return true;
}
