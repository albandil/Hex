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
const std::vector<std::string> CompleteCrossSection::Dependencies = {
    "ni", "li", "mi",
    "nf", "lf", "mf",
    "Ei"
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
    static const std::vector<std::string> cmd = {
        "CREATE TABLE IF NOT EXISTS '" + CompleteCrossSection::Id + "' "
        "("
            "ni INTEGER, "
            "li INTEGER, "
            "mi INTEGER, "
            "nf INTEGER, "
            "lf INTEGER, "
            "mf INTEGER, "
            "Ei DOUBLE PRECISION, "
            "sigma DOUBLE PRECISION, "
            "sigmaB DOUBLE PRECISION DEFAULT 0, "
            "PRIMARY KEY (ni,li,mi,nf,lf,mf,Ei)"
        ")"
    };
    
    return cmd;
}

std::vector<std::string> const & CompleteCrossSection::SQL_Update () const
{
    static std::vector<std::string> cmd = {
        
        //
        // Compute the complete cross section
        //                     ==
        //    S   kf 2S+1  1   \    LS  L'S*
        //   σ  = —— ———— ———  /   T   T
        //        ki   4  4π²  ==   ℓ   ℓ
        //                    ℓLL'
        
        "INSERT OR REPLACE INTO " + CompleteCrossSection::Id + " "
            "SELECT ni, li, mi, nf, lf, mf, Ei, "
                "sqrt(Ei-1./(ni*ni)+1./(nf*nf))/sqrt(Ei)*SUM((2*S+1)*sigma/157.91367), "
                "sqrt(Ei-1./(ni*ni)+1./(nf*nf))/sqrt(Ei)*SUM((2*S+1)*sigmaB/157.91367) "
            "FROM "
            "("
                "SELECT a.ni AS ni, a.li AS li, a.mi AS mi, a.nf AS nf, a.lf AS lf, a.mf AS mf, a.S AS S, a.Ei AS Ei, a.ell AS ell, "
                    "SUM(a.Re_T_ell    * b.Re_T_ell     + a.Im_T_ell     * b.Im_T_ell    ) AS sigma, "
                    "SUM(a.Re_TBorn_ell* b.Re_TBorn_ell + a.Im_TBorn_ell * b.Im_TBorn_ell) AS sigmaB "
                "FROM "
                "("
                    "SELECT ni,li,mi,nf,lf,mf,L,S,Ei,ell,Re_T_ell,Im_T_ell,Re_TBorn_ell,Im_TBorn_ell FROM " + TMatrix::Id + " "
                    "ORDER BY ell "
                ") AS a "
                "INNER JOIN "
                "("
                    "SELECT ni,li,mi,nf,lf,mf,L,S,Ei,ell,Re_T_ell,Im_T_ell,Re_TBorn_ell,Im_TBorn_ell FROM " + TMatrix::Id + " "
                    "ORDER BY ell "
                ") AS b "
                "ON a.ni = b.ni AND a.li = b.li AND a.mi = b.mi AND "
                "a.nf = b.nf AND a.lf = b.lf AND a.mf = b.mf AND "
                    "a.S = b.S  AND a.Ei = b.Ei AND a.ell = b.ell "
                "GROUP BY a.ni, a.li, a.mi, a.nf, a.lf, a.mf, a.S, a.Ei, a.ell "
            ") "
            "GROUP BY ni, li, mi, nf, lf, mf, Ei"
    };
    
    return cmd;
}

bool CompleteCrossSection::run
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
    
    // energies and cross sections
    double E, sigma, sigmab, sigmaB;
    rArray energies, E_arr, sigma_arr, sigmab_arr, EB_arr, sigmaB_arr;
    
    // get energy / energies
    try {
        
        // is there a single energy specified using command line ?
        energies.push_back(As<double>(sdata, "Ei", Id));
        
    } catch (std::exception e) {
        
        // are there more energies specified using the STDIN ?
        energies = readStandardInput<double>();
    }
    
    //
    // get contributions from the partial wave expansion
    //
    
    // compose query
    sqlitepp::statement st(db);
    st << "SELECT Ei, sigma, sigmaB FROM " + CompleteCrossSection::Id + " "
            "WHERE ni = :ni "
            "  AND li = :li "
            "  AND mi = :mi "
            "  AND nf = :nf "
            "  AND lf = :lf "
            "  AND mf = :mf "
            "GROUP BY Ei "
            "ORDER BY Ei ASC",
        sqlitepp::into(E), sqlitepp::into(sigma), sqlitepp::into(sigmab),
        sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
        sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf);
    
    // retrieve data
    while (st.exec())
    {
        E_arr.push_back(E);
        sigma_arr.push_back(sigma);
        sigmab_arr.push_back(sigmab);
    }
    
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
        sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
        sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf);
    
    // retrieve data
    while (stb.exec())
    {
        EB_arr.push_back(E);
        sigmaB_arr.push_back(sigmaB);
    }
    
    //
    // compute and write output
    //
    
    // write header
    std::cout << logo() <<
        "# Complete cross section in " << unit_name(Lunits) << " for\n" <<
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n" <<
        "# ordered by energy in " << unit_name(Eunits) << "\n" <<
        "# \n" <<
        "# E\t σ\n";
    
    // terminate if no data
    if (E_arr.empty())
        return true;
    
    // threshold for ionization
    double Eion = 1./(ni*ni);
    
    // negative energy indicates output of all available cross sections
    if (energies[0] < 0.)
    {
        // interpolate Born cross section to partial waves' energies
        rArray sigmaBorn = (efactor * energies.front() < Eion) ? 
            interpolate_real(EB_arr, sigmaB_arr, E_arr, gsl_interp_linear) :
            interpolate_real(EB_arr, sigmaB_arr, E_arr, gsl_interp_cspline);
        
        // output corrected cross section
        for (size_t i = 0; i < E_arr.size(); i++)
            std::cout << E_arr[i] / efactor << "\t" << (sigmaBorn[i] + (sigma_arr[i] - sigmab_arr[i])) * lfactor * lfactor << "\n";
    }
    else
    {
        // interpolate for given 'energies'
        rArray ccs = (efactor * energies.front() < Eion) ? 
            interpolate_real(E_arr, sigma_arr, energies * efactor, gsl_interp_linear) :
            interpolate_real(E_arr, sigma_arr, energies * efactor, gsl_interp_cspline);
        rArray ccsb = (efactor * energies.front() < Eion) ? 
            interpolate_real(E_arr, sigmab_arr, energies * efactor, gsl_interp_linear) :
            interpolate_real(E_arr, sigmab_arr, energies * efactor, gsl_interp_cspline);
        rArray ccsB = (efactor * energies.front() < Eion) ? 
            interpolate_real(EB_arr, sigmaB_arr, energies * efactor, gsl_interp_linear) :
            interpolate_real(EB_arr, sigmaB_arr, energies * efactor, gsl_interp_cspline);
        
        // corrected cross section
        rArray cs = ccsB + (ccs - ccsb);
        
        // output
        for (size_t i = 0; i < energies.size(); i++)
            std::cout << energies[i] << "\t" << (std::isfinite(cs[i]) ? cs[i] * lfactor * lfactor : 0.) << "\n";
    }
    
    return true;
}
