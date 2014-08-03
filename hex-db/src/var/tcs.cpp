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

#include "../interpolate.h"
#include "../variables.h"
#include "../version.h"

const std::string TotalCrossSection::Id = "tcs";
const std::string TotalCrossSection::Description = "Total cross section.";
const std::vector<std::string> TotalCrossSection::Dependencies = {
    "ni", "li", "mi", 
    "nf", "lf", "mf",
    "Ei"
};
const std::vector<std::string> TotalCrossSection::VecDependencies = { "Ei" };

bool TotalCrossSection::initialize(sqlitepp::session & db) const
{
    return true;
}

std::vector<std::string> const & TotalCrossSection::SQL_Update () const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

std::vector<std::string> const & TotalCrossSection::SQL_CreateTable () const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

bool TotalCrossSection::run
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
    
    // energies and cross sections
    double E, sigma, sigmab, sigmaB;
    rArray energies, E_arr, EB_arr, sigma_arr, sigmab_arr, sigmaB_arr;
    
    // get energy / energies
    try {
        
        // is there a single energy specified using command line ?
        energies.push_back(As<double>(sdata, "Ei", Id));
        
    } catch (std::exception e) {
        
        // are there more energies specified using the STDIN ?
        energies = readStandardInput<double>();
    }
    
    //
    // load partial wave constributions
    //
    
    // compose query
    sqlitepp::statement st(db);
    st << "SELECT Ei, sum(sigma), sum(sigmaB) FROM " + CompleteCrossSection::Id + " "
            "WHERE ni = :ni "
            "  AND li = :li "
            "  AND mi = :mi "
            "GROUP BY Ei "
            "ORDER BY Ei ASC",
        sqlitepp::into(E), sqlitepp::into(sigma), sqlitepp::into(sigmab),
        sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi);
    
    // retrieve data
    while (st.exec())
    {
        E_arr.push_back(E);
        sigma_arr.push_back(sigma);
        sigmab_arr.push_back(sigmab);
    }
    
    //
    // load Born cross section
    //
    
    // compose query
    sqlitepp::statement stb(db);
    stb << "SELECT Ei, SUM(borncs(cheb)) FROM " + BornFullTMatrix::Id + " "
            "WHERE ni = :ni "
            "  AND li = :li "
            "  AND mi = :mi "
            "GROUP BY Ei "
            "ORDER BY Ei ASC",
        sqlitepp::into(E), sqlitepp::into(sigmaB),
        sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi);
    
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
    std::cout << logo("#") <<
        "# Total cross section in " << unit_name(Lunits) << " for\n"
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << "\n" <<
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
        rArray tcs = (efactor * energies.front() < Eion) ? 
            interpolate_real(E_arr, sigma_arr, energies * efactor, gsl_interp_linear) :
            interpolate_real(E_arr, sigma_arr, energies * efactor, gsl_interp_cspline);
        rArray tcsb = (efactor * energies.front() < Eion) ? 
            interpolate_real(E_arr, sigmab_arr, energies * efactor, gsl_interp_linear) :
            interpolate_real(E_arr, sigmab_arr, energies * efactor, gsl_interp_cspline);
        rArray tcsB = (efactor * energies.front() < Eion) ? 
            interpolate_real(EB_arr, sigmaB_arr, energies * efactor, gsl_interp_linear) :
            interpolate_real(EB_arr, sigmaB_arr, energies * efactor, gsl_interp_cspline);
        
        // corrected cross section
        rArray cs = tcsB + (tcs - tcsb);
        
        // output
        for (size_t i = 0; i < energies.size(); i++)
            std::cout << energies[i] << "\t" << (std::isfinite(cs[i]) ? cs[i] * lfactor * lfactor : 0.) << "\n";
    }
    
    return true;
}
