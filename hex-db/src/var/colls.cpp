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

const std::string CollisionStrength::Id = "colls";
const std::string CollisionStrength::Description = "Collision strength (energy scaled integral cross section).";
const std::vector<std::string> CollisionStrength::Dependencies = {
    "ni", "li", "mi",
    "nf", "lf", "mf",
    "L", "S",
    "Ei"
};
const std::vector<std::string> CollisionStrength::VecDependencies = { "Ei" };

bool CollisionStrength::initialize(sqlitepp::session & db) const
{
    return true;
}

std::vector<std::string> const & CollisionStrength::SQL_CreateTable () const
{
    static std::vector<std::string> cmd;
    return cmd;
}

std::vector<std::string> const & CollisionStrength::SQL_Update () const
{
    static std::vector<std::string> cmd;
    return cmd;
}

bool CollisionStrength::run
(
    sqlitepp::session & db,
    std::map<std::string,std::string> const & sdata
) const
{
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
//     double lfactor = change_units(lUnit_au, Lunits);
    
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
    double E, sigma, sigmab;
    rArray energies, E_arr, sigma_arr, sigmab_arr;
    
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
    st << "SELECT Ei, sigma, sigmab FROM " + IntegralCrossSection::Id + " "
            "WHERE ni = :ni "
            "  AND li = :li "
            "  AND mi = :mi "
            "  AND nf = :nf "
            "  AND lf = :lf "
            "  AND mf = :mf "
            "  AND  L = :L  "
            "  AND  S = :S  "
            "ORDER BY Ei ASC",
        sqlitepp::into(E), sqlitepp::into(sigma), sqlitepp::into(sigmab),
        sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
        sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
        sqlitepp::use(L), sqlitepp::use(S);
    
    // retrieve data
    while (st.exec())
    {
        E_arr.push_back(E);
        sigma_arr.push_back(sigma);
        sigmab_arr.push_back(sigmab);
    }
    
    // write header
    std::cout << logo() <<
        "# Collision strength (dimensionless) for\n"
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n" <<
        "#     L = " << L << ", S = " << S << "\n" <<
        "# ordered by energy in " << unit_name(Eunits) << "\n" <<
        "# \n" <<
        "# E\tΩ\tΩBorn\n";
    
    // terminate if no data
    if (E_arr.size() == 0)
        return true;
    
    if (energies[0] < 0.)
    {
        // negative energy indicates full output
        for (size_t i = 0; i < E_arr.size(); i++)
        {
            double Omega = E_arr[i] * efactor * (2*L+1) * (2*S+1) * sigma_arr[i];
            double OmegaB = E_arr[i] * efactor * (2*L+1) * (2*S+1) * sigmab_arr[i];
            
            std::cout << E_arr[i] / efactor << "\t" << Omega << "\t" << OmegaB << "\n";
        }
    }
    else
    {
        // threshold for ionization
        double Eion = 1./(ni*ni);
        
        // interpolate
        rArray interp = (efactor * energies.front() < Eion) ? 
            interpolate_real(E_arr, sigma_arr, energies * efactor, gsl_interp_linear) :
            interpolate_real(E_arr, sigma_arr, energies * efactor, gsl_interp_cspline);
        rArray interpB = (efactor * energies.front() < Eion) ? 
            interpolate_real(E_arr, sigmab_arr, energies * efactor, gsl_interp_linear) :
            interpolate_real(E_arr, sigmab_arr, energies * efactor, gsl_interp_cspline);
        
        // compute collision strength
        rArray omegas = energies * (2*L+1) * (2*S+1) * interp * efactor;
        rArray omegasB = energies * (2*L+1) * (2*S+1) * interpB * efactor;
        
        // output
        for (size_t i = 0; i < energies.size(); i++)
            std::cout << energies[i] << "\t" << omegas[i] << "\t" << omegasB[i] << "\n";
    }
    
    return true;
}
