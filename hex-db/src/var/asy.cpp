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

#include <sqlitepp/sqlitepp.hpp>

#include "../interpolate.h"
#include "../special.h"
#include "../variables.h"
#include "../version.h"


/* * * * * * * * * * * * * * External prototypes * * * * * * * * * * * * * * */
rArray differential_cross_section (
    sqlitepp::session & db,
    int ni, int li, int mi,
    int nf, int lf, int mf,
    int S, double E,
    rArray const & angles
);
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


const std::string SpinAsymmetry::Id = "asy";
const std::string SpinAsymmetry::Description = "Spin asymetry.";
const std::vector<std::string> SpinAsymmetry::Dependencies = {
    "ni", "li", "mi", 
    "nf", "lf", "mf",
    "Ei", "theta"
};
const std::vector<std::string> SpinAsymmetry::VecDependencies = { "theta" };

bool SpinAsymmetry::initialize(sqlitepp::session & db) const
{
    return true;
}

std::vector<std::string> const & SpinAsymmetry::SQL_CreateTable() const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

std::vector<std::string> const & SpinAsymmetry::SQL_Update() const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

bool SpinAsymmetry::run (
    sqlitepp::session & db,
    std::map<std::string,std::string> const & sdata
) const {
    
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
    double afactor = change_units(Aunits, aUnit_rad);
    
    // scattering event parameters
    int ni = As<int>(sdata, "ni", Id);
    int li = As<int>(sdata, "li", Id);
    int mi = As<int>(sdata, "mi", Id);
    int nf = As<int>(sdata, "nf", Id);
    int lf = As<int>(sdata, "lf", Id);
    int mf = As<int>(sdata, "mf", Id);
    double E = As<double>(sdata, "Ei", Id) * efactor;
    
    // angles
    rArray angles;
    
    // get angle / angles
    try {
        
        // is there a single angle specified using command line ?
        angles.push_back(As<double>(sdata, "theta", Id));
        
    } catch (std::exception e) {
        
        // are there more angles specified using the STDIN ?
        angles = readStandardInput<double>();
    }
    
    // compute cross sections
    rArray dcs0 = differential_cross_section (db, ni,li,mi, nf,lf,mf, 0, E, angles * afactor);
    rArray dcs1 = differential_cross_section (db, ni,li,mi, nf,lf,mf, 1, E, angles * afactor);
    
    // compute spin asymetry
    rArray asy = (dcs0 - dcs1) / (dcs0 + 3. * dcs1);
    
    // substitute possible "nan"-s and "inf"-s by zeros
    for (double & x : asy)
    {
        if (not std::isfinite(x))
            x = 0.;
    }
    
    // write out
    std::cout << logo() <<
        "# Spin asymetry for \n" <<
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n" <<
        "#     E = " << E/efactor << " " << unit_name(Eunits)
                     << " ordered by angle in " << unit_name(Aunits) << "\n" <<
        "# θ\t dσ\n";
    for (size_t i = 0; i < angles.size(); i++)
        std::cout << angles[i] << "\t" << asy[i] << "\t" << dcs0[i] << "\t" << dcs1[i] << "\n";
    
    return true;
}
