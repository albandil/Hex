/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2013                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cmath>
#include <map>
#include <string>
#include <vector>

#include <sqlitepp/sqlitepp.hpp>

#include "../interpolate.h"
#include "../specf.h"
#include "../variables.h"
#include "../version.h"

const std::string DifferentialCrossSection::Id = "dcs";
const std::string DifferentialCrossSection::Description = "Differential cross section.";
const std::vector<std::string> DifferentialCrossSection::Dependencies = {
    "ni", "li", "mi", 
    "nf", "lf", "mf",
    "S",
    "Ei", "theta"
};
const std::vector<std::string> DifferentialCrossSection::VecDependencies = { "theta" };

bool DifferentialCrossSection::initialize(sqlitepp::session & db) const
{
    return true;
}

std::vector<std::string> const & DifferentialCrossSection::SQL_CreateTable() const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

std::vector<std::string> const & DifferentialCrossSection::SQL_Update() const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

rArray differential_cross_section(sqlitepp::session & db, int ni, int li, int mi, int nf, int lf, int mf, int S, double E, rArray const & angles)
{
    // the scattering amplitudes
    rArray dcs(angles.size());
    
    double ki = sqrt(E);
    double kf = sqrt(E - 1./(ni*ni) + 1./(nf*nf));
    
    // check if this is an allowed transition
    if (not finite(ki) or not finite(kf))
        return dcs;
    
    int ell;
    double Ei, sum_Re_T_ell, sum_Im_T_ell;
    rArrays E_ell;
    cArrays T_E_ell;
    
    // sum over L
    sqlitepp::statement st(db);
    st << "SELECT ell, Ei, SUM(Re_T_ell), SUM(Im_T_ell) "
          "FROM " + TMatrix::Id + " "
          "WHERE ni = :ni "
          "  AND li = :li "
          "  AND mi = :mi "
          "  AND nf = :nf "
          "  AND lf = :lf "
          "  AND mf = :mf "
          "  AND  S = :S  "
          "GROUP BY ell, Ei "
          "ORDER BY ell ASC, Ei ASC",
        sqlitepp::into(ell), sqlitepp::into(Ei),
        sqlitepp::into(sum_Re_T_ell), sqlitepp::into(sum_Im_T_ell),
        sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
        sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
        sqlitepp::use(S);
    
    // load data using the statement
    while (st.exec())
    {
        while ((int)E_ell.size() <= ell)
        {
            E_ell.push_back(rArray());
            T_E_ell.push_back(cArray());
        }
        
        E_ell[ell].push_back(Ei);
        T_E_ell[ell].push_back(Complex(sum_Re_T_ell,sum_Im_T_ell));
    }
    
    // terminate if no data
    if (E_ell.empty())
        return dcs;
    
    // for all angles
    for (int i = 0; i < (int)angles.size(); i++)
    {
        rArray e;    // energies
        cArray f;    // cross sections
        
        // for all projectile angular momenta sum arrays
        for (int l = 0; l < (int)E_ell.size(); l++)
            merge (e, f, E_ell[l], T_E_ell[l] * sphY(l,mi-mf,angles[i],0));
        
        // intepolate energies for unnormalized differential cross section
        if (not e.empty())
            dcs[i] = interpolate(e, sqrabs(f), {E})[0];
    }
    
    // normalize
    return dcs * kf * (2.*S + 1.) / (16 * M_PI * M_PI * ki);
}

bool DifferentialCrossSection::run (
    sqlitepp::session & db,
    std::map<std::string,std::string> const & sdata
) const {
    
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
    double lfactor = change_units(lUnit_au, Lunits);
    double afactor = change_units(Aunits, aUnit_rad);
    
    // scattering event parameters
    int ni = As<int>(sdata, "ni", Id);
    int li = As<int>(sdata, "li", Id);
    int mi = As<int>(sdata, "mi", Id);
    int nf = As<int>(sdata, "nf", Id);
    int lf = As<int>(sdata, "lf", Id);
    int mf = As<int>(sdata, "mf", Id);
    int  S = As<int>(sdata, "S", Id);
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
    
    // compute cross section
    rArray dcs = differential_cross_section(db, ni,li,mi, nf,lf,mf, S, E, angles * afactor);
    
    // write out
    std::cout << logo() <<
        "# Differential cross section in " << unit_name(Lunits) << " for \n" <<
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n" <<
        "#     S = " << S << ", E = " << E/efactor << " " << unit_name(Eunits)
                     << " ordered by angle in " << unit_name(Aunits) << "\n" <<
        "# θ\t dσ\n";
    for (size_t i = 0; i < angles.size(); i++)
        std::cout << angles[i] << "\t" << dcs[i]*lfactor*lfactor << "\n";
    
    return true;
}
