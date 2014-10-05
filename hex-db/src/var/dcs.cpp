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

#include <cmath>
#include <map>
#include <string>
#include <vector>

#include <sqlitepp/sqlitepp.hpp>

#include "../chebyshev.h"
#include "../interpolate.h"
#include "../special.h"
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

bool DifferentialCrossSection::initialize (sqlitepp::session & db) const
{
    return true;
}

std::vector<std::string> const & DifferentialCrossSection::SQL_CreateTable () const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

std::vector<std::string> const & DifferentialCrossSection::SQL_Update () const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

void hex_differential_cross_section_
(
    int * ni, int * li, int * mi,
    int * nf, int * lf, int * mf,
    int * S, double * E, int * N,
    double * angles, double * dcs
)
{
    double ki = std::sqrt((*E));
    double kf = std::sqrt((*E) - 1./((*ni)*(*ni)) + 1./((*nf)*(*nf)));
    
    // check if this is an allowed transition
    if (not std::isfinite(ki) or not std::isfinite(kf))
        return;
    
    int ell;
    double Ei, sum_Re_T_ell, sum_Im_T_ell, sum_Re_TBorn_ell, sum_Im_TBorn_ell;
    rArrays E_ell;
    cArrays T_E_ell, Tb_E_ell;
    
    //
    // load partial wave contributions
    //
    
    // sum over L
    sqlitepp::statement st(db);
    st << "SELECT ell, Ei, SUM(Re_T_ell), SUM(Im_T_ell), SUM(Re_TBorn_ell), SUM(Im_TBorn_ell) "
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
        sqlitepp::into(sum_Re_TBorn_ell), sqlitepp::into(sum_Im_TBorn_ell),
        sqlitepp::use(*ni), sqlitepp::use(*li), sqlitepp::use(*mi),
        sqlitepp::use(*nf), sqlitepp::use(*lf), sqlitepp::use(*mf),
        sqlitepp::use(*S);
    
    // load data using the statement
    while (st.exec())
    {
        while ((int)E_ell.size() <= ell)
        {
            E_ell.push_back(rArray());
            T_E_ell.push_back(cArray());
            Tb_E_ell.push_back(cArray());
        }
        
        E_ell[ell].push_back(Ei);
        T_E_ell[ell].push_back(Complex(sum_Re_T_ell,sum_Im_T_ell));
        Tb_E_ell[ell].push_back(Complex(sum_Re_TBorn_ell,sum_Im_TBorn_ell));
    }
    
    // terminate if no data
    if (E_ell.empty())
        return;
    
    //
    // load angle dependent Born T-matrices
    //
    
    std::string cheb;
    rArray E_arr;
    std::vector<Chebyshev<double,Complex>> bornT_arr;
    sqlitepp::statement stb(db);
    stb << "SELECT Ei, cheb FROM " + BornFullTMatrix::Id + " "
           "WHERE ni = :ni "
           "  AND li = :li "
           "  AND mi = :mi "
           "  AND nf = :nf "
           "  AND lf = :lf "
           "  AND mf = :mf "
           "ORDER BY Ei ASC",
        sqlitepp::into(Ei), sqlitepp::into(cheb),
        sqlitepp::use(*ni), sqlitepp::use(*li), sqlitepp::use(*mi),
        sqlitepp::use(*nf), sqlitepp::use(*lf), sqlitepp::use(*mf);
    while (stb.exec())
    {
        // add new energy
        E_arr.push_back(Ei);
        
        // decode Chebyshev expansion
        cArray coeffs;
        coeffs.fromBlob(cheb);
        bornT_arr.push_back(Chebyshev<double,Complex>(coeffs, -1, 1));
    }
    
    //
    // compute corrected cross sections
    //
    
    // for all angles
    for (int i = 0; i < (*N); i++)
    {
        rArray e;       // energies
        cArray f,fb,fB; // amplitudes
        
        // for all projectile angular momenta : sum partial wave arrays
        for (int l = 0; l < (int)E_ell.size(); l++)
        {
            Complex Y = special::sphY(l,mi-mf,angles[i],0);
            merge (e, f, E_ell[l], T_E_ell[l] * Y);
            merge (e, fb, E_ell[l], Tb_E_ell[l] * Y);
        }
        
        // skip empty cases
        if (e.empty())
            continue;
        
        // also, determine the angle-dependent Born T-matrix
        cArray fB0;
        double cosTheta = std::cos(angles[i]);
        
        // evaluate Chebyshev expansions for this angle
        for (unsigned j = 0; j < E_arr.size(); j++)
            fB0.push_back(bornT_arr[j].clenshaw(cosTheta, bornT_arr[j].tail(1e-8)));
        
        // interpolate angle-dependent Born T-matrices to partial wave energies ('e')
        fB = interpolate(E_arr, fB0, e);
        
        // intepolate corrected differential cross section
        dcs[i] = interpolate(e, sqrabs(fB + f - fb), { *E })[0] * kf * (2 * (*S) + 1) / (std::pow(4 * special::constant::pi, 2) * ki);
    }
}

bool DifferentialCrossSection::run (std::map<std::string,std::string> const & sdata) const
{
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
    rArray scaled_angles = angles * afactor, dcs(angles.size());
    hex_differential_cross_section(ni,li,mi, nf,lf,mf, S, E, angles.size(), scaled_angles.data(), dcs.data());
    
    // write out
    std::cout << logo("#") <<
        "# Differential cross section in " << unit_name(Lunits) << " for \n" <<
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n" <<
        "#     S = " << S << ", E = " << E/efactor << " " << unit_name(Eunits)
                     << " ordered by angle in " << unit_name(Aunits) << "\n" <<
        "# θ\t dσ\n";
    for (std::size_t i = 0; i < angles.size(); i++)
        std::cout << angles[i] << "\t" << dcs[i]*lfactor*lfactor << "\n";
    
    return true;
}
