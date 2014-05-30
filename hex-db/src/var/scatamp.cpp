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

#include "../interpolate.h"
#include "../chebyshev.h"
#include "../special.h"
#include "../variables.h"
#include "../version.h"

const std::string ScatteringAmplitude::Id = "scatamp";
const std::string ScatteringAmplitude::Description = "Scattering amplitude.";
const std::vector<std::string> ScatteringAmplitude::Dependencies = {
    "ni", "li", "mi", 
    "nf", "lf", "mf",
    "S", "Ei", "theta"
};
const std::vector<std::string> ScatteringAmplitude::VecDependencies = { "theta" };

bool ScatteringAmplitude::initialize(sqlitepp::session & db) const
{
    return true;
}

std::vector<std::string> const & ScatteringAmplitude::SQL_Update () const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

std::vector<std::string> const & ScatteringAmplitude::SQL_CreateTable () const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

cArray scattering_amplitude (sqlitepp::session & db, int ni, int li, int mi, int nf, int lf, int mf, int S, double E, rArray const & angles)
{
    cArray amplitudes(angles.size()), amplitudesb(angles.size()), bornf(angles.size());
    
    // total angular momentum projection (given by axis orientation)
    int M = mi;
    
    // get maximal partial wave angular momentum
    int max_L;
    sqlitepp::statement st3(db);
    st3 << "SELECT MAX(L) FROM " + TMatrix::Id + " "
           "WHERE ni = :ni "
           "  AND li = :li "
           "  AND mi = :mi "
           "  AND nf = :nf "
           "  AND lf = :lf "
           "  AND mf = :mf "
           "  AND  S = :S  ",
        sqlitepp::into(max_L),
        sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
        sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
        sqlitepp::use(S);
    st3.exec();
    
    //
    // sum the partial wave contributions for exact and Born partial waves
    //
    
    // for all total angular momenta
    for (int L = M; L <= max_L; L++)
    {
        // get maximal projectile partial wave angular momentum
        int max_ell;
        sqlitepp::statement st1(db);
        st1 << "SELECT MAX(ell) FROM " + TMatrix::Id + " "
               "WHERE ni = :ni "
               "  AND li = :li "
               "  AND mi = :mi "
               "  AND nf = :nf "
               "  AND lf = :lf "
               "  AND mf = :mf "
               "  AND  L = :L  "
               "  AND  S = :S  ",
            sqlitepp::into(max_ell),
            sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
            sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
            sqlitepp::use(L), sqlitepp::use(S);
        st1.exec();
        
        // for all outgoing partial waves
        for (int ell = std::abs(M - mf); ell <= max_ell; ell++)
        {
            // get all relevant lines from database
            double Ei, Re_T_ell, Im_T_ell, Re_TBorn_ell, Im_TBorn_ell;
            rArray db_Ei;
            cArray db_T_ell, db_TBorn_ell;
            std::string cheb;
            sqlitepp::statement st2(db);
            st2 << "SELECT Ei, Re_T_ell, Im_T_ell, Re_TBorn_ell, Im_TBorn_ell FROM " + TMatrix::Id + " "
                   "WHERE ni = :ni "
                   "  AND li = :li "
                   "  AND mi = :mi "
                   "  AND nf = :nf "
                   "  AND lf = :lf "
                   "  AND mf = :mf "
                   "  AND ell=:ell "
                   "  AND  L = :L  "
                   "  AND  S = :S  "
                   "ORDER BY Ei ASC",
                sqlitepp::into(Ei),
                sqlitepp::into(Re_T_ell), sqlitepp::into(Im_T_ell),
                sqlitepp::into(Re_TBorn_ell), sqlitepp::into(Im_TBorn_ell),
                sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
                sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
                sqlitepp::use(ell), sqlitepp::use(L), sqlitepp::use(S);
            while ( st2.exec() )
            {
                db_Ei.push_back(Ei);
                db_T_ell.push_back(Complex(Re_T_ell, Im_T_ell));
                db_TBorn_ell.push_back(Complex(Re_TBorn_ell, Im_TBorn_ell));
            }
            
            // skip empty queries
            if (db_Ei.size() == 0)
                continue;
            
            // update value of "f"
            Complex Tmatrix  = interpolate(db_Ei, db_T_ell, {E})[0];
            Complex Tmatrixb = (db_TBorn_ell.size() > 0 ? interpolate(db_Ei, db_TBorn_ell, {E})[0] : 0.);
            for (size_t i = 0; i < angles.size(); i++)
            {
                Complex Y = -0.5*special::constant::pi_inv * sphY(ell, std::abs(M-mf), angles[i], 0.);
                
                amplitudes[i]  += Tmatrix  * Y;
                amplitudesb[i] += Tmatrixb * Y;
            }
        }
    }
    
    //
    // load angle-dependent un-expanded Born T-matrix
    //
    
    rArray db_Ei;
    std::vector<Chebyshev<double,Complex>> db_bornf;
    double Ei;
    std::string cheb;
    sqlitepp::statement st4(db);
    st4 << "SELECT Ei, cheb FROM " + BornFullTMatrix::Id + " "
           "WHERE ni = :ni "
           "  AND li = :li "
           "  AND mi = :mi "
           "  AND nf = :nf "
           "  AND lf = :lf "
           "  AND mf = :mf ",
        sqlitepp::into(Ei), sqlitepp::into(cheb),
        sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
        sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf);
    while ( st4.exec() )
    {
        // add new energy
        db_Ei.push_back(Ei);
        
        // decode Chebyshev expansion
        cArray coeffs;
        coeffs.fromBlob(cheb);
        db_bornf.push_back(Chebyshev<double,Complex>(coeffs, -1., 1.));
    }
    for (size_t i = 0; i < angles.size(); i++)
    {
        // evaluate all energies for this angle
        cArray TE (db_Ei.size());
        for (unsigned j = 0; j < db_Ei.size(); j++)
            TE[i] = db_bornf[i].clenshaw(std::cos(angles[i]), db_bornf[i].tail(1e-10));
        
        // interpolate
        bornf[i] = -0.5*special::constant::pi_inv * (TE.size() > 0 ? interpolate(db_Ei, TE, { E })[0] : 0.);
    }
    
    //
    // return corrected expansion
    //
    
    return bornf + (amplitudes - amplitudesb);
}

bool ScatteringAmplitude::run
(
    sqlitepp::session & db,
    std::map<std::string,std::string> const & sdata
) const
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
    
    // the scattering amplitudes
    cArray amplitudes = scattering_amplitude(db, ni,li,mi, nf,lf,mf, S, E, angles * afactor);
    
    // write out
    std::cout << logo() <<
        "# Scattering amplitudes in " << unit_name(Lunits) << " for\n"
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n"
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n"
        "#     S = " << S << ", E = " << E/efactor << unit_name(Eunits) << "\n"
        "# ordered by angle in " << unit_name(Aunits) << "\n"
        "# \n"
        "# θ\t Re f\t Im f\n";
    for (size_t i = 0; i < angles.size(); i++)
    {
        std::cout << 
            angles[i] << "\t" << 
            amplitudes[i].real() * lfactor << "\t" <<
            amplitudes[i].imag() * lfactor << "\n";
    }
    
    return true;
}
