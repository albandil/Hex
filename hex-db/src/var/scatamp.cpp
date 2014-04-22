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

std::vector<std::string> const & ScatteringAmplitude::SQL_Update() const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

std::vector<std::string> const & ScatteringAmplitude::SQL_CreateTable() const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

cArray Born_direct_scattering_amplitude (int ni, int li, int mi, int nf, int lf, int mf, int S, double E, rArray const & angles)
{
    // form-factors for all given angles
    cArray ffactor (angles.size());
    
    // initial momentum
    if (E < 0)
        return ffactor;
    double ki = sqrt(E);
    
    // final momentum
    if (E < 1./(ni*ni) - 1./(nf*nf))
        return ffactor;
    double kf = sqrt(E - 1./(ni*ni) + 1./(nf*nf));
    
    // transferred momentum for all given angles
    rArray qsqr = ki*ki + kf*kf - 2.*ki*kf*cos(angles);
    rArray q = sqrt(qsqr);
    rArray zeta = asin(kf * sin(angles) / q);
    
    // for all transferred angular momenta
    int lmin = std::abs(li-lf), lmax = li+lf;
    for (int l = lmin; l <= lmax; l++)
    {
        // for all transferred linear momenta
        for (int iang = 0; iang < (int)angles.size(); iang++)
        {
            // radial integral
            double g = 0;
            
            // for all terms of the Laguerre polynomials
            for (int i = 0; i < ni - li - 1; i++)
            for (int f = 0; f < nf - lf - 1; f++)
            {
                double k = q[iang];       // transferred linear momentum
                double sign = (((i + f) % 2 == 0) ? 1 : -1);
                double factor = sqrt(0.5 * M_PI / k) * 0.25 * std::pow(2./ni,li+2+i) * std::pow(2./nf,lf+2+f) *
                    std::sqrt(gsl_sf_fact(ni-li-1)*gsl_sf_fact(ni+li)*gsl_sf_fact(nf-lf-1)*gsl_sf_fact(nf+lf)) /
                    (gsl_sf_fact(ni-li-1-i)*gsl_sf_fact(2*li+1+i)*gsl_sf_fact(nf-lf-1-f)*gsl_sf_fact(2*lf+1+f)*
                    gsl_sf_fact(i)*gsl_sf_fact(f));
                double a = -1./nf -1./ni; // exponent from exp(-ar)
                double b = li+lf+1.5+i+f; // exponent from r^b
                double n = l + 0.5;       // Bessel function degree
                double hyper = Hyper2F1(0.5*(b+n+1), 0.5*(b+n+2), n+1, -(k*k)/(a*a));
                double integral = std::pow(0.5*k,n) * std::exp(gsl_sf_lngamma(b+n+1) - gsl_sf_lngamma(n+1) - (b+n+1)*log(a)) * hyper;
                
                g += sign * factor * integral;
            }
            
            // add contribution to the whole form-factor
            ffactor[iang] += 4 * M_PI * std::pow(Complex(0.,1.),l) * Gaunt(li,mi,lf,mf,l,mi-mf) * sphY(l,mi-mf,zeta[iang],0) * g;
        }
    }
    
    // add "no scattering" contribution in elastic case
    if (ni == nf and li == lf and mi == mf)
        ffactor -= 1.;
    
    // return the result
    return -2. * ffactor / qsqr;
}

cArray scattering_amplitude (sqlitepp::session & db, int ni, int li, int mi, int nf, int lf, int mf, int S, double E, rArray const & angles)
{
    cArray amplitudes(angles.size());
    
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
        for (int ell = abs(M - mf); ell <= max_ell; ell++)
        {
            // get all relevant lines from database
            double __Ei, __Re_T_ell, __Im_T_ell;
            rArray db_Ei;
            cArray db_T_ell;
            sqlitepp::statement st2(db);
            st2 << "SELECT Ei, Re_T_ell, Im_T_ell FROM " + TMatrix::Id + " "
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
                sqlitepp::into(__Ei), sqlitepp::into(__Re_T_ell), sqlitepp::into(__Im_T_ell),
                sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
                sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
                sqlitepp::use(ell), sqlitepp::use(L), sqlitepp::use(S);
            while ( st2.exec() )
            {
                db_Ei.push_back(__Ei);
                db_T_ell.push_back(Complex(__Re_T_ell, __Im_T_ell));
            }
            
            if (db_Ei.size() == 0)
                continue;
            
            // update value of "f"
            Complex Tmatrix = interpolate(db_Ei, db_T_ell, {E})[0];
            
            for (size_t i = 0; i < angles.size(); i++)
                amplitudes[i] += -1./(2.*M_PI) * Tmatrix * sphY(ell, abs(M-mf), angles[i], 0.);
        }
    }
    
    return amplitudes;
}

bool ScatteringAmplitude::run (
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
