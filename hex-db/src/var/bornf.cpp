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

#include <gsl/gsl_sf.h>

#include "../arrays.h"
#include "../interpolate.h"
#include "../chebyshev.h"
#include "../special.h"
#include "../variables.h"
#include "../vec3d.h"
#include "../version.h"

const std::string BornFullTMatrix::Id = "bornf";
const std::string BornFullTMatrix::Description = "Full second Born T-matrix (angle dependent).";
const std::vector<std::string> BornFullTMatrix::Dependencies = {
    "ni", "li", "mi", 
    "nf", "lf", "mf", 
    "Ei", "theta"
};
const std::vector<std::string> BornFullTMatrix::VecDependencies = { "theta" };

bool BornFullTMatrix::initialize (sqlitepp::session & db) const
{
    return true;
}

std::vector<std::string> const & BornFullTMatrix::SQL_CreateTable () const
{
    static const std::vector<std::string> cmd = {
        "CREATE TABLE IF NOT EXISTS '" + BornFullTMatrix::Id + "' ("
            "ni INTEGER, "
            "li INTEGER, "
            "mi INTEGER, "
            "nf INTEGER, "
            "lf INTEGER, "
            "mf INTEGER, "
            "Ei DOUBLE PRECISION, "
            "cheb BLOB, "
            "PRIMARY KEY (ni,li,mi,nf,lf,mf,Ei)"
        ")"
    };
    
    return cmd;
}

std::vector<std::string> const & BornFullTMatrix::SQL_Update () const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

bool BornFullTMatrix::run
(
    sqlitepp::session & db,
    std::map<std::string,std::string> const & sdata
) const
{
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
    double lfactor = change_units(lUnit_au, Lunits);
    
    // atomic and projectile data
    int ni = As<int>(sdata, "ni", Id);
    int li = As<int>(sdata, "li", Id);
    int mi = As<int>(sdata, "mi", Id);
    int nf = As<int>(sdata, "nf", Id);
    int lf = As<int>(sdata, "lf", Id);
    int mf = As<int>(sdata, "mf", Id);
    double Ei = As<double>(sdata, "Ei", Id) * efactor;
    
    // read scattering angle (in user units)
    rArray angles;
    try {
        angles.push_back(As<double>(sdata, "theta", Id));
    } catch (std::exception e) {
        angles = readStandardInput<double>();
    }
    
    // energy and encoded Chebyshev approximation
    double E;
    std::string blob;
    
    // create query statement
    sqlitepp::statement st(db);
    st << "SELECT Ei, QUOTE(cheb) FROM " + BornFullTMatrix::Id + " "
          "WHERE ni = :ni "
          "  AND li = :li "
          "  AND mi = :mi "
          "  AND nf = :nf "
          "  AND lf = :lf "
          "  AND mf = :mf "
          "ORDER BY Ei ASC",
          
       sqlitepp::into(E), sqlitepp::into(blob),
       sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
       sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf);
    
    // get Chebyshev expansions
    rArray E_arr;
    cArray cb;
    std::vector<Chebyshev<double,Complex>> cheb_arr;
    while (st.exec())
    {
        // save energy
        E_arr.push_back(E);
        
        // decode Chebyshev expansion from hexadecimal format
        cb.fromBlob(blob);
        
        // save Chebyshev expansion
        cheb_arr.push_back
        (
            Chebyshev<double,Complex>
            (
                cb,                         // expansion coefficients
                0.,                         // smallest angle
                special::constant::pi_half  // largest angle
            )
        );
    }
    
    // terminate if no data
    if (E_arr.empty())
        return true;
    
    // for all angles
    cArray T_out(angles.size());
    for (std::size_t i = 0; i < angles.size(); i++)
    {
        // for all impact energies evaluate the T-matrix
        cArray Ti(E_arr.size());
        for (std::size_t ie = 0; ie < E_arr.size(); ie++)
            Ti[ie] = cheb_arr[ie].clenshaw(angles[i], cheb_arr[ie].tail(1e-8));
        
        // interpolate
        T_out[i] = interpolate(E_arr, Ti, { Ei })[0];
    }
    
    // write out
    std::cout << logo("#") <<
        "# Born T-matrix in " << unit_name(Lunits) << " for\n" <<
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n" <<
        "# and impact energy\n" <<
        "#     Ei = " << Ei << " in " << unit_name(Eunits) << "\n" <<
        "# ordered by angle in " << unit_name(Aunits) <<
        "# \n" <<
        "# angle\t Re TB\t Im TB\n";
    for (std::size_t i = 0; i < angles.size(); i++)
    {
        std::cout << 
            angles[i] << "\t" << 
            T_out[i].real()*lfactor << "\t" <<
            T_out[i].imag()*lfactor << "\n";
    }
    
    return true;
}
