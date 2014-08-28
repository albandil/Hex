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

#include "../arrays.h"
#include "../interpolate.h"
#include "../variables.h"
#include "../version.h"

const std::string TMatrix::Id = "tmat";
const std::string TMatrix::Description = "T-matrix.";
const std::vector<std::string> TMatrix::Dependencies = {
    "ni", "li", "mi", 
    "nf", "lf", "mf",
    "L", "S",
    "Ei", "ell"
};
const std::vector<std::string> TMatrix::VecDependencies = { "Ei" };

bool TMatrix::initialize (sqlitepp::session & db) const
{
    return true;
}

std::vector<std::string> const & TMatrix::SQL_CreateTable () const
{
    static const std::vector<std::string> cmd = {
        "CREATE TABLE IF NOT EXISTS '" + TMatrix::Id + "' ("
            "ni INTEGER, "
            "li INTEGER, "
            "mi INTEGER, "
            "nf INTEGER, "
            "lf INTEGER, "
            "mf INTEGER, "
            "L  INTEGER, "
            "S  INTEGER, "
            "Ei DOUBLE PRECISION, "
            "ell INTEGER, "
            "Re_T_ell DOUBLE PRECISION, "
            "Im_T_ell DOUBLE PRECISION, "
            "Re_TBorn_ell DOUBLE PRECISION DEFAULT 0, "
            "Im_TBorn_ell DOUBLE PRECISION DEFAULT 0, "
            "PRIMARY KEY (ni,li,mi,nf,lf,mf,L,S,Ei,ell)"
        ")"
    };
    return cmd;
}

std::vector<std::string> const & TMatrix::SQL_Update () const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

bool TMatrix::run
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
    int  L = As<int>(sdata,  "L", Id);
    int  S = As<int>(sdata,  "S", Id);
    int ell= As<int>(sdata, "ell",Id);
    
    // energies
    rArray energies;
    
    // get energy / energies
    try {
        
        // is there a single energy specified using command line ?
        energies.push_back(As<double>(sdata, "Ei", Id));
        
    } catch (std::exception e) {
        
        // are there more energies specified using the STDIN ?
        energies = readStandardInput<double>();
    }
    
    // energy and real and imarinary part of the T-matrix
    double E, Re_T_ell, Im_T_ell, Re_TBorn_ell, Im_TBorn_ell;
    
    // create query statement
    sqlitepp::statement st(db);
    st << "SELECT Ei, Re_T_ell, Im_T_ell, Re_TBorn_ell, Im_TBorn_ell, FROM " + TMatrix::Id + " "
          "WHERE ni = :ni "
          "  AND li = :li "
          "  AND mi = :mi "
          "  AND nf = :nf "
          "  AND lf = :lf "
          "  AND mf = :mf "
          "  AND  L = :L  "
          "  AND  S = :S  "
          "  AND ell=:ell "
          "ORDER BY Ei ASC",
       sqlitepp::into(E),
       sqlitepp::into(Re_T_ell), sqlitepp::into(Im_T_ell),
       sqlitepp::into(Re_TBorn_ell), sqlitepp::into(Im_TBorn_ell),
       sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
       sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
       sqlitepp::use(L), sqlitepp::use(S), sqlitepp::use(ell);
    
    // get T-matrices
    rArray E_arr;
    cArray T_arr, Tb_arr;
    while (st.exec())
    {
        E_arr.push_back(E);
        T_arr.push_back(Complex(Re_T_ell,Im_T_ell));
        Tb_arr.push_back(Complex(Re_TBorn_ell,Im_TBorn_ell));
    }
    
    // terminate if no data
    if (E_arr.empty())
        return true;
    
    // interpolate
    cArray T_out = interpolate(E_arr, T_arr, energies * efactor);
    cArray Tb_out = interpolate(E_arr, Tb_arr, energies * efactor);
    
    // write out
    std::cout << logo("#") <<
        "# T-matrices in " << unit_name(Lunits) << " for\n" <<
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n" <<
        "#     L = " << L << ", S = " << S << ", ℓ = " << ell << "\n" <<
        "# ordered by energy in " << unit_name(Eunits) << "\n" <<
        "# \n" <<
        "# E\tRe T\tIm T\tRe TBorn\tIm TBorn\n";
    for (std::size_t i = 0; i < energies.size(); i++)
    {
        std::cout << energies[i] << "\t" << 
            T_out[i].real()*lfactor << "\t" <<
            T_out[i].imag()*lfactor << "\t" <<
            Tb_out[i].real()*lfactor << "\t" <<
            Tb_out[i].imag()*lfactor << std::endl;
    }
    
    return true;
}
