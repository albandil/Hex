//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2014, Jakub Benda, Charles University in Prague                    //
//                                                                                   //
// MIT License:                                                                      //
//                                                                                   //
//  Permission is hereby granted, free of charge, to any person obtaining a          //
// copy of this software and associated documentation files (the "Software"),        //
// to deal in the Software without restriction, including without limitation         //
// the rights to use, copy, modify, merge, publish, distribute, sublicense,          //
// and/or sell copies of the Software, and to permit persons to whom the             //
// Software is furnished to do so, subject to the following conditions:              //
//                                                                                   //
//  The above copyright notice and this permission notice shall be included          //
// in all copies or substantial portions of the Software.                            //
//                                                                                   //
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS          //
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       //
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE       //
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, //
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF         //
// OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.  //
//                                                                                   //
//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //

#include <map>
#include <string>
#include <vector>

#include "../arrays.h"
#include "../interpolate.h"
#include "../variables.h"
#include "../version.h"

const std::string TMatrix::Id = "tmat";
const std::string TMatrix::Description = "T-matrix.";
const std::vector<std::pair<std::string,std::string>> TMatrix::Dependencies = {
    {"ni", "Initial atomic principal quantum number."},
    {"li", "Initial atomic orbital quantum number."},
    {"mi", "Initial atomic magnetic quantum number."},
    {"nf", "Final atomic principal quantum number."},
    {"lf", "Final atomic orbital quantum number."},
    {"mf", "Final atomic magnetic quantum number."},
    {"L", "Total orbital momentum of atomic + projectile electron."},
    {"S", "Total spin of atomic + projectile electron."},
    {"Ei", "Projectile impact energy (Rydberg)."},
    {"ell", "Outgoing projectile partial wave angular momentum."}
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

bool TMatrix::run (std::map<std::string,std::string> const & sdata) const
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
