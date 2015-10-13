//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2015, Jakub Benda, Charles University in Prague                    //
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

#include "hex-arrays.h"
#include "hex-interpolate.h"
#include "hex-version.h"

#include "variables.h"

const std::string TMatrix::Id = "tmat";
const std::string TMatrix::Description = "T-matrix.";
const std::vector<std::pair<std::string,std::string>> TMatrix::Dependencies = {
    {"ni", "Initial atomic principal quantum number."},
    {"li", "Initial atomic orbital quantum number."},
    {"mi", "Initial atomic magnetic quantum number."},
    {"nf", "Final atomic principal quantum number."},
    {"lf", "Final atomic orbital quantum number."},
    {"mf", "Final atomic magnetic quantum number."},
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
    int ni = Conv<int>(sdata, "ni", Id);
    int li = Conv<int>(sdata, "li", Id);
    int mi0= Conv<int>(sdata, "mi", Id);
    int nf = Conv<int>(sdata, "nf", Id);
    int lf = Conv<int>(sdata, "lf", Id);
    int mf0= Conv<int>(sdata, "mf", Id);
    int  S = Conv<int>(sdata,  "S", Id);
    int ell= Conv<int>(sdata, "ell",Id);
    
    // use mi >= 0; if mi < 0, flip both signs
    int mi = (mi0 < 0 ? -mi0 : mi0);
    int mf = (mi0 < 0 ? -mf0 : mf0);
    
    // energies
    rArray energies;
    
    // get energy / energies
    try {
        
        // is there a single energy specified using command line ?
        energies.push_back(Conv<double>(sdata, "Ei", Id));
        
    } catch (std::exception e) {
        
        // are there more energies specified using the STDIN ?
        energies = readStandardInput<double>();
    }
    
    // energy and real and imarinary part of the T-matrix
    double E, Re_T_ell, Im_T_ell, Re_TBorn_ell, Im_TBorn_ell;
    
    // create query statement
    sqlitepp::statement st(db);
    st << "SELECT Ei, SUM(Re_T_ell), SUM(Im_T_ell), SUM(Re_TBorn_ell), SUM(Im_TBorn_ell) FROM " + TMatrix::Id + " "
          "WHERE ni = :ni "
          "  AND li = :li "
          "  AND mi = :mi "
          "  AND nf = :nf "
          "  AND lf = :lf "
          "  AND mf = :mf "
          "  AND  S = :S  "
          "  AND ell=:ell "
          "GROUP BY L, Ei "
          "ORDER BY Ei ASC",
       sqlitepp::into(E),
       sqlitepp::into(Re_T_ell), sqlitepp::into(Im_T_ell),
       sqlitepp::into(Re_TBorn_ell), sqlitepp::into(Im_TBorn_ell),
       sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
       sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
       sqlitepp::use(S), sqlitepp::use(ell);
    
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
    {
        std::cout << "No data fot this selection." << std::endl;
        return true;
    }
    
    // get T-matrices
    cArray T_out, Tb_out;
    if (energies.size() > 0 and energies[0] < 0)
    {
        // use all
        energies = E_arr;
        T_out = T_arr;
        Tb_out = Tb_arr;
    }
    else
    {
        // interpolate
        T_out = interpolate(E_arr, T_arr, energies * efactor);
        Tb_out = interpolate(E_arr, Tb_arr, energies * efactor);
    }
    
    // write out
    std::cout << logo("#") <<
        "# T-matrices in " << unit_name(Lunits) << " for\n" <<
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n" <<
        "#     S = " << S << ", ℓ = " << ell << "\n" <<
        "# ordered by energy in " << unit_name(Eunits) << "\n" <<
        "# \n";
    OutputTable table;
    table.setWidth(15, 15, 15, 15);
    table.setAlignment(OutputTable::left);
    table.write("# E        ", "Re T     ", "Im T     ", "Re TBorn ", "Im TBorn ");
    table.write("# ---------", "---------", "---------", "---------", "---------");
    for (std::size_t i = 0; i < energies.size(); i++)
    {
        table.write
        (
            energies[i],
            T_out[i].real()*lfactor,  T_out[i].imag()*lfactor,
            Tb_out[i].real()*lfactor, Tb_out[i].imag()*lfactor
        );
    }
    
    return true;
}
