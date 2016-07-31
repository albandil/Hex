//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2016, Jakub Benda, Charles University in Prague                    //
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

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"
#include "hex-interpolate.h"
#include "hex-version.h"

// --------------------------------------------------------------------------------- //

#include "../quantities.h"
#include "../utils.h"

// --------------------------------------------------------------------------------- //

createNewScatteringQuantity(TMatrix);

// --------------------------------------------------------------------------------- //

std::string TMatrix::name ()
{
    return "tmat";
}

std::string TMatrix::description ()
{
    return "T-matrix.";
}

std::vector<std::string> TMatrix::dependencies ()
{
    return std::vector<std::string>();
}

std::vector<std::pair<std::string,std::string>> TMatrix::params ()
{
    return std::vector<std::pair<std::string,std::string>>
    {
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
}

std::vector<std::string> TMatrix::vparams ()
{
    return std::vector<std::string>
    {
        "Ei"
    };
}

// --------------------------------------------------------------------------------- //

bool TMatrix::initialize (sqlitepp::session & db)
{
    return ScatteringQuantity::initialize(db);
}

bool TMatrix::createTable ()
{
    sqlitepp::statement st (session());
    st <<
        "CREATE TABLE IF NOT EXISTS 'tmat' ("
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
            "PRIMARY KEY (ni,li,mi,nf,lf,mf,L,S,Ei,ell)"
        ")";
    
    try
    {
        st.exec();
    }
    catch (sqlitepp::exception & e)
    {
        std::cerr << "ERROR: Creation of table 'tmat' failed!" << std::endl;
        std::cerr << "       code = " << e.code() << " (\"" << e.what() << "\")" << std::endl;
        return false;
    }
    
    return ScatteringQuantity::createTable();
}

bool TMatrix::updateTable ()
{
    return ScatteringQuantity::updateTable();
}

// --------------------------------------------------------------------------------- //

bool TMatrix::run (std::map<std::string,std::string> const & sdata)
{
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
    double lfactor = change_units(lUnit_au, Lunits);
    
    // atomic and projectile data
    int ni = Conv<int>(sdata, "ni", name());
    int li = Conv<int>(sdata, "li", name());
    int mi0= Conv<int>(sdata, "mi", name());
    int nf = Conv<int>(sdata, "nf", name());
    int lf = Conv<int>(sdata, "lf", name());
    int mf0= Conv<int>(sdata, "mf", name());
    int  S = Conv<int>(sdata, "S",  name());
    int ell= Conv<int>(sdata, "ell",name());
    
    // use mi >= 0; if mi < 0, flip both signs
    int mi = (mi0 < 0 ? -mi0 : mi0);
    int mf = (mi0 < 0 ? -mf0 : mf0);
    
    // energies
    rArray energies;
    
    // get energy / energies
    try
    {
        // is there a single energy specified using command line ?
        energies.push_back(Conv<double>(sdata, "Ei", name()));
    }
    catch (std::exception e)
    {
        // are there more energies specified using the STDIN ?
        energies = readStandardInput<double>();
    }
    
    // energy and real and imarinary part of the T-matrix
    double E, Re_T_ell, Im_T_ell;
    
    // create query statement
    sqlitepp::statement st (session());
    st << "SELECT Ei, SUM(Re_T_ell), SUM(Im_T_ell) FROM 'tmat' "
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
       sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
       sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
       sqlitepp::use(S), sqlitepp::use(ell);
    
    // get T-matrices
    rArray E_arr;
    cArray T_arr;
    while (st.exec())
    {
        E_arr.push_back(E);
        T_arr.push_back(Complex(Re_T_ell,Im_T_ell));
    }
    
    // terminate if no data
    if (E_arr.empty())
    {
        std::cout << "No data fot this selection." << std::endl;
        return true;
    }
    
    // get T-matrices
    cArray T_out;
    if (energies.size() > 0 and energies[0] < 0)
    {
        // use all
        energies = E_arr;
        T_out = T_arr;
    }
    else
    {
        // interpolate
        T_out = interpolate(E_arr, T_arr, energies * efactor);
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
    table.write("# E        ", "Re T     ", "Im T     ");
    table.write("# ---------", "---------", "---------");
    for (std::size_t i = 0; i < energies.size(); i++)
    {
        table.write
        (
            energies[i],
            T_out[i].real()*lfactor,  T_out[i].imag()*lfactor
        );
    }
    
    return true;
}
