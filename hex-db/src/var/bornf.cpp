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
const std::vector<std::pair<std::string,std::string>> BornFullTMatrix::Dependencies = {
    {"ni", "Initial atomic principal quantum number."},
    {"li", "Initial atomic orbital quantum number."},
    {"mi", "Initial atomic magnetic quantum number."},
    {"nf", "Final atomic principal quantum number."},
    {"lf", "Final atomic orbital quantum number."},
    {"mf", "Final atomic magnetic quantum number."},
    {"Ei", "Projectile impact energy (Rydberg)."},
    {"theta", "Scattering angles for which to compute the amplitude."}
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

bool BornFullTMatrix::run (std::map<std::string,std::string> const & sdata) const
{
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
    double lfactor = change_units(lUnit_au, Lunits);
    
    // atomic and projectile data
    int ni = Conv<int>(sdata, "ni", Id);
    int li = Conv<int>(sdata, "li", Id);
    int mi = Conv<int>(sdata, "mi", Id);
    int nf = Conv<int>(sdata, "nf", Id);
    int lf = Conv<int>(sdata, "lf", Id);
    int mf = Conv<int>(sdata, "mf", Id);
    double Ei = Conv<double>(sdata, "Ei", Id) * efactor;
    
    // read scattering angle (in user units)
    rArray angles;
    try {
        angles.push_back(Conv<double>(sdata, "theta", Id));
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
        "# \n";
    OutputTable table;
    table.setWidth(15);
    table.setAlignment(OutputTable::left);
    table.write("# angle    ", "Re TBorn ", "Im TBorn ");
    table.write("# ---------", "---------", "---------");
    for (std::size_t i = 0; i < angles.size(); i++)
    {
        table.write
        (
            angles[i],
            T_out[i].real()*lfactor,
            T_out[i].imag()*lfactor
        );
    }
    
    return true;
}
