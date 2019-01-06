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

#include <gsl/gsl_sf.h>

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"
#include "hex-interpolate.h"
#include "hex-chebyshev.h"
#include "hex-vec3d.h"
#include "hex-version.h"

// --------------------------------------------------------------------------------- //

#include "../quantities.h"
#include "../utils.h"

// --------------------------------------------------------------------------------- //

createNewScatteringQuantity(IonizationF, "ionf")

// --------------------------------------------------------------------------------- //

std::string IonizationF::description ()
{
    return "Ionization amplitude radial part.";
}

std::vector<std::string> IonizationF::dependencies ()
{
    return std::vector<std::string>();
}

std::vector<std::pair<std::string,std::string>> IonizationF::params ()
{
    return std::vector<std::pair<std::string,std::string>>
    {
        {"ni", "Initial atomic principal quantum number."},
        {"li", "Initial atomic orbital quantum number."},
        {"mi", "Initial atomic magnetic quantum number."},
        {"L", "Total orbital momentum of atomic + projectile electron."},
        {"S", "Total spin of atomic + projectile electron."},
        {"Ei", "Projectile impact energy (Rydberg)."},
        {"l1", "Atomic electron orbital momentum in the final state."},
        {"l2", "Projectile orbital momentum in the final state."},
        {"Eshare", "Energy fraction (atomic vs projectile electron) in the final state."}
    };
}

std::vector<std::string> IonizationF::vparams ()
{
    return std::vector<std::string>
    {
        "Eshare"
    };
}

// --------------------------------------------------------------------------------- //

bool IonizationF::initialize (sqlitepp::session & db)
{
    return ScatteringQuantity::initialize(db);
}

bool IonizationF::createTable ()
{
    sqlitepp::statement st (session());
    st <<
        "CREATE TABLE IF NOT EXISTS 'ionf' ("
            "ni INTEGER, "
            "li INTEGER, "
            "mi INTEGER, "
            "L  INTEGER, "
            "S  INTEGER, "
            "Ei DOUBLE PRECISION, "
            "l1 INTEGER, "
            "l2 INTEGER, "
            "cheb BLOB, "
            "PRIMARY KEY (ni,li,mi,L,S,Ei,l1,l2)"
        ")";

    try
    {
        st.exec();
    }
    catch (sqlitepp::exception & e)
    {
        std::cerr << "ERROR: Creation of table 'ionf' failed!" << std::endl;
        std::cerr << "       code = " << e.code() << " (\"" << e.what() << "\")" << std::endl;
        return false;
    }

    return ScatteringQuantity::createTable();
}

bool IonizationF::updateTable ()
{
    return ScatteringQuantity::updateTable();
}

// --------------------------------------------------------------------------------- //

bool IonizationF::run (std::map<std::string,std::string> const & sdata)
{
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
    double lfactor = change_units(lUnit_au, Lunits);

    // atomic and projectile data
    int ni = Conv<int>(sdata, "ni", name());
    int li = Conv<int>(sdata, "li", name());
    int mi = Conv<int>(sdata, "mi", name());
    int  L = Conv<int>(sdata,  "L", name());
    int  S = Conv<int>(sdata,  "S", name());
    int l1 = Conv<int>(sdata, "l1", name());
    int l2 = Conv<int>(sdata, "l2", name());
    double Ei = Conv<double>(sdata, "Ei", name()) * efactor;

    if (mi < 0)
        mi = -mi;

    // read energy sharing (in user units)
    rArray Eshare;
    try
    {
        Eshare.push_back(Conv<double>(sdata, "Eshare", name()));
    }
    catch (std::exception e)
    {
        Eshare = readStandardInput<double>();
    }

    // energy and encoded Chebyshev approximation
    double E;
    std::string blob;

    // create query statement
    sqlitepp::statement st (session());
    st << "SELECT Ei, QUOTE(cheb) FROM 'ionf' "
          "WHERE ni = :ni "
          "  AND li = :li "
          "  AND mi = :mi "
          "  AND  L = :L  "
          "  AND  S = :S  "
          "  AND l1 = :l1 "
          "  AND l2 = :l2 "
          "ORDER BY Ei ASC",
       sqlitepp::into(E), sqlitepp::into(blob),
       sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
       sqlitepp::use(L), sqlitepp::use(S),
       sqlitepp::use(l1), sqlitepp::use(l2);

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
                cb,                   // expansion coefficients
                0.,                   // lowest energy
                sqrt(E - 1./(ni*ni))  // highest energy
            )
        );
    }

    // terminate if no data
    if (E_arr.empty())
        return true;

    // for all energy shares
    cArray f_out(Eshare.size());
    for (std::size_t i = 0; i < Eshare.size(); i++)
    {
        // initialize k₁ and k₂ so that
        //   1) (k₁)² + (k₂)² = Ei
        //   2) (k₁)² / (k₂)² = Eshare / (1 - Eshare)
        rArray k1 = sqrt((E_arr - 1./(ni*ni)) * Eshare[i]);
        rArray k2 = sqrt((E_arr - 1./(ni*ni)) * (1 - Eshare[i]));

        // for all impact energies evaluate the radial part
        cArray f0(E_arr.size());
        for (std::size_t ie = 0; ie < E_arr.size(); ie++)
            f0[ie] = cheb_arr[ie].clenshaw(k1[ie], cheb_arr[ie].tail(1e-8)) / gsl_sf_pow_int(k1[ie] * k2[ie], 2);

        // interpolate
        f_out[i] = interpolate(E_arr, f0, { Ei })[0];
    }

    // write out
    std::cout << logo("#") <<
        "# Ionization amplitudes (radial part) in " << unit_name(Lunits) << " for\n" <<
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
        "#     L = " << L << ", S = " << S << ", ℓ₁ = " << l1 << ", ℓ₂ = " << l2 << "\n" <<
        "# and impact energy\n" <<
        "#     Ei = " << Ei << " in " << unit_name(Eunits) << "\n" <<
        "# ordered by energy share " << 
        "# \n";
    OutputTable table;
    table.setWidth(15, 15);
    table.setAlignment(OutputTable::left);
    table.write("# Eshare   ", "Re f     ", "Im f     ");
    table.write("# ---------", "---------", "---------");

    for (std::size_t i = 0; i < Eshare.size(); i++)
    {
        table.write
        (
            Eshare[i],
            f_out[i].real()*lfactor,
            f_out[i].imag()*lfactor
        );
    }

    return true;
}
