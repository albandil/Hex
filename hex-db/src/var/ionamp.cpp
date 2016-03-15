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

#include "hex-arrays.h"
#include "hex-interpolate.h"
#include "hex-chebyshev.h"
#include "hex-special.h"
#include "hex-vec3d.h"
#include "hex-version.h"

#include "variables.h"

const std::string IonizationAmplitude::Id = "ionamp";
const std::string IonizationAmplitude::Description = "Ionization amplitude.";
const std::vector<std::pair<std::string,std::string>> IonizationAmplitude::Dependencies = {
    {"ni", "Initial atomic principal quantum number."},
    {"li", "Initial atomic orbital quantum number."},
    {"mi", "Initial atomic magnetic quantum number."},
    {"S", "Total spin of atomic + projectile electron."},
    {"Ei", "Projectile impact energy (Rydberg)."},
    {"dirs", "List of pairs of energy share and coordinate triplets in the, like this: '(ε₁,θ₁,φ₁) (ε₂,θ₂,φ₂)'."}
};
const std::vector<std::string> IonizationAmplitude::VecDependencies = { "dirs" };

bool IonizationAmplitude::initialize(sqlitepp::session & db) const
{
    return true;
}

std::vector<std::string> const & IonizationAmplitude::SQL_CreateTable () const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

std::vector<std::string> const & IonizationAmplitude::SQL_Update () const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

bool IonizationAmplitude::run (std::map<std::string,std::string> const & sdata) const
{
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
    double lfactor = change_units(lUnit_au, Lunits);
    double afactor = change_units(Aunits, aUnit_rad);
    
    // atomic and projectile data
    int ni = Conv<int>(sdata, "ni", Id);
    int li = Conv<int>(sdata, "li", Id);
    int mi = Conv<int>(sdata, "mi", Id);
    int  S = Conv<int>(sdata,  "S", Id);
    double Ei = Conv<double>(sdata, "Ei", Id) * efactor;
    
    // read directions
    //  dirs.first  = ( theta1, phi1, E1frac )
    //  dirs.second = ( theta2, phi2, E2frac )
    // NOTE: energy fractions will be normalized to become on-shell
    std::vector<std::pair<geom::vec3d,geom::vec3d>> dirs;
    try {
        dirs.push_back(Conv<std::pair<geom::vec3d,geom::vec3d>>(sdata, "dirs", Id));
    } catch (exception e) {
        dirs = readStandardInput<std::pair<geom::vec3d,geom::vec3d>>();
    }
    
    // energy and encoded Chebyshev approximation
    double E;
    std::string blob;
    
    // angular variables
    int L, l1, l2;
    
    // create query statement
    sqlitepp::statement st(db);
    st << "SELECT L, l1, l2, Ei, QUOTE(cheb) FROM " + IonizationF::Id + " "
          "WHERE ni = :ni "
          "  AND li = :li "
          "  AND mi = :mi "
          "  AND  S = :S  "
          "ORDER BY Ei ASC",
       sqlitepp::into(L), sqlitepp::into(l1), sqlitepp::into(l2),
       sqlitepp::into(E), sqlitepp::into(blob),
       sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
       sqlitepp::use(S);
    
    // get Chebyshev expansions
    rArray E_arr;
    std::vector<std::vector<std::tuple<int,int,int>>> Lll_arr;
    cArray cb;
    std::vector<std::vector<Chebyshev<double,Complex>>> cheb_arr;
    while (st.exec())
    {
        // on first iteration anf for every new energy
        if (E_arr.empty() or E != E_arr.back())
        {
            // save energy
            E_arr.push_back(E);
            
            // add new angular states for this energy
            Lll_arr.push_back(std::vector<std::tuple<int,int,int>>());
            
            // add new Chebyshev expansions for this energy
            cheb_arr.push_back(std::vector<Chebyshev<double,Complex>>());
        }
        
        // save angular states
        Lll_arr.back().push_back(std::make_tuple(L,l1,l2));
        
        // decode Chebyshev expansion from hexadecimal format
        cb.fromBlob(blob); 
        
        // save Chebyshev expansion
        cheb_arr.back().push_back
        (
            Chebyshev<double,Complex>
            (
                cb,                        // expansion coefficients
                0.,                        // lowest energy
                std::sqrt(E - 1./(ni*ni))  // highest energy
            )
        );
    }
    
    // terminate if no data
    if (E_arr.empty())
    {
        std::cout << "No data for this transition." << std::endl;
        return true;
    }
    
    // for all directions and energy shares evaluate the amplitude
    cArray ampls(dirs.size());
    for (std::size_t idir = 0; idir < dirs.size(); idir++)
    {
        // compute energy sharing factor
        double Eshare = 1. / (1 + dirs[idir].second.z / dirs[idir].first.z);
        
        // compute outgoing momenta so that
        //   1) (k₁)² + (k₂)² = Ei
        //   2) (k₂)² / (k₁)² = Eshare / (1 - Eshare)
        rArray k1 = sqrt((E_arr - 1./(ni*ni)) * Eshare);
        rArray k2 = sqrt((E_arr - 1./(ni*ni)) * (1 - Eshare));
        
        // evaluated amplitudes for all precomputed energies
        cArray ampls0(E_arr.size());
        
        // for all impact energies
        for (std::size_t ie = 0; ie < E_arr.size(); ie++)
        {
            // for all angular momenta
            for (std::size_t il = 0; il < Lll_arr[ie].size(); il++)
            {
                // get angular quantum numbers
                int  L = std::get<0>(Lll_arr[ie][il]);
                int l1 = std::get<1>(Lll_arr[ie][il]);
                int l2 = std::get<2>(Lll_arr[ie][il]);
                
                // evaluate the radial part for this angular & linear momenta
                Complex f = cheb_arr[ie][il].clenshaw(k1[ie],cheb_arr[ie][il].tail(1e-8)) / gsl_sf_pow_int(k1[ie] * k2[ie], 2);
                
                // evaluate bispherical function
                Complex YY = 0;
                for (int m = -l1; m <= l1; m++)
                {
                    YY += special::ClebschGordan(l1,m,l2,mi-m,L,mi)
                          * special::sphY(l1,m,dirs[idir].first.x*afactor,dirs[idir].first.y*afactor)
                          * special::sphY(l1,m,dirs[idir].second.x*afactor,dirs[idir].second.y*afactor);
                }
                
                // evaluate Coulomb phaseshifts
                Complex sig1 = special::coul_F_sigma(l1,k1[ie]);
                Complex sig2 = special::coul_F_sigma(l2,k2[ie]);
                
                // compute angular factors
                Complex angfact = std::pow(Complex(0.,1.),-l1-l2) * std::exp(Complex(0.,1.)*(sig1+sig2)) * YY;
                
                // sum the contribution
                ampls0[ie] += angfact * f;
            }
        }
        
        // interpolate
        if (not E_arr.empty())
            ampls[idir] = interpolate(E_arr, ampls0, { Ei })[0];
    }
    
    // write out
    std::cout << logo("#") <<
        "# Ionization amplitudes in " << unit_name(Lunits) << " for\n" <<
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
        "#     S = " << S << ", Ei = " << Ei << " in " << unit_name(Eunits) << "\n" <<
        "# ordered by direcion triplets (angles in " << unit_name(Aunits) << ")\n" <<
        "# \n";
    OutputTable table;
    table.setWidth(45, 45, 15, 15);
    table.setAlignment(OutputTable::left);
    table.write("# first electron directions", "second electron directions", "Re F     ", "Im F     ");
    table.write("# -------------------------", "--------------------------", "---------", "---------");
    
    for (std::size_t i = 0; i < dirs.size(); i++)
    {
        table.write
        (
            dirs[i].first,
            dirs[i].second,
            ampls[i].real()*lfactor,
            ampls[i].imag()*lfactor
        );
    }
    
    return true;
}
