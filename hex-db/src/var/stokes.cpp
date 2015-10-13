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
#include "hex-special.h"
#include "hex-version.h"

#include "interfaces.h"
#include "variables.h"

//
// StokesParameters members
//

const std::string StokesParameters::Id = "stokes";
const std::string StokesParameters::Description = "Reduced Stokes parameters for the ns->n'p transition.";
const std::vector<std::pair<std::string,std::string>> StokesParameters::Dependencies = {
    {"ni", "Initial atomic principal quantum number; the state (ni,0,0) is used."},
    {"nf", "Final atomic principal quantum number; the states (nf,1,*) are used."},
    {"Ei", "Projectile impact energy (Rydberg)."},
    {"beta", "Angle between impact direction and the quantization axis."},
    {"theta", "Scattering angles for which to compute the Stokes parameters."}
};
const std::vector<std::string> StokesParameters::VecDependencies = { "theta" };

bool StokesParameters::initialize (sqlitepp::session & db) const
{
    return true;
}

std::vector<std::string> const & StokesParameters::SQL_CreateTable () const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

std::vector<std::string> const & StokesParameters::SQL_Update () const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

bool StokesParameters::run (std::map<std::string,std::string> const & sdata) const
{
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
    double afactor = change_units(Aunits, aUnit_rad);
    
    // atomic and projectile data
    int ni = Conv<int>(sdata, "ni", Id);
    int nf = Conv<int>(sdata, "nf", Id);
    int Ei = Conv<double>(sdata, "Ei", Id) * efactor;
    double ki = sqrt(Ei);
    double kf = sqrt(Ei - 1./(ni*ni) + 1./(nf*nf));
    
    // angles
    rArray angles;
    
    // get angle / angles
    try {
        
        // is there a single angle specified using command line ?
        angles.push_back(Conv<double>(sdata, "theta", Id));
        
    } catch (std::exception e) {
        
        // are there more angles specified using the STDIN ?
        angles = readStandardInput<double>();
    }
    
    // Should return the following:
    //    P₁ = 2λ - 1
    //    P₂ = -2√2 R
    //    P₃ = 2√2 I
    // and
    //    Pl = √(P₁²+P₂²)
    //    γ  = arg (P₁ + iP₂) / 2
    //    P⁺ = √(P₁²+P₂² + P₃²)
    // where
    //    λ = <|f₀²|> / [ dσ/dΩ ]
    //    R = Re { <f₀* f₁> } / [ dσ/dΩ ]
    //    I = Im { <f₀* f₁> } / [ dσ/dΩ ]
    // where
    //    f₀ ... is the scattering amplitude to mf = 0
    //    f₁ ... is the scattering amplitude to mf = 1
    //    dσ/dΩ ... is the DCS summed over mf
    //    <.> ... stands for averaging over spin states, i.e.
    //        <a> = [ a(S=0) + 3a(S=1) ] / 4
    
    rArray scaled_angles = angles * afactor;
    
    // compute scattering amplitudes
    cArray f0_singlet(angles.size()), f0_triplet(angles.size());
    cArray f1_singlet(angles.size()), f1_triplet(angles.size());
    hex_scattering_amplitude(ni,0,0, nf,1,0, 0, Ei, angles.size(), scaled_angles.data(), reinterpret_cast<double*>(f0_singlet.data()));
    hex_scattering_amplitude(ni,0,0, nf,1,0, 1, Ei, angles.size(), scaled_angles.data(), reinterpret_cast<double*>(f0_triplet.data()));
    hex_scattering_amplitude(ni,0,0, nf,1,1, 0, Ei, angles.size(), scaled_angles.data(), reinterpret_cast<double*>(f1_singlet.data()));
    hex_scattering_amplitude(ni,0,0, nf,1,1, 1, Ei, angles.size(), scaled_angles.data(), reinterpret_cast<double*>(f1_triplet.data()));
    
    // compute differential cross sections
    rArray dcs_p0s(angles.size()), dcs_pps(angles.size()), dcs_pms(angles.size());
    rArray dcs_p0t(angles.size()), dcs_ppt(angles.size()), dcs_pmt(angles.size());
    hex_differential_cross_section(ni,0,0, nf,1, 0, 0, Ei, angles.size(), scaled_angles.data(), dcs_p0s.data());
    hex_differential_cross_section(ni,0,0, nf,1, 1, 0, Ei, angles.size(), scaled_angles.data(), dcs_pps.data());
    hex_differential_cross_section(ni,0,0, nf,1,-1, 0, Ei, angles.size(), scaled_angles.data(), dcs_pms.data());
    hex_differential_cross_section(ni,0,0, nf,1, 0, 1, Ei, angles.size(), scaled_angles.data(), dcs_p0t.data());
    hex_differential_cross_section(ni,0,0, nf,1, 1, 1, Ei, angles.size(), scaled_angles.data(), dcs_ppt.data());
    hex_differential_cross_section(ni,0,0, nf,1,-1, 1, Ei, angles.size(), scaled_angles.data(), dcs_pmt.data());
    rArray dcs = dcs_p0s + dcs_pps + dcs_pms + dcs_p0t + dcs_ppt + dcs_pmt;
    
    // compute basic parameters
    rArray lambda = (0.25 * abs(f0_singlet*f0_singlet) + 0.75 * abs(f0_triplet*f0_triplet)) / dcs * (kf/ki);
    rArray R = realpart(0.25 * f1_singlet*f0_singlet.conj() + 0.75 * f1_triplet*f0_triplet.conj()) / dcs * (kf/ki);
    rArray I = imagpart(0.25 * f1_singlet*f0_singlet.conj() + 0.75 * f1_triplet*f0_triplet.conj()) / dcs * (kf/ki);
    
    // compute derived parameters
    rArray P1 = 2. * lambda - 1.;
    rArray P2 = 2 * special::constant::sqrt_two * R;
    rArray P3 = -2 * special::constant::sqrt_two * I;
    rArray Pl = hypot(P1, P2);
    rArray gamma = 0.5 * atan2(P2, P1);
    rArray Pplus = hypot(Pl, P3);
    
    // write out
    std::cout << logo("#") <<
        "# Stokes parameters for\n"
        "#     ni = " << ni << ", li = " << 0 << ",\n"
        "#     nf = " << nf << ", lf = " << 1 << ",\n"
        "#     E = " << Ei/efactor << " " << unit_name(Eunits) << "\n"
        "# ordered by angle in " << unit_name(Aunits) << "\n"
        "# \n";
    OutputTable table;
    table.setWidth(15);
    table.setAlignment(OutputTable::left);
    table.write("# angle    ", "lambda   ", "R        ", "I        ", "P1       ", "P2       ", "P3       ", "Pl       ", "gamma    ", "P+       ");
    table.write("# ---------", "---------", "---------", "---------", "---------", "---------", "---------", "---------", "---------", "---------");
    
    for (std::size_t i = 0; i < angles.size(); i++)
    {
        table.write
        (
            angles[i],
            lambda[i], R[i], I[i],
            P1[i], P2[i], P3[i],
            Pl[i], gamma[i] / afactor, Pplus[i]
        );
    }
    
    return true;
}
