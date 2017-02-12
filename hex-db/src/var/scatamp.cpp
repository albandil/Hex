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

#include "hex-interpolate.h"
#include "hex-chebyshev.h"
#include "hex-special.h"
#include "hex-version.h"

// --------------------------------------------------------------------------------- //

#include "../quantities.h"
#include "../utils.h"

// --------------------------------------------------------------------------------- //

createNewScatteringQuantity(ScatteringAmplitude, "scatamp")

// --------------------------------------------------------------------------------- //

extern void hex_tmat_pw_transform
(
    int ni, int li, int mi,
    int nf, int lf, int mf,
    int nEnergies,
    double * energies,
    bool extra,
    std::function
    <
        void
        (
            int ell,
            iArrays const & converged,
            iArrays const & complete,
            cArrays const & tmatrices
        )
    > f
);

// --------------------------------------------------------------------------------- //

std::string ScatteringAmplitude::description ()
{
    return "Scattering amplitude.";
}

std::vector<std::string> ScatteringAmplitude::dependencies ()
{
    return std::vector<std::string>
    {
        "tmat"
    };
}

std::vector<std::pair<std::string,std::string>> ScatteringAmplitude::params ()
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
        {"theta", "Scattering angles for which to compute the amplitude."}
    };
}

std::vector<std::string> ScatteringAmplitude::vparams ()
{
    return std::vector<std::string>
    {
        "theta"
    };
}

// --------------------------------------------------------------------------------- //

bool ScatteringAmplitude::initialize (sqlitepp::session & db)
{
    return ScatteringQuantity::initialize(db);
}

bool ScatteringAmplitude::createTable ()
{
    return ScatteringQuantity::createTable();
}

bool ScatteringAmplitude::updateTable ()
{
    return ScatteringQuantity::updateTable();
}

// --------------------------------------------------------------------------------- //

void hex_scattering_amplitude
(
    int ni, int li, int mi,
    int nf, int lf, int mf,
    int S,
    int nEnergies, double * energies,
    int nAngles, double * angles,
    double * result, double * extra
)
{
    hex_scattering_amplitude_
    (
        &ni, &li, &mi,
        &nf, &lf, &mf,
        &S,
        &nEnergies, energies,
        &nAngles, angles,
        result, extra
    );
}

void hex_scattering_amplitude_
(
    int * ni, int * li, int * mi,
    int * nf, int * lf, int * mf,
    int * S,
    int * nEnergies, double * energies,
    int * nAngles, double * angles,
    double * result, double * extra
)
{
    // get views of the output arrays
    cArrayView results, extras;
    results.reset((*nEnergies) * (*nAngles), reinterpret_cast<Complex*>(result));
    if (extra)
    {
        extras.reset((*nEnergies) * (*nAngles), reinterpret_cast<Complex*>(extra));
    }
    
    // calculate cosines of the scattering angles
    rArray cos_angles (*nAngles);
    for (int i = 0; i < (*nAngles); i++)
    {
        cos_angles[i] = std::cos(angles[i]);
    }
    
    // This function will be called by "hex_tmat_pw_transform" for every partial wave
    //    -  "ell" is the angular momentum of the partial wave
    //    -  "converged" indicates whether the partial wave expansion is converged (0 = not yet, 1 = yes, 2 = convergence failed)
    //    -  "complete" indicates whether all T-matrix components were found in the database (0 = no, 1 = yes)
    //    -  "tmatrices" are singlet and triplet T-matrix arrays at given energies
    auto fun = [&]
    (
        int ell,
        iArrays const & converged,
        iArrays const & complete,
        cArrays const & tmatrices
    )
    {
        // for all angles
        for (int i = 0; i < (*nAngles); i++)
        {
            // evaluate the spherical harmonics
            double Y = gsl_sf_legendre_sphPlm(ell, std::abs((*mi) - (*mf)), cos_angles[i]);
            
            // all energies
            for (int j = 0; j < (*nEnergies); j++) if (converged[*S][j] == 0)
            {
                Complex contrib = -tmatrices[*S][j] * Y / special::constant::two_pi;
                
                // update scattering amplitude
                if (complete[*S][j])
                    results[j * (*nAngles) + i] += contrib;
                
                // update the extrapolated scattering amplitude
                if (extra)
                    extras[j * (*nAngles) + i] += contrib;
                
            }
        }
    };
    
    // call the T-matrix retrieval driver
    hex_tmat_pw_transform
    (
        *ni, *li, *mi, *nf, *lf, *mf,
        *nEnergies, energies,
        extra != nullptr,
        fun
    );
}

// --------------------------------------------------------------------------------- //

bool ScatteringAmplitude::run (std::map<std::string,std::string> const & sdata)
{
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
    double lfactor = change_units(lUnit_au, Lunits);
    double afactor = change_units(Aunits, aUnit_rad);
    
    // scattering event parameters
    int ni = Conv<int>(sdata, "ni", name());
    int li = Conv<int>(sdata, "li", name());
    int mi = Conv<int>(sdata, "mi", name());
    int nf = Conv<int>(sdata, "nf", name());
    int lf = Conv<int>(sdata, "lf", name());
    int mf = Conv<int>(sdata, "mf", name());
    int  S = Conv<int>(sdata, "S", name());
    double E = Conv<double>(sdata, "Ei", name()) * efactor;
    
    // angles
    rArray angles;
    
    // get angle / angles
    try
    {
        // is there a single angle specified using command line ?
        angles.push_back(Conv<double>(sdata, "theta", name()));
    }
    catch (std::exception e)
    {
        // are there more angles specified using the STDIN ?
        angles = readStandardInput<double>();
    }
    
    // the scattering amplitudes
    rArray scaled_angles = angles * afactor;
    cArray amplitudes(angles.size());
    cArray extra(angles.size());
    hex_scattering_amplitude
    (
        ni, li, mi,
        nf, lf, mf,
        S,
        1, &E,
        angles.size(), scaled_angles.data(),
        reinterpret_cast<double*>(amplitudes.data()),
        reinterpret_cast<double*>(extra.data())
    );
    
    // write out
    std::cout << logo("#") <<
        "# Scattering amplitudes in " << unit_name(Lunits) << " for\n"
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n"
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n"
        "#     S = " << S << ", E = " << E/efactor << " " << unit_name(Eunits) << "\n"
        "# ordered by angle in " << unit_name(Aunits) << "\n"
        "# \n";
    OutputTable table;
    table.setWidth(15, 15);
    table.setAlignment(OutputTable::left);
    table.write("# angle    ", "Re f     ", "Im f     ", "Re f [ex]", "Im f [ex]");
    table.write("# ---------", "---------", "---------", "---------", "---------");
    
    for (std::size_t i = 0; i < angles.size(); i++)
    {
        table.write
        (
            angles[i],
            amplitudes[i].real() * lfactor,
            amplitudes[i].imag() * lfactor,
            extra[i].real() * lfactor,
            extra[i].imag() * lfactor
        );
    }
    
    return true;
}
