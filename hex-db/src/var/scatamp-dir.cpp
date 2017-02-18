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

createNewScatteringQuantity(ScatteringAmplitudeDir, "scatamp-dir")

// --------------------------------------------------------------------------------- //

std::string ScatteringAmplitudeDir::description ()
{
    return "Scattering amplitude (arbitrary impact direction).";
}

std::vector<std::string> ScatteringAmplitudeDir::dependencies ()
{
    return std::vector<std::string>
    {
        "scatamp"
    };
}

std::vector<std::pair<std::string,std::string>> ScatteringAmplitudeDir::params ()
{
    return std::vector<std::pair<std::string,std::string>>
    {
        {"ni",    "Initial atomic principal quantum number."                  },
        {"li",    "Initial atomic orbital quantum number."                    },
        {"mi",    "Initial atomic magnetic quantum number."                   },
        {"nf",    "Final atomic principal quantum number."                    },
        {"lf",    "Final atomic orbital quantum number."                      },
        {"mf",    "Final atomic magnetic quantum number."                     },
        {"S",     "Total spin of atomic + projectile electron."               },
        {"Ei",    "Projectile impact energy (Rydberg)."                       },
        {"alpha", "Quantization axis orientation (1st Euler angle)."          },
        {"beta",  "Quantization axis orientation (2nd Euler angle)."          },
        {"gamma", "Quantization axis orientation (3rd Euler angle)."          },
        {"theta", "Scattering angles for which to compute the amplitude."     }
    };
}

std::vector<std::string> ScatteringAmplitudeDir::vparams ()
{
    return std::vector<std::string>
    {
        "theta"
    };
}

// --------------------------------------------------------------------------------- //

bool ScatteringAmplitudeDir::initialize (sqlitepp::session & db)
{
    return ScatteringQuantity::initialize(db);
}

bool ScatteringAmplitudeDir::createTable ()
{
    return ScatteringQuantity::createTable();
}

bool ScatteringAmplitudeDir::updateTable ()
{
    return ScatteringQuantity::updateTable();
}

// --------------------------------------------------------------------------------- //

void hex_scattering_amplitude_dir
(
    int ni, int li, int mi,
    int nf, int lf, int mf,
    int S, double E, int N,
    double alpha, double beta, double gamma,
    double * angles, double * result
)
{
    hex_scattering_amplitude_dir_
    (
        &ni, &li, &mi,
        &nf, &lf, &mf,
        &S, &E, &N,
        &alpha, &beta, &gamma,
        angles, result
    );
}

void hex_scattering_amplitude_dir_
(
    int * ni, int * li, int * mi,
    int * nf, int * lf, int * mf,
    int * S, double * E, int * N,
    double * alpha, double * beta, double * gamma,
    double * angles, double * result
)
{
    // set up and clean the output storage
    cArrayView resarr(*N,reinterpret_cast<Complex*>(result));
    resarr.fill(0.);
    
    // create the scattering amplitudes (for one term of the sum)
    cArray amplitudes(*N, 0.);
    
    // sum over all initial and final magnetic quantum numbers
    for (int mi_ = -(*li); mi_ <= (*li); mi_++)
    for (int mf_ = -(*lf); mf_ <= (*lf); mf_++)
    {
        int one = 1;
        
        // compute the scattering amplitude (ni,li,mi_)->(nf,lf,mf_)
        hex_scattering_amplitude_
        (
            ni,li,&mi_,
            nf,lf,&mf_,
            S,
            &one, E,
            N, angles,
            reinterpret_cast<double*>(amplitudes.data()),
            nullptr
        );
        
        // calculate the Wigner d-matrix factors
        double di = special::Wigner_d(2*(*li),2*(*mi),2*mi_,*beta);
        double df = special::Wigner_d(2*(*lf),2*(*mf),2*mf_,*beta);
        
        // calculate the whole transformation factor (D-matrices)
        Complex D = std::exp(Complex(0.,(*alpha)*((*mf) - (*mi)))) * di * df *
                    std::exp(Complex(0.,(*gamma)*(  mf_ -   mi_)));
        
        // contribute to the result
        resarr += D * amplitudes;
    }
}

// --------------------------------------------------------------------------------- //

bool ScatteringAmplitudeDir::run (std::map<std::string,std::string> const & sdata)
{
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
    double lfactor = change_units(lUnit_au, Lunits);
    double afactor = change_units(Aunits, aUnit_rad);
    
    // scattering event parameters
    int ni = Conv<int>(sdata, "ni", name());
    int li = Conv<int>(sdata, "li", name());
    int mi0= Conv<int>(sdata, "mi", name());
    int nf = Conv<int>(sdata, "nf", name());
    int lf = Conv<int>(sdata, "lf", name());
    int mf0= Conv<int>(sdata, "mf", name());
    int  S = Conv<int>(sdata, "S",  name());
    double E     = Conv<double>(sdata, "Ei",    name()) * efactor;
    double alpha = Conv<double>(sdata, "alpha", name()) * afactor;
    double beta  = Conv<double>(sdata, "beta",  name()) * afactor;
    double gamma = Conv<double>(sdata, "gamma", name()) * afactor;
    
    // use mi >= 0; if mi < 0, flip both signs
    int mi = (mi0 < 0 ? -mi0 : mi0);
    int mf = (mi0 < 0 ? -mf0 : mf0);
    
    // angles
    rArray angles;
    
    // get angle / angles
    try {
        
        // is there a single angle specified using command line ?
        angles.push_back(Conv<double>(sdata, "theta", name()));
        
    } catch (std::exception e) {
        
        // are there more angles specified using the STDIN ?
        angles = readStandardInput<double>();
    }
    
    // the scattering angles in radians
    rArray scaled_angles = angles * afactor;
    
    // the scattering amplitudes
    cArray amplitudes(angles.size(),0.0);
    hex_scattering_amplitude_dir
    (
        ni,li,mi,
        nf,lf,mf,
        S, E,
        angles.size(),
        alpha, beta, gamma,
        scaled_angles.data(),
        reinterpret_cast<double*>(amplitudes.data())
    );
    
    // write out
    std::cout << logo("#") <<
        "# Scattering amplitudes in " << unit_name(Lunits) << " for\n"
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi0 << ",\n"
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf0 << ",\n"
        "#     S = " << S << ", E = " << E/efactor << " " << unit_name(Eunits) << "\n"
        "#     impact direction = (" << beta << "," << gamma << ") " << unit_name(Aunits) << "\n"
        "# ordered by angle in " << unit_name(Aunits) << "\n"
        "# \n"
        "# θ\t Re f\t Im f\n";
    for (std::size_t i = 0; i < angles.size(); i++)
    {
        std::cout << 
            angles[i] << "\t" << 
            amplitudes[i].real() * lfactor << "\t" <<
            amplitudes[i].imag() * lfactor << "\n";
    }
    
    return true;
}
