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

#include "../interpolate.h"
#include "../chebyshev.h"
#include "../special.h"
#include "../variables.h"
#include "../version.h"

const std::string ScatteringAmplitudeDir::Id = "scatamp-dir";
const std::string ScatteringAmplitudeDir::Description = "Scattering amplitude (arbitrary impact direction).";
const std::vector<std::pair<std::string,std::string>> ScatteringAmplitudeDir::Dependencies = {
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
const std::vector<std::string> ScatteringAmplitudeDir::VecDependencies = { "theta" };

bool ScatteringAmplitudeDir::initialize(sqlitepp::session & db) const
{
    return true;
}

std::vector<std::string> const & ScatteringAmplitudeDir::SQL_Update () const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

std::vector<std::string> const & ScatteringAmplitudeDir::SQL_CreateTable () const
{
    static const std::vector<std::string> cmd;
    return cmd;
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
        // compute the scattering amplitude (ni,li,mi_)->(nf,lf,mf_)
        hex_scattering_amplitude_
        (
            ni,li,&mi_,
            nf,lf,&mf_,
            S, E,
            N, angles,
            reinterpret_cast<double*>(amplitudes.data())
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

bool ScatteringAmplitudeDir::run (std::map<std::string,std::string> const & sdata) const
{
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
    double lfactor = change_units(lUnit_au, Lunits);
    double afactor = change_units(Aunits, aUnit_rad);
    
    // scattering event parameters
    int ni = As<int>(sdata, "ni", Id);
    int li = As<int>(sdata, "li", Id);
    int mi = As<int>(sdata, "mi", Id);
    int nf = As<int>(sdata, "nf", Id);
    int lf = As<int>(sdata, "lf", Id);
    int mf = As<int>(sdata, "mf", Id);
    int  S = As<int>(sdata, "S",  Id);
    double E     = As<double>(sdata, "Ei",    Id) * efactor;
    double alpha = As<double>(sdata, "alpha", Id) * afactor;
    double beta  = As<double>(sdata, "beta",  Id) * afactor;
    double gamma = As<double>(sdata, "gamma", Id) * afactor;
    
    // angles
    rArray angles;
    
    // get angle / angles
    try {
        
        // is there a single angle specified using command line ?
        angles.push_back(As<double>(sdata, "theta", Id));
        
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
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n"
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n"
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
