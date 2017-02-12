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
    double Ei, ReT, ImT;
    rArray avail_energies;
    cArrays tmatrices;
    
    ScatteringQuantity * TMat = get_quantity("tmat");
    
    // get view of the output result array
    cArrayView results;
    results.reset((*nEnergies) * (*nAngles), reinterpret_cast<Complex*>(result));
    results.fill(0);
    
    // get view of the (optional) output extrapolated array
    cArrayView extras;
    if (extra)
    {
        extras.reset((*nEnergies) * (*nAngles), reinterpret_cast<Complex*>(extra));
        extras.fill(0);
    }
    
    // transform angles to cosines
    rArray cos_angles = rArray(*nAngles, angles).transform([](double x){ return std::cos(x); });
    
    // indicator of partial waves convergence (0 = not yet converged, 1 = converged, 2 = convergence not possible)
    iArray status (*nEnergies);
    status.fill(0);
    
    // determine highest partial wave
    int max_ell;
    sqlitepp::statement st0 (TMat->session());
    st0 << "SELECT MAX(ell) FROM '" + TMat->name() + "' "
        "WHERE ni = :ni "
        "  AND li = :li "
        "  AND mi = :mi "
        "  AND nf = :nf "
        "  AND lf = :lf "
        "  AND mf = :mf "
        "  AND  S = :S  ",
        sqlitepp::into(max_ell),
        sqlitepp::use(*ni), sqlitepp::use(*li),  sqlitepp::use(*mi),
        sqlitepp::use(*nf), sqlitepp::use(*lf),  sqlitepp::use(*mf),
        sqlitepp::use(*S);
    st0.exec();
    
    // loop over all partial waves (or until convergence)
    for (int ell = std::abs((*mi) - (*mf)); ell <= max_ell or extra; ell++)
    {
        cArray tmatrices_ell (*nEnergies);
        iArray complete (*nEnergies, true);
        
        //
        // check if all energies converged
        //
        
            if (std::find(status.begin(), status.end(), 0) == status.end())
                break;
        
        //
        // read existing data from the database
        //
        
            // for all total angular momenta contributing to this partial wave
            for (int L = std::abs((*lf) - ell); L <= (*lf) + ell; L++)
            {
                // skip forbidden contributions
                if (special::ClebschGordan(ell, (*mi) - (*mf), (*lf), (*mf), L, (*mi)) == 0)
                    continue;
                
                rArray energies_db, reT_db, imT_db;
                
                // prepare selection statement
                sqlitepp::statement st1 (TMat->session());
                st1 << "SELECT Ei, Re_T_ell, Im_T_ell FROM '" + TMat->name() + "' "
                    "WHERE ni = :ni "
                    "  AND li = :li "
                    "  AND mi = :mi "
                    "  AND nf = :nf "
                    "  AND lf = :lf "
                    "  AND mf = :mf "
                    "  AND  L = :L  "
                    "  AND  S = :S  "
                    "  AND ell = :ell "
                    "ORDER BY Ei ASC ",
                    sqlitepp::into(Ei), sqlitepp::into(ReT), sqlitepp::into(ImT),
                    sqlitepp::use(*ni), sqlitepp::use(*li),  sqlitepp::use(*mi),
                    sqlitepp::use(*nf), sqlitepp::use(*lf),  sqlitepp::use(*mf),
                    sqlitepp::use( L),  sqlitepp::use(*S),   sqlitepp::use(ell);
                
                // store the data produced by the statement
                while (st1.exec())
                {
                    energies_db.push_back(Ei);
                    reT_db.push_back(ReT * std::sqrt(Ei));
                    imT_db.push_back(ImT * std::sqrt(Ei));
                }
                
                // default to Akima spline interpolation, but fall back to linear if not enough points
                gsl_interp_type const * interp = gsl_interp_akima;
                if (gsl_interp_type_min_size(interp) <= energies_db.size() or gsl_interp_type_min_size(interp = gsl_interp_linear) <= energies_db.size())
                {
                    // create interpolators
                    gsl_interp* spline_Re = gsl_interp_alloc(interp, energies_db.size());
                    gsl_interp* spline_Im = gsl_interp_alloc(interp, energies_db.size());
                    gsl_interp_init(spline_Re, energies_db.data(), reT_db.data(), energies_db.size());
                    gsl_interp_init(spline_Im, energies_db.data(), imT_db.data(), energies_db.size());
                    gsl_interp_accel* accel_Re = gsl_interp_accel_alloc();
                    gsl_interp_accel* accel_Im = gsl_interp_accel_alloc();
                    
                    // interpolate retrieved data at requested energies
                    for (int i = 0; i < (*nEnergies); i++) if (status[i] == 0)
                    {
                        if (energies_db.empty() or energies[i] < energies_db.front() or energies_db.back() < energies[i])
                        {
                            // not available -> incomplete data for this partial wave
                            complete[i] = false;
                        }
                        else
                        {
                            int err_re = gsl_interp_eval_e(spline_Re, energies_db.data(), reT_db.data(), energies[i], accel_Re, &ReT);
                            int err_im = gsl_interp_eval_e(spline_Im, energies_db.data(), imT_db.data(), energies[i], accel_Im, &ImT);
                            
                            if (err_re != GSL_SUCCESS)
                            {
                                std::cout << "Warning: Failed to interpolate real data for transition (" << ni << "," << li << "," << mi << ") -> (" << nf << "," << lf << "," << mf << "), "
                                    "L = " << L << ", ell = " << ell << ", Ei = " << energies[i] << " (" << gsl_strerror(err_re) << ")" << std::endl;
                            }
                            
                            if (err_im != GSL_SUCCESS)
                            {
                                std::cout << "Warning: Failed to interpolate imag data for transition (" << ni << "," << li << "," << mi << ") -> (" << nf << "," << lf << "," << mf << "), "
                                    "L = " << L << ", ell = " << ell << ", Ei = " << energies[i] << " (" << gsl_strerror(err_re) << ")" << std::endl;
                            }
                            
                            tmatrices_ell[i] += Complex(ReT, ImT) / std::sqrt(energies[i]);
                        }
                    }
                    
                    // delete interpolators
                    gsl_interp_accel_free(accel_Re);
                    gsl_interp_accel_free(accel_Im);
                    gsl_interp_free(spline_Re);
                    gsl_interp_free(spline_Im);
                }
                else
                {
                    // not enough data for interpolation -> all T-matrices deemed incomplete
                    for (int i = 0; i < (*nEnergies); i++) if (status[i] == 0)
                        complete[i] = false;
                    break;
                }
            }
        
        //
        // extrapolate incomplete data
        //
        
            for (int i = 0; i < (*nEnergies); i++) if (not complete[i] and status[i] == 0)
            {
                // only extrapolate if we have two or more samples, that are decreasing at the end
                if
                (
                    tmatrices.size() >= 2 and
                    std::abs(tmatrices.back(0)[i].real()) < std::abs(tmatrices.back(1)[i].real()) and
                    std::abs(tmatrices.back(0)[i].imag()) < std::abs(tmatrices.back(1)[i].imag())
                )
                {
                    // degenerate transitions are notorious for slow convergence -> use finite range formula
                    if ((*ni) == (*nf))
                    {
                        // extrapolate T_ell for elastic scattering [T ~ L^(-2.5)]
                        tmatrices_ell[i] = Complex
                        (
                            tmatrices.back()[i].real() * std::pow((ell - 1.0) / ell, 2.5),
                            tmatrices.back()[i].imag() * std::pow((ell - 1.0) / ell, 1.5)
                        );
                    }
                    
                    // dipole transitions converge fast, no need to extrapolate
                    else if (std::abs((*li) - (*lf)) == 1)
                    {
                        // nothing
                    }
                    
                    // other transitions need simple geometric extrapolation
                    else
                    {
                        // extrapolate T_ell for inelastic scattering [use geometric series]
                        double extra_re = tmatrices.back(0)[i].real() * tmatrices.back(0)[i].real() / tmatrices.back(1)[i].real();
                        double extra_im = tmatrices.back(0)[i].imag() * tmatrices.back(0)[i].imag() / tmatrices.back(1)[i].imag();
                        tmatrices_ell[i] = Complex(extra_re,extra_im);
                    }
                    
                    // check convergence
                    Complex tmat_sum = 0;
                    for (unsigned n = 0; n < tmatrices.size(); n++)
                        tmat_sum += tmatrices[n][i];
                    if (std::abs(tmatrices_ell[i]) < 1e-5 * std::abs(tmat_sum))
                        status[i] = 1;
                }
                else
                {
                    status[i] = 2; // convergence not possible
                }
            }
            
            tmatrices.push_back(tmatrices_ell);
        
        //
        // update scattering amplitudes
        //
        
            // all angles
            for (int i = 0; i < (*nAngles); i++)
            {
                double Y = gsl_sf_legendre_sphPlm(ell, std::abs((*mi) - (*mf)), cos_angles[i]);
                
                // all energies
                for (int j = 0; j < (*nEnergies); j++) if (status[j] == 0)
                {
                    Complex contrib = -tmatrices_ell[j] * Y / special::constant::two_pi;
                    
                    // update scattering amplitude
                    if (complete[j])
                        results[j * (*nAngles) + i] += contrib;
                    
                    // update the extrapolated scattering amplitude
                    if (extra)
                        extras[j * (*nAngles) + i] += contrib;
                    
                }
            }
    }
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
