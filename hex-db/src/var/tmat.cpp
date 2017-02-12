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
#include "hex-special.h"
#include "hex-version.h"

// --------------------------------------------------------------------------------- //

#include "../quantities.h"
#include "../utils.h"

// --------------------------------------------------------------------------------- //

createNewScatteringQuantity(TMatrix, "tmat")

// --------------------------------------------------------------------------------- //

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

void hex_tmat_pw_transform
(
    int ni, int li, int mi,
    int nf, int lf, int mf,
    int nEnergies, double * energies,
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
    > callback
)
{
    double Ei, ReT, ImT;
    cArrays tmatrices_prev(2), tmatrices_prev_prev(2);
    rArrays tmatrices_sum = { rArray(nEnergies, 0.0), rArray(nEnergies, 0.0) };
    
    ScatteringQuantity * TMat = get_quantity("tmat");
    
    // indicator of partial waves convergence (0 = not yet converged, 1 = converged, 2 = convergence not possible)
    iArrays status;
    status.resize(2);
    status[0].resize(nEnergies);
    status[1].resize(nEnergies);
    
    // determine highest partial wave
    int max_ell;
    sqlitepp::statement st0 (TMat->session());
    st0 << "SELECT MAX(ell) FROM '" + TMat->name() + "' "
        "WHERE ni = :ni "
        "  AND li = :li "
        "  AND mi = :mi "
        "  AND nf = :nf "
        "  AND lf = :lf "
        "  AND mf = :mf ",
        sqlitepp::into(max_ell),
        sqlitepp::use(ni), sqlitepp::use(li),  sqlitepp::use(mi),
        sqlitepp::use(nf), sqlitepp::use(lf),  sqlitepp::use(mf);
    st0.exec();
    
    // loop over all partial waves (or until convergence)
    for (int ell = std::abs(mi - mf); ell <= max_ell or extra; ell++)
    {
        cArrays tmatrices = { cArray(nEnergies, 0.0_z), cArray(nEnergies, 0.0_z) };
        iArrays complete  = { iArray(nEnergies, true),  iArray(nEnergies, true)  };
        
        // check if all energies converged
        if
        (
            (
                std::find(status[0].begin(), status[0].end(), 0) == status[0].end() and
                std::find(status[1].begin(), status[1].end(), 0) == status[1].end()
            )
            or ell > 1000 // hard limit
        )
        {
            break;
        }
        
        // for both spins
        for (int S = 0; S < 2; S++)
        {
            //
            // read existing data from the database
            //
            
                // for all total angular momenta contributing to this partial wave
                for (int L = std::abs(lf - ell); L <= lf + ell; L++)
                {
                    // skip forbidden contributions
                    if (special::ClebschGordan(ell, mi - mf, lf, mf, L, mi) == 0)
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
                        sqlitepp::use(ni),  sqlitepp::use(li),   sqlitepp::use(mi),
                        sqlitepp::use(nf),  sqlitepp::use(lf),   sqlitepp::use(mf),
                        sqlitepp::use( L),  sqlitepp::use(S),    sqlitepp::use(ell);
                    
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
                        for (int i = 0; i < nEnergies; i++) if (status[S][i] == 0)
                        {
                            if (energies_db.empty() or energies[i] < energies_db.front() or energies_db.back() < energies[i])
                            {
                                // not available -> incomplete data for this partial wave
                                complete[S][i] = false;
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
                                
                                tmatrices[S][i] += Complex(ReT, ImT) / std::sqrt(energies[i]);
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
                        for (int i = 0; i < nEnergies; i++) if (status[S][i] == 0)
                            complete[S][i] = false;
                        break;
                    }
                }
            
            //
            // extrapolate incomplete data
            //
            
                for (int i = 0; i < nEnergies; i++) if (not complete[S][i] and status[S][i] == 0)
                {
                    // only extrapolate if this is at least the third contributing partial wave and the last three are decreasing
                    if
                    (
                        not tmatrices_prev[S].empty() and
                        std::abs(tmatrices[S][i].real()) < std::abs(tmatrices_prev[S][i].real()) and
                        std::abs(tmatrices[S][i].imag()) < std::abs(tmatrices_prev[S][i].imag()) and
                        not tmatrices_prev_prev[S].empty() and
                        std::abs(tmatrices_prev[S][i].real()) < std::abs(tmatrices_prev_prev[S][i].real()) and
                        std::abs(tmatrices_prev[S][i].imag()) < std::abs(tmatrices_prev_prev[S][i].imag())
                    )
                    {
                        // degenerate transitions are notorious for slow convergence -> use finite range formula
                        if (ni == nf)
                        {
                            // extrapolate T_ell for elastic scattering [T ~ L^(-2.5)]
                            tmatrices[S][i] = Complex
                            (
                                tmatrices_prev[S][i].real() * std::pow((ell - 1.0) / ell, 2.5),
                                tmatrices_prev[S][i].imag() * std::pow((ell - 1.0) / ell, 1.5)
                            );
                        }
                        
                        // dipole transitions converge fast, no need to extrapolate
                        else if (std::abs(li - lf) == 1)
                        {
                            // nothing
                        }
                        
                        // other transitions need simple geometric extrapolation
                        else
                        {
                            // extrapolate T_ell for inelastic scattering [use geometric series]
                            double extra_re = tmatrices_prev[S][i].real() * tmatrices_prev[S][i].real() / tmatrices_prev_prev[S][i].real();
                            double extra_im = tmatrices_prev[S][i].imag() * tmatrices_prev[S][i].imag() / tmatrices_prev_prev[S][i].imag();
                            tmatrices[S][i] = Complex(extra_re, extra_im);
                        }
                        
                        // update sum of (magnitudes of) the partial T-matrices
                        tmatrices_sum[S][i] += std::abs(tmatrices[S][i]);
                        
                        // check convergence
                        if (std::abs(tmatrices[S][i]) < 1e-5 * std::abs(tmatrices_sum[S][i]))
                            status[S][i] = 1;
                    }
                    else
                    {
                        status[S][i] = 2; // extrapolation not possible - too few data
                    }
                }
                
                
        }
        
        // update previous-partial-wave T-matrices
        tmatrices_prev_prev = tmatrices_prev;
        tmatrices_prev = tmatrices;
        
        // pass the partial T-matrices to the caller
        callback(ell, status, complete, tmatrices);
    }
    
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
