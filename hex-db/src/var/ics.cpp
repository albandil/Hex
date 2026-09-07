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

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

// --------------------------------------------------------------------------------- //

#include <sqlite3.h>

// --------------------------------------------------------------------------------- //

#include <gsl/gsl_interp.h>

// --------------------------------------------------------------------------------- //

#include "hex-chebyshev.h"
#include "hex-clenshawcurtis.h"
#include "hex-interpolate.h"
#include "hex-version.h"

// --------------------------------------------------------------------------------- //

#include "../quantities.h"
#include "../utils.h"

// --------------------------------------------------------------------------------- //

createNewScatteringQuantity(IntegralCrossSection, "ics")

// --------------------------------------------------------------------------------- //

//
// custom function for evaluation of square root within a SQL statement
//

void db_sqrt (sqlite3_context* pdb, int n, sqlite3_value** val)
{
    assert(n == 1);
    sqlite3_result_double(pdb, std::sqrt(sqlite3_value_double(*val)));
}

//
// custom function for integration of BLOB-represented Chebyshev
// expansion of the ionization amplitude
//

void db_ioncs (sqlite3_context* pdb, int n, sqlite3_value** val)
{
    assert(n == 1);

    // copy the coefficients out of the BLOB
    // NOTE: sqlite3_value_bytes has to be called after sqlite3_value_blob, never
    //       before it, and both return zero/null for a NULL column.
    cArray coeffs;
    void const * bytes = sqlite3_value_blob(*val);
    coeffs.fromBytes(bytes, sqlite3_value_bytes(*val));

    // some blobs can be empty
    if (coeffs.empty())
    {
        sqlite3_result_double(pdb, 0.0);
        return;
    }

    // reconstruct Chebyshev approximation object from the stored data
    Chebyshev<double,Complex> CB(coeffs, 0, 1);

    // integrate
    //
    // 1/√2                π/4
    //  ⌠                   ⌠
    //  ⎮            dκ     ⎮
    //  ⎮ |f(κ)|² ------- = ⎮ |f(sin β)|² dβ
    //  ⎮         √(1-κ²)   ⎮
    //  ⌡                   ⌡
    //  0                   0
    //
    // NOTE: The stored expansion holds Ξ / sqrt(k₁ k₂), so |Ξ|² is the square of the
    //       expansion times k₁ k₂ = E_tot sin β cos β. The energy is not available in a
    //       scalar SQL function, so only sin β cos β is taken here and the caller has to
    //       multiply the result by E_tot; see @ref IntegralCrossSection::updateTable.
    int tail = CB.tail(1e-10);
    auto fsqr = [&](double beta) -> double
    {
        return std::sin(beta) * std::cos(beta) * sqrabs(CB.clenshaw(std::sin(beta), tail));
    };
    ClenshawCurtis<decltype(fsqr),double> integrator(fsqr);
    double result = integrator.integrate(0, special::constant::pi_quart);

    // use result of the integration
    sqlite3_result_double(pdb, result);
}

//
// custom function for linear interpolation of two points
//

void db_interpolate (sqlite3_context* pdb, int n, sqlite3_value** val)
{
    assert(n == 5);

    // extract parameters
    double x1 = sqlite3_value_double(*(val + 0));
    double y1 = sqlite3_value_double(*(val + 1));
    double x2 = sqlite3_value_double(*(val + 2));
    double y2 = sqlite3_value_double(*(val + 3));
    double x  = sqlite3_value_double(*(val + 4));

    // handle degenerate case
    if (x1 == x2)
    {
        if (x == x1)
            sqlite3_result_double(pdb, y1);
        else
            sqlite3_result_double(pdb, special::constant::Nan);
    }

    // interpolate
    else
    {
        sqlite3_result_double(pdb, ((x - x1) * y2 + (x2 - x) * y1) / (x2 - x1));
    }
}

// --------------------------------------------------------------------------------- //

/**
 * @brief Sum the partial ionization cross sections over the total angular momentum.
 *
 * The partial cross sections decay geometrically in the total angular momentum. Behind
 * that is the adiabatic (Massey) cutoff of the energy transfer: a collision at the
 * impact parameter @f$ b @f$ can hand over the energy @f$ \Delta E @f$ only as long as
 * @f$ b @f$ stays below @f$ v/\Delta E @f$, which in terms of the partial wave reads
 * @f[
 *     \ell \lesssim k v / \Delta E = E_i / \Delta E \ \mathrm{(rydbergs)} ,
 * @f]
 * and past that limit the contributions fall off exponentially, the ratio of two
 * successive ones approaching @f$ q = \exp(-\Delta E / E_i) @f$. What is left of the
 * series behind the last computed partial wave is then a geometric progression,
 * @f[
 *     \sigma_{\mathrm{tail}} = \sigma_L \frac{q}{1-q} , \qquad
 *     q = \sigma_L / \sigma_{L-1} ,
 * @f]
 * whose addition is the same operation as Aitken's delta-squared process applied to the
 * partial sums. It saves roughly a third of the partial waves needed for a given
 * accuracy; measured on the runs of e-H ionization at 17.6, 30 and 54.4 eV, a part in
 * ten thousand of the total was reached at L = 7, 9 and 20 instead of 8, 16 and 29.
 * Being the estimate of what is left out, the remainder is returned as well.
 *
 * Two things spoil the ratio and both have to be kept out. Below the maximum of
 * @f$ \sigma_L @f$ the sequence is still growing, the ratio exceeds unity and the
 * formula would return a remainder of the wrong sign -- at 54.4 eV an extrapolation
 * from L = 4 is worse than the bare sum by a factor of four. Far above the maximum the
 * partial cross sections eventually reach the noise floor of the amplitudes, stop
 * decaying altogether and drive the ratio to one, where the factor @f$ q/(1-q) @f$ is
 * unbounded; in the 54.4 eV run this happens above L = 40 and the last available ratio
 * is 0.9998, which would multiply the noise by five thousand. The energy transfer of
 * ionization is at least the ionization threshold, so
 * @f[
 *     q \le \exp(-E_{\mathrm{ion}} / E_i)
 * @f]
 * is a bound that both regimes violate. Demanding it of two successive ratios singles
 * out the largest partial wave that is still within the geometric regime. Nothing more
 * elaborate is needed: requiring the two ratios to agree with each other as well was
 * tried and changed none of the three reference runs by more than one partial wave.
 *
 * @param sigma Partial cross sections, ordered by the partial wave.
 * @param Ei Impact energy (Rydberg).
 * @param Eion Ionization threshold (Rydberg).
 * @param tail Output estimate of the neglected remainder, NaN when there is no
 *             partial wave to extrapolate from.
 * @return The summed cross section, the remainder included.
 */
static double sum_partial_waves (rArray const & sigma, double Ei, double Eion, double & tail)
{
    // bound on the ratio implied by the adiabatic cutoff
    double qmax = (Ei > 0. ? std::exp(-Eion / Ei) : 0.);

    // look for the last partial wave that is still within the geometric regime
    std::size_t last = 0;
    double q = 0.;

    for (std::size_t n = 2; n < sigma.size(); n++)
    {
        if (sigma[n-2] <= 0. or sigma[n-1] <= 0. or sigma[n] <= 0.)
            continue;

        double qn = sigma[n] / sigma[n-1];
        double qp = sigma[n-1] / sigma[n-2];

        if (qn <= qmax and qp <= qmax)
        {
            last = n;
            q = qn;
        }
    }

    // there is none: hand back the bare sum and admit that nothing is known about the rest
    if (last == 0)
    {
        tail = special::constant::Nan;
        return sum(sigma);
    }

    // truncate there and replace everything beyond by the remainder of the progression
    // NOTE: The partial waves past 'last' are dropped on purpose. They are either the
    //       noise floor of the amplitudes or too few of them to judge, and the
    //       remainder that stands in for them is the better estimate of the two.
    double total = 0.;
    for (std::size_t n = 0; n <= last; n++)
        total += sigma[n];

    tail = sigma[last] * q / (1. - q);

    return total + tail;
}

// --------------------------------------------------------------------------------- //

std::string IntegralCrossSection::description ()
{
    return "Integral cross section.";
}

std::vector<std::string> IntegralCrossSection::dependencies ()
{
    return std::vector<std::string>
    {
        "tmat",
        "ionf"
    };
}

std::vector<std::pair<std::string,std::string>> IntegralCrossSection::params ()
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
        {"ell", "Partial wave. A negative value asks for the sum over all partial waves, "
                "which is available for ionization (nf = lf = mf = 0)."}
    };
}

std::vector<std::string> IntegralCrossSection::vparams ()
{
    return std::vector<std::string>
    {
        "Ei"
    };
}

// --------------------------------------------------------------------------------- //

bool IntegralCrossSection::initialize (sqlitepp::session & db)
{
    // define square root function
    sqlite3_create_function(db.impl(), "SQRT", 1, SQLITE_UTF8, nullptr, &db_sqrt, nullptr, nullptr);

    // define Gauss-Chebyshev integration of squared Chebyshev expansion
    sqlite3_create_function(db.impl(), "IONCS", 1, SQLITE_UTF8, nullptr, &db_ioncs, nullptr, nullptr);

    // define linear interpolation of two data points (5 arguments: x1, y1, x2, y2, x)
    sqlite3_create_function(db.impl(), "INTERPOLATE", 5, SQLITE_UTF8, nullptr, &db_interpolate, nullptr, nullptr);

    return ScatteringQuantity::initialize(db);
}

bool IntegralCrossSection::createTable ()
{
    sqlitepp::statement st (session());
    st <<
        "CREATE TABLE IF NOT EXISTS 'ics' "
        "("
            "ni  INTEGER, "
            "li  INTEGER, "
            "mi  INTEGER, "
            "nf  INTEGER, "
            "lf  INTEGER, "
            "mf  INTEGER, "
            "S   INTEGER, "
            "Ei  DOUBLE PRECISION, "
            "ell INTEGER, "
            "sigma DOUBLE PRECISION, "
            "PRIMARY KEY (ni,li,mi,nf,lf,mf,S,Ei,ell)"
        ")";

    try
    {
        st.exec();
    }
    catch (sqlitepp::exception & e)
    {
        std::cerr << "ERROR: Creation of table 'ics' failed!" << std::endl;
        std::cerr << "       code = " << e.code() << " (\"" << e.what() << "\")" << std::endl;
        return false;
    }

    return true;
}

bool IntegralCrossSection::updateTable ()
{
    // Get merged available energies from all total angular momenta for every discrete transition and partial wave.
    int ni, li, mi, nf, lf, mf, S, ell;
    sqlitepp::statement st1 (session());
    st1 << "SELECT DISTINCT ni,li,mi,nf,lf,mf,S,ell FROM 'tmat'",
        sqlitepp::into(ni), sqlitepp::into(li), sqlitepp::into(mi),
        sqlitepp::into(nf), sqlitepp::into(lf), sqlitepp::into(mf),
        sqlitepp::into(S),  sqlitepp::into(ell);

    // for all transitions and partial waves
    while (st1.exec())
    {
        // get all merged energies for this transition and partial wave
        double E = 0;
        sqlitepp::statement st2 (session());
        st2 << "SELECT DISTINCT Ei FROM 'tmat' "
               "WHERE ni = :ni AND li = :li AND mi = :mi "
               "  AND nf = :nf AND lf = :lf AND mf = :mf "
               "  AND S  = :S  AND ell = :ell ORDER BY Ei ASC",
               sqlitepp::into(E),
               sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
               sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
               sqlitepp::use(S),  sqlitepp::use(ell);
        rArray merged_energies;
        while (st2.exec())
        {
            merged_energies.push_back(E);
        }

        // get all total angular momenta and parities for this transition and partial wave
        int L;
        sqlitepp::statement st3 (session());
        st3 << "SELECT DISTINCT L FROM 'tmat' "
            "WHERE ni = :ni AND li = :li AND mi = :mi "
            "  AND nf = :nf AND lf = :lf AND mf = :mf "
            "  AND S  = :S  AND ell = :ell ORDER BY L ASC",
            sqlitepp::into(L),
            sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
            sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
            sqlitepp::use(S),  sqlitepp::use(ell);
        iArray angular_momenta;
        while (st3.exec())
        {
            angular_momenta.push_back(L);
        }

        // merge-interpolate all T-matrix data for this transition and partial wave
        cArray merged_T (merged_energies.size());
        for (int L : angular_momenta)
        {
            rArray energies, Re_T, Im_T, k;
            // retrieve data
            double ret, imt;
            sqlitepp::statement st4 (session());
            st4 << "SELECT Ei,Re_T_ell,Im_T_ell FROM 'tmat' "
                "WHERE ni = :ni AND li = :li AND mi = :mi "
                "  AND nf = :nf AND lf = :lf AND mf = :mf "
                "  AND S  = :S  AND ell = :ell AND L = :L "
                "ORDER BY Ei ASC",
                sqlitepp::into(E), sqlitepp::into(ret), sqlitepp::into(imt),
                sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
                sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
                sqlitepp::use(S),  sqlitepp::use(ell), sqlitepp::use(L);
            while (st4.exec())
            {
                energies.push_back(E);
                k.push_back(std::sqrt(E));
                Re_T.push_back(ret * k.back());
                Im_T.push_back(imt * k.back());
            }

            // default to Akima spline interpolation, but fall back to linear if not enough points
            gsl_interp_type const * interp = gsl_interp_akima;
            bool enoughPoints = energies.size() >= gsl_interp_type_min_size(interp) or energies.size() >= gsl_interp_type_min_size(interp = gsl_interp_linear);

            // interpolate data
            gsl_interp* spline_Re = enoughPoints ? gsl_interp_alloc(interp, energies.size()) : nullptr;
            gsl_interp* spline_Im = enoughPoints ? gsl_interp_alloc(interp, energies.size()) : nullptr;
            if (spline_Re) gsl_interp_init(spline_Re, energies.data(), Re_T.data(), energies.size());
            if (spline_Im) gsl_interp_init(spline_Im, energies.data(), Im_T.data(), energies.size());
            # pragma omp parallel
            {
                gsl_interp_accel* accel_Re = spline_Re ? gsl_interp_accel_alloc() : nullptr;
                gsl_interp_accel* accel_Im = spline_Im ? gsl_interp_accel_alloc() : nullptr;

                double ret, imt;

                # pragma omp for
                for (std::size_t i = 0; i < merged_energies.size(); i++)
                {
                    // try to find this merged energy within the data for this angular momentum
                    std::size_t j = std::lower_bound(energies.begin(), energies.end(), merged_energies[i]) - energies.begin();
                    if (j < energies.size() and energies[j] == merged_energies[i])
                    {
                        // great, no need to interpolate -- let's just use the value
                        merged_T[i] += Complex(Re_T[j], Im_T[j]) / k[j];
                    }

                    // interpolate if within dataset
                    else if (energies.front() <= merged_energies[i] and merged_energies[i] <= energies.back() and enoughPoints)
                    {
                        int err_re = gsl_interp_eval_e(spline_Re, energies.data(), Re_T.data(), merged_energies[i], accel_Re, &ret);
                        int err_im = gsl_interp_eval_e(spline_Im, energies.data(), Im_T.data(), merged_energies[i], accel_Im, &imt);

                        if (err_re != GSL_SUCCESS)
                        {
                            std::cout << "Warning: Failed to interpolate real data for transition (" << ni << "," << li << "," << mi << ") -> (" << nf << "," << lf << "," << mf << "), "
                                "L = " << L << ", ell = " << ell << ", Ei = " << merged_energies[i] << " (" << gsl_strerror(err_re) << ")" << std::endl;
                        }

                        if (err_im != GSL_SUCCESS)
                        {
                            std::cout << "Warning: Failed to interpolate imag data for transition (" << ni << "," << li << "," << mi << ") -> (" << nf << "," << lf << "," << mf << "), "
                                "L = " << L << ", ell = " << ell << ", Ei = " << merged_energies[i] << " (" << gsl_strerror(err_re) << ")" << std::endl;
                        }

                        merged_T[i] += Complex(ret, imt) / std::sqrt(merged_energies[i]);
                    }
                }

                if (accel_Re) gsl_interp_accel_free(accel_Re);
                if (accel_Im) gsl_interp_accel_free(accel_Im);
            }
            if (spline_Re) gsl_interp_free(spline_Re);
            if (spline_Im) gsl_interp_free(spline_Im);
        }

        // write interpolated data to the database
        sqlitepp::statement st5 (session());
        double ics, Ei, Ef;
        std::size_t rows = 0;
        st5 << "INSERT OR REPLACE INTO 'ics' VALUES ";
        for (std::size_t i = 0; i < merged_energies.size(); i++)
        {
            // initial and final energies
            Ei = merged_energies[i];
            Ef = Ei - 1./(ni*ni) + 1./(nf*nf);

            // skip forbidden channels
            if (Ef < 0)
                continue;

            // calculate partial integral cross section
            ics = std::sqrt(Ef/Ei) * (2 * S + 1) * sqrabs(merged_T[i]) / std::pow(4 * special::constant::pi, 2);

            // separate data n-tuples with comma
            if (rows++ > 0)
                st5.q() << ",";

            // add another data n-tuple
            st5.q() << " (" << ni << "," << li << "," << mi << "," << nf << "," << lf << "," << mf << "," << S << "," << Ei << "," << ell << "," << ics << ")";
        }
        if (rows > 0)
        {
            // add all data in one transaction
            st5.exec();
        }
    }

    // Insert ionization (no interpolation needed: L is used as the partial wave number here.

    // NOTE: The factor (Ei - 1/ni²) = E_tot is the one that IONCS leaves out; see the
    //       note at @ref db_ioncs.
    sqlitepp::statement st6 (session());
    st6 << "INSERT OR REPLACE INTO 'ics' "
            "SELECT ni, li, mi, "
                   "0,  0,  0,  "
                   "S,  Ei, L,  "
                   "SUM(0.25*(2*S+1)*(Ei-1.0/(ni*ni))*IONCS(cheb)/SQRT(Ei)) "
                "FROM 'ionf' "
                "GROUP BY ni, li, mi, S, Ei, L";
    st6.exec();

    return ScatteringQuantity::updateTable();
}

// --------------------------------------------------------------------------------- //

bool IntegralCrossSection::run (std::map<std::string,std::string> const & sdata)
{
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
    double lfactor = change_units(lUnit_au, Lunits);

    // scattering event parameters
    int ni = Conv<int>(sdata, "ni", name());
    int li = Conv<int>(sdata, "li", name());
    int mi0= Conv<int>(sdata, "mi", name());
    int nf = Conv<int>(sdata, "nf", name());
    int lf = Conv<int>(sdata, "lf", name());
    int mf0= Conv<int>(sdata, "mf", name());
    int ell= Conv<int>(sdata, "ell",name());
    int  S = Conv<int>(sdata, "S",  name());

    // use mi >= 0; if mi < 0, flip both signs
    int mi = (mi0 < 0 ? -mi0 : mi0);
    int mf = (mi0 < 0 ? -mf0 : mf0);

    // energies and cross sections
    double E, sigma;
    rArray energies, E_arr, sigma_arr;

    // get energy / energies
    try
    {
        // is there a single energy specified using command line ?
        energies.push_back(Conv<double>(sdata, "Ei", name()));
    }
    catch (std::exception const & e)
    {
        // are there more energies specified using the STDIN ?
        energies = readStandardInput<double>();
    }

    // Ionization summed over the partial waves.
    //   Asked for by a negative 'ell'. The partial cross sections are added up with the
    //   geometric extrapolation of what lies past the last of them; see the comment
    //   above sum_partial_waves. The estimate of that remainder is written out next to
    //   the cross section, being also the estimate of the error that is left in it.
    if (nf == 0 and lf == 0 and mf == 0 and ell < 0)
    {
        // ionization threshold
        double Eion = 1./(ni*ni);

        // number of partial waves available in the database
        int max_ell = -1;
        sqlitepp::statement ell_st (session());
        ell_st << "SELECT MAX(ell) FROM 'ics' "
                  "WHERE ni = :ni AND li = :li AND mi = :mi "
                  "  AND nf = 0 AND lf = 0 AND mf = 0 AND S = :S",
            sqlitepp::into(max_ell),
            sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi), sqlitepp::use(S);
        ell_st.exec();

        // retrieve the partial cross sections, one dataset per partial wave
        rArrays E_data, sigma_data;
        for (int l = 0; l <= max_ell; l++)
        {
            sqlitepp::statement st (session());
            st << "SELECT Ei, SUM(sigma) FROM 'ics' "
                  "WHERE ni = :ni AND li = :li AND mi = :mi "
                  "  AND nf = 0 AND lf = 0 AND mf = 0 "
                  "  AND  S = :S AND ell = :ell "
                  "GROUP BY Ei "
                  "ORDER BY Ei ASC",
                sqlitepp::into(E), sqlitepp::into(sigma),
                sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
                sqlitepp::use(S), sqlitepp::use(l);

            E_data.push_back(rArray());
            sigma_data.push_back(rArray());

            while (st.exec())
            {
                E_data.back().push_back(E);
                sigma_data.back().push_back(sigma);
            }
        }

        // merge the energy grids of all partial waves
        rArray E_all = join(E_data);
        std::sort(E_all.begin(), E_all.end());
        E_all.resize(std::unique(E_all.begin(), E_all.end()) - E_all.begin());

        // write header
        std::cout << logo("#") <<
            "# Integral ionization cross section in " << unit_name(Lunits) << " for\n" <<
            "#     ni = " << ni << ", li = " << li << ", mi = " << mi0 << ",\n" <<
            "#     S = " << S << ", summed over the partial waves 0 to " << max_ell << "\n" <<
            "# ordered by energy in " << unit_name(Eunits) << "\n" <<
            "#\n" <<
            "# The last column is the geometric extrapolation of the partial waves past the\n" <<
            "# last computed one. It is included in the cross section and is at the same time\n" <<
            "# the estimate of the truncation error left in it; 'nan' means that no partial\n" <<
            "# wave was found to extrapolate from and that the sum is bare.\n" <<
            "#\n";
        OutputTable table;
        table.setWidth(15, 15, 15);
        table.setAlignment(OutputTable::left, OutputTable::left, OutputTable::left);
        table.write("# E        ", "sigma    ", "tail     ");
        table.write("# ---------", "---------", "---------");

        // terminate if no data
        if (E_all.empty())
            return true;

        // energies to write out; a negative energy asks for all of them
        rArray Eout = (energies.front() < 0. ? E_all : energies * efactor);

        // put every partial wave on the requested energies
        // NOTE: A partial wave only contributes at the energies its own dataset covers.
        //       Without saying so a partial wave stored at a single energy would be
        //       carried over to every other energy by interpolate_real, which returns
        //       the lone value it was given whatever it is asked for; the sum would
        //       then quietly mix partial waves belonging to different impact energies.
        rArrays sig;
        for (int l = 0; l <= max_ell; l++)
        {
            rArray s = interpolate_real(E_data[l], E_data[l] * sigma_data[l], Eout, gsl_interp_akima) / Eout;

            for (std::size_t i = 0; i < Eout.size(); i++)
            {
                if (E_data[l].empty() or Eout[i] < E_data[l].front() or Eout[i] > E_data[l].back())
                    s[i] = 0.;
            }

            sig.push_back(s);
        }

        // sum over the partial waves, energy by energy
        for (std::size_t i = 0; i < Eout.size(); i++)
        {
            rArray sigma_ell (max_ell + 1);
            for (int l = 0; l <= max_ell; l++)
                sigma_ell[l] = sig[l][i];

            double tail = 0;
            double total = sum_partial_waves(sigma_ell, Eout[i], Eion, tail);

            table.write(Eout[i] / efactor, total * lfactor * lfactor, tail * lfactor * lfactor);
        }

        return true;
    }

    // compose query
    sqlitepp::statement st (session());
    st << "SELECT Ei, sigma FROM 'ics' "
            "WHERE ni = :ni "
            "  AND li = :li "
            "  AND mi = :mi "
            "  AND nf = :nf "
            "  AND lf = :lf "
            "  AND mf = :mf "
            "  AND ell= :ell"
            "  AND  S = :S  "
            "ORDER BY Ei ASC",
        sqlitepp::into(E), sqlitepp::into(sigma),
        sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
        sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
        sqlitepp::use(ell), sqlitepp::use(S);

    // retrieve data
    while (st.exec())
    {
        E_arr.push_back(E);
        sigma_arr.push_back(sigma);
    }

    // write header
    std::cout << logo("#") <<
        "# Integral cross section in " << unit_name(Lunits) << " for\n" <<
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi0 << ",\n" <<
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf0 << ",\n" <<
        "#     ell = " << ell << ", S = " << S << "\n" <<
        "# ordered by energy in " << unit_name(Eunits) << "\n" <<
        "#\n";
    OutputTable table;
    table.setWidth(15);
    table.setAlignment(OutputTable::left);
    table.write("# E        ", "sigma    ");
    table.write("# ---------", "---------");

    // terminate if no data
    if (E_arr.empty())
        return true;

    if (energies[0] < 0.)
    {
        // negative energy indicates full output
        for (std::size_t i = 0; i < E_arr.size(); i++)
            table.write(E_arr[i] / efactor, sigma_arr[i] * lfactor * lfactor);
    }
    else
    {
        // threshold for ionization
        double Eion = 1./(ni*ni);

        // interpolate (linear below ionization threshold, cspline above)
        rArray ics = (efactor * energies.front() < Eion) ? 
            interpolate_real(E_arr, E_arr * sigma_arr, energies * efactor, gsl_interp_linear) / (energies * efactor) :
            interpolate_real(E_arr, E_arr * sigma_arr, energies * efactor, gsl_interp_akima ) / (energies * efactor);

        // output
        for (std::size_t i = 0; i < energies.size(); i++)
            table.write(energies[i], ics[i] * lfactor * lfactor);
    }

    return true;
}
