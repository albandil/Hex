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

#include "../interpolate.h"
#include "../born.h"
#include "../chebyshev.h"
#include "../special.h"
#include "../variables.h"
#include "../version.h"

const std::string ScatteringAmplitude::Id = "scatamp";
const std::string ScatteringAmplitude::Description = "Scattering amplitude.";
const std::vector<std::pair<std::string,std::string>> ScatteringAmplitude::Dependencies = {
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
const std::vector<std::string> ScatteringAmplitude::VecDependencies = { "theta" };

bool ScatteringAmplitude::initialize (sqlitepp::session & db) const
{
    return true;
}

std::vector<std::string> const & ScatteringAmplitude::SQL_Update () const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

std::vector<std::string> const & ScatteringAmplitude::SQL_CreateTable () const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

cArray FirstBornFullTMatrix (int ni, int li, int mi, int nf, int lf, int mf, double E, const rArrayView angles)
{
    // get result of the symbolic integration
    auto Wb = Wb_symb_in(ni, li, mi, nf, lf, mf);
    
    // evaluate symbolic expression for all angles
    cArray tmat (angles.size());
    for (unsigned i = 0; i < angles.size(); i++)
    {
        // initial projectile linear momentum
        geom::vec3d ki = { 0., 0., std::sqrt(E) };
        
        // fianl projectile linear momentum
        double magkf = std::sqrt(E - 1./(ni*ni) + 1./(nf*nf));
        geom::vec3d kf = { magkf * std::sin(angles[i]), 0., magkf * std::cos(angles[i]) };
        
        // evaluate T-matrix
        tmat[i] = eval_Wb(Wb, 1./ni + 1./nf, kf - ki);
        
        // filter nan-s
        if (not Complex_finite(tmat[i]))
            tmat[i] = 0.;
    }
    
    return tmat;
}

void hex_scattering_amplitude_Born_
(
    int * ni, int * li, int * mi,
    int * nf, int * lf, int * mf,
    int * S, double * E, int * N,
    double * angles, double * result
)
{
    cArrayView results(*N,reinterpret_cast<Complex*>(result));
    results.fill(0.);
    
    // we need cosines, not angles
    rArray cos_angles = rArray(*N,angles).transform([](double x){ return std::cos(x); });
    
    // get lowest and highest partial wave
    int min_ell, max_ell;
    sqlitepp::statement st(db);
    st << "SELECT MIN(ell), MAX(ell) FROM " + TMatrix::Id + " "
           "WHERE ni = :ni "
           "  AND li = :li "
           "  AND mi = :mi "
           "  AND nf = :nf "
           "  AND lf = :lf "
           "  AND mf = :mf "
           "  AND  S = :S  ",
        sqlitepp::into(min_ell), sqlitepp::into(max_ell),
        sqlitepp::use(*ni), sqlitepp::use(*li), sqlitepp::use(*mi),
        sqlitepp::use(*nf), sqlitepp::use(*lf), sqlitepp::use(*mf),
        sqlitepp::use(*S);
    st.exec();
    
    // resulting T-matrices
    cArray tmatrices;
    
    // loop over all partial waves
    for (int ell = min_ell; ell <= max_ell; ell++)
    {
        Complex Tmat;
        double Ei, ReT, ImT;
        
        // prepare selection statement
        sqlitepp::statement st1(db);
        st1 << "SELECT Ei, SUM(Re_TBorn_ell), SUM(Im_TBorn_ell) FROM " + TMatrix::Id + " "
            "WHERE ni = :ni "
            "  AND li = :li "
            "  AND mi = :mi "
            "  AND nf = :nf "
            "  AND lf = :lf "
            "  AND mf = :mf "
            "  AND  S = :S  "
            "  AND ell = :ell "
            "GROUP BY Ei "
            "ORDER BY Ei ASC ",
            sqlitepp::into(Ei), sqlitepp::into(ReT), sqlitepp::into(ImT),
            sqlitepp::use(*ni), sqlitepp::use(*li), sqlitepp::use(*mi),
            sqlitepp::use(*nf), sqlitepp::use(*lf), sqlitepp::use(*mf),
            sqlitepp::use(*S), sqlitepp::use(ell);
        
        // load the data corresponding to the statement
        rArray energies_ell;
        cArray tmatrices_ell;
        while (st1.exec())
        {
            energies_ell.push_back(Ei);
            tmatrices_ell.push_back(Complex(ReT,ImT));
        }
        
        // is the requested energy in the available energy interval?
        if (energies_ell.size() != 0 and energies_ell.front() <= *E and *E <= energies_ell.back())
        {
            // interpolate requested energy
            tmatrices.push_back(interpolate(energies_ell, tmatrices_ell, { *E })[0]);
            
            std::cout << "# " << ell << " " << tmatrices.back() << std::endl;
        }
        else
            break;
        
        // update scattering amplitudes¨
        for (int i = 0; i < *N; i++)
        {
            double Y = gsl_sf_legendre_sphPlm(ell,std::abs((*mi)-(*mf)),cos_angles[i]);
            results[i] += -tmatrices.back() * Y / special::constant::two_pi;
        }
    }
}

void hex_scattering_amplitude_
(
    int * ni, int * li, int * mi,
    int * nf, int * lf, int * mf,
    int * S, double * E, int * N,
    double * angles, double * result
)
{
    cArrayView results(*N,reinterpret_cast<Complex*>(result));
    results.fill(0.);
    
    // we need cosines, not angles
    rArray cos_angles = rArray(*N,angles).transform([](double x){ return std::cos(x); });
    
    // get lowest and highest partial wave
    int min_ell, max_ell;
    sqlitepp::statement st(db);
    st << "SELECT MIN(ell), MAX(L)-lf FROM " + TMatrix::Id + " "
           "WHERE ni = :ni "
           "  AND li = :li "
           "  AND mi = :mi "
           "  AND nf = :nf "
           "  AND lf = :lf "
           "  AND mf = :mf "
           "  AND  S = :S  ",
        sqlitepp::into(min_ell), sqlitepp::into(max_ell),
        sqlitepp::use(*ni), sqlitepp::use(*li), sqlitepp::use(*mi),
        sqlitepp::use(*nf), sqlitepp::use(*lf), sqlitepp::use(*mf),
        sqlitepp::use(*S);
    st.exec();
    
    // resulting T-matrices
    cArray tmatrices, tmatricesBorn;
    
    // loop over all partial waves
    for (int ell = min_ell; ; ell++)
    {
        Complex Tmat;
        double Ei, ReT, ImT;
        bool missing = true;
        
        // check if there are data in the database
        if (ell <= max_ell)
        {
            // prepare selection statement
            sqlitepp::statement st1(db);
            // WARNING Sign convention changed.
//             st1 << "SELECT Ei, SUM(-Re_T_ell - Re_TBorn_ell), SUM(-Im_T_ell - Im_TBorn_ell) FROM " + TMatrix::Id + " "
            st1 << "SELECT Ei, SUM(Re_T_ell), SUM(Im_T_ell) FROM " + TMatrix::Id + " "
                "WHERE ni = :ni "
                "  AND li = :li "
                "  AND mi = :mi "
                "  AND nf = :nf "
                "  AND lf = :lf "
                "  AND mf = :mf "
                "  AND  S = :S  "
                "  AND ell = :ell "
                "GROUP BY Ei "
                "ORDER BY Ei ASC ",
                sqlitepp::into(Ei), sqlitepp::into(ReT), sqlitepp::into(ImT),
                sqlitepp::use(*ni), sqlitepp::use(*li), sqlitepp::use(*mi),
                sqlitepp::use(*nf), sqlitepp::use(*lf), sqlitepp::use(*mf),
                sqlitepp::use(*S), sqlitepp::use(ell);
            
            // load the data corresponding to the statement
            rArray energies_ell;
            cArray tmatrices_ell;
            while (st1.exec())
            {
                energies_ell.push_back(Ei);
                tmatrices_ell.push_back(Complex(ReT,ImT));
            }
            
            // is the requested energy in the available energy interval?
            if (energies_ell.size() != 0 and energies_ell.front() <= *E and *E <= energies_ell.back())
            {
                // interpolate requested energy
                tmatrices.push_back(interpolate(energies_ell, tmatrices_ell, { *E })[0]);
                missing = false;
                
                if (not missing)
                    std::cout << "# retrieve " << ell << " " << tmatrices.back() << std::endl;
            }
        }
        
//         if (missing) break;
        
        // if there are no more partial waves in the database, try to extrapolate
        if (ell > max_ell or missing)
        {
            // only extrapolate if we have two or more samples, that are decreasing at the end
            if
            (
                tmatrices.size() >= 2 and
                std::abs(tmatrices.back(0).real()) < std::abs(tmatrices.back(1).real()) and
                std::abs(tmatrices.back(0).imag()) < std::abs(tmatrices.back(1).imag())
            )
            {
                // choose the asymptotics
                if (*ni == *nf and *li == *lf and *mi == *mf)
                {
                    // extrapolate T_ell for elastic scattering [~ L⁻³]
                    tmatrices.push_back(tmatrices.back() * gsl_sf_pow_int((ell-1.)/ell,3));
                }
                else /*if (std::abs((*ni) - (*nf)) == 1 and std::abs((*li) - (*lf)) == 1)
                {
                    // dipole transitions are special, they will need to be Born-subtracted
                    break;
                }
                else*/
                {
                    // extrapolate T_ell for inelastic scattering [use geometric series]
                    double extra_re = tmatrices.back(0).real() * tmatrices.back(0).real() / tmatrices.back(1).real();
                    double extra_im = tmatrices.back(0).imag() * tmatrices.back(0).imag() / tmatrices.back(1).imag();
                    tmatrices.push_back(Complex(extra_re,extra_im));
                    std::cout << "# extrapolate " << ell << " " << tmatrices.back() << std::endl;
                }
            }
            else
            {
                // unable to extrapolate => stop
                break;
            }
        }
        
        // update scattering amplitudes¨
        for (int i = 0; i < *N; i++)
        {
            double Y = gsl_sf_legendre_sphPlm(ell,std::abs((*mi)-(*mf)),cos_angles[i]);
            results[i] += -tmatrices.back() * Y / special::constant::two_pi;
        }
        
        // stop if the T-matrix is small
        if (std::abs(tmatrices.back()) < 1e-5 * std::abs(sum(tmatrices)))
            break;
    }
    
    // finally, Born-subtract dipole transition
    /*if (std::abs((*ni) - (*nf)) == 1 and std::abs((*li) - (*lf)) == 1)
        results += -FirstBornFullTMatrix(*ni, *li, *mi, *nf, *lf, *mf, *E, rArrayView(*N,angles)) / special::constant::two_pi;*/
}

bool ScatteringAmplitude::run (std::map<std::string,std::string> const & sdata) const
{
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
    double lfactor = change_units(lUnit_au, Lunits);
    double afactor = change_units(Aunits, aUnit_rad);
    
    // scattering event parameters
    int ni = Conv<int>(sdata, "ni", Id);
    int li = Conv<int>(sdata, "li", Id);
    int mi = Conv<int>(sdata, "mi", Id);
    int nf = Conv<int>(sdata, "nf", Id);
    int lf = Conv<int>(sdata, "lf", Id);
    int mf = Conv<int>(sdata, "mf", Id);
    int  S = Conv<int>(sdata, "S", Id);
    double E = Conv<double>(sdata, "Ei", Id) * efactor;
    
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
    
    // the scattering amplitudes
    rArray scaled_angles = angles * afactor;
    cArray amplitudes(angles.size());
    hex_scattering_amplitude
    (
        ni,li,mi,
        nf,lf,mf,
        S, E,
        angles.size(),
        scaled_angles.data(),
        reinterpret_cast<double*>(amplitudes.data())
    );
    
    int N = angles.size();
    cArray amplitudesBorn(angles.size());
    hex_scattering_amplitude_Born_
    (
        &ni,&li,&mi,
        &nf,&lf,&mf,
        &S, &E,
        &N,
        scaled_angles.data(),
        reinterpret_cast<double*>(amplitudesBorn.data())
    );
    
    cArray bampl = -FirstBornFullTMatrix(ni, li, mi, nf, lf, mf, E, scaled_angles) / special::constant::two_pi;
    
    // write out
    std::cout << logo("#") <<
        "# Scattering amplitudes in " << unit_name(Lunits) << " for\n"
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n"
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n"
        "#     S = " << S << ", E = " << E/efactor << unit_name(Eunits) << "\n"
        "# ordered by angle in " << unit_name(Aunits) << "\n"
        "# \n";
    OutputTable table;
    table.setWidth(15, 15);
    table.setAlignment(OutputTable::left);
    table.write("# angle    ", "Re f     ", "Im f     ", "Re f [DB]", "Im f [DB]", "Re f [1B]", "Im f [1B]");
    table.write("# ---------", "---------", "---------", "---------", "---------", "---------", "---------");
    
    for (std::size_t i = 0; i < angles.size(); i++)
    {
        table.write
        (
            angles[i],
            amplitudes[i].real() * lfactor,
            amplitudes[i].imag() * lfactor,
            amplitudesBorn[i].real() * lfactor,
            amplitudesBorn[i].imag() * lfactor,
            bampl[i].real() * lfactor,
            bampl[i].imag() * lfactor
        );
    }
    
    return true;
}
