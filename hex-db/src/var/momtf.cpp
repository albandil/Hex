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

#include <deque>
#include <map>
#include <string>
#include <vector>

// --------------------------------------------------------------------------------- //

#include "hex-special.h"
#include "hex-interpolate.h"
#include "hex-version.h"

// --------------------------------------------------------------------------------- //

#include "../quantities.h"
#include "../utils.h"

// --------------------------------------------------------------------------------- //

createNewScatteringQuantity(MomentumTransfer, "momtf")

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

std::string MomentumTransfer::description ()
{
    return "Momentum transfer.";
}

std::vector<std::string> MomentumTransfer::dependencies ()
{
    return std::vector<std::string>
    {
        "tmat"
    };
}

std::vector<std::pair<std::string,std::string>> MomentumTransfer::params ()
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
        {"Ei", "Projectile impact energy (Rydberg)."}
    };
}

std::vector<std::string> MomentumTransfer::vparams ()
{
    return std::vector<std::string>
    {
        "Ei"
    };
}

// --------------------------------------------------------------------------------- //

bool MomentumTransfer::initialize(sqlitepp::session & db)
{
    return ScatteringQuantity::initialize(db);
}

bool MomentumTransfer::createTable ()
{
    return ScatteringQuantity::createTable();
}

bool MomentumTransfer::updateTable ()
{
    return ScatteringQuantity::updateTable();
}

// --------------------------------------------------------------------------------- //

bool MomentumTransfer::run (std::map<std::string,std::string> const & sdata)
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
    int  S = Conv<int>(sdata,  "S", name());
    
    // use mi >= 0; if mi < 0, flip both signs
    int mi = (mi0 < 0 ? -mi0 : mi0);
    int mf = (mi0 < 0 ? -mf0 : mf0);
    
    // energies
    rArray energies;
    
    // get energy / energies
    try {
        
        // is there a single energy specified using command line ?
        energies.push_back(Conv<double>(sdata, "Ei", name()));
        
    } catch (std::exception e) {
        
        // are there more energies specified using the STDIN ?
        energies = readStandardInput<double>();
    }
    
    // convert energies to Ry
    rArray scaled_energies = energies * efactor;
    
    // get available energies if requested
    if (not scaled_energies.empty() and scaled_energies.front() == -1)
    {
        scaled_energies.clear();
        double E;
        sqlitepp::statement st (session());
        st << "SELECT DISTINCT Ei FROM 'tmat' "
              "WHERE ni = :ni "
              "  AND li = :li "
              "  AND mi = :mi "
              "  AND nf = :nf "
              "  AND lf = :lf "
              "  AND mf = :mf "
              "ORDER BY Ei ASC",
              sqlitepp::into(E),
              sqlitepp::use(ni), sqlitepp::use(li),  sqlitepp::use(mi),
              sqlitepp::use(nf), sqlitepp::use(lf),  sqlitepp::use(mf);
        while (st.exec())
        {
            scaled_energies.push_back(E);
        }
    }
    
    // output data
    rArray eta (scaled_energies.size()), eta_ex (scaled_energies.size());
    
    // T-matrices for the last three partial waves
    std::deque<int> ells;
    std::deque<iArrays> conv;
    std::deque<iArrays> cmpl;
    std::deque<cArrays> tmat;
    
    // callback for collecting T-matrices (called for every partial wave starting with ell = |mi - mf|)
    auto callback = [&]
    (
        int ell,
        iArrays const & converged,
        iArrays const & complete,
        cArrays const & tmatrices
    )
    {
        ells.push_back(ell);
        conv.push_back(converged);
        cmpl.push_back(complete);
        tmat.push_back(tmatrices);
        
        while (ells.size() > 3)
        {
            ells.pop_front();
            conv.pop_front();
            cmpl.pop_front();
            tmat.pop_front();
        }
        
        if (ells.size() >= 2)
        {
            // sum over
            //   l = |mi-mf|, |mi-mf|+1, ...
            // and
            //   lp = l - 1, l, l + 1
            
            int il = ells.size() - 2;
            int l = ells[il];
            
            for (unsigned ilp = 0; ilp < ells.size(); ilp++)
            {
                int lp = ells[ilp];
                
                // angular integral
                double angintg = (
                    l == lp ?
                    1.0 :
                    std::sqrt((2*l + 1.) / (2*lp + 1.)) * special::ClebschGordan(l,mi-mf,1,0,lp,mi-mf) * special::ClebschGordan(l,0,1,0,lp,0)
                );
                
                for (unsigned i = 0; i < scaled_energies.size(); i++)
                {
                    // update momentum transfer with this partial wave if contained in database (and not yet converged)
                    if (not conv[il][S][i] and cmpl[il][S][i])
                        eta[i] += angintg * ( tmat[il][S][i] * std::conj(tmat[ilp][S][i]) ).real();
                    
                    // update extrapolated momentum transfer always
                    if (true)
                        eta_ex[i] += angintg * ( tmat[il][S][i] * std::conj(tmat[ilp][S][i]) ).real();
                }
            }
        }
    };
    
    // call the T-matrix retrieval driver
    hex_tmat_pw_transform
    (
        ni, li, mi, nf, lf, mf,
        scaled_energies.size(), scaled_energies.data(),
        true,
        callback
    );
    
    // momenta
    rArray ki = sqrt(scaled_energies);
    rArray kf = sqrt(scaled_energies - 1.0 / (ni * ni) + 1.0 / (nf * nf));
    
    // add prefactor
    eta *= kf / ki * (2 * S + 1) / 157.9136704174297;
    eta_ex *= kf / ki * (2 * S + 1) / 157.9136704174297;
    
    // write header
    std::cout << logo("#") <<
        "# Momentum transfer in " << unit_name(Lunits) << " for\n" <<
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi0 << ",\n" <<
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf0 << ",\n" <<
        "#     S = " << S << "\n" <<
        "# ordered by energy in " << unit_name(Eunits) << "\n" <<
        "#\n";
    OutputTable table;
    table.setWidth(20, 20, 20);
    table.setAlignment(OutputTable::left);
    table.write("# E        ", "momtransfer", "momtransfer [ex]");
    table.write("# ---------", "-----------", "----------------");
    
    for (unsigned i = 0; i < scaled_energies.size(); i++)
    {
        table.write
        (
            scaled_energies[i] / efactor,
            std::isfinite(eta[i]) ? eta[i] * lfactor * lfactor : 0,
            std::isfinite(eta_ex[i]) ? eta_ex[i] * lfactor * lfactor : 0
        );
    }
    
    return true;
}
