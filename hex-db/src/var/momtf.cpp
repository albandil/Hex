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

#include "hex-special.h"
#include "hex-interpolate.h"
#include "hex-version.h"

#include "variables.h"

const std::string MomentumTransfer::Id = "momtf";
const std::string MomentumTransfer::Description = "Momentum transfer.";
const std::vector<std::pair<std::string,std::string>> MomentumTransfer::Dependencies = {
    {"ni", "Initial atomic principal quantum number."},
    {"li", "Initial atomic orbital quantum number."},
    {"mi", "Initial atomic magnetic quantum number."},
    {"nf", "Final atomic principal quantum number."},
    {"lf", "Final atomic orbital quantum number."},
    {"mf", "Final atomic magnetic quantum number."},
    {"S", "Total spin of atomic + projectile electron."},
    {"Ei", "Projectile impact energy (Rydberg)."}
};
const std::vector<std::string> MomentumTransfer::VecDependencies = { "Ei" };

bool MomentumTransfer::initialize(sqlitepp::session & db) const
{
    return true;
}

std::vector<std::string> const & MomentumTransfer::SQL_CreateTable () const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

std::vector<std::string> const & MomentumTransfer::SQL_Update () const
{
    static const std::vector<std::string> cmd;
    return cmd;
}

bool MomentumTransfer::run (std::map<std::string,std::string> const & sdata) const
{
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
    double lfactor = change_units(lUnit_au, Lunits);
    
    // atomic and projectile data
    int ni = Conv<int>(sdata, "ni", Id);
    int li = Conv<int>(sdata, "li", Id);
    int mi0= Conv<int>(sdata, "mi", Id);
    int nf = Conv<int>(sdata, "nf", Id);
    int lf = Conv<int>(sdata, "lf", Id);
    int mf0= Conv<int>(sdata, "mf", Id);
    int  S = Conv<int>(sdata,  "S", Id);
    
    // use mi >= 0; if mi < 0, flip both signs
    int mi = (mi0 < 0 ? -mi0 : mi0);
    int mf = (mi0 < 0 ? -mf0 : mf0);
    
    // energies
    rArray energies;
    
    // get energy / energies
    try {
        
        // is there a single energy specified using command line ?
        energies.push_back(Conv<double>(sdata, "Ei", Id));
        
    } catch (std::exception e) {
        
        // are there more energies specified using the STDIN ?
        energies = readStandardInput<double>();
    }
    
    // SQL interface variables
    double E, Re_T_ell, Im_T_ell;
    int L, ell;
    
    // compose the statement
    sqlitepp::statement st(db);
    st << "SELECT Ei, L, ell, Re_T_ell, Im_T_ell FROM " + TMatrix::Id + " "
          "WHERE ni = :ni "
          "  AND li = :li "
          "  AND mi = :mi "
          "  AND nf = :nf "
          "  AND lf = :lf "
          "  AND mf = :mf "
          "  AND  S = :S  "
          "ORDER BY L, ell, Ei ASC",
        sqlitepp::into(E), sqlitepp::into(L), sqlitepp::into(ell), 
        sqlitepp::into(Re_T_ell), sqlitepp::into(Im_T_ell),
        sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
        sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
        sqlitepp::use(S);
    
    // auxiliary variables
    rArray *pE = nullptr;    // pointer to current energy array
    cArray *pT = nullptr;    // pointer to current T-matrix array
    std::pair<int,int> key(-1,-1);    // current (L,ell) pair
    int maxL = -1;    // maximal L
    int maxl = -1;    // maximal ell
        
    // load T-matrices
    std::map<std::pair<int,int>, std::pair<rArray*, cArray*>> Tmatrices;
    while (st.exec())
    {
        // create new arrays if necessary
        if (key.first != L or key.second != ell)
        {
            pE = new rArray();
            pT = new cArray();
            key = std::make_pair(L,ell);
            Tmatrices[key] = std::make_pair(pE,pT);
            maxL = (L > maxL) ? L : maxL;
            maxl = (ell > maxl) ? ell : maxl;
        }
        
        // insert new data
        pE->push_back(E);
        pT->push_back(Complex(Re_T_ell,Im_T_ell));
    }
    
    // terminate if no data
    if (Tmatrices.empty())
        return true;
    
    // compute momentum transfer
    rArray eta_energies;
    rArray eta;
    for (int L = 0; L <= maxL; L++)
    for (int Lp = 0; Lp <= maxL; Lp++)
    for (int l = 0; l <= maxl; l++)
    for (int lp = 0; lp <= maxl; lp++)
    {
        std::pair<int,int> key(L,l), keyp(Lp,lp);
        
        // if no data, skip
        if (Tmatrices.find(key) == Tmatrices.end() or Tmatrices.find(keyp) == Tmatrices.end())
            continue;
        
        // get energies and T-matrices pointers
        rArray *pE = Tmatrices[key].first;
        cArray *pT = Tmatrices[key].second;
        rArray *pEp = Tmatrices[keyp].first;
        cArray *pTp = Tmatrices[keyp].second;
        
        // normalize
        cArray emptyp(pEp->size());    // zero-filled ghost
        merge (*pE, *pT, *pEp, emptyp);    // inflate *pE and *pT to accomodate *pEp and emptyp
        cArray empty(pE->size()); // zero-filled ghost
        merge (*pEp, *pTp, *pE, empty);
        
        // compute factor
        double factor = 0;
        for (int m = -l; m <= l; m++)
            factor += special::ClebschGordan(l,m,1,0,lp,m);
        factor *= -sqrt((2.*l+1.)/(2.*lp+1.)) * special::ClebschGordan(l,0,1,0,lp,0);
        if (l == lp)
            factor += 1.;
        
        // update momentum transfer
        merge (eta_energies, eta, *pE, factor * abs((*pT) * (*pTp).conj()));
    }
    eta *= 1. / std::pow(2 * special::constant::pi, 2);
    
    // write header
    std::cout << logo("#") <<
        "# Momentum transfer in " << unit_name(Lunits) << " for\n" <<
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
        "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n" <<
        "#     S = " << S << "\n" <<
        "# ordered by energy in " << unit_name(Eunits) << "\n" <<
        "#\n";
    OutputTable table;
    table.setWidth(15, 15);
    table.setAlignment(OutputTable::left);
    table.write("# E        ", "momtrans.");
    table.write("# ---------", "---------");
    
    if (energies[0] < 0.)
    {
        // negative energy indicates full output
        for (std::size_t i = 0; i < eta_energies.size(); i++)
            table.write(eta_energies[i] / efactor, eta[i] * lfactor * lfactor);
    }
    else
    {
        // interpolate
        eta = interpolate (eta_energies, eta, energies * efactor);
        
        // output
        for (std::size_t i = 0; i < energies.size(); i++)
            table.write(energies[i], eta[i] * lfactor * lfactor);
    }
    
    return true;
}
