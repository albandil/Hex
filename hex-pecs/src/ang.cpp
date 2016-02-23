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

#include "ang.h"

AngularBasis::AngularBasis (CommandLine const & cmd, InputFile const & inp)
    : L_(inp.L), S_(0), Pi_(inp.Pi), nL_(inp.nL), maxlambda_(inp.L + 2 * inp.nL), maxell_(nL_ + L_ + Pi_)
{
    std::cout << "Setting up the angular states..." << std::endl;
    
    // for given L, Π and levels list all available (ℓ₁ℓ₂) pairs
    for (int ell = 0; ell <= inp.nL; ell++)
    {
        std::cout << "\t- [" << ell << "] ";
        
        // get sum of the angular momenta for this angular level
        int sum = 2 * ell + inp.L + inp.Pi;
        
        // for all angular momentum pairs that do compose L
        for (int l1 = ell; l1 <= sum - ell; l1++)
        {
            int l2 = sum - l1;
            if (std::abs(l1 - l2) <= inp.L and inp.L <= l1 + l2)
            {
                std::cout << "(" << l1 << "," << l2 << ") ";
                states_.push_back(std::make_pair(l1, l2));
            }
        }
        std::cout << std::endl;
    }
    
    // set coupled groups
    if (cmd.fully_coupled)
    {
        std::cout << "\t- using fully coupled system" << std::endl;
        groups_.push_back(linspace<unsigned>(0, states_.size() - 1, states_.size()));
    }
    else if (cmd.group_coupled)
    {
        std::cout << "\t- using " << inp.nL + 1 << " independently coupled groups (" << L_ + Pi_ + 1 << " states each)" << std::endl;
        for (int igroup = 0; igroup <= inp.nL; igroup++)
            groups_.push_back(linspace(igroup * (L_ + Pi_ + 1), (igroup + 1) * (L_ + Pi_ + 1) - 1, L_ + Pi_ + 1));
    }
    else
    {
        std::cout << "\t- using fully decoupled system" << std::endl;
        for (unsigned iang = 0; iang < states_.size(); iang++)
            groups_.push_back(NumberArray<unsigned>(1, iang));
    }
    
    // precompute angular integrals
    std::cout << "\t- calculating angular integrals ... " << std::flush;
    for (unsigned lambda = 0; lambda <= maxlambda_; lambda++)
    for (unsigned ill = 0; ill < states_.size(); ill++)
    for (unsigned illp = 0; illp < states_.size(); illp++)
    {
        f_.push_back
        (
            special::computef
            (
                lambda,
                states_[ill].first,
                states_[ill].second,
                states_[illp].first,
                states_[illp].second,
                L_
            )
        );
        
        if (not std::isfinite(f_.back()))
        {
            HexException
            (
                "\nFailed to evaluate the angular integral f[%d](%d,%d,%d,%d;%d).",
                lambda,
                states_[ill].first,
                states_[ill].second,
                states_[illp].first,
                states_[illp].second,
                L_
            );
        }
    }
    std::cout << "ok" << std::endl << std::endl;;
}

double AngularBasis::f (unsigned lambda, unsigned ill, unsigned illp) const
{
    return std::max(ill,illp) < states_.size() ? f_[(lambda * states_.size() + ill) * states_.size() + illp] : 0.;
}

unsigned AngularBasis::index (unsigned l1, unsigned l2) const
{
    return std::find(states_.begin(), states_.end(), std::make_pair(l1,l2)) - states_.begin();
}

double AngularBasis::f (unsigned lambda, unsigned l1, unsigned l2, unsigned l1p, unsigned l2p) const
{
    // find the matching angular states
    unsigned m = index(l1,l2);
    unsigned n = index(l1p,l2p);
    
    // call the other getter
    return f(lambda, m, n);
}
