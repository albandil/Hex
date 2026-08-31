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

#include "hex-special.h"

#include "ang.h"

AngularBasis::AngularBasis (InputFile const & inp, bool fold)
    : L_(inp.L), S_(0), Pi_(inp.Pi), maxell_(inp.maxell), folded_(false)
{
    std::cout << "Setting up the coupled angular states..." << std::endl;

    // for given L, Π and levels list all available (ℓ₁ℓ₂) pairs
    int maxl1 = 0, maxl2 = 0;
    for (int ell = 0; ell <= inp.levels; ell++)
    {
        std::cout << "\t- [" << ell << "] ";

        // get sum of the angular momenta for this angular level
        int sum = 2 * ell + inp.L + inp.Pi;

        // for all angular momentum pairs that do compose L
        for (int l1 = ell; l1 <= sum - ell; l1++)
        {
            int l2 = sum - l1;

            if ((l1 <= l2 or inp.exchange) and std::abs(l1 - l2) <= inp.L and inp.L <= l1 + l2 and (inp.limit < 0 or std::min(l1, l2) <= inp.limit))
            {
                std::cout << "(" << l1 << "," << l2 << ") ";

                maxl1 = std::max(maxl1, l1);
                maxl2 = std::max(maxl2, l2);

                states_full_.push_back(std::make_pair(l1, l2));
            }
        }
        std::cout << std::endl;
    }

    // get maximal angular momentum transfer
    maxlambda_ = 2 * std::min(maxl1, maxl2);

    // Fold the exchange symmetry into the equations, if asked to. The l1 > l2 states are
    // moved past the end of the solved list; they stay in 'states_full_', because the
    // matrix blocks that couple them to the solved states are still needed.
    if (fold)
    {
        if (not inp.exchange)
            HexException("The option --fold-exchange needs the full angular basis; set 'exchange = 1' in the input file.");

        if (inp.Zp != -1)
            HexException("The option --fold-exchange is only valid for two identical particles, but the projectile charge is %g.", inp.Zp);

        // sort the solved states to the front, keeping the relative order of both halves
        std::stable_partition
        (
            states_full_.begin(),
            states_full_.end(),
            [](std::pair<int,int> const & ll) { return ll.first <= ll.second; }
        );

        // the solved states are those that precede the first state with l1 > l2
        std::size_t nsolved = std::find_if
        (
            states_full_.begin(),
            states_full_.end(),
            [](std::pair<int,int> const & ll) { return ll.first > ll.second; }
        ) - states_full_.begin();

        states_.assign(states_full_.begin(), states_full_.begin() + nsolved);
        folded_ = true;
    }
    else
    {
        states_ = states_full_;
    }

    // find the mirror state (l2,l1) of every state (l1,l2)
    for (std::pair<int,int> const & ll : states_full_)
    {
        std::pair<int,int> mirror = std::make_pair(ll.second, ll.first);

        mirror_.push_back
        (
            std::find(states_full_.begin(), states_full_.end(), mirror) - states_full_.begin()
        );

        // this cannot happen for a basis built with 'exchange = 1', which is the only
        // case where the mirror blocks are ever asked for
        if (mirror_.back() == (int)states_full_.size())
            mirror_.back() = -1;
    }

    if (folded_)
    {
        // every solved state needs its mirror, which is where its counterpart is taken from
        for (std::size_t ill = 0; ill < states_full_.size(); ill++)
        if (mirror_[ill] < 0)
        {
            HexException
            (
                "The angular state (%d,%d) has no counterpart (%d,%d) in the basis, so the "
                "exchange symmetry cannot be folded in.",
                states_full_[ill].first, states_full_[ill].second,
                states_full_[ill].second, states_full_[ill].first
            );
        }

        std::cout << "\t- exchange symmetry folded in: " << states_.size() << " of "
                  << states_full_.size() << " blocks are solved for" << std::endl;
    }

    // precompute angular integrals
    std::cout << "\t- calculating angular integrals ... " << std::flush;
    for (int lambda = 0; lambda <= maxlambda_; lambda++)
    for (unsigned ill = 0; ill < states_full_.size(); ill++)
    for (unsigned illp = 0; illp < states_full_.size(); illp++)
    {
        f_.push_back
        (
            special::computef
            (
                lambda,
                states_full_[ill].first,
                states_full_[ill].second,
                states_full_[illp].first,
                states_full_[illp].second,
                L_
            )
        );

        if (not std::isfinite(f_.back()))
        {
            HexException
            (
                "\nFailed to evaluate the angular integral f[%d](%d,%d,%d,%d;%d).",
                lambda,
                states_full_[ill].first,
                states_full_[ill].second,
                states_full_[illp].first,
                states_full_[illp].second,
                L_
            );
        }
    }
    std::cout << "ok" << std::endl;
}

