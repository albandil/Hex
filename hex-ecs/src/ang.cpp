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

#include "hex-special.h"

#include "ang.h"

AngularBasis::AngularBasis (InputFile const & inp)
    : L_(inp.L), S_(0), Pi_(inp.Pi), maxell_(inp.maxell), maxlambda_(inp.L + 2 * inp.levels)
{
    std::cout << "Setting up the coupled angular states..." << std::endl;
    
    // for given L, Π and levels list all available (ℓ₁ℓ₂) pairs
    for (int ell = 0; ell <= inp.levels; ell++)
    {
        std::cout << "\t- [" << ell << "] ";
        
        // get sum of the angular momenta for this angular level
        int sum = 2 * ell + inp.L + inp.Pi;
        
        // for all angular momentum pairs that do compose L
        for (int l1 = ell; l1 <= sum - ell; l1++)
        {
            int l2 = sum - l1;
            if (std::abs(l1 - l2) <= inp.L and inp.L <= l1 + l2 and (inp.limit < 0 or std::min(l1, l2) <= inp.limit))
            {
                std::cout << "(" << l1 << "," << l2 << ") ";
                states_.push_back(std::make_pair(l1, l2));
            }
        }
        std::cout << std::endl;
    }
    
    // precompute angular integrals
    std::cout << "\t- calculating angular integrals ... " << std::flush;
    for (int lambda = 0; lambda <= maxlambda_; lambda++)
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
    std::cout << "ok" << std::endl;
}
