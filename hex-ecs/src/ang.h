//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2017, Jakub Benda, Charles University in Prague                    //
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

#ifndef HEX_ANG_H
#define HEX_ANG_H

#include "hex-arrays.h"

#include "inout.h"

/**
 * @brief Angular basis.
 * 
 * List of coupled angular states included in the calculation.
 */
class AngularBasis
{
    public:

        /// Constructor.
        AngularBasis (InputFile const & inp);
        AngularBasis (AngularBasis const & ang)
            : L_(ang.L_), S_(ang.S_), Pi_(ang.Pi_), maxell_(ang.maxell_),
            maxlambda_(ang.maxlambda_), states_(ang.states_), f_(ang.f_) {}

        /// List of coupled angular states.
        std::vector<std::pair<int,int>> const & states () const { return states_; }

        /// Highest orbital number.
        int maxell () const { return maxell_; }

        /// Highest multipole.
        int maxlambda () const { return maxlambda_; }

        /// Angular integrals.
        Real f (int ill, int illp, int lambda) const
        {
            return f_[(lambda * states_.size() + ill) * states_.size() + illp];
        }

        /// Angular integrals.
        rArray const & f () const { return f_; }

        /// Total angular momentum.
        int L () const { return L_; }

        /// Total spin.
        int S () const { return S_; }
        int & S () { return S_; }

        /// Total parity.
        int Pi () const { return Pi_; }

    private:

        // Quantum numbers.
        int L_, S_, Pi_, nL_;

        // Highest orbital number.
        int maxell_;

        // Highest multipole.
        int maxlambda_;

        // List of coupled angular states.
        std::vector<std::pair<int,int>> states_;

        // Angular integrals.
        rArray f_;
};

#endif
