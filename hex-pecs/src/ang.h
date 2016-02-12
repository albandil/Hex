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

#ifndef HEX_PECS_ANG_H
#define HEX_PECS_ANG_H

#include <iostream>
#include <vector>

#include "hex-arrays.h"
#include "hex-special.h"

#include "io.h"

class AngularBasis
{
    public:
        
        AngularBasis (CommandLine const & cmd, InputFile const & inp);
        
        /// Get particular state.
        std::pair<unsigned,unsigned> const & state (unsigned ill) const { return states_[ill]; }
        
        /// Angular integrals.
        double f (unsigned lambda, unsigned ill, unsigned illp) const;
        double f (unsigned lambda, unsigned l1, unsigned l2, unsigned l1p, unsigned l2p) const;
        
        /// Get index of a specific angular momentum state.
        unsigned index (unsigned l1, unsigned l2) const;
        
        /// Highest multipole.
        unsigned maxlambda () const { return maxlambda_; }
        
        /// Larges angular momentum.
        unsigned maxell () const { return maxell_; }
        
        /// Number of angular momentum states.
        unsigned size () const { return states_.size(); }
        
        /// Number of coupled angular groups.
        unsigned Ngroups () const { return groups_.size(); }
        
        /// Number of angular momentum states in a group.
        unsigned groupsize () const { return states_.size() / groups_.size(); }
        
    private:
        
        // Quantum numbers.
        unsigned L_, S_, Pi_, nL_;
        
        // Highest multipole.
        unsigned maxlambda_;
        
        // Largest angular momentum.
        unsigned maxell_;
        
        // List of angular states.
        std::vector<std::pair<unsigned,unsigned>> states_;
        
        // Coupled groups.
        std::vector<NumberArray<unsigned>> groups_;
        
        // Angular integrals.
        rArray f_;
        
};

#endif
