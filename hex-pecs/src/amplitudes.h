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

#ifndef HEX_PECS_AMPLITUDES_H
#define HEX_PECS_AMPLITUDES_H

#include <vector>

#include "ang.h"
#include "io.h"
#include "radial.h"

class Amplitudes
{
    public:
        
        /// Constructor - allocates the arrays.
        Amplitudes
        (
            CommandLine const & cmd,
            InputFile const & inp,
            AngularBasis const & ang,
            RadialBasis const & rad
        );
        
        /// Extract amplitudes from a given wave function.
        void extract
        (
            unsigned istate,
            unsigned iEpert,
            Complex const * psi
        );
        
        /// Write all extracted data to the data files.
        void write ();
    
    private:
        
        /// Command line parameters.
        CommandLine const & cmd_;
        
        /// Input file contents.
        InputFile const & inp_;
        
        /// Angular basis states.
        AngularBasis const & ang_;
        
        /// Radial grid data.
        RadialBasis const & rad_;
        
        /// Partial T-matrices for every transition, shape: Nenergy × istates × fstates × (maxell + 1).
        Array<Array<Array<cArray>>> T_matrices;
};

#endif
