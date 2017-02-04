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

#ifndef HEX_ECS_HLDATA_H
#define HEX_ECS_HLDATA_H

// --------------------------------------------------------------------------------- //

#include "hex-densematrix.h"

// --------------------------------------------------------------------------------- //

/**
 * @brief One-electron diagonalization data.
 * 
 * This data structure is used to manage the one-electron eigenstates and some other data.
 * It is used by @ref NoPreconditioner to construct extended matrix of the system
 * and to calculate the right-hand side, by @ref KPAPreconditioner to precondition the
 * system by solution of the independent electrons and by the @ref Amplitudes class
 * to extract the scattering T-matrix.
 */
class HlData
{
    public:
        
        /// Link the structure to a disk file.
        void hdflink (const char * file);
        
        /// Check that the file exists and can be opened for reading.
        bool hdfcheck (const char * file = nullptr) const;
        
        /// Try to load data from a disk file.
        bool hdfload (const char * file = nullptr);
        
        /// Save data to disk.
        bool hdfsave (const char * file = nullptr) const;
        
        /// Read a pseudo bound state.
        cArray readPseudoState (unsigned l, unsigned ichan) const;
        
        /// Release memory.
        void drop ();
        
        /// One-electron hamiltonian eigenvalues.
        cArray Dl;
        
        /// One-electron hamiltonian eigenvectors.
        ColMatrix<Complex> Cl;
        
        /// Other combinations, used by @ref KPAPreconditioner only.
        RowMatrix<Complex> invCl_invsqrtS, invsqrtS_Cl;
        
        /// Filename.
        std::string filename;
};

// --------------------------------------------------------------------------------- //

#endif
