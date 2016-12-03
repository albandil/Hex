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

#ifndef HEX_DOM_PRECONDITIONER_H
#define HEX_DOM_PRECONDITIONER_H

// --------------------------------------------------------------------------------- //

#include <array>
#include <vector>

// --------------------------------------------------------------------------------- //

#include "NoPreconditioner.h"

// --------------------------------------------------------------------------------- //

/**
 * @brief Multiplicative Schwarz domain decomposition preconditioner.
 * 
 * This preconditioner splits the simulated domain into several overlapping subdomains
 * and solves the smaller per-subdomain problems in a sequence.
 */
class DOMPreconditioner : public NoPreconditioner
{
    public:
        
        // run-time selection mechanism
        preconditionerRunTimeSelectionDefinitions(DOMPreconditioner, "DOM")
        
        // default constructor needed by the RTS mechanism
        DOMPreconditioner () {}
        
        // constructor
        DOMPreconditioner
        (
            Parallel const & par,
            InputFile const & inp,
            AngularBasis const & ll,
            Bspline const & bspline_inner,
            Bspline const & bspline_full,
            CommandLine const & cmd
        );
        
        // description of the preconditioner
        virtual std::string description () const;
        
        // reuse parent definitions
        using NoPreconditioner::rhs;
        using NoPreconditioner::multiply;
        
        // declare own definitions
        virtual void setup ();
        virtual void update (Real E);
        virtual void precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const;
        virtual void finish ();
    
    protected:
        
        // solutions on the sub-domains
        typedef struct
        {
            cArray r;               // original source
            cArray z;               // solution
            std::array<cArray,4> ssrc;  // surrogate sources from neighbour panels
            std::array<cArray,4> outf;  // outgoing field to neighbour panels
        }
        PanelSolution;
        
        // find solution on a sub-domain
        void solvePanel (int n, std::vector<PanelSolution> & p, int i, int j) const;
};

// --------------------------------------------------------------------------------- //

#endif
