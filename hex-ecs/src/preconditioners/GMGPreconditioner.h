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

#ifndef HEX_GMGPRECONDITIONER_H
#define HEX_GMGPRECONDITIONER_H

// --------------------------------------------------------------------------------- //

#include "NoPreconditioner.h"

// --------------------------------------------------------------------------------- //

class GMGPreconditioner : public NoPreconditioner
{
    public:
        
        // run-time selection mechanism
        preconditionerRunTimeSelectionDefinitions(GMGPreconditioner, "GMG")
        
        // constructor
        GMGPreconditioner
        (
            Parallel const & par,
            InputFile const & inp,
            AngularBasis const & ll,
            Bspline const & bspline_inner,
            Bspline const & bspline_full,
            CommandLine const & cmd
        );
        
        // sub-grid contructor
        GMGPreconditioner
        (
            Parallel const & par,
            InputFile const & inp,
            AngularBasis const & ll,
            Bspline const & bspline_inner,
            Bspline const & bspline_full,
            CommandLine const & cmd,
            int level
        );
        
        // destructor
        virtual ~GMGPreconditioner ();
        
        // preconditioner description
        virtual std::string description () const;
        
        // reuse parent definitions
        using NoPreconditioner::rhs;
        using NoPreconditioner::multiply;
        using NoPreconditioner::finish;
        
        // declare own definitions
        virtual void setup ();
        virtual void update (Real E);
        virtual void precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const;
    
    private:
        
        int level_;
        
        Bspline const & bspline_inner_fine_;
        Bspline const & bspline_full_fine_;
        
        Bspline bspline_inner_coarse_;
        Bspline bspline_full_coarse_;
        
        RowMatrix<Complex> restrictor_inner_, restrictor_outer_;
        RowMatrix<Complex> prolongator_inner_, prolongator_outer_;
        
        cArrays D;
        
        PreconditionerBase * subgrid_;
};

// --------------------------------------------------------------------------------- //

#endif
