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

#ifndef HEX_ECS_COUPLED_PRECONDITIONER_H
#define HEX_ECS_COUPLED_PRECONDITIONER_H

// --------------------------------------------------------------------------------- //

#include "KPAPreconditioner.h"

// --------------------------------------------------------------------------------- //

class CoupledPreconditioner : public virtual KPACGPreconditioner
{
    public:

        // run-time selection mechanism
        preconditionerRunTimeSelectionDefinitions(CoupledPreconditioner, "coupled")

        // default constructor needed by the RTS mechanism
        CoupledPreconditioner () {}

        // constructor
        CoupledPreconditioner
        (
            CommandLine  const & cmd,
            InputFile    const & inp,
            Parallel     const & par,
            AngularBasis const & ang,
            Bspline const & bspline_inner,
            Bspline const & bspline_full,
            Bspline const & bspline_panel_x,
            Bspline const & bspline_panel_y
        ) : CGPreconditioner
            (
                cmd, inp, par, ang,
                bspline_inner,
                bspline_full,
                bspline_panel_x,
                bspline_panel_y
            ),
            KPACGPreconditioner
            (
                cmd, inp, par, ang,
                bspline_inner,
                bspline_full,
                bspline_panel_x,
                bspline_panel_y
            )
        {
            // nothing to do
        }

        // preconditioner description
        virtual std::string description () const;

        // reuse parent definitions
        using KPACGPreconditioner::setup;
        using KPACGPreconditioner::rhs;
        using KPACGPreconditioner::multiply;

        // override segregated preconditioner routine
        virtual int solve_block (int ill, const cArrayView r, cArrayView z) const;

        // declare own definitions
        virtual void update (Real E);
        virtual void precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const;
        virtual void finish ();

    protected:

        // LU factorization.
        std::shared_ptr<LUft> lu_;

        // Workspace used for the solution.
        mutable cArray X, Y;
};

// --------------------------------------------------------------------------------- //

#endif // HEX_ECS_COUPLED_PRECONDITIONER_H
