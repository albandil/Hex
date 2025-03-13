//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2025, Jakub Benda, Charles University in Prague                    //
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

#ifndef HEX_ECS_DIAG_PRECONDITIONER_H
#define HEX_ECS_DIAG_PRECONDITIONER_H

// --------------------------------------------------------------------------------- //

#include "NoPreconditioner.h"

// --------------------------------------------------------------------------------- //

/**
 * @brief Plain Jacobi preconditioner.
 *
 * Extends NoPreconditioner with diagonal (Jacobi) preconditioning.
 */
class DiagPreconditioner : public NoPreconditioner
{
    public:

        // run-time selection mechanism
        preconditionerRunTimeSelectionDefinitions(DiagPreconditioner, "diag")

        // default constructor needed by the RTS mechanism
        DiagPreconditioner () {}

        // constructor
        DiagPreconditioner
        (
            CommandLine  const & cmd,
            InputFile    const & inp,
            Parallel     const & par,
            AngularBasis const & ang,
            Bspline const & bspline_inner,
            Bspline const & bspline_full,
            Bspline const & bspline_panel_x,
            Bspline const & bspline_panel_y
        ) : NoPreconditioner
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

        // description of the preconditioner
        std::string description () const override;

        // reuse parent definitions
        using NoPreconditioner::rhs;
        using NoPreconditioner::multiply;
        using NoPreconditioner::update;
        using NoPreconditioner::finish;

        // declare own definitions
        void setup () override;
        void precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const override;

    protected:

        mutable cBlockArray diagonal;
};

// --------------------------------------------------------------------------------- //

#endif // HEX_ECS_DIAG_PRECONDITIONER_H
