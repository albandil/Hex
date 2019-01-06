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

#ifndef HEX_HYBPRECONDITIONER_H
#define HEX_HYBPRECONDITIONER_H

// --------------------------------------------------------------------------------- //

#include <set>
#include <string>
#include <vector>

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"
#include "hex-matrix.h"

// --------------------------------------------------------------------------------- //

#include "ILUPreconditioner.h"
#include "KPAPreconditioner.h"

// --------------------------------------------------------------------------------- //

/**
 * @brief Hybrid preconditioner.
 * 
 * Combination of ILU and KPA:
 * - KPA is used for angular blocks with no asymptotic channels.
 * - ILU is used for angular blocks with asymptotic channels.
 */
class HybCGPreconditioner : public ILUCGPreconditioner, public KPACGPreconditioner
{
    public:

        // run-time selection mechanism
        preconditionerRunTimeSelectionDefinitions(HybCGPreconditioner, "HYB")

        // default constructor needed by the RTS mechanism
        HybCGPreconditioner () {}

        // constructor
        HybCGPreconditioner
        (
            CommandLine  const & cmd,
            InputFile    const & inp,
            Parallel     const & par,
            AngularBasis const & ang,
            Bspline const & bspline_x_inner,
            Bspline const & bspline_x_full,
            Bspline const & bspline_y_inner,
            Bspline const & bspline_y_full
        ) : CGPreconditioner
            (
                cmd, inp, par, ang,
                bspline_x_inner, bspline_x_full,
                bspline_y_inner, bspline_y_full
            ),
            ILUCGPreconditioner
            (
                cmd, inp, par, ang,
                bspline_x_inner, bspline_x_full,
                bspline_y_inner, bspline_y_full
            ),
            KPACGPreconditioner
            (
                cmd, inp, par, ang,
                bspline_x_inner, bspline_x_full,
                bspline_y_inner, bspline_y_full
            )
        {
            // nothing more to do
        }

        // preconditioner description
        virtual std::string description () const;

        // reuse parent definitions
        using CGPreconditioner::multiply;
        using CGPreconditioner::rhs;
        using CGPreconditioner::precondition;

        // declare own definitions
        virtual void setup ();
        virtual void update (Real E);
        virtual void finish ();

        // inner CG callback (needed by parent)
        virtual void CG_init (int iblock) const;
        virtual void CG_prec (int iblock, const cArrayView r, cArrayView z) const;
        virtual void CG_mmul (int iblock, const cArrayView r, cArrayView z) const;
        virtual void CG_exit (int iblock) const;

        // decide whether to use the ILU preconditioner
        bool ilu_needed (int iblock) const;
};

// --------------------------------------------------------------------------------- //

#endif
