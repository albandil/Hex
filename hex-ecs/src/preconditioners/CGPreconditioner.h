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

#ifndef HEX_CGPRECONDITIONER_H
#define HEX_CGPRECONDITIONER_H

// --------------------------------------------------------------------------------- //

#include "NoPreconditioner.h"

// --------------------------------------------------------------------------------- //

/**
 * @brief CG iteration-based preconditioner.
 * 
 * This class adds some preconditioning capabilities to its base class
 * NoPreconditioner. The preconditioning is done by diagonal block solution
 * using the conjugate gradients solver (which itself is non-preconditioned).
 */
class CGPreconditioner : public NoPreconditioner
{
    public:
        
        // run-time selection mechanism
        preconditionerRunTimeSelectionDefinitions(CGPreconditioner, "cg")
        
        // default constructor needed by the RTS mechanism
        CGPreconditioner () {}
        
        // constructor
        CGPreconditioner
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
            ),
            n_(ang.states().size(), -1) {}
        
        // description of the preconditioner
        virtual std::string description () const;
        
        // reuse parent definitions
        using NoPreconditioner::setup;
        using NoPreconditioner::update;
        using NoPreconditioner::rhs;
        using NoPreconditioner::multiply;
        
        // declare own definitions
        virtual void precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const;
        virtual void finish ();
        
        // solve diagonal block
        int solve_block (int ill, const cArrayView r, cArrayView z) const;
        
        // inner CG driver
        virtual void CG_init (int iblock) const;
        virtual void CG_mmul (int iblock, const cArrayView p, cArrayView q) const;
        virtual void CG_prec (int iblock, const cArrayView r, cArrayView z) const;
        virtual void CG_exit (int iblock) const;
        
        // inner CG callback routines
        virtual Real CG_compute_norm (const cArrayView a) const;
        virtual Complex CG_scalar_product (const cArrayView a, const cArrayView b) const;
        virtual void CG_axby_operation (Complex a, cArrayView x, Complex b, const cArrayView y) const;
        virtual void CG_constrain (cArrayView r) const;
    
    protected:
        
        // last iterations
        mutable iArray n_;
        
        // timing
        mutable std::size_t us_axby_, us_mmul_, us_norm_, us_prec_, us_spro_;
};

// --------------------------------------------------------------------------------- //

#endif
