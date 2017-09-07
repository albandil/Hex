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

#ifndef HEX_ILUPRECONDITIONER_H
#define HEX_ILUPRECONDITIONER_H

// --------------------------------------------------------------------------------- //

#include "luft.h"

// --------------------------------------------------------------------------------- //

#include "CGPreconditioner.h"

// --------------------------------------------------------------------------------- //

#ifdef _OPENMP
    #include <omp.h>
#endif

// --------------------------------------------------------------------------------- //

#ifdef WITH_SUPERLU_DIST
    #include <superlu_zdefs.h>
#endif

// --------------------------------------------------------------------------------- //

/**
 * @brief ILU-preconditioned CG-based preconditioner.
 * 
 * Enhances CGPreconditioner conjugate gradients solver by incomplete LU factorization
 * preconditioning. This is done by redefining virtual function CG_prec. The factorization
 * is drop tolerance based and is computed by UMFPACK.
 */
class ILUCGPreconditioner : public virtual CGPreconditioner
{
    public:
        
        // run-time selection mechanism
        preconditionerRunTimeSelectionDefinitions(ILUCGPreconditioner, "ILU")
        
        // default constructor needed by the RTS mechanism
        ILUCGPreconditioner () {}
        
        // constructor
        ILUCGPreconditioner
        (
            CommandLine  const & cmd,
            InputFile    const & inp,
            Parallel     const & par,
            AngularBasis const & ang,
            Bspline const & bspline_x_inner,
            Bspline const & bspline_x_full,
            Bspline const & bspline_y_inner,
            Bspline const & bspline_y_full
        );
        
        // destructor
        ~ILUCGPreconditioner ();
        
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
        virtual void CG_exit (int iblock) const;
        
    protected:
        
        // LU decompositions of the diagonal blocks
        mutable std::vector<std::shared_ptr<LUft>> lu_;
        
        // prepare data structures for LU factorizations
        void reset_lu ();
        
#ifdef _OPENMP
        // factorization lock
        mutable omp_lock_t lu_lock_;
#endif
        
#ifdef WITH_SUPERLU_DIST
        // process grid
        gridinfo_t grid_;
#endif
};

// --------------------------------------------------------------------------------- //

#endif
