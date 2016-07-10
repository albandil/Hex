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

#include <set>
#include <string>
#include <vector>

#include "hex-arrays.h"
#include "hex-matrix.h"

#include "preconditioners.h"

#ifndef NO_LAPACK

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
        
        static const std::string prec_name;
        static const std::string prec_description;
        
        virtual std::string const & name () const { return prec_name; }
        virtual std::string const & description () const { return prec_description; }
        
        HybCGPreconditioner
        (
            Parallel const & par,
            InputFile const & inp,
            AngularBasis const & ll,
            Bspline const & bspline_inner,
            Bspline const & bspline_full,
            CommandLine const & cmd
        ) : CGPreconditioner(par, inp, ll, bspline_inner, bspline_full, cmd),
            ILUCGPreconditioner(par, inp, ll, bspline_inner, bspline_full, cmd),
            KPACGPreconditioner(par, inp, ll, bspline_inner, bspline_full, cmd)
        {
            // nothing more to do
        }
        
        // reuse parent definitions
        virtual void multiply (BlockArray<Complex> const & p, BlockArray<Complex> & q) const { CGPreconditioner::multiply(p,q); }
        virtual void rhs (BlockArray<Complex> & chi, int ienergy, int instate) const { CGPreconditioner::rhs(chi,ienergy,instate); }
        virtual void precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const { CGPreconditioner::precondition(r,z); }
        
        // declare own definitions
        virtual void setup ();
        virtual void update (Real E);
        virtual void finish ();
        
        // inner CG callback (needed by parent)
        virtual void CG_init (int iblock) const;
        virtual void CG_prec (int iblock, const cArrayView r, cArrayView z) const;
        virtual void CG_mmul (int iblock, const cArrayView r, cArrayView z) const;
        virtual void CG_exit (int iblock) const;
        
    protected:
        
        // decide whether to use the ILU preconditioner
        bool ilu_needed (int iblock) const;
};

#endif

#endif
