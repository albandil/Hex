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

#ifndef HEX_ILUPRECONDITIONER_H
#define HEX_ILUPRECONDITIONER_H

#include "preconditioners.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef WITH_SUPERLU_DIST
#include <superlu_zdefs.h>
#endif

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
        
        static const std::string prec_name;
        static const std::string prec_description;
        
        virtual std::string const & name () const { return prec_name; }
        virtual std::string const & description () const { return prec_description; }
        
        ILUCGPreconditioner
        (
            Parallel const & par,
            InputFile const & inp,
            AngularBasis const & ll,
            Bspline const & bspline_atom,
            Bspline const & bspline_proj,
            Bspline const & bspline_proj_full,
            CommandLine const & cmd
        ) : CGPreconditioner(par, inp, ll, bspline_atom, bspline_proj, bspline_proj_full, cmd),
            csr_blocks_(ll.states().size()), lu_(ll.states().size())
        {
#ifdef _OPENMP
            omp_init_lock(&lu_lock_);
#endif
        }
        
        ~ILUCGPreconditioner ()
        {
#ifdef _OPENMP
            omp_destroy_lock(&lu_lock_);
#endif
        }
        
        // reuse parent definitions
        virtual void multiply (BlockArray<Complex> const & p, BlockArray<Complex> & q) const { CGPreconditioner::multiply(p,q); }
        virtual void rhs (BlockArray<Complex> & chi, int ienergy, int instate) const { CGPreconditioner::rhs(chi,ienergy,instate); }
        virtual void precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const { CGPreconditioner::precondition(r,z); }
        
        // declare own definitions
        virtual void setup ();
        virtual void update (double E);
        virtual void finish ();
        
        // inner CG callback (needed by parent)
        virtual void CG_init (int iblock) const;
        virtual void CG_prec (int iblock, const cArrayView r, cArrayView z) const;
        virtual void CG_exit (int iblock) const;
        
    protected:
        
        // diagonal CSR block for every coupled state
        mutable std::vector<CsrMatrix<LU_int_t,Complex>> csr_blocks_;
        
        // LU decompositions of the CSR blocks
        mutable std::vector<std::shared_ptr<LUft<LU_int_t,Complex>>> lu_;
        
#ifdef _OPENMP
        // factorization lock
        mutable omp_lock_t lu_lock_;
#endif
        
#ifdef WITH_SUPERLU_DIST
        // process grid
        gridinfo_t grid_;
#endif
};

#endif
