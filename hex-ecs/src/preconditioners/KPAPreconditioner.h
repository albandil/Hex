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

#ifndef HEX_KPAPRECONDITIONER_H
#define HEX_KPAPRECONDITIONER_H

// --------------------------------------------------------------------------------- //

#include <deque>
#include <set>
#include <string>
#include <vector>

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"
#include "hex-matrix.h"

// --------------------------------------------------------------------------------- //

#include "CGPreconditioner.h"

// --------------------------------------------------------------------------------- //

/**
 * @brief KPA-preconditioned CG-preconditioner.
 * 
 * This nested preconditioner simplifies the Hamiltonian matrix by omitting electron-electron
 * interaction. The block of the matrix can then be written as a sum of simple Kronecker products
 * @f[
 *     \mathsf{A} = E\mathsf{S}\otimes\mathsf{S}
 *      - \mathsf{H}_1\otimes\mathsf{S}
 *      - \mathsf{S}\otimes\mathsf{H}_2 \,,
 * @f]
 * which can be easily diagonalized (only diagonalization of the small matrices are needed)
 * and so also inverted, which we need for solution of the equations.
 */
class KPACGPreconditioner : public virtual CGPreconditioner
{
    public:
        
        // run-time selection mechanism
        preconditionerRunTimeSelectionDefinitions(KPACGPreconditioner, "KPA")
        
        typedef struct sData
        {
            /// Link the structure to a disk file.
            void hdflink (const char * file);
            
            /// Check that the file exists and can be opened for reading.
            bool hdfcheck (const char * file = nullptr) const;
            
            /// Try to load data from a disk file.
            bool hdfload (const char * file = nullptr);
            
            /// Save data to disk.
            bool hdfsave (const char * file = nullptr) const;
            
            /// Release memory.
            void drop ();
            
            /// One-electron hamiltonian diagonalization.
            RowMatrix<Complex> invCl_invsqrtS, invsqrtS_Cl;
            
            /// Diagonal parts of one-electron hamiltonians.
            cArray Dl;
            
            /// Filename.
            std::string filename;
        } Data;
        
        // default constructor needed by the RTS mechanism
        KPACGPreconditioner () {}
        
        // constructor
        KPACGPreconditioner
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
            prec_atom_(inp.maxell + 1),
            prec_proj_(inp.maxell + 1),
            refcount_atom_(inp.maxell + 1, 0),
            refcount_proj_(inp.maxell + 1, 0)
        {
#ifdef _OPENMP
            omp_init_lock(&lck_);
#endif
        }
        
        virtual ~KPACGPreconditioner ()
        {
#ifdef _OPENMP
            omp_destroy_lock(&lck_);
#endif
        }
        
        // preconditioner description
        virtual std::string description () const;
        
        // reuse parent definitions
        using CGPreconditioner::update;
        using CGPreconditioner::rhs;
        using CGPreconditioner::multiply;
        using CGPreconditioner::precondition;
        
        // declare own definitions
        virtual void setup ();
        virtual void finish ();
        
        // inner CG callback (needed by parent)
        virtual void CG_init (int iblock) const;
        virtual void CG_mmul (int iblock, const cArrayView r, cArrayView z) const;
        virtual void CG_prec (int iblock, const cArrayView r, cArrayView z) const;
        virtual void CG_exit (int iblock) const;
        virtual void CG_constrain (cArrayView r) const;
        
    protected:
        
        // internal setup routine
        void prepare
        (
            std::vector<Data> & prec,
            std::size_t Nspline,
            SymBandMatrix<Complex> const & mS,
            SymBandMatrix<Complex> const & mD,
            SymBandMatrix<Complex> const & mMm1_tr,
            SymBandMatrix<Complex> const & mMm2,
            Array<bool> done,
            std::set<int> comp_l,
            std::set<int> needed_l
        );
        
        // internal concurrent access lock
        void lock_kpa_access () const;
        void unlock_kpa_access () const;
        
        // preconditioner data
        mutable std::vector<Data> prec_atom_, prec_proj_;
        mutable std::vector<unsigned> refcount_atom_, refcount_proj_;
        
        // memory access lock
#ifdef _OPENMP
        mutable omp_lock_t lck_;
#endif
        
        // workspace
        mutable cArrays workspace_;
};

// --------------------------------------------------------------------------------- //

#endif
