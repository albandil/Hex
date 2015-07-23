//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2015, Jakub Benda, Charles University in Prague                    //
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

#include <set>
#include <string>
#include <vector>

#include "../arrays.h"
#include "../matrix.h"
#include "../preconditioners.h"

#ifndef NO_LAPACK

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
class KPACGPreconditioner : public CGPreconditioner
{
    public:
        
        static const std::string prec_name;
        static const std::string prec_description;
        
        virtual std::string const & name () const { return prec_name; }
        virtual std::string const & description () const { return prec_description; }
        
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
        
        KPACGPreconditioner
        (
            Parallel const & par,
            InputFile const & inp,
            std::vector<std::pair<int,int>> const & ll,
            Bspline const & bspline_atom,
            Bspline const & bspline_proj,
            CommandLine const & cmd
        ) : CGPreconditioner(par, inp, ll, bspline_atom, bspline_proj, cmd),
            prec_atom_(inp.maxell+1), prec_proj_(inp.maxell+1)
        {
            // nothing more to do
        }
        
        // reuse parent definitions
        virtual void multiply (BlockArray<Complex> const & p, BlockArray<Complex> & q) const { CGPreconditioner::multiply(p,q); }
        virtual void rhs (BlockArray<Complex> & chi, int ienergy, int instate, int Spin) const { CGPreconditioner::rhs(chi,ienergy,instate,Spin); }
        virtual void precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const { CGPreconditioner::precondition(r,z); }
        virtual void update (double E) { CGPreconditioner::update(E); }
        
        // declare own definitions
        virtual void setup ();
        
        // inner CG callback (needed by parent)
        virtual void CG_init (int iblock) const;
        virtual void CG_prec (int iblock, const cArrayView r, cArrayView z) const;
        virtual void CG_mmul (int iblock, const cArrayView r, cArrayView z) const;
        virtual void CG_exit (int iblock) const;
        
    protected:
        
        // internal setup routine
        void prepare
        (
            std::vector<Data> & prec,
            int Nspline,
            SymBandMatrix<Complex> const & mS,
            SymBandMatrix<Complex> const & mD,
            SymBandMatrix<Complex> const & mMm1_tr,
            SymBandMatrix<Complex> const & mMm2,
            Array<bool> done,
            std::set<int> comp_l,
            std::set<int> needed_l
        );
        
        // preconditioner data
        mutable std::vector<Data> prec_atom_;
        mutable std::vector<Data> prec_proj_;
};

#endif

#endif
