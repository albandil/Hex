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

#ifndef HEX_GPUPRECONDITIONER_H
#define HEX_GPUPRECONDITIONER_H

#include <string>

#include "hex-arrays.h"
#include "hex-matrix.h"

#include "preconditioners.h"

#ifdef WITH_OPENCL

#include "clarrays.h"

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
class GPUCGPreconditioner : public virtual KPACGPreconditioner
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
            
            /// Try to load data from a disk file.
            bool hdfload (const char * file = nullptr);
            
            /// Save data to disk.
            bool hdfsave (const char * file = nullptr) const;
            
            /// One-electron hamiltonian diagonalization.
            RowMatrix<Complex> invCl_invsqrtS, invsqrtS_Cl;
            
            /// Diagonal parts of one-electron hamiltonians.
            cArray Dl;
            
            /// Filename.
            std::string filename;
        } Data;
        
        GPUCGPreconditioner
        (
            Parallel const & par,
            InputFile const & inp,
            AngularBasis const & ll,
            Bspline const & bspline_inner,
            Bspline const & bspline_outer,
            Bspline const & bspline_full,
            CommandLine const & cmd
        ) : CGPreconditioner(par, inp, ll, bspline_inner, bspline_outer, bspline_full, cmd),
            KPACGPreconditioner(par, inp, ll, bspline_inner, bspline_outer, bspline_full, cmd)
        {
            // nothing more to do
        }
        
        // reuse parent definitions
        virtual void rhs (BlockArray<Complex> & chi, int ienergy, int instate) const { KPACGPreconditioner::rhs(chi,ienergy,instate); }
        virtual void update (Real E) { KPACGPreconditioner::update(E); }
        
        // declare own definitions
        virtual void setup ();
        virtual void finish ();
        virtual void multiply (BlockArray<Complex> const & p, BlockArray<Complex> & q) const;
        virtual void precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const;
        
    protected:
        
        // OpenCL environment
        cl_platform_id platform_;
        cl_device_id device_;
        cl_context context_;
        cl_command_queue queue_;
        cl_program program_;
        
        // size of a workgroup
        std::size_t Nlocal_;
        std::size_t blocksize_;
        
        // memory allocation flags
        cl_mem_flags smallDataFlags_;
        cl_mem_flags largeDataFlags_;
        
        // computational kernels
        cl_kernel mabt_;
        cl_kernel mml1_;
        cl_kernel mml2_;
        cl_kernel mml2_dcpl_, mml2_cpld_;
        cl_kernel mml2_dcpl_offset_, mml2_cpld_offset_;
        cl_kernel axby_;
        cl_kernel norm_;
        cl_kernel spro_;
        cl_kernel krdv_;
        
        // device data connections
        mutable clArray<Complex> tmp_, tmA_;
        mutable clArray<Real> nrm_;
        clArrayView<Complex> t_atom_, t_proj_;
        clArrayView<Complex> S_atom_p_, D_atom_p_, Mm1_tr_atom_p_, Mm2_atom_p_;
        clArrayView<Complex> S_proj_p_, D_proj_p_, Mm1_tr_proj_p_, Mm2_proj_p_;
        std::vector<clArrayView<Complex>> Mi_L_atom_, Mi_mLm1_atom_, M_L_atom_, M_mLm1_atom_;
        std::vector<clArrayView<Complex>> Mi_L_proj_, Mi_mLm1_proj_, M_L_proj_, M_mLm1_proj_;
        std::vector<clArrayView<Complex>> Rdia_;
        
        cl_short nsrcseg_, ndstseg_;
};

#endif

#endif
