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

#if (!defined(HEX_GPUPRECONDITIONER_H) && defined(WITH_OPENCL))
#define HEX_GPUPRECONDITIONER_H

// --------------------------------------------------------------------------------- //

#include <string>

// --------------------------------------------------------------------------------- //

#include <CL/cl.h>

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"
#include "hex-matrix.h"

// --------------------------------------------------------------------------------- //

#include "KPAPreconditioner.h"

// --------------------------------------------------------------------------------- //

#include "clarrays.h"

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
class GPUCGPreconditioner : public virtual KPACGPreconditioner
{
    public:
        
        // run-time selection mechanism
        preconditionerRunTimeSelectionDefinitions(GPUCGPreconditioner, "GPU")
        
        // default constructor needed by the RTS mechanism
        GPUCGPreconditioner ();
        
        // constructor
        GPUCGPreconditioner
        (
            CommandLine  const & cmd,
            InputFile    const & inp,
            Parallel     const & par,
            AngularBasis const & ang,
            Bspline const & bspline_inner,
            Bspline const & bspline_full,
            Bspline const & bspline_panel_x,
            Bspline const & bspline_panel_y
        );
        
        // preconditioner description
        virtual std::string description () const;
        
        // reuse parent definitions
        using KPACGPreconditioner::rhs;
        using KPACGPreconditioner::update;
        
        // declare own definitions
        virtual void setup ();
        virtual void finish ();
        virtual void multiply (BlockArray<Complex> const & p, BlockArray<Complex> & q, MatrixSelection::Selection tri) const;
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
        cl_kernel mmls_, mabt_, matbt_;
        cl_kernel mms1_, mms2d_, mms2c_;
        cl_kernel mml1_, mml2d_, mml2c_;
        cl_kernel axby_, norm_, spro_, krdv_;
        
        // device data connections
        mutable clArray<Complex> tmp_, tmA_;
        mutable clArray<Real> nrm_;
        clArrayView<Complex> t_inner_a_, t_inner_p_;
        clArrayView<Complex> S_inner_a_, D_inner_a_, Mm1_inner_a_, Mm2_inner_a_;
        clArrayView<Complex> S_inner_p_, D_inner_p_, Mm1_inner_p_, Mm2_inner_p_;
        std::vector<clArrayView<Complex>> M_L_inner_a_, M_mLm1_inner_a_;
        std::vector<clArrayView<Complex>> M_L_inner_p_, M_mLm1_inner_p_;
        std::vector<clArrayView<LU_int_t>> R_coupled_p_, R_coupled_i_;
        std::vector<clArrayView<Complex>> R_coupled_x_;
};

// --------------------------------------------------------------------------------- //

#endif
