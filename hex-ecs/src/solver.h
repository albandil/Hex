//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2018, Jakub Benda, Charles University in Prague                    //
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

#ifndef HEX_ECS_SOLVER_H
#define HEX_ECS_SOLVER_H

// --------------------------------------------------------------------------------- //

#include <vector>

// --------------------------------------------------------------------------------- //

#include "hex-itersolve.h"

// --------------------------------------------------------------------------------- //

#include "ang.h"
#include "bspline.h"
#include "inout.h"
#include "parallel.h"
#include "preconditioners.h"

// --------------------------------------------------------------------------------- //

class Solver
{
    public:

        /// Constructor.
        Solver
        (
            CommandLine        & cmd,
            InputFile    const & inp,
            Parallel     const & par,
            AngularBasis const & ang,
            Bspline const & bspline_inner,
            Bspline const & bspline_full
        );

        /// Pick preconditioner for a given panel.
        void choose_preconditioner ();

        /// Precompute (or load) the radial data.
        void setup_preconditioner ();

        /// Find solution of the scattering equations on chosen panel.
        void solve ();

        /// Release resources.
        void finish ();

    protected:

        // CG preconditioner callback
        void apply_preconditioner_ (BlockArray<Complex> const & r, BlockArray<Complex> & z) const;

        // CG matrix multiplication callback
        void matrix_multiply_ (BlockArray<Complex> const & p, BlockArray<Complex> & q) const;

        // CG scalar product function callback
        Complex scalar_product_ (BlockArray<Complex> const & x, BlockArray<Complex> const & y) const;

        // CG norm function that broadcasts master's result to all nodes
        Real compute_norm_ (BlockArray<Complex> const & r) const;

        // CG linear combination
        void axby_operation_ (Complex a, BlockArray<Complex> & x, Complex b, BlockArray<Complex> const & y) const;

        // CG new array
        BlockArray<Complex> new_array_ (std::size_t N, std::string name) const;

        // CG optional write current solution
        void process_solution_ (unsigned iteration, BlockArray<Complex> const & x) const;

        // concatenate previous-panel full solution and new single-panel solution
        void concatenate_panels_ (cArray & psi, cArray const & psip) const;

        // save array to disk
        void checkpoint_array_ (BlockArray<Complex> const & psi) const;

        // read array from disk
        void recover_array_ (BlockArray<Complex> & psi) const;

        // constrain the residual
        void constrain_ (BlockArray<Complex> & r) const;

        // filter data in the Amplitudes::TmatArray to contain only transition for the current initial state, spin, energy
        cArray filter_Tmat_data_ (Amplitudes::TmatArray const & T) const;

        // compare two T-matrix arrays
        bool check_convergence_ (cArray const & T1, cArray const & T2) const;

    private:

        /// Command line parameters.
        CommandLine & cmd_;

        /// Input file settings.
        InputFile const & inp_;

        /// Parallel environment.
        Parallel const & par_;

        /// Angular basis.
        AngularBasis ang_;

        /// Radial bases.
        Bspline const & bspline_inner_;
        Bspline const & bspline_full_;

        /// Solver preconditioner.
        PreconditionerBase * prec_;

        /// Linear solver (conjugate gradients).
        ConjugateGradients <Complex, cBlockArray, cBlockArray&> CG_;

        /// States currently being solved.
        iArray instates_;

        /// Runtime information.
        Real E_;
        int iE_;
        mutable Real progress_;
        mutable Real bnorm_;
        mutable int autostop_;

        /// Asymptotic bound channels for every angular momentum state (l₁,l₂) and r₁- or r₂-asymptotics.
        std::vector<std::pair<iArray,iArray>> bstates_;

        /// Number of channels for individual angular states.
        std::vector<std::pair<int,int>> channels_;
};

// --------------------------------------------------------------------------------- //

#endif // HEX_ECS_SOLVER_H
