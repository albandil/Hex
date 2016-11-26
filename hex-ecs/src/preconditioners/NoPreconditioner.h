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

#ifndef HEX_ECS_NOPRECONDITIONER_H
#define HEX_ECS_NOPRECONDITIONER_H

// --------------------------------------------------------------------------------- //

#include "preconditioners.h"

// --------------------------------------------------------------------------------- //

/**
 * @brief Solution driver without actual preconditioner.
 * 
 * This class "preconditions" by identity matrix, but implements all other important
 * routines, that can be used by derived classes, namely:
 * - setup : Loads / computed radial integrals for construction of the matrix of the set
 *           and for the construction of the right-hand side.
 * - update : Creates the diagonal blocks.
 * - finish : Cleanup of memory. Un-do the setup.
 * - rhs : Composes the right-hand side.
 * - multiply : Multiplies a vector by the matrix of the set of equations.
 * - precondition : This preconditioner uses identity as its matrix, but overload do more.
 */
class NoPreconditioner : public PreconditionerBase
{
    public:
        
        // run-time selection mechanism
        preconditionerRunTimeSelectionDefinitions(NoPreconditioner, "none")
        
        // default constructor required by RTS system
        NoPreconditioner ();
        
        // constructor
        NoPreconditioner
        (
            Parallel const & par,
            InputFile const & inp,
            AngularBasis const & ll,
            Bspline const & bspline_inner,
            Bspline const & bspline_full,
            CommandLine const & cmd
        );
        
        // destructor
        ~NoPreconditioner ();
        
        // description of the preconditioner
        virtual std::string description () const;
        
        // member functions
        virtual void setup ();
        virtual void update (Real E);
        virtual void finish ();
        virtual void rhs (BlockArray<Complex> & chi, int ienergy, int instate) const;
        virtual void multiply (BlockArray<Complex> const & p, BlockArray<Complex> & q, MatrixSelection::Selection tri = MatrixSelection::Both) const;
        virtual void precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const;
        
        // internal routines
        BlockSymBandMatrix<Complex> calc_A_block (int ill, int illp, bool twoel = true) const;
        RadialIntegrals const & rad () const { return *rad_; }
        
    protected:
        
        // energy
        Real E_;
        
        // command line switches
        CommandLine const * cmd_;
        
        // parallel environment
        Parallel const * par_;
        
        // input parameters
        InputFile const * inp_;
        
        // coupled states
        AngularBasis const * ang_;
        
        // Sub-blocks composing the angular blocks of the full matrix:
        //  ┏━━━━━━┯━━━━━━━━┓
        //  ┃ A    │  Cu    ┃
        //  ┃    A │        ┃
        //  ┠──────┼───┬────┨
        //  ┃      │B1 │  0 ┃
        //  ┃ Cl   ├───┼────┨
        //  ┃      │ 0 │ B2 ┃
        //  ┗━━━━━━┷━━━┷━━━━┛
        // The off-diagonal blocks Cu and Cl are actually stored as COO matrices with the dimension
        // of the whole matrix, but only the elements of the respective blocks are non-zero.
        std::vector<BlockSymBandMatrix<Complex>> A_blocks_;
        std::vector<std::vector<SymBandMatrix<Complex>>> B1_blocks_;
        std::vector<std::vector<SymBandMatrix<Complex>>> B2_blocks_;
        std::vector<CooMatrix<LU_int_t,Complex>> Cu_blocks_;
        std::vector<CooMatrix<LU_int_t,Complex>> Cl_blocks_;
        
        // maximal bound state principal quantum number for given energy
        int max_n_;
        
        // number of channels when r1 -> inf and r2 -> inf, respectively
        std::vector<std::pair<int,int>> Nchan_;
        
        // radial integrals for the solution
        RadialIntegrals * rad_;
        
        // hydrogen orbitals B-spline overlaps and expansions (on inner basis)
        std::vector<cArrays> Spx, Xpx, Spy, Xpy;
};

// --------------------------------------------------------------------------------- //

#endif // HEX_ECS_NOPRECONDITIONER_H
