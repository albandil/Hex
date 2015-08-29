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

#ifndef HEX_NOPRECONDITIONER_H
#define HEX_NOPRECONDITIONER_H

#include "../preconditioners.h"

/**
 * @brief Solution driver without actual preconditioner.
 * 
 * This class "preconditions" by identity matrix, but implements all other important
 * routines, that can be used by derived classes, namely:
 * - setup : Loads / computed radial integrals for construction of the matrix of the set
 *           and for the construction of the right-hand side.
 * - update : Creates the diagonal blocks.
 * - rhs : Composes the right-hand side.
 * - multiply : Multiplies a vector by the matrix of the set of equations.
 */
class NoPreconditioner : public PreconditionerBase
{
    public:
        
        static const std::string prec_name;
        static const std::string prec_description;
        
        virtual std::string const & name () const { return prec_name; }
        virtual std::string const & description () const { return prec_description; }
        
        NoPreconditioner
        (
            Parallel const & par,
            InputFile const & inp,
            std::vector<std::pair<int,int>> const & ll,
            Bspline const & bspline_atom,
            Bspline const & bspline_proj,
            CommandLine const & cmd
        ) : PreconditionerBase(), cmd_(cmd), par_(par), inp_(inp), l1_l2_(ll),
            dia_blocks_(l1_l2_.size()), bspline_atom_(bspline_atom), bspline_proj_(bspline_proj),
            rad_(bspline_atom, bspline_proj, inp.L + 2 * inp.levels + 1)
        {
            // nothing to do
        }
        
        virtual void setup ();
        virtual void update (double E);
        virtual void finish ();
        virtual void rhs (BlockArray<Complex> & chi, int ienergy, int instate, int Spin, Bspline const & bfull) const;
        virtual void multiply (BlockArray<Complex> const & p, BlockArray<Complex> & q) const;
        virtual void precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const { z = r; }
        
    protected:
        
        // energy
        double E_;
        
        // command line switches
        CommandLine const & cmd_;
        
        // parallel environment
        Parallel const & par_;
        
        // input parameters
        InputFile const & inp_;
        
        // coupled states
        std::vector<std::pair<int,int>> const & l1_l2_;
        
        // diagonal blocks in DIA format (these will be used in matrix multiplication)
        mutable std::vector<BlockSymBandMatrix<Complex>> dia_blocks_;
        
        // B-spline environment for the solution
        Bspline const & bspline_atom_;
        Bspline const & bspline_proj_;
            
        // radial integrals for the solution
        RadialIntegrals rad_;
};

#endif
