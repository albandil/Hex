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

#include "GMGPreconditioner.h"

const std::string GMGPreconditioner::prec_name = "GMG";
const std::string GMGPreconditioner::prec_description = "Geometric multigrid.";

GMGPreconditioner::GMGPreconditioner
(
    Parallel const & par,
    InputFile const & inp,
    AngularBasis const & ll,
    Bspline const & bspline_inner,
    Bspline const & bspline_full,
    CommandLine const & cmd
) : GMGPreconditioner(par, inp, ll, bspline_inner, bspline_full, cmd, cmd.multigrid_depth)
{
    // nothing
}

GMGPreconditioner::GMGPreconditioner
(
    Parallel const & par,
    InputFile const & inp,
    AngularBasis const & ll,
    Bspline const & bspline_inner,
    Bspline const & bspline_full,
    CommandLine const & cmd,
    int level
) : NoPreconditioner(par, inp, ll, bspline_inner, bspline_full, cmd), level_(level)
{
    // setup subgrid
    if (level_ > 0)
    {
        subgrid_ = new GMGPreconditioner
        (
            par, inp, ll,
            bspline_inner, bspline_full,
            cmd, level_ - 1
        );
    }
    else
    {
        subgrid_ = Preconditioners::choose
        (
            cmd.multigrid_coarse_prec,
            par, inp, ll,
            bspline_inner, bspline_full,
            cmd
        );
    }
}

void GMGPreconditioner::setup ()
{
    // setup self
    if (level_ > 0)
        NoPreconditioner::setup();
    
    // setup subgrid
    subgrid_->setup();
}

void GMGPreconditioner::update (Real E)
{
    // update self
    if (level_ > 0)
        NoPreconditioner::update(E);
    
    // update subgrid
    subgrid_->update(E);
}

void GMGPreconditioner::precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const
{
    BlockArray<Complex> rn (r.size()), zn (z.size());
    
    // restriction: down-sample the residual
    if (level_ > 0)
    {
        
    }
    
    // solution: precondition by sub-grid
    subgrid_->precondition(rn, zn);
    
    // prolongation: interpolate the solution
    if (level_ > 0)
    {
        
    }
    
    // correct high-frequency error using Gauss-Seidel iterations
    int nSmoothCycles = 5;
    for (int i = 0; i < nSmoothCycles; i++)
    {
        
    }
}
