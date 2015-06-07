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

#include <iostream>
#include <cstdio>

#include "../arrays.h"
#include "../itersolve.h"
#include "../misc.h"
#include "../preconditioners.h"

const std::string CGPreconditioner::prec_name = "CG";
const std::string CGPreconditioner::prec_description = 
    "Block inversion using plain conjugate gradients. "
    "Use --tolerance option to set the termination tolerance.";

void CGPreconditioner::precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const
{
    // shorthands
    int Nspline = s_rad_.bspline().Nspline();
    
    // iterations
    iArray n (l1_l2_.size());
    
    # pragma omp parallel for schedule (dynamic, 1) if (cmd_.parallel_block && cmd_.groupsize == 1)
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyGroupWork(ill))
    {
        try
        {
            // load blocks, if necessary
            if (cmd_.outofcore)
            {
                const_cast<BlockArray<Complex>&>(r).hdfload(ill);
                z.hdfload(ill);
            }
            
            // wrappers around the callbacks
            auto inner_mmul = [&](const cArrayView a, cArrayView b) { this->CG_mmul(ill, a, b); };
            auto inner_prec = [&](const cArrayView a, cArrayView b) { this->CG_prec(ill, a, b); };
            
            // prepare the block-preconditioner for run
            this->CG_init(ill);
            
            // solve using the CG solver
            ConjugateGradients < cArray, cArrayView > CG;
            CG.reset();
            n[ill] = CG.solve
            (
                r[ill],                 // rhs
                z[ill],                 // solution
                cmd_.prec_itertol,      // preconditioner tolerance
                0,                      // min. iterations
                Nspline * Nspline,      // max. iteration
                inner_prec,             // preconditioner
                inner_mmul,             // matrix multiplication
                false                   // verbose output
            );
            
            // release block-preconditioner block-specific data
            this->CG_exit(ill);
            
            // unload blocks
            if (cmd_.outofcore)
            {
                const_cast<BlockArray<Complex>&>(r)[ill].drop();
                
                z.hdfsave(ill);
                z[ill].drop();
            }
        }
        catch (std::exception const & e)
        {
            HexException("%s", e.what());
        }
    }
    
    // broadcast inner preconditioner iterations
    par_.sync(n.data(), 1, l1_l2_.size());
    
    // inner preconditioner info (max and avg number of iterations)
    std::cout << " | ";
    std::cout << std::setw(5) << (*std::min_element(n.begin(), n.end()));
    std::cout << std::setw(5) << (*std::max_element(n.begin(), n.end()));
    std::cout << std::setw(5) << format("%g", std::accumulate(n.begin(), n.end(), 0) / float(n.size()));
}

void CGPreconditioner::CG_init (int iblock) const
{
    if (cmd_.outofcore and cmd_.wholematrix)
        dia_blocks_[iblock].hdfload();
}

void CGPreconditioner::CG_mmul (int iblock, const cArrayView p, cArrayView q) const
{
    if (cmd_.lightweight_full)
        HexException("Preconditioner %s is not compatible with the option --lightweight-full.", this->name().c_str());
        
    q = dia_blocks_[iblock].dot(p, true);
}

void CGPreconditioner::CG_prec (int iblock, const cArrayView r, cArrayView z) const
{
    z = r;
}

void CGPreconditioner::CG_exit (int iblock) const
{
    if (cmd_.outofcore and cmd_.wholematrix)
        dia_blocks_[iblock].drop();
}
