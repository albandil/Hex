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

#include <algorithm>

#include "../preconditioners.h"

const std::string HybCGPreconditioner::prec_name = "HYB";
const std::string HybCGPreconditioner::prec_description = "Combination of ILU and KPA.";

bool HybCGPreconditioner::ilu_needed (int iblock) const
{
    // decide which preconditioner to use
    if (CGPreconditioner::n_[iblock] >= 0 and
        NoPreconditioner::cmd_.kpa_max_iter >= 0 and
        CGPreconditioner::n_[iblock] > NoPreconditioner::cmd_.kpa_max_iter)
    {
        // count ILU blocks
        int nILU = std::count_if
        (
            prec_.begin(),
            prec_.end(),
            [](int i) -> bool { return i == UseILU; }
        );
        
        // only allow next factorization if fitting in restriction
        if (NoPreconditioner::cmd_.ilu_max_blocks < 0 or
            NoPreconditioner::cmd_.ilu_max_blocks > nILU)
            prec_[iblock] = UseILU;
    }
    
    return prec_[iblock] == UseILU;
}

void HybCGPreconditioner::setup ()
{
    ILUCGPreconditioner::setup();
    KPACGPreconditioner::setup();
}

void HybCGPreconditioner::update (double E)
{
    // update ILU
    if (E != CGPreconditioner::E_)
    {
        // release outdated LU factorizations
        for (auto & lu : ILUCGPreconditioner::lu_)
        {
            lu->drop();
            lu->unlink();
        }
        
        // release outdated CSR diagonal blocks
        for (auto & csr : ILUCGPreconditioner::csr_blocks_)
        {
            csr.drop();
            csr.unlink();
        }
    }
    
    // update common ancestor
    CGPreconditioner::update(E);
    
    // reset counters
    CGPreconditioner::n_.fill(0);
    prec_.fill(Undecided);
}

void HybCGPreconditioner::CG_init (int iblock) const
{
    if (ilu_needed(iblock))
        ILUCGPreconditioner::CG_init(iblock);
    else
        KPACGPreconditioner::CG_init(iblock);
}

void HybCGPreconditioner::CG_mmul (int iblock, const cArrayView r, cArrayView z) const
{
    if (ilu_needed(iblock))
        ILUCGPreconditioner::CG_mmul(iblock, r, z);
    else
        KPACGPreconditioner::CG_mmul(iblock, r, z);
}

void HybCGPreconditioner::CG_prec (int iblock, const cArrayView r, cArrayView z) const
{
    if (ilu_needed(iblock))
        ILUCGPreconditioner::CG_prec(iblock, r, z);
    else
        KPACGPreconditioner::CG_prec(iblock, r, z);
}

void HybCGPreconditioner::CG_exit (int iblock) const
{
    ILUCGPreconditioner::CG_exit(iblock);
    KPACGPreconditioner::CG_exit(iblock);
}

void HybCGPreconditioner::finish()
{
    ILUCGPreconditioner::finish();
    KPACGPreconditioner::finish();
}
