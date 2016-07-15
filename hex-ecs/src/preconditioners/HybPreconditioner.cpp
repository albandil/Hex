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

#include <algorithm>

#include "preconditioners.h"

const std::string HybCGPreconditioner::prec_name = "HYB";
const std::string HybCGPreconditioner::prec_description = "Combination of ILU and KPA.";

bool HybCGPreconditioner::ilu_needed (int iblock) const
{
//     std::cout << "  1 channels = " << Nchan_[iblock].first << std::endl;
//     std::cout << "  2 channels = " << Nchan_[iblock].second << std::endl;
    return Nchan_[iblock].first > 0 or Nchan_[iblock].second > 0;
}

void HybCGPreconditioner::setup ()
{
    ILUCGPreconditioner::setup();
    KPACGPreconditioner::setup();
}

void HybCGPreconditioner::update (Real E)
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
    }
    
    // update common ancestor
    CGPreconditioner::update(E);
    
    // reset counters
    CGPreconditioner::n_.fill(-1);
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
//     std::cout << "HYB:CG_mmul:" << std::endl;
//     std::cout << "  iblock = " << iblock << std::endl;
//     std::cout << "  l1 = " << ang_.states()[iblock].first << std::endl;
//     std::cout << "  l2 = " << ang_.states()[iblock].second << std::endl;
//     std::cout << "  r.size() = " << r.size() << std::endl;
//     std::cout << "  z.size() = " << z.size() << std::endl;
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
