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

#include "../arrays.h"
#include "../misc.h"
#include "../preconditioners.h"

const std::string ILUCGPreconditioner::prec_name = "ILU";
const std::string ILUCGPreconditioner::prec_description = 
    "Block inversion using conjugate gradients preconditioned by Incomplete LU. "
    "The drop tolerance can be given as the --droptol parameter.";

void ILUCGPreconditioner::setup ()
{
    // setup parent
    CGPreconditioner::setup();
    
    // prepare initial (empty) factorization data
    for (auto & lu : lu_)
        lu.reset(new LUft<LU_int_t,Complex>());
    
#ifdef WITH_SUPERLU_DIST
    // create process grid for SuperLU-dist
    if (cmd_.factorizer == LUFT_SUPERLU_DIST)
    {
        // list processes
        for (int igroup = 0; igroup < par_.Nproc() / cmd_.groupsize; igroup++)
        {
            // list member processes
            NumberArray<int_t> usermap;
            for (int iproc = 0; iproc < cmd_.groupsize; iproc++)
                usermap.push_back(igroup * cmd_.groupsize + iproc);
            
            // create the grid
            gridinfo_t grid;
            superlu_gridmap(MPI_COMM_WORLD, cmd_.groupsize, 1, &usermap[0], cmd_.groupsize, &grid);
            
            // assign this process to the grid, if it is member of current group
            if (par_.isMyGroupWork(igroup))
                grid_ = grid;
        }
    }
#endif // WITH_SUPERLU_DIST
}

void ILUCGPreconditioner::update (double E)
{
    // reset data on energy change
    if (E != E_)
    {
        // release outdated LU factorizations
        for (auto & lu : lu_)
        {
            lu->drop();
            lu->unlink();
        }
        
        // release outdated CSR diagonal blocks
        for (auto & csr : csr_blocks_)
        {
            csr.drop();
            csr.unlink();
        }
    }
    
    // update parent
    CGPreconditioner::update(E);
}

void ILUCGPreconditioner::CG_prec (int iblock, const cArrayView r, cArrayView z) const
{
    // load data from linked disk files
    if (cmd_.outofcore)
    {
        csr_blocks_[iblock].hdfload();
        # pragma omp critical
        lu_[iblock]->silent_load();
    }
    
    // check that the factorization is loaded
    if (lu_[iblock]->size() == 0)
    {
        // create CSR block
        csr_blocks_[iblock] = dia_blocks_[iblock].tocoo<LU_int_t>().tocsr();
        
        // start timer
        Timer timer;
        
        // factorize the block and store it
        lu_[iblock] = csr_blocks_[iblock].factorize
        (
                cmd_.droptol,
                cmd_.factorizer,
#ifdef WITH_SUPERLU_DIST
                cmd_.factorizer == LUFT_SUPERLU_DIST ?
                    const_cast<gridinfo_t*>(&grid_) :
#endif
                    nullptr
        );
        
        // print time and memory info for this block (one thread at a time)
        # pragma omp critical
        std::cout << std::endl << std::setw(37) << format
        (
            "\tLU #%d (%d,%d) in %d:%02d (%d MiB)",
            iblock, l1_l2_[iblock].first, l1_l2_[iblock].second,    // block identification (id, ℓ₁, ℓ₂)
            timer.seconds() / 60, timer.seconds() % 60,             // factorization time
            lu_[iblock]->size() / 1048576                           // final memory size
        );
        
        // save the diagonal block
        if (cmd_.outofcore)
        {
            csr_blocks_[iblock].hdflink(format("csr-%d.ooc", iblock));
            csr_blocks_[iblock].hdfsave();
        }
    }
    
    // precondition by LU
    lu_[iblock]->solve(r, z, 1);
    
    // release memory
    if (cmd_.outofcore)
    {
        // link to a disk file and save (if not already done)
        if (lu_[iblock]->name().size() == 0)
        {
            lu_[iblock]->link(format("lu-%d.ooc", iblock));
            # pragma omp critical
            lu_[iblock]->save();
        }
        
        // release memory objects
        lu_[iblock]->drop();
        csr_blocks_[iblock].drop();
    }
}
