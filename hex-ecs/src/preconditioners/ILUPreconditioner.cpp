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

#include <iostream>

#include "hex-arrays.h"
#include "hex-luft.h"
#include "hex-misc.h"

#include "preconditioners.h"

const std::string ILUCGPreconditioner::prec_name = "ILU";
const std::string ILUCGPreconditioner::prec_description = 
    "Block inversion using conjugate gradients preconditioned by Incomplete LU. "
    "The drop tolerance can be given as the --droptol parameter.";

ILUCGPreconditioner::ILUCGPreconditioner
(
    Parallel const & par,
    InputFile const & inp,
    AngularBasis const & ll,
    Bspline const & bspline_inner,
    Bspline const & bspline_outer,
    Bspline const & bspline_full,
    CommandLine const & cmd
) : CGPreconditioner(par, inp, ll, bspline_inner, bspline_outer, bspline_full, cmd),
    csr_blocks_(ll.states().size()),
    lu_(ll.states().size())
{
#ifdef _OPENMP
    omp_init_lock(&lu_lock_);
#endif
}

ILUCGPreconditioner::~ILUCGPreconditioner ()
{
#ifdef _OPENMP
    omp_destroy_lock(&lu_lock_);
#endif
}

void ILUCGPreconditioner::setup ()
{
    // setup parent
    CGPreconditioner::setup();
    
    // check compatibility with the command line setup
    if (cmd_.lightweight_full)
        HexException("ILU preconditioner is not compatible with --lightweight-full. But you can still try --lightweight-radial-cache.");
    
    // setup attributes
    for (unsigned iblock = 0; iblock < ang_.states().size(); iblock++)
    {
        // prepare initial (empty) factorization data
        lu_[iblock].reset(LUft<LU_int_t,Complex>::New(cmd_.factorizer));
        
#ifdef WITH_UMFPACK
        if (cmd_.factorizer == LUFT_UMFPACK)
        {
            // associate the matrix pointer
            LUft_UMFPACK<LU_int_t,Complex> * ptr = dynamic_cast<LUft_UMFPACK<LU_int_t,Complex>*>(lu_[iblock].get());
            if (ptr)
                ptr->matrix(&csr_blocks_[iblock]);
        }
#endif
        
        // associate existing disk files
        lu_[iblock]->link(format("lu-%d.bin", iblock));
    }
    
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

void ILUCGPreconditioner::update (Real E)
{
    // reset data on energy change
    if (E != E_ and not cmd_.noluupdate)
    {
        // release outdated LU factorizations
        for (auto & lu : lu_)
            lu->drop();
        
        // release outdated CSR diagonal blocks
        for (auto & csr : csr_blocks_)
            csr.drop();
    }
    
    // update parent
    CGPreconditioner::update(E);
}

void ILUCGPreconditioner::CG_init (int iblock) const
{
    // update parent
    CGPreconditioner::CG_init(iblock);
    
    // load data from linked disk files
    if (lu_[iblock]->size() == 0)
        lu_[iblock]->silent_load();
    
    // check that the factorization is loaded
    if (lu_[iblock]->size() == 0)
    {
#ifdef _OPENMP
        // allow only one factorization at a time when not using SuperLU DIST
        if (cmd_.factorizer != LUFT_SUPERLU_DIST and not cmd_.parallel_factorization)
            omp_set_lock(&lu_lock_);
#endif
        
        // start timer
        Timer timer;
        
        // number of asymptotic channels
        int Nchan1 = Nchan_[iblock].first;
        int Nchan2 = Nchan_[iblock].second;
        
        // number of B-splines
        LU_int_t Nspline_inner = rad_.bspline_inner().Nspline();
        LU_int_t Nspline_outer = rad_.bspline_outer().Nspline();
        
        // angular block
        int iang = iblock * ang_.states().size() + iblock;
        
        // convert inner region matrix block to COO matrix
        CooMatrix<LU_int_t,Complex> A_coo = A_blocks_[iang].tocoo<LU_int_t>();
        A_coo.tocsr().tocoo().write("A-coo.txt"); // DEBUG
        
        // add the A-block
        CooMatrix<LU_int_t,Complex> coo_block
        (
            Nspline_inner * Nspline_inner + (Nchan1 + Nchan2) * Nspline_outer,
            Nspline_inner * Nspline_inner + (Nchan1 + Nchan2) * Nspline_outer,
            A_coo.i(), A_coo.j(), A_coo.v()
        );
        
        if (not inp_.inner_only)
        {
            // add the outer region C-blocks
            coo_block += Cu_blocks_[iang];
            coo_block += Cl_blocks_[iang];
            Cu_blocks_[iang].write("Cu-coo.txt"); // DEBUG
            Cl_blocks_[iang].write("Cl-coo.txt"); // DEBUG
            
            // add the B-blocks
            for (int m = 0; m < Nchan1; m++)
            for (int n = 0; n < Nchan1; n++)
            {
                CooMatrix<LU_int_t,Complex> B_coo_small = B1_blocks_[iang][m * Nchan1 + n].tocoo<LU_int_t>();
                CooMatrix<LU_int_t,Complex> B_coo_large
                (
                    coo_block.rows(), coo_block.cols(),
                    B_coo_small.i() + Nspline_inner * Nspline_inner + m * Nspline_outer,
                    B_coo_small.j() + Nspline_inner * Nspline_inner + n * Nspline_outer,
                    B_coo_small.v()
                );
                coo_block += B_coo_large;
                B_coo_large.tocsr().tocoo().write(format("B1-%d-%d-coo.txt", m, n)); // DEBUG
            }
            for (int m = 0; m < Nchan2; m++)
            for (int n = 0; n < Nchan2; n++)
            {
                CooMatrix<LU_int_t,Complex> B_coo_small = B2_blocks_[iang][m * Nchan2 + n].tocoo<LU_int_t>();
                CooMatrix<LU_int_t,Complex> B_coo_large
                (
                    coo_block.rows(), coo_block.cols(),
                    B_coo_small.i() + Nspline_inner * Nspline_inner + (Nchan1 + m) * Nspline_outer,
                    B_coo_small.j() + Nspline_inner * Nspline_inner + (Nchan1 + n) * Nspline_outer,
                    B_coo_small.v()
                );
                coo_block += B_coo_large;
                B_coo_large.tocsr().tocoo().write(format("B2-%d-%d-coo.txt", m, n)); // DEBUG
            }
        }
        
        // create the CSR block that will be factorized
        csr_blocks_[iblock] = coo_block.tocsr();
        csr_blocks_[iblock].tocoo().write("full-coo.txt"); // DEBUG
        csr_blocks_[iblock].plot(format("csr-%d.png", iblock)); // DEBUG
        
        // factorize the block and store it
        lu_[iblock] = csr_blocks_[iblock].factorize
        (
            cmd_.droptol,
            cmd_.factorizer,
#ifdef WITH_SUPERLU_DIST
            // pass pointer to grid info if needed
            cmd_.factorizer == LUFT_SUPERLU_DIST ? const_cast<gridinfo_t*>(&grid_) : nullptr
#else
            // do not pass any additional data
            nullptr
#endif
        );
        
        // print time and memory info for this block (one thread at a time)
        # pragma omp critical
        std::cout << std::endl << std::setw(37) << format
        (
            "\tLU #%d (%d,%d) in %d:%02d (%s, cond %1.0e)",
            iblock, ang_.states()[iblock].first, ang_.states()[iblock].second,      // block identification (id, ℓ₁, ℓ₂)
            timer.seconds() / 60, timer.seconds() % 60,                             // factorization time
            nice_size(lu_[iblock]->size()).c_str(),                                 // final memory size
            lu_[iblock]->cond()                                                     // estimation of the condition number
        );
        
        // save the diagonal block's CSR representation and its factorization
        lu_[iblock]->link(format("lu-%d.bin", iblock));
//         lu_[iblock]->save();
        
#ifdef _OPENMP
        // release lock
        if (cmd_.factorizer != LUFT_SUPERLU_DIST and not cmd_.parallel_factorization)
            omp_unset_lock(&lu_lock_);
#endif
    }
}

void ILUCGPreconditioner::CG_prec (int iblock, const cArrayView r, cArrayView z) const
{
    // precondition by LU
    lu_[iblock]->solve(r, z, 1);
}

void ILUCGPreconditioner::CG_exit (int iblock) const
{
    // release memory
    if (cmd_.outofcore)
    {
        csr_blocks_[iblock].drop();
        lu_[iblock]->drop();
    }
    
    // exit parent
    CGPreconditioner::CG_exit(iblock);
}

void ILUCGPreconditioner::finish ()
{
    lu_.resize(0);
    CGPreconditioner::finish();
}
