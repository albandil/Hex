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
#include "NoPreconditioner.h"

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
    Bspline const & bspline_full,
    CommandLine const & cmd
) : CGPreconditioner(par, inp, ll, bspline_inner, bspline_full, cmd),
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

void ILUCGPreconditioner::reset_lu ()
{
    for (unsigned iblock = 0; iblock < ang_.states().size(); iblock++)
    {
        // prepare initial (empty) factorization data
        lu_[iblock].reset(LUft<LU_int_t,Complex>::New(cmd_.factorizer));
        
        // associate existing disk files
        lu_[iblock]->link(format("lu-%d.ooc", iblock));
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

void ILUCGPreconditioner::setup ()
{
    // setup parent
    CGPreconditioner::setup();
    
    // prepare data structures for LU factorizations
    reset_lu();
}

void ILUCGPreconditioner::update (Real E)
{
    // reset data on energy change
    if (E != E_ and not cmd_.noluupdate)
    {
        // release outdated LU factorizations
        reset_lu();
    }
    
    // update parent
    CGPreconditioner::update(E);
}

void ILUCGPreconditioner::CG_init (int iblock) const
{
    // update parent
    CGPreconditioner::CG_init(iblock);
    
    // load data from linked disk files
    if (not lu_[iblock]->valid())
        lu_[iblock]->silent_load();
    
    // check that the factorization is loaded
    if (not lu_[iblock]->valid())
    {
#ifdef _OPENMP
        // allow only one factorization at a time when not using SuperLU DIST or MUMPS
        if (not cmd_.parallel_factorization)
            omp_set_lock(&lu_lock_);
#endif
        
        // start timer
        Timer timer;
        
        // number of asymptotic channels
        int Nchan1 = Nchan_[iblock].first;
        int Nchan2 = Nchan_[iblock].second;
        
        // number of B-splines
        LU_int_t Nspline_inner = rad_.bspline_inner().Nspline();
        LU_int_t Nspline_outer = rad_.bspline_full().Nspline() - Nspline_inner;
        
        // angular block
        int iang = iblock * ang_.states().size() + iblock;
        
        // convert inner region matrix block to COO matrix
        CooMatrix<LU_int_t,Complex> coo_block;
        if (cmd_.lightweight_full)
            coo_block = std::move(dynamic_cast<NoPreconditioner const*>(this)->calc_A_block(iblock, iblock).tocoo<LU_int_t>());
        else
            coo_block = std::move(A_blocks_[iang].tocoo<LU_int_t>());
        
        // add the A-block
        coo_block.resize
        (
            Nspline_inner * Nspline_inner + (Nchan1 + Nchan2) * Nspline_outer,
            Nspline_inner * Nspline_inner + (Nchan1 + Nchan2) * Nspline_outer
        );
        
        if (not inp_.inner_only)
        {
            // add the outer region C-blocks
            coo_block += Cu_blocks_[iang];
            coo_block += Cl_blocks_[iang];
            
            // add the B-blocks
            for (int m = 0; m < Nchan1; m++)
            for (int n = 0; n < Nchan1; n++)
            {
                if (cmd_.outofcore) const_cast<SymBandMatrix<Complex>&>(B1_blocks_[iang][m * Nchan1 + n]).hdfload();
                CooMatrix<LU_int_t,Complex> B_coo_small = B1_blocks_[iang][m * Nchan1 + n].tocoo<LU_int_t>();
                if (cmd_.outofcore) const_cast<SymBandMatrix<Complex>&>(B1_blocks_[iang][m * Nchan1 + n]).drop();
                CooMatrix<LU_int_t,Complex> B_coo_large
                (
                    coo_block.rows(), coo_block.cols(),
                    B_coo_small.i() + Nspline_inner * Nspline_inner + m * Nspline_outer,
                    B_coo_small.j() + Nspline_inner * Nspline_inner + n * Nspline_outer,
                    B_coo_small.v()
                );
                coo_block += B_coo_large;
            }
            for (int m = 0; m < Nchan2; m++)
            for (int n = 0; n < Nchan2; n++)
            {
                if (cmd_.outofcore) const_cast<SymBandMatrix<Complex>&>(B2_blocks_[iang][m * Nchan2 + n]).hdfload();
                CooMatrix<LU_int_t,Complex> B_coo_small = B2_blocks_[iang][m * Nchan2 + n].tocoo<LU_int_t>();
                if (cmd_.outofcore) const_cast<SymBandMatrix<Complex>&>(B2_blocks_[iang][m * Nchan2 + n]).drop();
                CooMatrix<LU_int_t,Complex> B_coo_large
                (
                    coo_block.rows(), coo_block.cols(),
                    B_coo_small.i() + Nspline_inner * Nspline_inner + (Nchan1 + m) * Nspline_outer,
                    B_coo_small.j() + Nspline_inner * Nspline_inner + (Nchan1 + n) * Nspline_outer,
                    B_coo_small.v()
                );
                coo_block += B_coo_large;
            }
        }
        
        // create the CSR block that will be factorized
        CsrMatrix<LU_int_t,Complex> csr = std::move(coo_block.tocsr());
        coo_block = CooMatrix<LU_int_t,Complex>();
        
        // set up factorization data
        void * data = nullptr;
#ifdef WITH_SUPERLU_DIST
        if (cmd_.factorizer == LUFT_SUPERLU_DIST)
            data = const_cast<gridinfo_t*>(&grid_);
#endif
#ifdef WITH_MUMPS
        if (cmd_.factorizer == LUFT_MUMPS)
            data = (void*)std::intptr_t(cmd_.mumps_outofcore + 2 * cmd_.mumps_verbose);
#endif
        
        // factorize the block and store it
        lu_[iblock] = csr.factorize(cmd_.droptol, cmd_.factorizer, data);
        
        // print time and memory info for this block (one thread at a time)
        # pragma omp critical
        {
            if (lu_[iblock]->cond() != 0)
            {
                std::cout << std::endl << std::setw(37) << format
                (
                    "\tLU #%d (%d,%d) in %d:%02d (%s, cond %1.0e)",
                    iblock, ang_.states()[iblock].first, ang_.states()[iblock].second,      // block identification (id, ℓ₁, ℓ₂)
                    timer.seconds() / 60, timer.seconds() % 60,                             // factorization time
                    nice_size(lu_[iblock]->size()).c_str(),                                 // final memory size
                    lu_[iblock]->cond()                                                     // estimation of the condition number
                );
            }
            else
            {
                std::cout << std::endl << std::setw(37) << format
                (
                    "\tLU #%d (%d,%d) in %d:%02d (%s)",
                    iblock, ang_.states()[iblock].first, ang_.states()[iblock].second,      // block identification (id, ℓ₁, ℓ₂)
                    timer.seconds() / 60, timer.seconds() % 60,                             // factorization time
                    nice_size(lu_[iblock]->size()).c_str()                                  // final memory size
                );
            }
        }
        
        // save the diagonal block's CSR representation and its factorization
        lu_[iblock]->link(format("lu-%d.bin", iblock));
        if (cmd_.outofcore)
        {
            lu_[iblock]->save();
        }
        
#ifdef _OPENMP
        // release lock
        if (not cmd_.parallel_factorization)
            omp_unset_lock(&lu_lock_);
#endif
    }
}

void ILUCGPreconditioner::CG_prec (int iblock, const cArrayView r, cArrayView z) const
{
#ifdef _OPENMP
    if (cmd_.factorizer == LUFT_MUMPS)
        omp_set_lock(&lu_lock_);
#endif
    
    // precondition by LU
    lu_[iblock]->solve(r, z, 1);
    
#ifdef _OPENMP
    if (cmd_.factorizer == LUFT_MUMPS)
        omp_unset_lock(&lu_lock_);
#endif    
}

void ILUCGPreconditioner::CG_exit (int iblock) const
{
#ifdef _OPENMP
    if (cmd_.factorizer == LUFT_MUMPS)
        omp_set_lock(&lu_lock_);
#endif
    
    // release memory
    if (cmd_.outofcore)
    {
        lu_[iblock]->drop();
    }
    
#ifdef _OPENMP
    if (cmd_.factorizer == LUFT_MUMPS)
        omp_unset_lock(&lu_lock_);
#endif
    
    // exit parent
    CGPreconditioner::CG_exit(iblock);
}

void ILUCGPreconditioner::finish ()
{
#ifdef _OPENMP
    if (cmd_.factorizer == LUFT_MUMPS)
        omp_set_lock(&lu_lock_);
#endif
    
    lu_.resize(0);
    
#ifdef _OPENMP
    if (cmd_.factorizer == LUFT_MUMPS)
        omp_unset_lock(&lu_lock_);
#endif
    
    CGPreconditioner::finish();
}
