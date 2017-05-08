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

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"
#include "hex-csrmatrix.h"
#include "hex-misc.h"

// --------------------------------------------------------------------------------- //

#include "luft.h"

// --------------------------------------------------------------------------------- //

#include "NoPreconditioner.h"
#include "ILUPreconditioner.h"

// --------------------------------------------------------------------------------- //

std::string ILUCGPreconditioner::description () const
{
    return "Block inversion using conjugate gradients preconditioned by Incomplete LU. "
    "The drop tolerance can be given as the --droptol parameter.";
}

ILUCGPreconditioner::ILUCGPreconditioner
(
    CommandLine  const & cmd,
    InputFile    const & inp,
    Parallel     const & par,
    AngularBasis const & ang,
    Bspline const & bspline_x_inner,
    Bspline const & bspline_x_full,
    Bspline const & bspline_y_inner,
    Bspline const & bspline_y_full
) : CGPreconditioner
    (
        cmd, inp, par, ang,
        bspline_x_inner, bspline_x_full,
        bspline_y_inner, bspline_y_full
    ),
    data_(ang.states().size()),
    lu_(ang.states().size())
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
    
#ifdef WITH_SUPERLU_DIST
    if (cmd_->factorizer == "superlu_dist")
        superlu_gridexit(&grid_);
#endif
}

void ILUCGPreconditioner::reset_lu ()
{
    for (unsigned iblock = 0; iblock < ang_->states().size(); iblock++)
    {
        // prepare initial (empty) factorization data
        lu_[iblock].reset(LUft::Choose(cmd_->factorizer));
        
        // feedback to the user
        if (iblock == 0 and cmd_->factorizer == "any")
            std::cout << "\tUsing factorizer \"" << lu_[iblock]->name() << "\" (use --lu <factorizer> to change this)" << std::endl;
        
        // associate existing disk files
        lu_[iblock]->link(format("lu-%d.ooc", iblock));
    }
    
#ifdef WITH_SUPERLU_DIST
    // create process grid for SuperLU-dist
    if (cmd_->factorizer == "superlu_dist")
    {
        int_t nprow = std::sqrt(par_->groupsize());
        int_t npcol = par_->groupsize() / nprow;
        
        while (nprow * npcol != par_->groupsize())
        {
            nprow--;
            npcol = par_->groupsize() / nprow;
        }
        
        superlu_gridinit((MPI_Comm)par_->groupcomm(), nprow, npcol, &grid_);
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
    if (E != E_ and not cmd_->noluupdate)
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
        if (not cmd_->parallel_factorization)
            omp_set_lock(&lu_lock_);
#endif
        
        // start timer
        Timer timer;
        
        // number of asymptotic channels
        int Nchan1 = Nchan_[iblock].first;
        int Nchan2 = Nchan_[iblock].second;
        
        // panel x basis
        LU_int_t Nspline_x = rad_panel_->bspline_x().Nspline();
        LU_int_t Nspline_inner_x = rad_panel_->bspline_x().knot(rad_inner_->bspline_x().R2()) - inp_->order;
        LU_int_t Nspline_outer_x = Nspline_x - Nspline_inner_x;
        
        // panel y basis
        LU_int_t Nspline_y = rad_panel_->bspline_y().Nspline();
        LU_int_t Nspline_inner_y = rad_panel_->bspline_y().knot(rad_inner_->bspline_y().R2()) - inp_->order;
        LU_int_t Nspline_outer_y = Nspline_y - Nspline_inner_y;
        
        // angular block
        int iang = iblock * ang_->states().size() + iblock;
        
        // convert inner region matrix block to COO matrix
        CooMatrix<LU_int_t,Complex> coo_block;
        if (cmd_->lightweight_full)
            coo_block = std::move(dynamic_cast<NoPreconditioner const*>(this)->calc_A_block(iblock, iblock).tocoo<LU_int_t>());
        else
            coo_block = std::move(A_blocks_[iang].tocoo<LU_int_t>());
        
        // add the A-block
        coo_block.resize
        (
            Nspline_inner_x * Nspline_inner_y + Nchan1 * Nspline_outer_x + Nchan2 * Nspline_outer_y,
            Nspline_inner_x * Nspline_inner_y + Nchan1 * Nspline_outer_x + Nchan2 * Nspline_outer_y
        );
        
        if (not inp_->inner_only)
        {
            // add the outer region C-blocks
            coo_block += Cu_blocks_[iang];
            coo_block += Cl_blocks_[iang];
            
            // add the B-blocks
            for (int m = 0; m < Nchan1; m++)
            for (int n = 0; n < Nchan1; n++)
            {
                if (cmd_->outofcore) const_cast<SymBandMatrix<Complex>&>(B1_blocks_[iang][m * Nchan1 + n]).hdfload();
                CooMatrix<LU_int_t,Complex> B_coo_small = B1_blocks_[iang][m * Nchan1 + n].tocoo<LU_int_t>();
                if (cmd_->outofcore) const_cast<SymBandMatrix<Complex>&>(B1_blocks_[iang][m * Nchan1 + n]).drop();
                CooMatrix<LU_int_t,Complex> B_coo_large
                (
                    coo_block.rows(), coo_block.cols(),
                    B_coo_small.i() + Nspline_inner_x * Nspline_inner_y + m * Nspline_outer_x,
                    B_coo_small.j() + Nspline_inner_x * Nspline_inner_y + n * Nspline_outer_x,
                    B_coo_small.v()
                );
                coo_block += B_coo_large;
            }
            for (int m = 0; m < Nchan2; m++)
            for (int n = 0; n < Nchan2; n++)
            {
                if (cmd_->outofcore) const_cast<SymBandMatrix<Complex>&>(B2_blocks_[iang][m * Nchan2 + n]).hdfload();
                CooMatrix<LU_int_t,Complex> B_coo_small = B2_blocks_[iang][m * Nchan2 + n].tocoo<LU_int_t>();
                if (cmd_->outofcore) const_cast<SymBandMatrix<Complex>&>(B2_blocks_[iang][m * Nchan2 + n]).drop();
                CooMatrix<LU_int_t,Complex> B_coo_large
                (
                    coo_block.rows(), coo_block.cols(),
                    B_coo_small.i() + Nspline_inner_x * Nspline_inner_y + Nchan1 * Nspline_outer_x + m * Nspline_outer_y,
                    B_coo_small.j() + Nspline_inner_x * Nspline_inner_y + Nchan1 * Nspline_outer_x + n * Nspline_outer_y,
                    B_coo_small.v()
                );
                coo_block += B_coo_large;
            }
        }
        
        // create the CSR block that will be factorized
        CsrMatrix<LU_int_t,Complex> csr = coo_block.tocsr();
        coo_block = CooMatrix<LU_int_t,Complex>();
        
        // set up factorization data
        data_[iblock].drop_tolerance = cmd_->droptol;
        data_[iblock].groupsize = cmd_->groupsize;
        data_[iblock].ooc_dir = cmd_->scratch.c_str();
#ifdef WITH_SUPERLU_DIST
        if (cmd_->factorizer == "superlu_dist")
        {
            data_[iblock].superlu_dist_grid = const_cast<gridinfo_t*>(&grid_);
        }
#endif
#ifdef WITH_MUMPS
        if (cmd_->factorizer == "mumps")
        {
            data_[iblock].out_of_core = cmd_->mumps_outofcore;
            data_[iblock].verbosity = cmd_->mumps_verbose;
    #ifdef WITH_MPI
        #ifdef _WIN32
            data_[iblock].fortran_comm = MPI_Comm_c2f((MPI_Fint)(std::intptr_t) par_->groupcomm());
        #else
            data_[iblock].fortran_comm = MPI_Comm_c2f((ompi_communicator_t*) par_->groupcomm());
        #endif
    #else
            data_[iblock].fortran_comm = 0;
    #endif
        }
#endif
        
        // factorize the block and store it
        lu_[iblock]->factorize(csr, data_[iblock]);
        
        // print time and memory info for this block (one thread at a time)
        # pragma omp critical
        {
            if (lu_[iblock]->cond() != 0)
            {
                std::cout << std::endl << std::setw(37) << format
                (
                    "\tLU #%d (%d,%d) in %d:%02d (%s, cond %1.0e)",
                    iblock, ang_->states()[iblock].first, ang_->states()[iblock].second,      // block identification (id, ℓ₁, ℓ₂)
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
                    iblock, ang_->states()[iblock].first, ang_->states()[iblock].second,      // block identification (id, ℓ₁, ℓ₂)
                    timer.seconds() / 60, timer.seconds() % 60,                             // factorization time
                    nice_size(lu_[iblock]->size()).c_str()                                  // final memory size
                );
            }
        }
        
        // save the diagonal block's factorization
        lu_[iblock]->link(format("lu-%d.bin", iblock));
        if (cmd_->outofcore)
            lu_[iblock]->save();
        
#ifdef _OPENMP
        // release lock
        if (not cmd_->parallel_factorization)
            omp_unset_lock(&lu_lock_);
#endif
    }
}

void ILUCGPreconditioner::CG_prec (int iblock, const cArrayView r, cArrayView z) const
{
#ifdef _OPENMP
    if (cmd_->factorizer == "mumps")
        omp_set_lock(&lu_lock_);
#endif
    
    // precondition by LU
    lu_[iblock]->solve(r, z, 1);
    
    if (cmd_->factorizer == "mumps")
    {
#ifdef _OPENMP
        omp_unset_lock(&lu_lock_);
#endif
        par_->bcast_g(par_->igroup(), 0, z.data(), z.size());
    }
}

void ILUCGPreconditioner::CG_exit (int iblock) const
{
#ifdef _OPENMP
    if (cmd_->factorizer == "mumps")
        omp_set_lock(&lu_lock_);
#endif
    
    // release memory
    if (cmd_->outofcore)
        lu_[iblock]->drop();
    
#ifdef _OPENMP
    if (cmd_->factorizer == "mumps")
        omp_unset_lock(&lu_lock_);
#endif
    
    // exit parent
    CGPreconditioner::CG_exit(iblock);
}

void ILUCGPreconditioner::finish ()
{
#ifdef _OPENMP
    if (cmd_->factorizer == "mumps")
        omp_set_lock(&lu_lock_);
#endif
    
    lu_.resize(0);
    
#ifdef _OPENMP
    if (cmd_->factorizer == "mumps")
        omp_unset_lock(&lu_lock_);
#endif
    
    CGPreconditioner::finish();
}

// --------------------------------------------------------------------------------- //

addClassToParentRunTimeSelectionTable(PreconditionerBase, ILUCGPreconditioner)

// --------------------------------------------------------------------------------- //
