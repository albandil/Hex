//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2017, Jakub Benda, Charles University in Prague                    //
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

#include "hex-csrmatrix.h"

// --------------------------------------------------------------------------------- //

#include "CoupledPreconditioner.h"

// --------------------------------------------------------------------------------- //

std::string CoupledPreconditioner::description () const
{
    return "Coupled preconditioner that uses LU decomposition "
    "of a matrix composed of all diagonal and selected off-diagonal blocks. This is just a testing feature "
    "and is likely to severely exceed your RAM.";
}

void CoupledPreconditioner::update (Real E)
{
    if (par_->Nproc() != par_->groupsize())
        HexException("Coupled preconditioner must be executed with full MPI groupsize.");
    
    KPACGPreconditioner::update(E);
    
    // do not update if running in energy-perturbation mode
    if (lu_ != nullptr and lu_->valid() and cmd_->noluupdate)
        return;
    
    // shorthands
    LU_int_t order = inp_->order;
    LU_int_t Nang = ang_->states().size();
    
    // B-spline bases
    Bspline const & bspline_full  = rad_full_ ->bspline();
    Bspline const & bspline_inner = rad_inner_->bspline();
    Bspline const & bspline_x     = rad_panel_->bspline_x();
    Bspline const & bspline_y     = rad_panel_->bspline_y();
    
    // global basis
    int Nspline_full  = bspline_full.Nspline();
    int Nspline_inner = bspline_inner.Nspline();
    int Nspline_outer = Nspline_full - Nspline_inner;
    
    // panel x basis
    int Nspline_x_full  = bspline_x.Nspline();
    int Nspline_x_inner = bspline_x.hash() == bspline_full.hash() ? bspline_inner.Nspline() : bspline_x.Nspline();
    int Nspline_x_outer = Nspline_x_full - Nspline_x_inner;
    
    // panel y basis
    int Nspline_y_full  = bspline_y.Nspline();
    int Nspline_y_inner = bspline_y.hash() == bspline_full.hash() ? bspline_inner.Nspline() : bspline_y.Nspline();
    int Nspline_y_outer = Nspline_y_full - Nspline_y_inner;
    
    // concatenate matrix blocks as a COO matrix
    NumberArray<LU_int_t> I, J;
    cArray V;
    
    // calculate full matrix rank
    std::size_t N = 0;
    for (std::size_t n : block_rank_) N += n;
    
    // calculate approximate number of non-zero elements and pre-allocate memory
    std::size_t NZ = 0;
    for (LU_int_t ill  = 0; ill  < Nang; ill ++)
    for (LU_int_t illp = 0; illp < Nang; illp++)
    {
        // number of asymptotic bound channels
        int Nchan1 = Nchan_[ill].first;     // # r1 -> inf, l2 bound
        int Nchan2 = Nchan_[ill].second;    // # r2 -> inf, l1 bound
        int Nchan1p = Nchan_[illp].first;   // # r1 -> inf, l2p bound
        int Nchan2p = Nchan_[illp].second;  // # r2 -> inf, l1p bound
        
        // update non-zero count (NOTE: Needs checking.)
        NZ += Nspline_x_inner * Nspline_y_inner * (2*order + 1) * (2*order + 1);    // A-block
        NZ += Nchan1 * Nchan1p * Nspline_x_outer * (2*order + 1);                   // B1-blocks
        NZ += Nchan2 * Nchan2p * Nspline_y_outer * (2*order + 1);                   // B2-blocks
        NZ += Nchan1 * (2*order + 1) * (2*order + 1) * Nspline_outer;               // Cu-blocks
        NZ += Nchan2 * (2*order + 1) * (2*order + 1) * Nspline_outer;               // Cu-blocks
    }
    I.reserve(NZ);
    J.reserve(NZ);
    V.reserve(NZ);
    
    std::vector<std::pair<int,int>> segregated;
    std::cout << "\tAssemble full matrix of the system ... " << std::flush;
    Timer timer;
    
    for (LU_int_t ill  = 0, offset  = 0; ill  < Nang; offset  += block_rank_[ill ++])
    for (LU_int_t illp = 0, offsetp = 0; illp < Nang; offsetp += block_rank_[illp++])
    if ((ill * Nang + illp) % par_->groupsize() == par_->igroupproc())
    {
        int l1 = ang_->states()[ill].first;
        int l2 = ang_->states()[ill].second;
        
        // when this angular state is only weakly coupled, do not include it in the matrix
        if (not cmd_->couple_all and Xp_[0][l1].empty() and Xp_[1][l2].empty())
        {
            // ignore angular coupling blocks
            if (ill == illp)
            {
                // this block will be handled by segregated preconditioner
                segregated.push_back(std::make_pair(l1,l2));
                
                // put identity on diagonal
                I.append(identity<LU_int_t>(block_rank_[ill]) + offset);
                J.append(identity<LU_int_t>(block_rank_[ill]) + offsetp);
                V.append(cArray(block_rank_[ill], 1.0_z));
            }
        }
        else
        {
            // convert matrix block to COO format
            CooMatrix<LU_int_t,Complex> coo = calc_full_block(ill, illp);
            
            // copy elements to the assembled storage
            I.append(coo.i() + offset);
            J.append(coo.j() + offsetp);
            V.append(coo.v());
        }
    }
    
    std::cout << "done after " << timer.nice_time() << std::endl;
    std::cout << "\tUncoupled blocks:";
    if (segregated.empty())
        std::cout << " none";
    else for (std::pair<int,int> const & p : segregated)
        std::cout << " (" << p.first << "," << p.second << ")";
    std::cout << std::endl;
    
    // convert the blocks to CSR
    CsrMatrix<LU_int_t, Complex> csr = CooMatrix<LU_int_t, Complex>
    (
        N, N,
        std::move(I), std::move(J), std::move(V)
    ).tocsr();
    
    // set up factorization data
    LUftData data = defaultLUftData;
    data.drop_tolerance = cmd_->droptol;
#ifdef WITH_SUPERLU_DIST
    if (cmd_->factorizer == "superlu_dist")
    {
        HexException("Coupled preconditioner cannot be used with SuperLU_DIST yet - try MUMPS.");
    }
#endif
#ifdef WITH_MUMPS
    if (cmd_->factorizer == "mumps")
    {
        data.out_of_core = cmd_->mumps_outofcore;
        data.verbosity = cmd_->mumps_verbose;
        data.centralized_matrix = false;
    #ifdef WITH_MPI
        #ifdef _WIN32
        data.fortran_comm = MPI_Comm_c2f((MPI_Fint)(std::intptr_t) par_->groupcomm());
        #else
        data.fortran_comm = MPI_Comm_c2f((ompi_communicator_t*) par_->groupcomm());
        #endif
    #else
        data.fortran_comm = 0;
    #endif
    }
#endif
    
    // factorize
    timer.reset();
    lu_.reset(LUft::Choose(cmd_->factorizer));
    lu_->factorize(csr, data);
    
    // print statistics
    if (lu_->cond() != 0)
    {
        std::cout << format
        (
            "\tCoupled preconditioner: LU in %d:%02d (%s, cond %1.0e)",
            timer.seconds() / 60, timer.seconds() % 60,   // factorization time
            nice_size(lu_->size()).c_str(),               // final memory size
            lu_->cond()                                   // estimation of the condition number
        ) << std::endl;
    }
    else
    {
        std::cout << format
        (
            "\tCoupled preconditioner: LU in %d:%02d (%s)",
            timer.seconds() / 60, timer.seconds() % 60,   // factorization time
            nice_size(lu_->size()).c_str()                // final memory size
        ) << std::endl;
    }
}

int CoupledPreconditioner::solve_block (int ill, const cArrayView r, cArrayView z) const
{
    int l1 = ang_->states()[ill].first;
    int l2 = ang_->states()[ill].second;
    
    if (not cmd_->couple_all and Xp_[0][l1].empty() and Xp_[1][l2].empty())
    {
        // uncoupled block - solve using KPA
        return CGPreconditioner::solve_block(ill, r, z);
    }
    else
    {
        // coupled block - already solved
        return 1;
    }
}

void CoupledPreconditioner::precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const
{
    // number of angular blocks
    std::size_t Nang = r.size();
    
    // number of initial states (right-hand sides) solved at once
    std::size_t Nini = r[0].size() / block_rank_[0];
    
    // in parallel, only the master process and its fellows from the MPI group have the complete segment r[0]
    // so to get correct number of right-hand sides we need to let the process broadcast the number;
    par_->bcast(0, &Nini, 1);
    
    // total solution size (all blocks, all initial states)
    std::size_t N = 0;
    for (cArray const & a : r) N += a.size();
    X.resize(N);
    Y.resize(N);
    
    // convert block array to monolithic array
    for (std::size_t ini = 0, Xoffset = 0; ini < Nini; ini++)
    for (std::size_t ill = 0; ill < Nang; ill++)
    {
        std::size_t chunk = r[ill].size() / Nini;
        std::size_t roffset = chunk * ini;
        
        cArrayView(X, Xoffset, chunk) = cArrayView(r[ill], roffset, chunk);
        
        Xoffset += r[ill].size() / Nini;
    }
    
    // solve
    lu_->solve(X, Y, Nini);
    
    // broadcast solution from master to everyone
    par_->bcast(0, Y);
    
    // copy solution to result
    for (std::size_t ini = 0, Yoffset = 0; ini < Nini; ini++)
    for (std::size_t ill = 0; ill < Nang; ill++)
    {
        std::size_t chunk = z[ill].size() / Nini;
        std::size_t zoffset = chunk * ini;
        
        cArrayView(z[ill], zoffset, chunk) = cArrayView(Y, Yoffset, chunk);
        
        Yoffset += z[ill].size() / Nini;
    }
    
    // use segregated KPA preconditioner to solve uncoupled blocks
    KPACGPreconditioner::precondition(r, z);
}

void CoupledPreconditioner::finish ()
{
    // delete the factorization object
    lu_.reset();
    
    // finish parent class
    KPACGPreconditioner::finish();
}

// --------------------------------------------------------------------------------- //

addClassToParentRunTimeSelectionTable(PreconditionerBase, CoupledPreconditioner)

// --------------------------------------------------------------------------------- //
