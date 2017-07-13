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
    if (not inp_->inner_only)
        HexException("The coupled preconditioner does not yet support channel reduction");
    
    if (par_->Nproc() != par_->groupsize())
        HexException("Coupled preconditioner must be executed with full MPI groupsize.");
    
    NoPreconditioner::update(E);
    
    // shorthands
    unsigned order = inp_->order;
    std::size_t Nxspline = rad_inner().bspline_x().Nspline();
    std::size_t Nyspline = rad_inner().bspline_y().Nspline();
    std::size_t Nang = ang_->states().size();
    std::size_t N = Nang * Nxspline * Nyspline;
    std::size_t NZ = Nang * Nang * Nxspline * (2*order + 1) * Nyspline * (2*order + 1);
    
    // concatenate matrix blocks as a COO matrix
    NumberArray<LU_int_t> I, J;
    cArray V;
    
    // allocate memory
    X.resize(N);
    Y.resize(N);
    if (par_->IamMaster())
    {
        I.reserve(NZ);
        J.reserve(NZ);
        V.reserve(NZ);
    }
    
    SymBandMatrix<Complex> subblock (Nyspline, order + 1);
    
    for (unsigned ill = 0; ill < Nang; ill++) if (par_->IamMaster()) /*if (par_->isMyWork(ill))*/
    for (unsigned illp = 0; illp < Nang; illp++)
    {
        int l1 = ang_->states()[ill].first;
        int l2 = ang_->states()[ill].second;
        
        for (std::size_t i = 0; i < Nxspline; i++)
        for (std::size_t k = (i > order ? i - order : 0); k < Nxspline and k <= i + order; k++)
        {
            subblock.data().fill(0);
            
            // one-electron part
            if (ill == illp)
            {
                subblock.populate
                (
                    [&](int j, int l)
                    {
                        Complex Sx = rad_inner().S_x()(i,k), Dx = rad_inner().D_x()(i,k), Mm1x = rad_inner().Mm1_x()(i,k), Mm2x = rad_inner().Mm2_x()(i,k);
                        Complex Sy = rad_inner().S_y()(j,l), Dy = rad_inner().D_y()(j,l), Mm1y = rad_inner().Mm1_y()(j,l), Mm2y = rad_inner().Mm2_y()(j,l);
                        
                        return E_ * Sx * Sy - 0.5_r * Dx * Sy + inp_->Za *            Mm1x * Sy - 0.5_r * l1 * (l1 + 1) * Mm2x * Sy
                                            - 0.5_r * Dy * Sx - inp_->Za * inp_->Zp * Mm1y * Sx - 0.5_r * l2 * (l2 + 1) * Mm2y * Sx;
                    }
                );
            }
            
            // two electron part
            for (int lambda = 0; lambda <= std::min(rad_inner().maxlambda(), cmd_->coupling_limit); lambda++)
            {
                double f = ang_->f(ill, illp, lambda);
                if (f != 0)
                    subblock += Complex(inp_->Zp * f) * rad_inner().calc_R_tr_dia_block(lambda, i, k);
            }
            
            // convert block to COO matrix
            CooMatrix<LU_int_t, Complex> coo = subblock.tocoo<LU_int_t>();
            
            // block index shift
            LU_int_t ishift = (ill  * Nxspline + i) * Nyspline;
            LU_int_t jshift = (illp * Nxspline + k) * Nyspline;
            
            // update matrix
            I.append(coo.i() + ishift);
            J.append(coo.j() + jshift);
            V.append(coo.v());
        }
    }
    
    // construct COO matrix
    CooMatrix<LU_int_t, Complex> coo (N, N, I, J, V);
    
    I.drop();
    J.drop();
    V.drop();
    
    // convert the blocks to CSR
    CsrMatrix<LU_int_t, Complex> csr = coo.tocsr();
    
    // set up factorization data
    LUftData data;
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
        data.centralized_matrix = true;
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
    Timer timer;
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

void CoupledPreconditioner::precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const
{
    // some useful constants
    std::size_t Nang = r.size(), Nchunk = r[0].size();
    
    // convert block array to monolithic array
    for (std::size_t ill = 0; ill < Nang; ill++)
        cArrayView(X, ill * Nchunk, Nchunk) = r[ill];
    
    // solve
    lu_->solve(X, Y, 1);
    
    // broadcast solution from master to everyone
    par_->bcast(0, Y);
    
    // copy solution to result
    for (std::size_t ill = 0; ill < Nang; ill++)
        z[ill] = cArrayView(Y, ill * Nchunk, Nchunk);
}

void CoupledPreconditioner::finish ()
{
    // delete the factorization object
    lu_.reset();
    
    // finish parent class
    NoPreconditioner::finish();
}

// --------------------------------------------------------------------------------- //

addClassToParentRunTimeSelectionTable(PreconditionerBase, CoupledPreconditioner)

// --------------------------------------------------------------------------------- //
