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
    HexException("The coupled preconditioner is broken in this version of the program!");
    
    // concatenate all matrix blocks
    NumberArray<LU_int_t> I, J;
    cArray A;
    // TODO
    
    // convert the blocks to CSR
    CsrMatrix<LU_int_t, Complex> csr;
    // TODO
    
    // set up factorization data
    LUftData data;
    data.drop_tolerance = cmd_->droptol;
#ifdef WITH_SUPERLU_DIST
    if (cmd_->factorizer == "superlu_dist")
    {
        // TODO : setup grid
    }
#endif
#ifdef WITH_MUMPS
    if (cmd_->factorizer == "mumps")
    {
        data.out_of_core = cmd_->mumps_outofcore;
        data.verbosity = cmd_->mumps_verbose;
    #ifdef WITH_MPI
        data.fortran_comm = MPI_Comm_c2f((ompi_communicator_t*) par_->groupcomm());
    #else
        data.fortran_comm = 0;
    #endif
    }
#endif
    
    // factorize
    lu_->factorize(csr, data);
}

void CoupledPreconditioner::precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const
{
    // some useful constants
    std::size_t Nang = r.size(), Nchunk = r[0].size();
    
    // convert block array to monolithic array
    for (unsigned ill = 0; ill < Nang; ill++)
    for (unsigned i = 0; i < Nchunk; i++)
        X[ill * Nchunk + i] = r[ill][i];
    
    // solve
    lu_->solve(X, X, 1);
    
    // copy solution to result
    for (unsigned ill = 0; ill < Nang; ill++)
    for (unsigned i = 0; i < Nchunk; i++)
        z[ill][i] = X[ill * Nchunk + i];
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
