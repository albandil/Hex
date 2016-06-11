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

#ifdef WITH_MUMPS

#include "preconditioners.h"

const std::string CoupledPreconditioner::prec_name = "coupled";
const std::string CoupledPreconditioner::prec_description = "Coupled solver that uses MUMPS out of core; nevertheless, it still consumes a huge amount of RAM.";

void CoupledPreconditioner::update (Real E)
{
    // concatenate all matrix blocks
    NumberArray<LU_int_t> I, J;
    cArray A;
    // TODO
    
    // convert the blocks to CSR
    CsrMatrix<LU_int_t, Complex> csr;
    // TODO
    
    // factorize
    lu_ = csr.factorize_mumps(0, (void*)std::intptr_t(cmd_.mumps_outofcore + 2 * cmd_.mumps_verbose));
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

#endif // WITH_MUMPS
