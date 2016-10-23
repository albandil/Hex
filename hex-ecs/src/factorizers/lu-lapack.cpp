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

// --------------------------------------------------------------------------------- //

#include "hex-blas.h"
#include "lu-lapack.h"

// --------------------------------------------------------------------------------- //

template<>
LUft_LAPACK<LU_int_t,Complex>::LUft_LAPACK ()
    : LUft<LU_int_t,Complex>()
{
    
}

template<>
void LUft_LAPACK<LU_int_t,Complex>::drop ()
{
    LU_.drop();
    ipiv_.drop();
}

template<>
LUft_LAPACK<LU_int_t,Complex>::~LUft_LAPACK ()
{
    
}

template<>
void LUft_LAPACK<LU_int_t,Complex>::factorize (CsrMatrix<LU_int_t,Complex> const & matrix, LUftData data)
{
    // determine bandwidth
    k_ = 0;
    for (blas::Int irow = 0; irow < (blas::Int)matrix.p().size(); irow++)
    {
        for (blas::Int idx = matrix.p()[irow]; idx < (blas::Int)matrix.p()[irow + 1]; idx++)
        {
            blas::Int icol = matrix.i()[idx];
            if (matrix.x()[idx] != 0.0_z)
                k_ = std::max(k_, std::abs(irow - icol));
        }
    }
    
    // allocate pivot array and the working matrix
    n_ = matrix.rows();
    ipiv_.resize(n_);
    LU_.resize((3*k_ + 1) * n_);
    LU_.fill(0);
    
    // populate the input matrix
    ColMatrixView<Complex> LU (3*k_ + 1, n_, LU_);
    for (blas::Int irow = 0; irow < (blas::Int)matrix.p().size(); irow++)
    {
        for (blas::Int idx = matrix.p()[irow]; idx < (blas::Int)matrix.p()[irow + 1]; idx++)
        {
            blas::Int icol = matrix.i()[idx];
            if (matrix.x()[idx] != 0.0_z)
                LU(2*k_ + irow - icol, icol) = matrix.x()[idx];
        }
    }
    
    // factorize the symmetric banded matrix
    blas::Int info = blas::sbtrf(n_, k_, LU_, ipiv_);
    
    // check success indicator
    if (info < 0)
        HexException("LAPACK: The argument %d has illegal value.", -info);
    if (info > 0)
        HexException("LAPACK: Matrix is singular.");
}

template<>
void LUft_LAPACK<LU_int_t,Complex>::solve (const cArrayView b, cArrayView x, int eqs) const
{
    // solve the system
    blas::Int info = blas::sbtrs(n_, k_, LU_, ipiv_, b, x);
    
    // check success indicator
    if (info < 0)
        HexException("LAPACK: The argument %d has illegal value.", -info);
}

template<>
std::size_t LUft_LAPACK<LU_int_t,Complex>::size () const
{
    return LU_.size() * sizeof(Complex) + ipiv_.size() * sizeof(blas::Int);
}

template<>
void LUft_LAPACK<LU_int_t,Complex>::save (std::string name) const
{ 
    HexException("LAPACK factorizer does not yet support --out-of-core option.");
}

template<>
void LUft_LAPACK<LU_int_t,Complex>::load (std::string name, bool throw_on_io_failure)
{
    if (throw_on_io_failure)
        HexException("LAPACK factorizer does not yet support --out-of-core option.");
}

// --------------------------------------------------------------------------------- //

addFactorizerToRuntimeSelectionTable(LAPACK, LU_int_t, Complex)

// --------------------------------------------------------------------------------- //
