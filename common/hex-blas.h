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

#ifndef HEX_BLAS_H
#define HEX_BLAS_H

#include "hex-arrays.h"
#include "hex-numbers.h"
#include "hex-matrix.h"

/// BLAS wrapper
namespace blas
{
    // integer type used by BLAS (32 or 64 bits)
#ifdef BLAS64
    typedef int64_t Int;
#else
    typedef int Int;
#endif
    
    void gemm
    (
        Complex a,
            DenseMatrixView<Complex> const & A,
            DenseMatrixView<Complex> const & B,
        Complex b,
            ColMatrixView<Complex> C
    );
    
    void gemm
    (
        Complex a,
            DenseMatrixView<Complex> const & A,
            DenseMatrixView<Complex> const & B,
        Complex b,
            RowMatrixView<Complex> C
    );
    
    void gemv
    (
        Real a,
            DenseMatrixView<Real> const & A,
            const rArrayView v,
        Real b,
            rArrayView w
    );
    
    void gemv
    (
        Complex a,
            DenseMatrixView<Complex> const & A,
            const cArrayView v,
        Complex b,
            cArrayView w
    );
    
    void xpby
    (
        cArrayView x,
        Complex b,
        const cArrayView y
    );
}

#endif
