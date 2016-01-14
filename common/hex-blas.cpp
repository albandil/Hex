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

#include "hex-blas.h"
#include "hex-densematrix.h"

// ---------------------------------------------------------------------------------- //

// some used Blas prototypes (Fortran convention!)
extern "C" void dgemv_ (char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
extern "C" void zgemv_ (char*, int*, int*, Complex*, Complex*, int*, Complex*, int*, Complex*, Complex*, int*);
extern "C" void dgemm_ (char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
extern "C" void zgemm_ (char*, char*, int*, int*, int*, Complex*, Complex*, int*, Complex*, int*, Complex*, Complex*, int*);

// ---------------------------------------------------------------------------------- //

// some used Lapack prototypes (Fortran convention!)
extern "C" void zgetrf_ (int*, int*, Complex*, int*, int*, int*);
extern "C" void zgetri_ (int*, Complex*, int*, int*, Complex*, int*, int*);
extern "C" void zgeev_ (char*, char*, int*, Complex*, int*, Complex*, Complex*, int*, Complex*, int*, Complex*, int*, double*, int*);

// ---------------------------------------------------------------------------------- //

void blas::gemm (Complex a, DenseMatrixView<Complex> const & A, DenseMatrixView<Complex> const & B, Complex b, ColMatrixView<Complex> C)
{
    // get matrix layouts
    char layout_A = A.layout(), layout_B = B.layout();
    
    // get matrix leading dimensions
    Int ld_A = A.ld(), ld_B = B.ld(), ld_C = C.ld();
    
    // get matrix dimensions
    Int m = A.rows(), n = B.cols(), k = A.cols();
    
    // get data pointers
    Complex * pA = const_cast<Complex*>(A.data().data());
    Complex * pB = const_cast<Complex*>(B.data().data());
    Complex * pC = C.data().data();
    
    /*std::cout << "ZGEMM" << std::endl;
    std::cout << "layout_A = " << layout_A << std::endl;
    std::cout << "layout_B = " << layout_B << std::endl;
    std::cout << "ld_A = " << ld_A << std::endl;
    std::cout << "ld_B = " << ld_B << std::endl;
    std::cout << "ld_C = " << ld_C << std::endl;
    std::cout << "m = " << m << std::endl;
    std::cout << "n = " << n << std::endl;
    std::cout << "k = " << k << std::endl;
    std::cout << "pA = " << pA << std::endl;
    std::cout << "pB = " << pB << std::endl;
    std::cout << "pC = " << pC << std::endl;*/
    
    // call the dedicated BLAS routine
    zgemm_(&layout_A, &layout_B, &m, &n, &k, &a, pA, &ld_A, pB, &ld_B, &b, pC, &ld_C);
}

void blas::gemm (Complex a, DenseMatrixView<Complex> const & A, DenseMatrixView<Complex> const & B, Complex b, RowMatrixView<Complex> C)
{
    // NOTE: To obtain row-major C, we need to multiply transpose: C' = B' A'
    
    // get transposed matrix layouts
    char layout_A = (A.layout() == 'N' ? 'T' : 'N');
    char layout_B = (B.layout() == 'N' ? 'T' : 'N');
    
    // get matrix leading dimensions
    Int ld_A = A.ld(), ld_B = B.ld(), ld_C = C.ld();
    
    // get matrix dimensions
    Int m = B.cols(), n = A.rows(), k = B.rows();
    
    // get data pointers
    Complex * pA = const_cast<Complex*>(A.data().data());
    Complex * pB = const_cast<Complex*>(B.data().data());
    Complex * pC = C.data().data();
    
    // call the dedicated BLAS routine
    zgemm_(&layout_A, &layout_B, &m, &n, &k, &a, pB, &ld_B, pA, &ld_A, &b, pC, &ld_C);
}

void blas::gemv (double a, DenseMatrixView<double> const & A, const rArrayView v, double b, rArrayView w)
{
    // get matrix data
    char layout_A = A.layout();
    Int ld_A = A.ld();
    Int m = A.rows();
    Int n = A.cols();
    double * pA = const_cast<double*>(A.data().data());
    
    // get array data
    Int incv = 1, incw = 1;
    double * pv = const_cast<double*>(v.data());
    double * pw = w.data();
    
    // call the dedicated BLAS routine
    dgemv_(&layout_A, &m, &n, &a, pA, &ld_A, pv, &incv, &b, pw, &incw);
}

void blas::gemv (Complex a, DenseMatrixView<Complex> const & A, const cArrayView v, Complex b, cArrayView w)
{
    // get matrix data
    char layout_A = A.layout();
    Int ld_A = A.ld();
    Int m = A.rows();
    Int n = A.cols();
    Complex * pA = const_cast<Complex*>(A.data().data());
    
    // get array data
    Int incv = 1, incw = 1;
    Complex * pv = const_cast<Complex*>(v.data());
    Complex * pw = w.data();
    
    // call the dedicated BLAS routine
    zgemv_(&layout_A, &m, &n, &a, pA, &ld_A, pv, &incv, &b, pw, &incw);
}
