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

#ifndef HEX_MATRIX
#define HEX_MATRIX

#include "hex-numbers.h"

// Forward declaration of sparse matrix classes in order to enable them as return types
// of other classes before proper definition.
template <class IdxT, class DataT> class CooMatrix;
template <class IdxT, class DataT> class CscMatrix;
template <class IdxT, class DataT> class CsrMatrix;
template <class DataT> class SymBandMatrix;

// Forward declaration of dense matrix classes in order to enable them as return types
// of other classes before proper definition. Note that RowMatrix and ColMatrix
// are implicitly declared as "deep", i.e. allocating their own storage.
// However, they can be declared using DenseMatrixView<T> to make them shallow
// (i.e. carry only pointer to the array of matrix elements).
template <class T> class DenseMatrixView;
template <class T> class DenseMatrix;
template <class T, class Base = DenseMatrix<T>> class RowMatrix;
template <class T, class Base = DenseMatrix<T>> class ColMatrix;

// Shorthand (a C++11 alias) for shallow dense matrices.
template <typename T> using RowMatrixView = RowMatrix<T,DenseMatrixView<T>>;
template <typename T> using ColMatrixView = ColMatrix<T,DenseMatrixView<T>>;

// some used Blas prototypes (Fortran convention!)
extern "C" void dgemv_ (char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
extern "C" void zgemv_ (char*, int*, int*, Complex*, Complex*, int*, Complex*, int*, Complex*, Complex*, int*);
extern "C" void dgemm_ (char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
extern "C" void zgemm_ (char*, char*, int*, int*, int*, Complex*, Complex*, int*, Complex*, int*, Complex*, Complex*, int*);

// some used Lapack prototypes (Fortran convention!)
extern "C" void zgetrf_ (int*, int*, Complex*, int*, int*, int*);
extern "C" void zgetri_ (int*, Complex*, int*, int*, Complex*, int*, int*);
extern "C" void zgeev_ (char*, char*, int*, Complex*, int*, Complex*, Complex*, int*, Complex*, int*, Complex*, int*, double*, int*);

// the classes are in separate files
#include "hex-densematrix.h"
#include "hex-cscmatrix.h"
#include "hex-csrmatrix.h"
#include "hex-coomatrix.h"
#include "hex-symbandmatrix.h"

#endif // HEX_MATRIX
