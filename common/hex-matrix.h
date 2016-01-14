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

#endif // HEX_MATRIX
