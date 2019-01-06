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

extern "C" void sgemv_
(
    char*, blas::Int*, blas::Int*, float*, float*,
    blas::Int*, float*, blas::Int*, float*, float*, blas::Int*
);
extern "C" void dgemv_
(
    char*, blas::Int*, blas::Int*, double*, double*,
    blas::Int*, double*, blas::Int*, double*, double*, blas::Int*
);
extern "C" void cgemv_
(
    char*, blas::Int*, blas::Int*, std::complex<float>*, std::complex<float>*,
    blas::Int*, std::complex<float>*, blas::Int*, std::complex<float>*,
    std::complex<float>*, blas::Int*
);
extern "C" void zgemv_
(
    char*, blas::Int*, blas::Int*, std::complex<double>*, std::complex<double>*,
    blas::Int*, std::complex<double>*, blas::Int*, std::complex<double>*,
    std::complex<double>*, blas::Int*
);
extern "C" void sgemm_
(
    char*, char*, blas::Int*, blas::Int*, blas::Int*, float*, float*,
    blas::Int*, float*, blas::Int*, float*, float*, blas::Int*
);
extern "C" void dgemm_
(
    char*, char*, blas::Int*, blas::Int*, blas::Int*, double*, double*,
    blas::Int*, double*, blas::Int*, double*, double*, blas::Int*
);
extern "C" void cgemm_
(
    char*, char*, blas::Int*, blas::Int*, blas::Int*, std::complex<float>*, std::complex<float>*,
    blas::Int*, std::complex<float>*, blas::Int*, std::complex<float>*, std::complex<float>*, blas::Int*
);
extern "C" void zgemm_
(
    char*, char*, blas::Int*, blas::Int*, blas::Int*, std::complex<double>*, std::complex<double>*, blas::Int*,
    std::complex<double>*, blas::Int*, std::complex<double>*, std::complex<double>*, blas::Int*
);
extern "C" void sgbtrf_
(
    blas::Int*, blas::Int*, blas::Int*, blas::Int*, float*,
    blas::Int*, blas::Int*, blas::Int*
);
extern "C" void dgbtrf_
(
    blas::Int*, blas::Int*, blas::Int*, blas::Int*, double*,
    blas::Int*, blas::Int*, blas::Int*
);
extern "C" void cgbtrf_
(
    blas::Int*, blas::Int*, blas::Int*, blas::Int*, std::complex<float>*,
    blas::Int*, blas::Int*, blas::Int*
);
extern "C" void zgbtrf_
(
    blas::Int*, blas::Int*, blas::Int*, blas::Int*, std::complex<double>*,
    blas::Int*, blas::Int*, blas::Int*
);
extern "C" void sgbtrs_
(
    char*, blas::Int*, blas::Int*, blas::Int*, blas::Int*, float*,
    blas::Int*, blas::Int*, float*, blas::Int*, blas::Int*
);
extern "C" void dgbtrs_
(
    char*, blas::Int*, blas::Int*, blas::Int*, blas::Int*, double*,
    blas::Int*, blas::Int*, double*, blas::Int*, blas::Int*
);
extern "C" void cgbtrs_
(
    char*, blas::Int*, blas::Int*, blas::Int*, blas::Int*, std::complex<float>*,
    blas::Int*, blas::Int*, std::complex<float>*, blas::Int*, blas::Int*
);
extern "C" void zgbtrs_
(
    char*, blas::Int*, blas::Int*, blas::Int*, blas::Int*, std::complex<double>*,
    blas::Int*, blas::Int*, std::complex<double>*, blas::Int*, blas::Int*
);

// ---------------------------------------------------------------------------------- //

// some used Lapack prototypes (Fortran convention!)
extern "C" void cgetrf_ (int*, int*, std::complex<float>*, int*, int*, int*);
extern "C" void zgetrf_ (int*, int*, std::complex<double>*, int*, int*, int*);
extern "C" void cgetri_ (int*, std::complex<float>*, int*, int*, std::complex<float>*, int*, int*);
extern "C" void zgetri_ (int*, std::complex<double>*, int*, int*, std::complex<double>*, int*, int*);
extern "C" void cgeev_ (char*, char*, int*, std::complex<float>*, int*, std::complex<float>*, std::complex<float>*, int*, std::complex<float>*, int*, std::complex<float>*, int*, float*, int*);
extern "C" void zgeev_ (char*, char*, int*, std::complex<double>*, int*, std::complex<double>*, std::complex<double>*, int*, std::complex<double>*, int*, std::complex<double>*, int*, double*, int*);

// ---------------------------------------------------------------------------------- //

void blas::gemm (Complex a, DenseMatrixView<Complex> const & A, DenseMatrixView<Complex> const & B, Complex b, ColMatrixView<Complex> C)
{
#ifndef NDEBUG
    // verify compatibility of dimensions
    if (A.cols() != B.rows() or A.rows() != C.rows() or B.cols() != C.cols())
        HexException("Incompatible dimensions: [%d x %d] [%d x %d] != [%d x %d].", A.rows(), A.cols(), B.rows(), B.cols(), C.rows(), C.cols());
#endif

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

    // call the dedicated BLAS routine
#ifdef SINGLE
    cgemm_(&layout_A, &layout_B, &m, &n, &k, &a, pA, &ld_A, pB, &ld_B, &b, pC, &ld_C);
#else
    zgemm_(&layout_A, &layout_B, &m, &n, &k, &a, pA, &ld_A, pB, &ld_B, &b, pC, &ld_C);
#endif
}

void blas::gemm (Complex a, DenseMatrixView<Complex> const & A, DenseMatrixView<Complex> const & B, Complex b, RowMatrixView<Complex> C)
{
#ifndef NDEBUG
    // verify compatibility of dimensions
    if (A.cols() != B.rows() or A.rows() != C.rows() or B.cols() != C.cols())
        HexException("Incompatible dimensions: [%d x %d] [%d x %d] != [%d x %d].", A.rows(), A.cols(), B.rows(), B.cols(), C.rows(), C.cols());
#endif

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
#ifdef SINGLE
    cgemm_(&layout_B, &layout_A, &m, &n, &k, &a, pB, &ld_B, pA, &ld_A, &b, pC, &ld_C);
#else
    zgemm_(&layout_B, &layout_A, &m, &n, &k, &a, pB, &ld_B, pA, &ld_A, &b, pC, &ld_C);
#endif
}

void blas::gemv (Real a, DenseMatrixView<Real> const & A, const rArrayView v, Real b, rArrayView w)
{
#ifndef NDEBUG
    // verify compatibility of dimensions
    if (A.cols() != (int)v.size() or A.rows() != (int)w.size())
        HexException("Incompatible dimensions: [%d x %d] [%d x 1] != [%d x 1].", A.rows(), A.cols(), v.size(), w.size());
#endif

    // get matrix data
    char layout_A = A.layout();
    Int ld_A = A.ld();
    Int m = A.rows();
    Int n = A.cols();
    Real * pA = const_cast<Real*>(A.data().data());

    // get array data
    Int incv = 1, incw = 1;
    Real * pv = const_cast<Real*>(v.data());
    Real * pw = w.data();

    // call the dedicated BLAS routine
#ifdef SINGLE
    sgemv_(&layout_A, &m, &n, &a, pA, &ld_A, pv, &incv, &b, pw, &incw);
#else
    dgemv_(&layout_A, &m, &n, &a, pA, &ld_A, pv, &incv, &b, pw, &incw);
#endif
}

void blas::gemv (Complex a, DenseMatrixView<Complex> const & A, const cArrayView v, Complex b, cArrayView w)
{
#ifndef NDEBUG
    // verify compatibility of dimensions
    if (A.cols() != (int)v.size() or A.rows() != (int)w.size())
        HexException("Incompatible dimensions: [%d x %d] [%d x 1] != [%d x 1].", A.rows(), A.cols(), v.size(), w.size());
#endif

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
#ifdef SINGLE
    cgemv_(&layout_A, &m, &n, &a, pA, &ld_A, pv, &incv, &b, pw, &incw);
#else
    zgemv_(&layout_A, &m, &n, &a, pA, &ld_A, pv, &incv, &b, pw, &incw);
#endif
}

void blas::xpby (cArrayView x, Complex b, const cArrayView y)
{
#ifndef NDEBUG
    // verify compatibility of dimensions
    if (x.size() != y.size())
        HexException("Incompatible dimensions: %d != %d", x.size(), y.size());
#endif

    // get length of the arrays
    std::size_t N = x.size();

    // get pointers to the data
    Complex       * const restrict px = x.data();
    Complex const * const restrict py = y.data();

    // run the vectorized loop
    # pragma omp simd
    for (std::size_t i = 0; i < N; i++)
        px[i] += b * py[i];
}

blas::Int blas::sbtrf (blas::Int n, blas::Int k, rArrayView ab, ArrayView<blas::Int> ipiv)
{
    blas::Int info;
    blas::Int ldab = 3*k + 1;

#ifdef SINGLE
    sgbtrf_(&n, &n, &k, &k, &ab[0], &ldab, &ipiv[0], &info);
#else
    dgbtrf_(&n, &n, &k, &k, &ab[0], &ldab, &ipiv[0], &info);
#endif

    return info;
}


blas::Int blas::sbtrf (blas::Int n, blas::Int k, cArrayView ab, ArrayView<blas::Int> ipiv)
{
    blas::Int info;
    blas::Int ldab = 3*k + 1;

#ifdef SINGLE
    cgbtrf_(&n, &n, &k, &k, &ab[0], &ldab, &ipiv[0], &info);
#else
    zgbtrf_(&n, &n, &k, &k, &ab[0], &ldab, &ipiv[0], &info);
#endif

    return info;
}

blas::Int blas::sbtrs (blas::Int n, blas::Int k, rArrayView ab, ArrayView<blas::Int> ipiv, rArrayView bx)
{
    blas::Int info;
    blas::Int ldab = 3*k + 1;
    blas::Int nrhs = bx.size() / n;
    char trans = 'N';

#ifdef SINGLE
    sgbtrs_(&trans, &n, &k, &k, &nrhs, &ab[0], &ldab, &ipiv[0], &bx[0], &n, &info);
#else
    dgbtrs_(&trans, &n, &k, &k, &nrhs, &ab[0], &ldab, &ipiv[0], &bx[0], &n, &info);
#endif

    return info;
}

blas::Int blas::sbtrs (blas::Int n, blas::Int k, cArrayView ab, ArrayView<blas::Int> ipiv, cArrayView bx)
{
    blas::Int info;
    blas::Int ldab = 3*k + 1;
    blas::Int nrhs = bx.size() / n;
    char trans = 'N';

#ifdef SINGLE
    cgbtrs_(&trans, &n, &k, &k, &nrhs, &ab[0], &ldab, &ipiv[0], &bx[0], &n, &info);
#else
    zgbtrs_(&trans, &n, &k, &k, &nrhs, &ab[0], &ldab, &ipiv[0], &bx[0], &n, &info);
#endif

    return info;
}
