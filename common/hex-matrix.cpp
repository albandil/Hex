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

#include <algorithm>
#include <complex>
#include <cstdio>
#include <cstring>
#include <memory>
#include <set>
#include <vector>

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"
#include "hex-csrmatrix.h"
#include "hex-densematrix.h"
#include "hex-hdffile.h"
#include "hex-matrix.h"
#include "hex-misc.h"
#include "hex-symbandmatrix.h"

// --------------------------------------------------------------------------------- //

#ifdef SINGLE
extern "C" void cgetrf_ (int*, int*, Complex*, int*, int*, int*);
extern "C" void cgetri_ (int*, Complex*, int*, int*, Complex*, int*, int*);
extern "C" void cgetrs_ (char*, int*, int*, Complex*, int*, int*, Complex*, int*, int*);
extern "C" void cgeev_ (char*, char*, int*, Complex*, int*, Complex*, Complex*, int*, Complex*, int*, Complex*, int*, Real*, int*);
extern "C" void cggev_ (char*, char*, int*, Complex*, int*, Complex*, int*, Complex*, Complex*, Complex*, int*, Complex*, int*, Complex*, int*, Real*, int*);
#else
extern "C" void zgetrf_ (int*, int*, Complex*, int*, int*, int*);
extern "C" void zgetrs_ (char*, int*, int*, Complex*, int*, int*, Complex*, int*, int*);
extern "C" void zgetri_ (int*, Complex*, int*, int*, Complex*, int*, int*);
extern "C" void zgeev_ (char*, char*, int*, Complex*, int*, Complex*, Complex*, int*, Complex*, int*, Complex*, int*, Real*, int*);
extern "C" void zggev_ (char*, char*, int*, Complex*, int*, Complex*, int*, Complex*, Complex*, Complex*, int*, Complex*, int*, Complex*, int*, Real*, int*);
#endif

// -------------------------------------------------------------------------------------
// Dense matrix routines
//

template<>
void ColMatrix<Complex>::invert (ColMatrix<Complex> & inv) const
{
#ifndef NO_LAPACK
    if (rows() != cols())
        HexException("Only square matrices can be diagonalized.");

    int N = rows();
    iArray IPIV(N);
    int LWORK = -1;
    cArray WORK(1);
    int INFO;

    // create copy of current matrix (will be overwritten)
    inv = *this;

    // compute the LU factorization
#ifdef SINGLE
    cgetrf_
#else
    zgetrf_
#endif
    (
        &N,             // number of rows in the matrix
        &N,             // number of columns of the matrix
        inv.begin(),    // matrix data (on return the LU factors)
        &N,             // leading dimension of A
        &IPIV[0],       // pivot indices
        &INFO           // diagnostic information
    );

    // query for optimal work size for inversion
#ifdef SINGLE
    cgetri_
#else
    zgetri_
#endif
    (
        &N,             // order of the matrix
        inv.begin(),    // LU-factors (on return the inverse matrix)
        &N,             // leading dimension of A
        &IPIV[0],       // pivot indices
        &WORK[0],       // work array
        &LWORK,         // work array length
        &INFO           // diagnostic information
    );

    // resize work array
    LWORK = WORK.front().real();
    WORK.resize(LWORK);

    // compute the inverse
#ifdef SINGLE
    cgetri_
#else
    zgetri_
#endif
    (
        &N,             // order of the matrix
        inv.begin(),    // LU-factors (on return the inverse matrix)
        &N,             // leading dimension of A
        &IPIV[0],       // pivot indices
        &WORK[0],       // work array
        &LWORK,         // work array length
        &INFO           // diagnostic information
    );

#else
    HexException("Cannot invert dense matrix without LAPACK support.");
#endif
}

template<>
void ColMatrix<Complex>::diagonalize_g
(
    ColMatrix<Complex> const & S,
    NumberArray<Complex> & W,
    ColMatrix<Complex> * VL,
    ColMatrix<Complex> * VR
) const
{
#ifndef NO_LAPACK
    if (rows() != cols())
        HexException("Only square matrices can be diagonalized.");

    Complex* B = const_cast<Complex*>(S.data().data());

    int N = rows();
    cArray WORK(1);
    rArray RWORK(8*N);
    int LWORK = -1, INFO;
    char JOBL = 'N', JOBR = 'N';
    Complex *pVL = nullptr, *pVR = nullptr;

    cArray U (N);
    W.resize(N);

    if (VL != nullptr)
    {
        // compute left eigenvectors
        JOBL = 'V';

        // resize output data
        if (VL->rows() != N or VL->cols() != N)
            *VL = std::move(ColMatrix<Complex>(N,N));
        pVL = VL->begin();
    }

    if (VR != nullptr)
    {
        // compute right eigenvectors
        JOBR = 'V';

        // resize output data
        if (VR->rows() != N or VR->cols() != N)
            *VR = std::move(ColMatrix<Complex>(N,N));
        pVR = VR->begin();
    }

    // copy matrix data (it will be overwritten)
    cArray A = data();

    // get work size
#ifdef SINGLE
    cggev_
#else
    zggev_
#endif
    (
        &JOBL,      // compute left eigenvectors
        &JOBR,      // compute right eigenvectors
        &N,         // order of the matrix
        &A[0],      // matrix A elements
        &N,         // leading dimension of the matrix A
        &B[0],      // matrix S elements
        &N,         // leading dimension of the matrix S
        &W[0],      // eigenvalue numerators
        &U[0],      // eigenvalue denominators
        pVL,        // left eigenvectors elements
        &N,         // leading dimension of the matrix VL
        pVR,        // right eigenvectors elements
        &N,         // leading dimension of the matrix VR
        &WORK[0],   // work array
        &LWORK,     // length of the work array
        &RWORK[0],  // work array
        &INFO       // diagnostic information
    );

    // resize work array
    LWORK = WORK[0].real();
    WORK.resize(LWORK);

    // run the diagonalization
#ifdef SINGLE
    cggev_
#else
    zggev_
#endif
    (
        &JOBL,      // compute left eigenvectors
        &JOBR,      // compute right eigenvectors
        &N,         // order of the matrix
        &A[0],      // matrix A elements
        &N,         // leading dimension of the matrix A
        &B[0],      // matrix S elements
        &N,         // leading dimension of the matrix S
        &W[0],      // eigenvalue numerators
        &U[0],      // eigenvalue denominators
        pVL,        // left eigenvectors elements
        &N,         // leading dimension of the matrix VL
        pVR,        // right eigenvectors elements
        &N,         // leading dimension of the matrix VR
        &WORK[0],   // work array
        &LWORK,     // length of the work array
        &RWORK[0],  // work array
        &INFO       // diagnostic information
    );

    if (INFO != 0)
    {
        std::cerr << "Diagonalization error: " << INFO << std::endl;
    }

    W /= U;

#else
    HexException("Cannot diagonalize dense matrix without LAPACK support.");
#endif
}

// specialization of diagonalization of complex dense column-major ("Fortran-like") matrices
template<>
void ColMatrix<Complex>::diagonalize
(
    NumberArray<Complex> & W,
    ColMatrix<Complex> * VL,
    ColMatrix<Complex> * VR
) const
{
#ifndef NO_LAPACK
    if (rows() != cols())
        HexException("Only square matrices can be diagonalized.");

    int N = rows();
    cArray WORK(1);
    rArray RWORK(2*N);
    int LWORK = -1, INFO;
    char JOBL = 'N', JOBR = 'N';
    Complex *pVL = nullptr, *pVR = nullptr;

    W.resize(N);

    if (VL != nullptr)
    {
        // compute left eigenvectors
        JOBL = 'V';

        // resize output data
        if (VL->rows() != N or VL->cols() != N)
            *VL = std::move(ColMatrix<Complex>(N,N));
        pVL = VL->begin();
    }

    if (VR != nullptr)
    {
        // compute right eigenvectors
        JOBR = 'V';

        // resize output data
        if (VR->rows() != N or VR->cols() != N)
            *VR = std::move(ColMatrix<Complex>(N,N));
        pVR = VR->begin();
    }

    // copy matrix data (it will be overwritten)
    cArray A = data();

    // get work size
#ifdef SINGLE
    cgeev_
#else
    zgeev_
#endif
    (
        &JOBL,      // compute left eigenvectors
        &JOBR,      // compute right eigenvectors
        &N,         // order of the matrix
        &A[0],      // matrix A elements
        &N,         // leading dimension of the matrix A
        &W[0],      // eigenvalues
        pVL,        // left eigenvectors elements
        &N,         // leading dimension of the matrix VL
        pVR,        // right eigenvectors elements
        &N,         // leading dimension of the matrix VR
        &WORK[0],   // work array
        &LWORK,     // length of the work array
        &RWORK[0],  // work array
        &INFO       // diagnostic information
    );

    // resize work array
    LWORK = WORK[0].real();
    WORK.resize(LWORK);

    // run the diagonalization
#ifdef SINGLE
    cgeev_
#else
    zgeev_
#endif
    (
        &JOBL,      // compute left eigenvectors
        &JOBR,      // compute right eigenvectors
        &N,         // order of the matrix
        &A[0],      // matrix A elements
        &N,         // leading dimension of the matrix A
        &W[0],      // eigenvalues
        pVL,        // left eigenvectors elements
        &N,         // leading dimension of the matrix VL
        pVR,        // right eigenvectors elements
        &N,         // leading dimension of the matrix VR
        &WORK[0],   // work array
        &LWORK,     // length of the work array
        &RWORK[0],  // work array
        &INFO       // diagnostic information
    );

    if (INFO != 0)
    {
        std::cerr << "Diagonalization error: " << INFO << std::endl;
    }
#else
    HexException("Cannot diagonalize dense matrix without LAPACK support.");
#endif
}

// -------------------------------------------------------------------------------------
// Sparse matrix routines.
//

template<>
void SymBandMatrix<Complex>::sym_band_dot (int n, int d, const cArrayView M, Complex alpha, const cArrayView A, Complex beta, cArrayView B, MatrixSelection::Selection tri)
{
    // check dimensions
    if (A.size() % n != 0)
        HexException("Incompatible dimensions: %d (mat) × %ld (vec). You are probably mixing radial data for different grids.", n, A.size());

    // get number of source vectors
    int Nvec = A.size() / n;

    // skip if there are no data
    if (n == 0 or d == 0)
        return;

    // matrix data pointer
    Complex const * const restrict pM = M.data();

    // source and destination vector pointers
    Complex const * restrict pA = A.data();
    Complex       * restrict pB = B.data();

    // for all source vectors (columns)
    for (int j = 0; j < Nvec; j++)
    {
        // scale the destination vector
        for (int k = 0; k < n; k++)
            pB[k] *= beta;

        // for all rows/cols of the matrix (and of the target vector)
        for (int k = 0; k < n; k++)
        {
            // Diagonal element

                if (tri & MatrixSelection::Diagonal)
                {
                    pB[k] += pM[k * d] * alpha * pA[k];
                }

            // Direct pass (~ strict upper triangle)

                if (tri & MatrixSelection::StrictUpper)
                {
/*#if !defined(SINGLE) && defined(__AVX512F__)
                    if (d == 5 and k + 4 < n)
                    {
                        // AVX 512 - manual optimization for order = 4
                        //
                        //   pB[k] += alpha * (pA[k + 1] * pM[k * d + 1]);
                        //   pB[k] += alpha * (pA[k + 2] * pM[k * d + 2]);
                        //   pB[k] += alpha * (pA[k + 3] * pM[k * d + 3]);
                        //   pB[k] += alpha * (pA[k + 4] * pM[k * d + 4]);

                        Complex X[4] = { 0., 0., 0., 0. };

                        cmul4xd(X, pA + k + 1, pM + k * d + 1); // AVX 512

                        pB[k] += alpha * (X[0] + X[1] + X[2] + X[3]);
                    }
                    else
#elif !defined(SINGLE) && defined(__AVX2__) && defined(__FMA__)
                    if (d == 5 and k + 4 < n)
                    {
                        // AVX 256 - manual optimization for order = 4
                        //
                        //   pB[k] += alpha * pA[k + 1] * pM[k * d + 1];
                        //   pB[k] += alpha * pA[k + 2] * pM[k * d + 2];
                        //
                        //   pB[k] += alpha * pA[k + 3] * pM[k * d + 3];
                        //   pB[k] += alpha * pA[k + 4] * pM[k * d + 4];

                        Complex X[2] = { 0., 0. };

                        cmul2xd(X, pA + k + 1, pM + k * d + 1); // AVX2
                        cmul2xd(X, pA + k + 3, pM + k * d + 3); // AVX2

                        pB[k] += alpha * (X[0] + X[1]);
                    }
                    else
#endif*/
                    {
                        for (int l = 1; l < d and k + l < n; l++)
                            pB[k] += pM[k * d + l] * alpha * pA[k + l];
                    }
                }

            // Mirror pass (~ strict lower triangle)

                if (tri & MatrixSelection::StrictLower)
                {
/*#if !defined(SINGLE) && defined(__AVX512F__)
                    if (d == 5 and k + 4 < n)
                    {
                        // AVX 512 - manual optimization for order = 4
                        //
                        //   pB[k + 1] += (alpha * pA[k]) * pM[k * d + 1];
                        //   pB[k + 2] += (alpha * pA[k]) * pM[k * d + 2];
                        //   pB[k + 3] += (alpha * pA[k]) * pM[k * d + 3];
                        //   pB[k + 4] += (alpha * pA[k]) * pM[k * d + 4];

                        Complex X[4];

                        X[0] = X[1] = X[2] = X[3] = alpha * pA[k];

                        cmul4xd(pB + k + 1, X, pM + k * d + 1); // AVX 512
                    }
                    else
#elif !defined(SINGLE) && defined(__AVX2__) && defined(__FMA__)
                    if (d == 5 and k + 4 < n)
                    {
                        // AVX 256 - manual optimization for order = 4
                        //
                        //   pB[k + 1] += (alpha * pA[k]) * pM[k * d + 1];
                        //   pB[k + 2] += (alpha * pA[k]) * pM[k * d + 2];
                        //
                        //   pB[k + 3] += (alpha * pA[k]) * pM[k * d + 3];
                        //   pB[k + 4] += (alpha * pA[k]) * pM[k * d + 4];

                        Complex X[2];

                        X[0] = X[1] = alpha * pA[k];

                        cmul2xd(pB + k + 1, X, pM + k * d + 1); // AVX2
                        cmul2xd(pB + k + 3, X, pM + k * d + 3); // AVX2
                    }
                    else
#endif*/
                    {
                        for (int l = 1; l < d and k + l < n; l++)
                            pB[k + l] += pM[k * d + l] * alpha * pA[k];
                    }
                }
        }

        // move on to the next source vector
        pA += n;
        pB += n;
    }

    return;
}

template<>
CsrMatrix<LU_int_t,Complex> CooMatrix<LU_int_t,Complex>::tocsr () const
{
    // get number of structurally non-zero elements
    LU_int_t nz = x_.size();
    assert((LU_int_t)i_.size() == nz);

    // get row lengths
    std::vector<LU_int_t> len (m_, 0);
    for (LU_int_t n = 0; n < nz; n++) 
    {
        if (x_[n] != 0.0_r)
        {
            assert(i_[n] < m_);
            len[i_[n]]++;
        }
    }

    // create element pointer array for each matrix row
    std::vector<std::vector<LU_int_t>> elem_ptrs (m_);

    // reserve memory for the row data
    for (LU_int_t n = 0; n < m_; n++)
        elem_ptrs[n].reserve(len[n]);

    // store index of each element into the appropriate row
    for (LU_int_t n = 0; n < nz; n++) if (x_[n] != 0.0_r)
        elem_ptrs[i_[n]].push_back(n);

    // sort element pointers by their column index
    # pragma omp parallel for
    for (LU_int_t n = 0; n < m_; n++)
        std::sort(elem_ptrs[n].begin(), elem_ptrs[n].end(), [&](LU_int_t a, LU_int_t b) { return j_[a] < j_[b]; });

    // allocate output arrays
    NumberArray<LU_int_t> Ap(m_ + 1), Ai(nz);
    cArray Ax(nz);

    // copy & sum entries
    LU_int_t pos = 0;
    for (LU_int_t m = 0; m < m_; m++)
    {
        // for all non-zero elements
        for (LU_int_t n = 0; n < (LU_int_t)elem_ptrs[m].size(); n++)
        {
            // get global index
            LU_int_t gid = elem_ptrs[m][n];

            // update elements
            Ai[pos]  = j_[gid];
            Ax[pos] += x_[gid];

            // if the next row entry is in the same column, do not advance position counter,
            // so that the element will be added to current position
            if (n < (LU_int_t)elem_ptrs[m].size() - 1 and j_[gid] == j_[elem_ptrs[m][n + 1]])
                continue;

            // otherwise move on to the next position
            pos++;
        }

        // insert next row pointer
        Ap[m + 1] = pos;
    }

    return CsrMatrix<LU_int_t,Complex> (m_, n_, std::move(Ap), std::move(Ai), std::move(Ax));
}
