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

#include "hex-arrays.h"
#include "hex-hdffile.h"

#include "matops.h"

// --------------------------------------------------------------------------------- //

// Standard LAPACK routine for dense LU factorization (real variant).
extern "C" void dgetrf_
(
    int * m, int * n,
    double * A, int * lda,
    int * pivots, int * Info
);

// Standard LAPACK routine for dense LU backsubstitution (real variant).
extern "C" void dgetrs_
(
    char * trans,
    int * N, int * nrhs,
    double * A, int * lda, int * pivots,
    double * B, int * ldb, int *info
);

// Standard LAPACK routine for dense LU factorization (complex variant).
extern "C" void zgetrf_
(
    int * m, int * n,
    Complex * A, int * lda,
    int * pivots, int * Info
);

// Standard LAPACK routine for dense LU backsubstitution (complex variant).
extern "C" void zgetrs_
(
    char * trans,
    int * N, int * nrhs,
    Complex * A, int * lda, int * pivots,
    Complex * B, int * ldb, int *info
);

// --------------------------------------------------------------------------------- //

namespace matops
{
    // definition of 'save' template
    template <class T> bool save (T const * data, std::size_t N, std::string filename)
    {
        HDFFile hdf (filename, HDFFile::overwrite);
        return hdf.valid() and hdf.write("size", &N, 1) and hdf.write("data", data, N);
    }

    // explicit instantiation of 'save' for needed types
    template bool save<int> (int const * data, std::size_t N, std::string filename);
    template bool save<double> (double const * data, std::size_t N, std::string filename);
    template bool save<Complex> (Complex const * data, std::size_t N, std::string filename);

    // definition of 'load' template
    template <class T> bool load (T * data, std::size_t N, std::string filename)
    {
        HDFFile hdf (filename, HDFFile::readonly);
        std::size_t M;
        return hdf.valid() and hdf.read("size", &M, 1) and M == N and hdf.read("data", data, N);
    }

    // explicit instantiation of 'load' for needed types
    template bool load<int> (int * data, std::size_t N, std::string filename);
    template bool load<double> (double * data, std::size_t N, std::string filename);
    template bool load<Complex> (Complex * data, std::size_t N, std::string filename);
}

// --------------------------------------------------------------------------------- //

namespace matops
{
    // definition of 'sum' template
    template <class T> void sum (std::size_t N, T const * restrict x, T const * restrict y, T * restrict z)
    {
        for (std::size_t i = 0; i < N; i++)
            z[i] = x[i] + y[i];
    }

    // explicit instantiation of 'sum' for needed types
    template void sum<double> (std::size_t N, double const * restrict x, double const * restrict y, double * restrict z);
    template void sum<Complex> (std::size_t N, Complex const * restrict x, Complex const * restrict y, Complex * restrict z);

    // definition of 'subtract' template
    template <class T> void subtract (std::size_t N, T const * restrict x, T const * restrict y, T * restrict z)
    {
        for (std::size_t i = 0; i < N; i++)
            z[i] = x[i] - y[i];
    }

    // explicit instantiation of 'subtract' for needed types
    template void subtract<double> (std::size_t N, double const * restrict x, double const * restrict y, double * restrict z);
    template void subtract<Complex> (std::size_t N, Complex const * restrict x, Complex const * restrict y, Complex * restrict z);

    // definition of 'flip_sign' template
    template <class T> void flip_sign (std::size_t N, T * restrict x)
    {
        for (std::size_t i = 0; i < N; i++)
            x[i] = -x[i];
    }

    // explicit instantiation of 'flip_sign' for needed types
    template void flip_sign<double> (std::size_t N, double * restrict x);
    template void flip_sign<Complex> (std::size_t N, Complex * restrict x);
}

// --------------------------------------------------------------------------------- //

namespace matops
{
    // definition of 'dense_mul_vector' template
    template <class T> void dense_mul_vector (std::size_t M, std::size_t N, T const * restrict A, T const * restrict v, T * restrict w)
    {
        for (std::size_t i = 0; i < M; i++)
        {
            w[i] = 0;
            T const * restrict pA = A + i * N;
            
            for (std::size_t j = 0; j < N; j++)
                w[i] += pA[j] * v[j];
        }
    }

    // explicit instantiation of 'dense_mul_vector' for needed types
    template void dense_mul_vector<double> (std::size_t M, std::size_t N, double const * restrict A, double const * restrict v, double * restrict w);
    template void dense_mul_vector<Complex> (std::size_t M, std::size_t N, Complex const * restrict A, Complex const * restrict v, Complex * restrict w);
}

// --------------------------------------------------------------------------------- //

namespace matops
{
    template <class T>
    void dense_add_blockband (std::size_t Nblocks, std::size_t N, std::size_t Ndiag, T * restrict D, T const * restrict B)
    {
        // for all rows
        for (std::size_t iblock = 0; iblock < Nblocks; iblock++)
        for (std::size_t i = 0; i < N; i++)
        {
            // calculate row index
            std::size_t irow = iblock * N + i;
            
            // for all structurally nonzero elements of B in this row
            for (std::size_t jblock = 0; jblock < Nblocks; jblock++)
            for (std::size_t j = (i < Ndiag ? 0 : i - Ndiag); j <= i + Ndiag and j < N; j++)
            {
                // dense column
                std::size_t jcol = jblock * N + j;
                
                // update element of the dense matrix
                D[irow * Nblocks * N + jcol] += B[((iblock * Nblocks + jblock) * N + i) * (2 * Ndiag + 1) + j + Ndiag - i];
//                 std::cout << "B's contribution from " << i << " " << j << " (@ " << ((iblock * Nblocks + jblock) * N + i) * (2 * Ndiag + 1) + j + Ndiag - i << "): " << B[((iblock * Nblocks + jblock) * N + i) * (2 * Ndiag + 1) + j + Ndiag - i] << std::endl;
//                 std::cout << "new D's value (@ " << irow * Nblocks * N + jcol << "): " << D[irow * Nblocks * N + jcol] << std::endl;
            }
        }
    }

    template void dense_add_blockband<double> (std::size_t Nblocks, std::size_t N, std::size_t Ndiag, double * restrict D, double const * restrict B);
    template void dense_add_blockband<Complex> (std::size_t Nblocks, std::size_t N, std::size_t Ndiag, Complex * restrict D, Complex const * restrict B);
}

// --------------------------------------------------------------------------------- //

namespace matops
{
    template<class T> void dense_LU_factor (std::size_t N, T * A, int * pivots)
    {
        HexException("Unsupported type.");
    }
    
    template<> void dense_LU_factor<double> (std::size_t N, double * A, int * pivots)
    {
        int m = N, n = N, lda = N, info;
        dgetrf_(&m, &n, A, &lda, pivots, &info);
    }
    
    template<> void dense_LU_factor<Complex> (std::size_t N, Complex * A, int * pivots)
    {
        int m = N, n = N, lda = N, info;
        zgetrf_(&m, &n, A, &lda, pivots, &info);
    }
}

// --------------------------------------------------------------------------------- //

namespace matops
{
    template<> void dense_LU_solve_vector<double> (std::size_t N, double const * restrict M, int const * restrict pivots, double const * restrict v, double * restrict w)
    {
        // copy v to w
        for (std::size_t i = 0; i < N; i++)
            w[i] = v[i];
        
        char trans = 'T';
        int n = N, nrhs = 1, lda = N, ldb = N, info;
        dgetrs_(&trans, &n, &nrhs, const_cast<double*>(M), &lda, const_cast<int*>(pivots), w, &ldb, &info);
    }

    template<> void dense_LU_solve_vector<Complex> (std::size_t N, Complex const * restrict M, int const * restrict pivots, Complex const * restrict v, Complex * restrict w)
    {
        // copy v to w
        for (std::size_t i = 0; i < N; i++)
            w[i] = v[i];
        
        char trans = 'T';
        int n = N, nrhs = 1, lda = N, ldb = N, info;
        zgetrs_(&trans, &n, &nrhs, const_cast<Complex*>(M), &lda, const_cast<int*>(pivots), w, &ldb, &info);
    }
}

// --------------------------------------------------------------------------------- //

namespace matops
{
    template <class T>
    void dense_LU_solve_blockband (std::size_t Nblocks, std::size_t M, std::size_t N, std::size_t Ndiag, T const * LU, const int* pivots, T const * B, T * C, T * w)
    {
        // The LU is for transposed matrix, so instead of x = U⁻¹L⁻¹Py we need x = P⁻¹L'⁻¹U'⁻¹y.
        
        // FIXME : Avoid calls to dense_LU_solve_vector, which is O(N^2), making this routine O(N^3).
        
//         std::cout << "LU.norm = " << ArrayView<T>(Nblocks * M * Nblocks * M, const_cast<T*>(LU)).norm() << std::endl;
        for (std::size_t n = 0; n < Nblocks; n++)
        for (std::size_t icol = 0; icol < N; icol++)
        {
            for (std::size_t m = 0; m < Nblocks; m++)
            for (std::size_t irow = 0; irow < M; irow++)
            {
                if (std::max(irow,icol) - std::min(irow,icol) <= Ndiag)
                    w[m * M + irow] = B[((m * Nblocks + n) * M + irow) * (2*Ndiag + 1) + icol + Ndiag - irow];
                else
                    w[m * M + irow] = 0.;
            }
            
//             std::cout << ArrayView<T>(Nblocks * M, w) << " ";
            dense_LU_solve_vector(Nblocks * M, LU, pivots, w, C + (n * Nblocks + icol) * M);
//             std::cout << ArrayView<T>(Nblocks * M, C + (n * Nblocks + icol) * M).norm() << " | ";
        }
//         std::cout << std::endl;
    }
    
    template void dense_LU_solve_blockband<double> (std::size_t Nblocks, std::size_t M, std::size_t N, std::size_t Ndiag, double const * LU, int const * pivots, double const * B, double * C, double * workspace);
    template void dense_LU_solve_blockband<Complex> (std::size_t Nblocks, std::size_t M, std::size_t N, std::size_t Ndiag, Complex const * LU, int const * pivots, Complex const * B, Complex * C, Complex * workspace);
}

// --------------------------------------------------------------------------------- //

namespace matops
{
    template <class T> void blockband_mul_vector (std::size_t Nblocks, std::size_t M, std::size_t N, std::size_t Ndiag, T const * restrict A, T const * restrict v, T * restrict w)
    {
        // for all row blocks
        for (std::size_t m = 0; m < Nblocks; m++)
        {
            // for all column blocks
            for (std::size_t n = 0; n < Nblocks; n++)
            {
                // Multiply n-the segment of the source vector by current block,
                // store to m-th segment of destination vector.
                
                // for all rows of the block
                for (std::size_t irow = 0; irow < M; irow++)
                {
                    // erase the output element
                    w[m * M + irow] = 0.;
                    
                    // for all structurally non-zero column entries in this row
                    for (std::size_t icol = (irow < Ndiag ? 0 : irow - Ndiag); icol <= irow + Ndiag and icol < N; icol++)
                        w[m * M + irow] += A[((m * Nblocks + n) * M + irow) * (2*Ndiag + 1) + (icol + Ndiag - irow)] * v[n * N + icol];
                }
            }
        }
    }
    
    template void blockband_mul_vector<double> (std::size_t Nblocks, std::size_t M, std::size_t N, std::size_t Ndiag, double const * restrict A, double const * restrict v, double * restrict w);
    template void blockband_mul_vector<Complex> (std::size_t Nblocks, std::size_t M, std::size_t N, std::size_t Ndiag, Complex const * restrict A, Complex const * restrict v, Complex * restrict w);
    
    template <class T> void blockband_mul_dense (std::size_t Nblocks, std::size_t M, std::size_t K, std::size_t N, std::size_t Ndiag, T const * restrict A, T const * restrict D, T * restrict B)
    {
        // for all blocks of the block-band matrix
        for (std::size_t m = 0; m < Nblocks; m++)
        for (std::size_t n = 0; n < Nblocks; n++)
        {
            // for all columns of the dense matrix
            for (std::size_t dcol = 0; dcol < Nblocks * N; dcol++)
            {
                // for all rows of the block
                for (std::size_t irow = 0; irow < M; irow++)
                {
                    // erase the output element
                    B[(m * M + irow) * (Nblocks * N) + dcol] = 0.;
                    
                    // for all structurally non-zero column entries in this block's row
                    for (std::size_t icol = (irow < Ndiag ? 0 : irow - Ndiag); icol <= irow + Ndiag and icol < K; icol++)
                    {
//                         std::cout << "dcol = " << dcol << ", irow = " << irow << ", icol = " << icol << std::endl;
//                         std::cout << "A[" << ((m * Nblocks + n) * M + irow) * (2*Ndiag + 1) + (icol + Ndiag - irow) << "] = " << A[((m * Nblocks + n) * M + irow) * (2*Ndiag + 1) + (icol + Ndiag - irow)] << std::endl;
//                         std::cout << "D[" << (dcol * Nblocks + n) * K + icol << "] = " << D[(dcol * Nblocks + n) * K + icol] << std::endl;
                        B[(m * M + irow) * (Nblocks * N) + dcol] += A[((m * Nblocks + n) * M + irow) * (2*Ndiag + 1) + (icol + Ndiag - irow)] * D[(dcol * Nblocks + n) * K + icol];
                    }
                }
            }
        }
    }
    
    template void blockband_mul_dense<double> (std::size_t Nblocks, std::size_t M, std::size_t K, std::size_t N, std::size_t Ndiag, double const * restrict A, double const * restrict D, double * restrict B);
    template void blockband_mul_dense<Complex> (std::size_t Nblocks, std::size_t M, std::size_t K, std::size_t N, std::size_t Ndiag, Complex const * restrict A, Complex const * restrict D, Complex * restrict B);
}
