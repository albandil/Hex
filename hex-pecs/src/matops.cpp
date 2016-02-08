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

// Standard LAPACK routine for dense inversion (real variant).
extern "C" void dgetri_
(
    int * n,
    double * A, int * lda, int * ipiv,
    double * work, int * lwork, int * info
);

// Standard LAPACK routine for dense LU factorization (complex variant).
extern "C" void zgetrf_
(
    int * m, int * n,
    Complex * A, int * lda,
    int * pivots, int * Info
);

// Standard LAPACK routine for dense inversion (complex variant).
extern "C" void zgetri_
(
    int * n,
    Complex * A, int * lda, int * ipiv,
    Complex * work, int * lwork, int * info
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
        // erase destination vector
        for (std::size_t i = 0; i < M; i++)
            w[i] = 0;
        
        // for all columns of the dense matrix
        for (std::size_t j = 0; j < N; j++)
        {
            // get column pointer
            T const * restrict pA = A + j * M;
            
            // for all elements of the column: update elements of the destination vector
            for (std::size_t i = 0; i < M; i++)
                w[i] += pA[i] * v[j];
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
                D[jcol * Nblocks * N + irow] += B[((iblock * Nblocks + jblock) * N + i) * (2 * Ndiag + 1) + j + Ndiag - i];
            }
        }
    }
    
    template void dense_add_blockband<double> (std::size_t Nblocks, std::size_t N, std::size_t Ndiag, double * restrict D, double const * restrict B);
    template void dense_add_blockband<Complex> (std::size_t Nblocks, std::size_t N, std::size_t Ndiag, Complex * restrict D, Complex const * restrict B);
}

// --------------------------------------------------------------------------------- //

namespace matops
{
    template<class T> void dense_invert (std::size_t N, T * A, int * pivots, T * work)
    {
        HexException("Unsupported type.");
    }
    
    template<> void dense_invert<double> (std::size_t N, double * A, int * pivots, double * work)
    {
        int n = N, lda = N, lwork = N * N, info;
        
        dgetrf_(&n, &n, A, &lda, pivots, &info);
        if (info != 0)
            HexException("DGETRF failed with error code %d.", info);
        
        dgetri_(&n, A, &lda, pivots, work, &lwork, &info);
        if (info != 0)
            HexException("DGETRI failed with error code %d.", info);
    }
    
    template<> void dense_invert<Complex> (std::size_t N, Complex * A, int * pivots, Complex * work)
    {
        int n = N, lda = N, lwork = N * N, info;
        
        zgetrf_(&n, &n, A, &lda, pivots, &info);
        if (info != 0)
            HexException("ZGETRF failed with error code %d.", info);
        
        zgetri_(&n, A, &lda, pivots, work, &lwork, &info);
        if (info != 0)
            HexException("ZGETRI failed with error code %d.", info);
    }
}

// --------------------------------------------------------------------------------- //

namespace matops
{
    template <class T>
    void dense_mul_blockband (std::size_t Nblocks, std::size_t M, std::size_t K, std::size_t N, std::size_t Ndiag, T const * D, T const * B, T * C)
    {
        // for all rows of the dense matrix
        for (std::size_t drow = 0; drow < Nblocks * M; drow++)
        {
            // for all columns of the block-band matrix
            for (std::size_t col_block = 0; col_block < Nblocks; col_block++)
            for (std::size_t icol = 0; icol < N; icol++)
            {
                // erase the output element
                C[(col_block * N + icol) * Nblocks * M + drow] = 0.;
                
                // for all rows of the block-band matrix (and columns of the matrix 'D')
                for (std::size_t row_block = 0; row_block < Nblocks; row_block++)
                for (std::size_t irow = (icol < Ndiag ? 0 : icol - Ndiag); irow <= icol + 1 and irow < K; irow++)
                {
                    C[(col_block * N + icol) * Nblocks * M + drow] += D[(row_block * K + irow) * Nblocks * M + drow] * B[((row_block * Nblocks + col_block) * K + irow) * (2*Ndiag + 1) + (icol + Ndiag - irow)];
                }
            }
        }
    }
    
    template void dense_mul_blockband<double> (std::size_t Nblocks, std::size_t M, std::size_t K, std::size_t N, std::size_t Ndiag, double const * D, double const * B, double * C);
    template void dense_mul_blockband<Complex> (std::size_t Nblocks, std::size_t M, std::size_t K, std::size_t N, std::size_t Ndiag, Complex const * D, Complex const * B, Complex * C);
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
        // multiply all columns independently
        for (std::size_t dcol = 0; dcol < Nblocks * N; dcol++)
            blockband_mul_vector(Nblocks, M, K, Ndiag, A, D + dcol * Nblocks * K, B + dcol * Nblocks * M);
    }
    
    template void blockband_mul_dense<double> (std::size_t Nblocks, std::size_t M, std::size_t K, std::size_t N, std::size_t Ndiag, double const * restrict A, double const * restrict D, double * restrict B);
    template void blockband_mul_dense<Complex> (std::size_t Nblocks, std::size_t M, std::size_t K, std::size_t N, std::size_t Ndiag, Complex const * restrict A, Complex const * restrict D, Complex * restrict B);
}
