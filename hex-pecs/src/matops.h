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

#include "hex-memory.h"

namespace matops
{
    /// Save dense object (matrix, vector) to disk.
    template <class T>
    bool save (T const * data, std::size_t N, std::string filename);
    
    /// Load dense object (matrix, vector) from disk.
    template <class T>
    bool load (T * data, std::size_t N, std::string filename);
    
    /// Sum two dense arrays.
    template <class T>
    void sum (std::size_t N, T const * restrict x, T const * restrict y, T * restrict z);
    
    /// Subtract two dense arrays.
    template <class T>
    void subtract (std::size_t N, T const * restrict x, T const * restrict y, T * restrict z);
    
    /// Flip sign of a dense array.
    template <class T>
    void flip_sign (std::size_t N, T * restrict x);
    
    /// Multiply vector by a dense column-matrix.
    template <class T>
    void dense_mul_vector (std::size_t M, std::size_t N, T const * restrict A, T const * restrict v, T * restrict w);
    
    /**
     * @brief Multiply banded matrix by a dense column-matrix.
     * 
     * The dense matrix 'D' is Nblocks*M-by-Nblocks*K, the block-band matrix 'B' has Nblocks blocks on the diagonal
     * and each of them is K-by-N. The block-band matrix has Ndiag upper and lower diagonals. The resulting
     * matrix is Nblocks*M-by-Nblocks*N.
     */
    template <class T>
    void dense_mul_blockband
    (
        std::size_t Nblocks, std::size_t M, std::size_t K, std::size_t N, std::size_t Ndiag,
        T const * restrict D, T const * restrict B, T * restrict R
    );
    
    /**
     * @brief Multiply block-band matrix by a dense matrix.
     * 
     * Row-major storage assumed for the dense matrix.
     */
    template <class T>
    void dense_mul_blockband
    (
        std::size_t Nblocks, std::size_t M, std::size_t K, std::size_t N, std::size_t Ndiag,
        T const * restrict D, T const * restrict B, T * restrict C
    );
    
    /// Add a block-band matrix to dense matrix; 'Ndiag' is the number of upper diagonals.
    template <class T>
    void dense_add_blockband (std::size_t Nblocks, std::size_t N, std::size_t Ndiag, T * restrict D, T const * restrict B);
    
    /**
     * @brief A simple Wrapper around xGETRF + xGETRI.
     * 
     * Inverts the dense matrix, overwriting the input array. The parameter 'pivots' is a workspace of N integers.
     * The parameter 'work' is another workspace of N*N matrix elements.
     */
    template <class T>
    void dense_invert (std::size_t N, T * A, int * pivots, T * work);
    
    /// Multiply vector by a block-band matrix.
    template <class T>
    void blockband_mul_vector (std::size_t Nblocks, std::size_t M, std::size_t N,  std::size_t Ndiag, T const * restrict A, T const * restrict v, T * restrict w);
    
    /**
     * @brief Multiply dense matrix by a block-band matrix.
     * 
     * Multiply AD -> B, where A is M'-by-K' and D is K'-by-N' and B is M'-by-N' matrix,
     * where M' = Nblocks * M, K' = Nblocks * K and N' = Nblocks = N.
     * The matrix A is block-banded. It has Nblock x Nblock blocks, every one of them
     * is M-by-K and has one main diagonal and Ndiag upper diagonals and Ndiag lower diagonals.
     * The storage scheme of D has to be column-major, as will be the storage scheme of B.
     */
    template <class T>
    void blockband_mul_dense
    (
        std::size_t Nblocks,
        std::size_t M, std::size_t K, std::size_t N,
        std::size_t Ndiag,
        T const * restrict A, T const * restrict D, T * restrict B
    );
}

