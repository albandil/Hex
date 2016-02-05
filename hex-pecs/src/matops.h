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
    
    /// Multiply vector by a dense matrix.
    template <class T>
    void dense_mul_vector (std::size_t M, std::size_t N, T const * restrict A, T const * restrict v, T * restrict w);
    
    /**
     * @brief Multiply block-band matrix by an inverse matrix (given the LU decomposition of the dense matrix).
     * 
     * The LU decomposition is assumed to come from @ref dense_LU_factor, i.e. it is expected to be of a transposed matrix.
     * The solution then needs to be rephrased from common LUx = Py to (P⁻¹LU)'x = y, i.e.
     * @f[
     *         x = P^{-1} L^{-\top} U^{-\top} y \,.
     * @f]
     * 
     * @param Ndiag Number of upper diagonals.
     * @param workspace Workspave of size Nblocks * N.
     */
    template <class T>
    void dense_LU_solve_blockband
    (
        std::size_t Nblocks, std::size_t M, std::size_t N, std::size_t Ndiag,
        T const * restrict LU, int const * restrict pivots, T const * restrict B, T * restrict C, T * restrict workspace
    );
    
    /// Add a block-band matrix to dense matrix; 'Ndiag' is the number of upper diagonals.
    template <class T>
    void dense_add_blockband (std::size_t Nblocks, std::size_t N, std::size_t Ndiag, T * restrict D, T const * restrict B);
    
    /**
     * @brief A simple Wrapper around xGETRF.
     * 
     * The LAPACK routine for LU factorization assumes the column-major storage. However, we are passing in a row-major stored matrix A.
     * In such a case it will be
     * @f[
     *         PA^T = LU \,.
     * @f]
     */
    template <class T>
    void dense_LU_factor (std::size_t N, T * A, int * pivots);
    
    /**
     * @brief A simple Wrapper around xGETRS.
     * 
     * This routine expects to receive LU factorization of a transposed column-matrix, or of a non-transposed row-matrix,
     * which is exactly what the routine @ref dense_LU_factor does.
     */
    template <class T>
    void dense_LU_solve_vector (std::size_t N, T const * restrict M, int const * restrict pivots, T const * restrict v, T * restrict w);
    
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
     * The storage scheme of D has to be column-major. The storage of B will be row-major.
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

