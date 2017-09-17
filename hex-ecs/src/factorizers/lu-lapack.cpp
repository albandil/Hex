//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2017, Jakub Benda, Charles University in Prague                    //
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

// --------------------------------------------------------------------------------- //

#include "hex-blas.h"
#include "lu-lapack.h"

// --------------------------------------------------------------------------------- //

/**
 * @brief Compute the Cuthill-McKee ordering.
 * 
 * Returns a row/column permutation array that will minimize the bandwidth of a symmetric
 * matrix using the Cuthill-McKee algorithm. The matrix structure is given using the
 * CSR index arrays (row pointers and column indices). Only upper triangle positions 
 * are considered.
 * 
 * The algorithm works as follows:
 * - A row with the lowest number of non-zeros is taken as the starting row and inserted
 *   into the resulting row sequence R.
 * - While there are some rows not in R:
 * -- Construct the adjacency set A of R, i.e. all rows that are coupled by the matrix
 *    rows already in R, excluding those already in R. If there is none such row, the
 *    reordered matrix has a block structure and it is necessary to restart the algorithm
 *    from the first step on the remaining rows.
 * -- Sort A by number of non-zeros on the rows.
 * -- Append A to R.
 */
template <class IdxT> NumberArray<IdxT> cuthill_mckee
(
    NumberArray<IdxT> const & P,    // row pointers
    NumberArray<IdxT> const & I     // column indices
)
{
    // number of rows
    IdxT n = P.size() - 1;
    
    // result sequence
    NumberArray<IdxT> R;
    R.reserve(n);
    
    // adjacency set
    NumberArray<IdxT> A;
    A.reserve(n);
    
    // number of non-zeros (in upper triangle)
    NumberArray<IdxT> O (n);
    for (IdxT irow = 0; irow < n; irow++)
        O[irow] = P[irow + 1] - P[irow];
    
    // flags indicating that the row has been already used
    NumberArray<bool> used (n, false);
    
    // number of elements in R whose adjacents were already added
    IdxT adjdone = 0;
    
    // while some rows are undecided
    while ((IdxT)R.size() != n)
    {
        // find an unused row with the lowest degree
        IdxT pos = -1, deg = -1;
        for (IdxT i = 0; i < n; i++)
        {
            if (not used[i] and O[i] > deg)
            {
                pos = i;
                deg = O[i];
            }
        }
        
        // add the row as a new seed
        R.push_back(pos);
        used[pos] = true;
        
        // while some rows are undecided
        while ((IdxT)R.size() != n)
        {
            // check that we added some elements in the previous pass
            if (adjdone == (IdxT)R.size())
                break;
            
            // construct an adjacency set of R
            A.resize(0);
            for (IdxT iadj = adjdone; iadj < (IdxT)R.size(); iadj++)
            {
                IdxT irow = R[iadj];
                
                for (IdxT idx = P[irow]; idx < P[irow + 1]; idx++)
                {
                    IdxT icol = I[idx];
                    
                    if (icol > irow and not used[icol])
                    {
                        A.push_back(icol);
                        used[icol] = true;
                    }
                }
            }
            
            // sort the adjacency sequence by non-zeros
            std::sort(A.begin(), A.end(), [&O](IdxT i, IdxT j) { return O[i] < O[j]; });
            
            // append the adjacency sequence to the result set
            adjdone = R.size();
            R.append(A);
        }
    }
    
    return R;
}

LUft_LAPACK::LUft_LAPACK () : LUft()
{
    
}

void LUft_LAPACK::drop ()
{
    R_.drop();
    LU_.drop();
    ipiv_.drop();
}

LUft_LAPACK::~LUft_LAPACK ()
{
    
}

void LUft_LAPACK::factorize (CsrMatrix<LU_int_t,Complex> const & matrix)
{
    // reorder the matrix to make the bandwidth as small as possible
    //R_ = cuthill_mckee(matrix.p(), matrix.i());
    R_.resize(matrix.rows());
    
    // reverse the order to minimize fill-in
    //std::reverse(R_.begin(), R_.end());
    std::iota(R_.begin(), R_.end(), 0);
    
    // determine bandwidth
    k_ = 0;
    for (blas::Int irow = 0; irow < (blas::Int)matrix.p().size() - 1; irow++)
    {
        for (blas::Int idx = matrix.p()[R_[irow]]; idx < (blas::Int)matrix.p()[R_[irow] + 1]; idx++)
        {
            blas::Int icol = matrix.i()[idx];
            if (matrix.x()[idx] != 0.0_z)
                k_ = std::max<blas::Int>(k_, std::abs(R_[irow] - R_[icol]));
        }
    }
    
    // allocate pivot array and the working matrix
    n_ = matrix.rows();
    ipiv_.resize(n_);
    LU_.resize((3*k_ + 1ULL) * n_);
    LU_.fill(0);
    
    // populate the input matrix
    ColMatrixView<Complex> LU (3*k_ + 1ULL, n_, LU_);
    for (blas::Int irow = 0; irow < (blas::Int)matrix.p().size() - 1; irow++)
    {
        for (blas::Int idx = matrix.p()[R_[irow]]; idx < (blas::Int)matrix.p()[R_[irow] + 1]; idx++)
        {
            blas::Int icol = matrix.i()[idx];
            if (matrix.x()[idx] != 0.0_z)
                LU(2ULL*k_ + R_[irow] - R_[icol], R_[icol]) = matrix.x()[idx];
        }
    }
    
    // factorize the symmetric banded matrix
    blas::Int info = blas::sbtrf(n_, k_, LU_, ipiv_);
    
    // check success indicator
    if (info < 0)
        HexException("LAPACK: The argument %d has illegal value.", -info);
    if (info > 0)
        HexException("LAPACK: Matrix is singular.");
}

void LUft_LAPACK::solve (const cArrayView b, cArrayView x, int eqs) const
{
    // number of equation in a set
    std::size_t neq = b.size() / eqs;
    
    // work array
    cArray xp (x.size());
    
    // reorder all righ-hand sides
    for (int s = 0; s < eqs; s++)
    for (std::size_t i = 0; i < neq; i++)
        xp[s * neq + i] = b[s * neq + R_[i]];
    
    // solve the reordered system
    blas::Int info = blas::sbtrs(n_, k_, LU_, ipiv_, xp);
    
    // check success indicator
    if (info < 0)
        HexException("LAPACK: The argument %d has illegal value.", -info);
    
    // reorder all solutions
    for (int s = 0; s < eqs; s++)
    for (std::size_t i = 0; i < neq; i++)
        x[s * neq + R_[i]] = xp[s * neq + i];
}

std::size_t LUft_LAPACK::size () const
{
    return R_.size()    * sizeof(LU_int_t)
         + LU_.size()   * sizeof(Complex)
         + ipiv_.size() * sizeof(blas::Int);
}

void LUft_LAPACK::save (std::string name) const
{ 
    // TODO
}

void LUft_LAPACK::load (std::string name, bool throw_on_io_failure)
{
    // TODO
    
    if (throw_on_io_failure)
        HexException("Failed to load LAPACK LU decomposition from disk.");
}

// --------------------------------------------------------------------------------- //

addClassToParentRunTimeSelectionTable(LUft, LUft_LAPACK)

// --------------------------------------------------------------------------------- //
