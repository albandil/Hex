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

#ifndef HEX_CSCMATRIX_H
#define HEX_CSCMATRIX_H

#ifndef NO_UMFPACK
    #include <umfpack.h>
#endif

#include "arrays.h"
#include "coomatrix.h"

/**
 * @brief Complex CSC matrix.
 * 
 * Complex compressed-storage-column-major-ordered sparse matrix.
 * The data are stored in three arrays,
 * @code
 * p, i, x
 * @endcode
 * The format is that used in UMFPACK with interleaved real and comples parts;
 * for explanation see UMFPACK Uses's Guide.
 */
class CscMatrix
{
    
public:
    
    // Constructors
    
    CscMatrix ()
        : m_(0), n_(0) {}
    CscMatrix (std::size_t m, std::size_t n)
        : m_(m), n_(n) {}
    CscMatrix (CscMatrix const & A)
        : m_(A.m_), n_(A.n_), p_(A.p_), i_(A.i_), x_(A.x_) {}
    CscMatrix (std::size_t m, std::size_t n, const lArrayView p, const lArrayView i, const cArrayView x)
        : m_(m), n_(n), p_(p), i_(i), x_(x) {}
    
    // Destructor
    
    ~CscMatrix () {}
    
    /**
     * Convert to COO-matrix.
     */
    CooMatrix tocoo () const;
    
    /**
     * Matrix-vector product using a transposed matrix, \f$ A^T \cdot b \f$.
     * For ordinary matrix vector product convert this matrix first to CSR
     * format.
     * @param b Vector to multiply with.
     */
    cArray dotT (const cArrayView b) const;
    
    // Getters
    
    std::size_t size () const { return i_.size(); }
    std::size_t rows () const { return m_; }
    std::size_t cols () const { return n_; }
    lArray const & p () const { return p_; }
    lArray const & i () const { return i_; }
    cArray const & x () const { return x_; }
    
    /**
     * Multiplication by a real number.
     */
    CscMatrix & operator *= (double r);
    
    /**
     * Addition of another CSC matrix.
     * The matrices MUST HAVE THE SAME SPARSE STRUCTURE, as no indices are
     * checked.
     */
    CscMatrix & operator &= (const CscMatrix&  B);
    
    /**
     * Subtraction of another CSC matrix.
     * The matrices MUST HAVE THE SAME SPARSE STRUCTURE, as no indices are
     * checked.
     */
    CscMatrix& operator ^= (const CscMatrix&  B);
    
    /**
     * Save matrix to HDF file.
     * @param name Filename.
     */
    bool hdfsave (const char* name) const;
    
    /**
     * Load matrix from HDF file.
     * @param name Filename.
     */
    bool hdfload (const char* name);
    
private:
    
    // dimensions
    std::int64_t m_;
    std::int64_t n_;
    
    // representation
    lArray p_;
    lArray i_;
    cArray x_;
};

/**
 * Computes a sum of two csc-matrices OF THE SAME SPARSE STRUCTURE.
 */
inline CscMatrix operator & (CscMatrix const & A, CscMatrix const & B)
{
    CscMatrix C = A;
    return C &= B;
}

/**
 * Computes a difference of two csc-matrices OF THE SAME SPARSE STRUCTURE.
 */
inline CscMatrix operator ^ (CscMatrix const & A, CscMatrix const & B)
{
    CscMatrix C = A;
    return C ^= B;
}

/**
 * Multiplication of csc-matrix by a number.
 */
inline CscMatrix operator * (double r, CscMatrix const & B)
{
    CscMatrix C = B;
    return C *= r;
}

/**
 * Multiplication of csc-matrix by a number.
 */
inline CscMatrix operator * (CscMatrix const & A, double r)
{
    CscMatrix C = A;
    return C *= r;
}

#endif // HEX_CSCMATRIX_H
