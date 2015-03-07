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

#ifndef HEX_COOMATRIX_H
#define HEX_COOMATRIX_H

#include "arrays.h"
#include "cscmatrix.h"
#include "csrmatrix.h"
#include "densematrix.h"

/**
 * @brief Complex COO matrix.
 * 
 * Coordinate (ijv) storage format sparse matrix.
 * The data are stored in four arrays,
 * @code
 * i, j, x, z
 * @endcode
 * first two of which contain coordinates of nonzero element and the rest
 * contain its real and imaginary part. Same coordinates amy be used several
 * times; resulting element is then sum of these entries.
 */
class CooMatrix
{
    
public:
    
    // Empty constructors
    
    CooMatrix ()
        : m_(0), n_(0), sorted_(true) {}
    CooMatrix (size_t m, size_t n)
        : m_(m), n_(n), sorted_(true) {}
    CooMatrix (CooMatrix const & A)
        : m_(A.m_), n_(A.n_), i_(A.i_), j_(A.j_), x_(A.x_), sorted_(false) {}
    CooMatrix (std::size_t m, std::size_t n, lArrayView i, lArrayView j, cArrayView x)
        : m_(m), n_(n), i_(i), j_(j), x_(x), sorted_(false) {}
    CooMatrix (std::size_t m, std::size_t n, lArray && i, lArray && j, cArray && x)
        : m_(m), n_(n), i_(i), j_(j), x_(x), sorted_(false) {}
    
    /**
     * Copy constructor initialized from dense array.
     * @param m Row counz of new matrix.
     * @param n Column count of new matrix.
     * @param a Column-major ordered dense array with matrix elements.
     *          Only nonzero elements are copied into internal storage.
     */
    template <class T> CooMatrix (std::size_t m, std::size_t n, T a) : m_(m), n_(n), sorted_(false)
    {
        // initialize from column-major formatted input
        std::size_t i = 0;
        for (std::size_t col = 0; col < n; col++)
        {
            for (std::size_t row = 0; row < m; row++)
            {
                // get element from array
                Complex val = *(a + i);
                
                // if nonzero, store
                if (val != 0.)
                {
                    i_.push_back(row);
                    j_.push_back(col);
                    x_.push_back(val);
                }
                
                // move to next element
                i++;
            }
        }
    }
    
    // Destructor
    ~CooMatrix () {}
    
    /// Convert 1×1 matrix to a complex number.
    operator Complex () const
    {
        if (m_ == 1 and n_ == 1)
        {
            // get number of nonzero matrix elements
            std::size_t elems = this->shake().x_.size();
            if (elems > 1)
                HexException("[CooMatrix::operator Complex] more elements than nominal volume!");
            else if (elems == 1)
                return this->shake().x_[0];
            else
                return 0.;
        }
        else
            HexException("[CooMatrix::operator Complex] matrix is not 1×1!");
    }
    
    // Getters
    
    std::size_t rows () const { return m_; }
    std::size_t cols () const { return n_; }
    std::size_t size () const { return i_.size(); }
    lArray const & i () const { return i_; }
    lArray const & j () const { return j_; }
    cArray const & v () const { return x_; }
    
    /// Index operator. Returns the existing value or zero.
    Complex operator() (std::size_t ix, std::size_t iy) const
    {
        for (std::size_t n = 0; n < i_.size(); n++)
             if ((std::size_t)i_[n] == ix and (std::size_t)j_[n] == iy)
                return x_[n];
        return 0.;
    }
    
    /**
     * @brief Symmetrical band populator.
     * 
     * Sets values in the specified band. If there are already some values,
     * these are added as independent (ijv)-triplets and thus summed together
     * with existing entries.
     * @param d Band halfwidth (full width is 2 * d + 1).
     * @param f Object compatible with signature:
     * @code
     *  Complex (*) (long, long)
     * @endcode
     */
    template <class Functor> CooMatrix& symm_populate_band (std::size_t d, Functor f)
    {
        Complex val;
        
        for (std::size_t row = 0; row < m_; row++)
        {
            for (std::size_t col = row; col < n_ and col - row <= d; col++)
            {
                val = f(row,col);
                
                if (val != 0.)
                {
                    i_.push_back(row);
                    j_.push_back(col);
                    x_.push_back(val);
                    
                    if (row != col)
                    {
                        i_.push_back(col);
                        j_.push_back(row);
                        x_.push_back(val);
                    }
                }
            }
        }
        
        sorted_ = false;
        
        return *this;
    }
    
    /**
     * @brief Full populator.
     * 
     * Sets all values. The matrix can become easily dense!
     * @param f Object compatible with signature:
     * @code
     *  Complex (*) (long, long)
     * @endcode
     */
    template <class Functor> CooMatrix& populate (Functor f)
    {
        Complex val;
        
        for (std::size_t row = 0; row < m_; row++)
        {
            for (std::size_t col = 0; col < n_; col++)
            {
                val = f(row,col);
                
                if (val != 0.)
                {
                    i_.push_back(row);
                    j_.push_back(col);
                    x_.push_back(val);
                }
            }
        }
        
        sorted_ = false;
        
        return *this;
    }
    
    // Assignment
    CooMatrix & operator = (CooMatrix const & A)
    {
        m_ = A.m_;
        n_ = A.n_;
        i_ = A.i_;
        j_ = A.j_;
        x_ = A.x_;
        
        sorted_ = A.sorted_;
        
        return *this;
    }
    
    /**
     * @brief Addition of an element to matrix.
     * 
     * Adds a new element to the matrix. If an existing coordinates are used,
     * the numbers will be summed.
     */
    void add (std::int64_t i, std::int64_t j, Complex v)
    {
        i_.push_back(i);
        j_.push_back(j);
        x_.push_back(v);
        
        sorted_ = false;
    }
    
    /// Transposition, implemented as an interchange of "i" and "j" data.
    CooMatrix transpose () const
    {
        CooMatrix tr;
        
        tr.m_ = n_;
        tr.n_ = m_;
        tr.i_ = j_;
        tr.j_ = i_;
        tr.x_ = x_;
        
        tr.sorted_ = false;
        
        return tr;
    }
    
    /// Addition.
    CooMatrix& operator += (CooMatrix const & A)
    {
        assert(m_ == A.m_);
        assert(n_ == A.n_);
        
        i_.append(A.i_.begin(), A.i_.end());
        j_.append(A.j_.begin(), A.j_.end());
        x_.append(A.x_.begin(), A.x_.end());
        
        sorted_ = false;
        
        return *this;
    }
    
    /// Subtraction.
    CooMatrix& operator -= (CooMatrix const & A)
    {
        assert(m_ == A.m_);
        assert(n_ == A.n_);
        
        size_t prev_size = x_.size();
        
        i_.append(A.i_.begin(), A.i_.end());
        j_.append(A.j_.begin(), A.j_.end());
        x_.append(A.x_.begin(), A.x_.end());
        
        // negate the newly added elements
        for (std::size_t i = prev_size; i < x_.size(); i++)
            x_[i] = -x_[i];
        
        sorted_ = false;
        
        return *this;
    }
    
    /// Element-wise multiplication by complex number.
    CooMatrix& operator *= (Complex c)
    {
        std::size_t nz = i_.size();
        for (std::size_t i = 0; i < nz; i++)
            x_[i] *= c;
        
        return *this;
    }
    
    /// SpMV multiplication.
    CooMatrix dot (const cArrayView B) const;
    
    /**
     * @brief Double inner matrix-matrix product.
     * 
     * Double inner matrix-matrix product, \f$ A : B \f$.
     * @note Works only on sorted data.
     * @param B Other matrix.
     */
    Complex ddot (CooMatrix const & B) const;
    
    /**
     * @brief SpMV multiplication.
     * @param B Dense column-major ordered 1D-array. It's row count
     *          is assumed to be equal to the column count of *this
     *          and the column count is computed from the size of the
     *          array, which must be integer multiple of *this's column
     *          count.
     */
    CooMatrix& operator *= (const cArrayView B);
    
    /// Element-wise divide by a complex number.
    CooMatrix& operator /= (Complex c)
    {
        return *this *= 1./c;
    }
    
    /**
     * @brief Change dimension of the matrix.
     * @warning No row/column index range checking.
     */
    void resize (std::size_t m, std::size_t n)
    {
        m_ = m;
        n_ = n;
    }
    
    /**
     * @brief Change matrix shape.
     * 
     * Conserves volume, i.e. it holds
     * @f[
     *                m \cdot n = m_0 \cdot n_0
     * @f]
     * @param m New row count.
     * @param n New column coount.
     */
    CooMatrix reshape (std::size_t m, std::size_t n) const;
    
    /// Convert matrix to dense column-major ordered 1D-array.
    cArray todense () const;
    
    /// Sort indices (by i_, then by j_)
    void sort ();
    bool sorted () const { return sorted_; }
    
    /// Convert to CSC matrix.
    CscMatrix tocsc () const;
    
    /// Convert to CSR matrix.
    CsrMatrix tocsr () const;
    
    /// Convert to dense matrix of a given underlying type.
    template <typename DenseMatrixType> DenseMatrixType todense () const
    {
        DenseMatrixType M (rows(), cols());
        for (std::size_t idx = 0; idx < x_.size(); idx++)
            M (i_[idx], j_[idx]) = x_[idx];
        return M;
    }
    
    /// Convert to dense matrix (row-ordered).
    RowMatrix<Complex> torow () const
    {
        return todense<RowMatrix<Complex>>();
    }
    
    /// Convert to dense matrix (column-ordered).
    ColMatrix<Complex> tocol () const
    {
        return todense<ColMatrix<Complex>>();
    }
    
    /**
     * @brief Solve matrix equation.
     * 
     * Solve the Ax = b problem, where "b" can be a matrix.
     * @param b Complex vector containing column-major ordered data; it may be
     *          a flattened matrix.
     * @param eqs Number of columns.
     * @return Array of roots in the same shape as "b".
     */
    cArray solve (const cArrayView b, std::size_t eqs = 1) const
    {
        // COO format is not optimal for solving -> covert to CSC
        return tocsr().solve(b, eqs);
    }
    
    /**
     * @brief Write the matrix data to a file.
     * 
     * @param filename Filename.
     * 
     * Format of the fields i,j,x,z is:
     * @code
     * %d\t%d\t%g\t%g
     * @endcode
     */
    void write (const char* filename) const;
    
    /// Shake the content, i.e. sum same element entries.
    CooMatrix shake () const;
    
    /**
     * @brief Save matrix to HDF file.
     * @param name Filename.
     */
    bool hdfsave (const char* name) const;
    
    /**
     * @brief Load matrix from HDF file.
     * @param name Filename.
     */
    bool hdfload (const char* name);

private:

    // dimensions
    std::size_t m_, n_;

    // ijv-representation
    lArray i_, j_;
    cArray x_;
    
    bool sorted_;
};

/**
 * Computes a sum of two coo-matrices.
 * @param A COO-matrix.
 * @param B COO-matrix.
 */
inline CooMatrix operator + (CooMatrix const & A, CooMatrix const & B)
{
    CooMatrix C = A;
    return C += B;
}

/**
 * Computes a difference of two coo-matrices.
 * @param A COO-matrix.
 * @param B COO-matrix.
 */
inline CooMatrix operator - (CooMatrix const & A, CooMatrix const & B)
{
    CooMatrix C = A;
    return C -= B;
}

/**
 * Computes product of a matrix and a complex number.
 * @param z Number.
 * @param B Matrix.
 */
inline CooMatrix operator * (Complex z, CooMatrix const & B)
{
    CooMatrix C = B;
    return C *= z;
}

/**
 * Computes product of a sparse matrix and a dense matrix.
 * @param A COO-Matrix.
 * @param B Dense matrix (as a column-major ordered array).
 */
inline CooMatrix operator * (CooMatrix const & A, const cArrayView B)
{
    CooMatrix C = A;
    return C *= B;
}

/**
 * Computes product of a matrix and a complex number.
 * @param A COO-Matrix.
 * @param z Number.
 */
inline CooMatrix operator * (CooMatrix const & A, Complex const & z)
{
    CooMatrix C = A;
    return C *= z;
}

/**
 * Identity matrix.
 * @param N Dimension.
 */
CooMatrix eye (std::size_t N);

/**
 * Diagonal matrix with diagonal consisting of 0, 1, 2, 3, ..., N - 1.
 * @param N Dimension.
 */
CooMatrix stairs (std::size_t N);

/**
 * @brief Kronecker product.
 * 
 * This routine computes the Kronecker product (“flattened tensor product”)
 * @f[
 *     C = A \otimes B \ ,
 * @f]
 * where the matrix @f$ C @f$ is of size @f$ mn \times rs @f$ (for matrix @f$ A @f$
 * of size @f$ m \times r @f$ and matrix @f$ B @f$ of size @f$ n \times s @f$)
 * and defined by
 * @f[
 *     c_{i,j} = a_{i \,\mathrm{div}\, m, j \,\mathrm{div}\, r} \cdot b_{i \,\mathrm{mod}\, m, j \,\mathrm{mod}\, r}
 * @f]
 * 
 * @param A First matrix
 * @param B Second matrix.
 */
CooMatrix kron (CooMatrix const & A, CooMatrix const & B);

#endif // HEX_COOMATRIX_H
