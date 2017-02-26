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

#ifndef HEX_COOMATRIX_H
#define HEX_COOMATRIX_H

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"
#include "hex-csrmatrix.h"
#include "hex-matrix.h"

// --------------------------------------------------------------------------------- //

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
template <class IdxT, class DataT> class CooMatrix
{
    
public:
    
    // Empty constructors
    
    CooMatrix ()
        : m_(0), n_(0), sorted_(true) {}
    CooMatrix (IdxT m, IdxT n)
        : m_(m), n_(n), sorted_(true) {}
    CooMatrix (CooMatrix const & A)
        : m_(A.m_), n_(A.n_), i_(A.i_), j_(A.j_), x_(A.x_), sorted_(false) {}
    CooMatrix (IdxT m, IdxT n, ArrayView<IdxT> i, ArrayView<IdxT> j, ArrayView<DataT> x)
        : m_(m), n_(n), i_(i), j_(j), x_(x), sorted_(false) {}
    CooMatrix (std::size_t m, std::size_t n, NumberArray<IdxT> && i, NumberArray<IdxT> && j, NumberArray<DataT> && x)
        : m_(m), n_(n), i_(std::move(i)), j_(std::move(j)), x_(std::move(x)), sorted_(false) {}
    
    /**
     * Copy constructor initialized from dense array.
     * @param m Row counz of new matrix.
     * @param n Column count of new matrix.
     * @param a Column-major ordered dense array with matrix elements.
     *          Only nonzero elements are copied into internal storage.
     */
    template <class T> CooMatrix (IdxT m, IdxT n, T a) : m_(m), n_(n), sorted_(false)
    {
        // initialize from column-major formatted input
        IdxT i = 0;
        for (IdxT col = 0; col < n; col++)
        {
            for (IdxT row = 0; row < m; row++)
            {
                // get element from array
                DataT val = *(a + i);
                
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
    
    // Getters
    
    std::size_t rows () const { return m_; }
    std::size_t cols () const { return n_; }
    std::size_t size () const { return i_.size(); }
    NumberArray<IdxT> const & i () const { return i_; }
    NumberArray<IdxT> const & j () const { return j_; }
    NumberArray<DataT> const & v () const { return x_; }
    
    /// Index operator. Returns the existing value or zero.
    DataT operator() (IdxT ix, IdxT iy) const
    {
        for (std::size_t n = 0; n < i_.size(); n++)
             if (i_[n] == ix and j_[n] == iy)
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
    template <class Functor> CooMatrix& symm_populate_band (IdxT d, Functor f)
    {
        Complex val;
        
        for (IdxT row = 0; row < m_; row++)
        {
            for (IdxT col = row; col < n_ and col - row <= d; col++)
            {
                val = f(row,col);
                
                if (val != 0.0_r)
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
        
        for (IdxT row = 0; row < m_; row++)
        {
            for (IdxT col = 0; col < n_; col++)
            {
                val = f(row,col);
                
                if (val != 0.0_r)
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
    void add (IdxT i, IdxT j, DataT v)
    {
        i_.push_back(i);
        j_.push_back(j);
        x_.push_back(v);
        
        sorted_ = false;
    }
    
    /// Transposition, implemented as an interchange of "i" and "j" data.
    CooMatrix<IdxT,DataT> transpose () const
    {
        CooMatrix<IdxT,DataT> tr;
        
        tr.m_ = n_;
        tr.n_ = m_;
        tr.i_ = j_;
        tr.j_ = i_;
        tr.x_ = x_;
        
        tr.sorted_ = false;
        
        return tr;
    }
    
    /// Addition.
    CooMatrix<IdxT,DataT>& operator += (CooMatrix<IdxT,DataT> const & A)
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
    CooMatrix<IdxT,DataT>& operator -= (CooMatrix<IdxT,DataT> const & A)
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
    CooMatrix<IdxT,DataT>& operator *= (DataT c)
    {
        std::size_t nz = i_.size();
        for (std::size_t i = 0; i < nz; i++)
            x_[i] *= c;
        
        return *this;
    }
    
    /// Inplace SpMV multiplication.
    void dot (DataT a, const ArrayView<DataT> x, DataT b, ArrayView<DataT> y) const
    {
        assert(x.size() == (unsigned)n_);
        assert(y.size() == (unsigned)m_);
        
        for (IdxT i = 0; i < (IdxT)x_.size(); i++)
        {
            IdxT row = i_[i];
            IdxT col = j_[i];
            
            y[row] = a * x_[i] * x[col] + b * y[row];
        }
    }
    
    /**
     * @brief Double inner matrix-matrix product.
     * 
     * Double inner matrix-matrix product, \f$ A : B \f$.
     * @note Works only on sorted data.
     * @param B Other matrix.
     */
    DataT ddot (CooMatrix<IdxT,DataT> const & B) const
    {
        assert(m_ == B.m_);
        assert(n_ == B.n_);
        
        // sort by i_ and j_
        if (not sorted() or not B.sorted())
            HexException("[CooMatrix] Sort matrices before ddot!");
            
        DataT result = 0;
        
        auto Ai = i_.begin();
        auto Aj = j_.begin();
        auto Av = x_.begin();
        
        auto Bi = B.i_.begin();
        auto Bj = B.j_.begin();
        auto Bv = x_.begin();
        
        while (Av != x_.end() and Bv != B.x_.end())
        {
            if (*Ai < *Bi or (*Ai == *Bi and *Aj < *Bj))
            {
                Ai++; Aj++; Av++;
            }
            else if (*Ai > *Bi or (*Ai == *Bi and *Aj > *Bj))
            {
                Bi++; Bj++; Bv++;
            }
            else // (*Ai == *Bi and *Aj == *Bj)
            {
                result += (*Av) * (*Bv);
                Ai++; Aj++; Av++;
                Bi++; Bj++; Bv++;
            }
        }
        
        return result;
    }
    
    /**
     * @brief SpMV multiplication.
     * @param B Dense column-major ordered 1D-array. It's row count
     *          is assumed to be equal to the column count of *this
     *          and the column count is computed from the size of the
     *          array, which must be integer multiple of *this's column
     *          count.
     */
    CooMatrix<IdxT,DataT>& operator *= (const ArrayView<DataT> B)
    {
        return *this = this->dot(B);
    }
    
    /// Element-wise divide by a complex number.
    CooMatrix<IdxT,DataT>& operator /= (DataT c)
    {
        return *this *= 1./c;
    }
    
    /**
     * @brief Change dimension of the matrix.
     * @warning No row/column index range checking.
     */
    void resize (IdxT m, IdxT n)
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
    CooMatrix<IdxT,DataT> reshape (IdxT m, IdxT n) const
    {
        CooMatrix<IdxT,DataT> C = *this;
        
        // conserved dimensions
        IdxT N = C.i_.size();
        IdxT H = C.m_;
        
        // check dimensions
        assert(m * n == C.m_ * C.n_);
        
        // reshape
        for (IdxT i = 0; i < N; i++)
        {
            // conserved position in column-ordered array
            IdxT idx = C.i_[i] + C.j_[i] * H;
            
            // new coordinates
            IdxT row = idx % m;
            IdxT col = idx / m;
            
            // update values
            C.i_[i] = row;
            C.j_[i] = col;
        }
        
        C.m_ = m;
        C.n_ = n;
        return C;
    }
    
    /// Convert matrix to dense column-major ordered 1D-array.
    NumberArray<DataT> todense () const
    {
        NumberArray<DataT> v (m_ * n_);
        
        IdxT N = i_.size();
        for (IdxT i = 0; i < N; i++)
            v[i_[i] + j_[i] * m_] += x_[i];
        
        return v;
    }
    
    /// Sort indices (by i_, then by j_)
    void sort ();
    bool sorted () const { return sorted_; }
    
    /// Convert to CSC matrix.
    CscMatrix<IdxT,DataT> tocsc () const;
    
    /// Convert to CSR matrix.
    CsrMatrix<IdxT,DataT> tocsr () const;
    
    /// Convert to dense matrix of a given underlying type.
    template <typename DenseMatrixType> DenseMatrixType todense () const
    {
        DenseMatrixType M (rows(), cols());
        for (std::size_t idx = 0; idx < x_.size(); idx++)
            M (i_[idx], j_[idx]) = x_[idx];
        return M;
    }
    
    /// Convert to dense matrix (row-ordered).
    RowMatrix<DataT> torow () const
    {
        return todense<RowMatrix<DataT>>();
    }
    
    /// Convert to dense matrix (column-ordered).
    ColMatrix<DataT> tocol () const
    {
        return todense<ColMatrix<DataT>>();
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
    NumberArray<DataT> solve (const ArrayView<DataT> b, std::size_t eqs = 1) const
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
    void write (std::string const & filename) const
    {
        std::ofstream f (filename);
        
        f << "# Matrix " << m_ << " × " << n_ << " with " << x_.size() << " nonzero elements:\n\n";
        
        for (IdxT i = 0; i < (IdxT)i_.size(); i++)
        {
            f << i_[i] << "\t" << j_[i];
            for (unsigned icmp = 0; icmp < typeinfo<DataT>::ncmpt; icmp++)
                f << "\t" << typeinfo<DataT>::cmpt(icmp, x_[i]);
            f << "\n";
        }
    }
    
    /// Shake the content, i.e. sum same element entries.
    CooMatrix<IdxT,DataT> shake () const
    {
        // ugly and memory inefficient method... FIXME
        return tocsr().tocoo();
    }
    
    /**
     * @brief Save matrix to HDF file.
     * @param name Filename.
     */
    bool hdfsave (std::string const & name) const
    {
        HDFFile hdf (name, HDFFile::overwrite);
        
        // write dimensions
        hdf.write("m", &m_, 1);
        hdf.write("n", &n_, 1);
        
        // write indices
        if (not i_.empty())
            hdf.write("i", &(i_[0]), i_.size());
        if (not j_.empty())
            hdf.write("j", &(j_[0]), j_.size());
        
        // write data
        if (not x_.empty())
        {
            hdf.write
            (
                "x",
                reinterpret_cast<typename typeinfo<DataT>::cmpttype const*>(&x_),
                x_.size() * typeinfo<DataT>::ncmpt
            );
        }

        return true;
    }
    
    /**
     * @brief Load matrix from HDF file.
     * @param name Filename.
     */
    bool hdfload (std::string const & name)
    {
        sorted_ = false;
        
        HDFFile hdf (name, HDFFile::readonly);
        
        // read dimensions
        hdf.read("m", &m_, 1);
        hdf.read("n", &n_, 1);
        
        // read indices
        if (i_.resize(hdf.size("i")))
            hdf.read("i", &(i_[0]), i_.size());
        if (j_.resize(hdf.size("j")))
            hdf.read("j", &(j_[0]), j_.size());
        
        // read data
        if (x_.resize(hdf.size("x") / 2))
        {
            hdf.read
            (
                "x",
                reinterpret_cast<typename typeinfo<DataT>::cmpttype*>(&(x_[0])),
                x_.size() * typeinfo<DataT>::ncmpt
            );
        }
        
        return true;
    }

private:

    // dimensions
    IdxT m_, n_;
    
    // ijv-representation
    NumberArray<IdxT> i_, j_;
    NumberArray<DataT> x_;
    
    bool sorted_;
};

/**
 * Computes a sum of two coo-matrices.
 * @param A COO-matrix.
 * @param B COO-matrix.
 */
template <class IdxT, class DataT>
CooMatrix<IdxT,DataT> operator + (CooMatrix<IdxT,DataT> const & A, CooMatrix<IdxT,DataT> const & B)
{
    CooMatrix<IdxT,DataT> C = A;
    return C += B;
}

/**
 * Computes a difference of two coo-matrices.
 * @param A COO-matrix.
 * @param B COO-matrix.
 */
template <class IdxT, class DataT>
CooMatrix<IdxT,DataT> operator - (CooMatrix<IdxT,DataT> const & A, CooMatrix<IdxT,DataT> const & B)
{
    CooMatrix<IdxT,DataT> C = A;
    return C -= B;
}

/**
 * Computes product of a matrix and a complex number.
 * @param z Number.
 * @param B Matrix.
 */
template <class IdxT, class DataT>
CooMatrix<IdxT,DataT> operator * (DataT z, CooMatrix<IdxT,DataT> const & B)
{
    CooMatrix<IdxT,DataT> C = B;
    return C *= z;
}

/**
 * Computes product of a sparse matrix and a dense matrix.
 * @param A COO-Matrix.
 * @param B Dense matrix (as a column-major ordered array).
 */
template <class IdxT, class DataT>
CooMatrix<IdxT,DataT> operator * (CooMatrix<IdxT,DataT> const & A, const ArrayView<DataT> B)
{
    CooMatrix<IdxT,DataT> C = A;
    return C *= B;
}

/**
 * Computes product of a matrix and a complex number.
 * @param A COO-Matrix.
 * @param z Number.
 */
template <class IdxT, class DataT>
CooMatrix<IdxT,DataT> operator * (CooMatrix<IdxT,DataT> const & A, DataT z)
{
    CooMatrix<IdxT,DataT> C = A;
    return C *= z;
}

/**
 * Identity matrix.
 * @param N Dimension.
 */
template <class IdxT, class DataT>
CooMatrix<IdxT,DataT> eye (IdxT N)
{
    return CooMatrix<IdxT,DataT>(N,N).symm_populate_band
    (
        0,
        [](IdxT i, IdxT j) -> DataT { return 1.; }
    );
}

/**
 * Diagonal matrix with diagonal consisting of 0, 1, 2, 3, ..., N - 1.
 * @param N Dimension.
 */
template <class IdxT, class DataT>
CooMatrix<IdxT,DataT> stairs (IdxT N)
{
    return CooMatrix<IdxT,DataT>(N,N).symm_populate_band
    (
        0,
        [](IdxT i, IdxT j) -> DataT { return i; }
    );
}

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
template <class IdxT, class DataT>
CooMatrix<IdxT,DataT> kron (CooMatrix<IdxT,DataT> const & A, CooMatrix<IdxT,DataT> const & B);

// --------------------------------------------------------------------------------- //

#endif // HEX_COOMATRIX_H
