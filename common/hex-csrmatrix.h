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

#ifndef HEX_CSRMATRIX_H
#define HEX_CSRMATRIX_H

// --------------------------------------------------------------------------------- //

#include <memory>

// --------------------------------------------------------------------------------- //

#ifdef WITH_PNG
#include <png++/png.hpp>
#endif

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"
#include "hex-coomatrix.h"
#include "hex-densematrix.h"

// --------------------------------------------------------------------------------- //

/**
 * @brief Complex CSR matrix.
 * 
 * Complex compressed-storage-row-major-ordered sparse matrix.
 * The data are stored in four arrays,
 * @code
 * p, i, x, z
 * @endcode
 * Here, \c Ap is an array of N+1 elements, where N is the number of rows.
 * It contains indices of beginnings of single rows in all remaining arrays.
 * The last element contains one-past-end pointer, i.e. it should be aqual to
 * the element count of the matrix (zero and nonzero). The array \c Ai
 * holds column indices of elements, which are stored in \c Ax and \c Az, real
 * ang imaginary part being separated.
 */
template <class IdxT, class DataT> class CsrMatrix
{

public:

    // Constructors

    CsrMatrix ()
        : m_(0), n_(0), name_() {}
    CsrMatrix (std::size_t m, std::size_t n)
        : m_(m), n_(n), name_() {}
    CsrMatrix (CsrMatrix const & A)
        : m_(A.m_), n_(A.n_), p_(A.p_), i_(A.i_), x_(A.x_), name_() {}
    CsrMatrix (std::size_t m, std::size_t n, const ArrayView<IdxT> p, const ArrayView<IdxT> i, const ArrayView<DataT> x)
        : m_(m), n_(n), p_(p), i_(i), x_(x), name_() {}
    CsrMatrix (std::size_t m, std::size_t n, NumberArray<IdxT> && p, NumberArray<IdxT> && i, NumberArray<DataT> && x)
        : m_(m), n_(n), p_(std::move(p)), i_(std::move(i)), x_(std::move(x)), name_() {}

    /// Destructor
    ~CsrMatrix () {}

    /// Clear data.
    void drop ()
    {
        m_ = n_ = 0;
        p_.drop();
        i_.drop();
        x_.drop();
    }

    /**
     * Ordinary matrix-vector product, \f$ A\cdot b \f$.
     * @param b Vector to multiply with.
     */
    NumberArray<DataT> dot (const ArrayView<DataT> b) const
    {
        // create output array
        NumberArray<DataT> c(m_);

        for (IdxT irow = 0; irow < m_; irow++)
        {
            IdxT idx1 = p_[irow];
            IdxT idx2 = p_[irow+1];

            // for all nonzero elements in this row
            for (IdxT idx = idx1; idx < idx2; idx++)
            {
                // get column number
                IdxT icol = i_[idx];

                // store product
                c[irow] += x_[idx] * b[icol];
            }
        }

        return c;
    }

    /// Get number of structurally non-zero elements.
    std::size_t size () const { return i_.size(); }

    /// Get row count.
    std::size_t rows () const { return m_; }

    /// Get column count.
    std::size_t cols () const { return n_; }

    /// Get row pointers.
    NumberArray<IdxT> const & p () const { return p_; }
    NumberArray<IdxT> & p () { return p_; }

    /// Get column indices.
    NumberArray<IdxT> const & i () const { return i_; }
    NumberArray<IdxT> & i () { return i_; }

    /// Get data array.
    NumberArray<DataT> const & x () const { return x_; }
    NumberArray<DataT> & x () { return x_; }

    /// Return absolute value of the (in absolute value) largest element.
    double norm () const
    {
        std::size_t N = i_.size();
        double res = 0.;

        // return the abs(largest element)
        for (std::size_t i = 0; i < N; i++)
        {
            // compute the absolute value
            double val = std::abs(x_[i]);

            // update the winner
            if (val > res)
                res = val;
        }

        return res;
    }

    /**
     * Write the matrix data to a file.
     * 
     * @param filename Filename.
     * 
     * Format of the fields i,j,x,z is:
     * @code
     * %d\t%d\t%g\t%g
     * @endcode
     */
    void write (std::string filename) const
    {
        std::ofstream out (filename);
        out << "# Matrix " << m_ << " × " << n_ << " with " << x_.size() << " nonzero elements:\n\n";
        for (unsigned irow = 0; irow < m_; irow++)
        {
            IdxT idx1 = p_[irow];
            IdxT idx2 = p_[irow + 1];

            for (IdxT idx = idx1; idx < idx2; idx++)
            {
                out << irow << '\t' << i_[idx];
                for (unsigned icomp = 0; icomp < typeinfo<DataT>::ncmpt; icomp++)
                    out << '\t' << typeinfo<DataT>::cmpt(icomp, x_[idx]);
                out << '\n';
            }
        }
        out.close();
    }

#ifdef WITH_PNG
    /**
     * PNG row data generator for use in \ref plot function.
     * 
     * @note This class requires PNG++.
     */
    class PngGenerator : public png::generator<png::gray_pixel_1,PngGenerator>
    {
        typedef png::generator<png::gray_pixel_1,PngGenerator> base_t;
        typedef png::packed_pixel_row<png::gray_pixel_1> row;
        typedef png::row_traits<row> row_traits;

        public:

            PngGenerator (CsrMatrix<IdxT,DataT> const * mat, double threshold)
                : base_t(mat->cols(), mat->rows()), M_(mat), buff_(mat->cols()), threshold_(threshold)
            {
                // do nothing
            }

            ~PngGenerator ()
            {
                // do nothing
            }

            png::byte * get_next_row (std::size_t irow)
            {
                // get column indices
                IdxT idx_min = M_->p_[irow];
                IdxT idx_max = M_->p_[irow + 1];

                // clear memory
                for (std::size_t icol = 0; icol < M_->cols(); icol++)
                    buff_[icol] = 1;

                // for all nonzero columns
                for (IdxT idx = idx_min; idx < idx_max; idx++)
                    if (std::abs(M_->x_[idx]) > threshold_)
                        buff_[M_->i_[idx]] = 0;

                // pass the buffer
                return reinterpret_cast<png::byte*>(row_traits::get_data(buff_));
            }

        private:

            CsrMatrix<IdxT,DataT> const * M_;
            row buff_;
            double threshold_;
    };
#endif // WITH_PNG

    /**
     * Save matrix structure as a black-and-white image.
     * @param filename File name.
     * @param threshold Largest absolute value represented by white colour.
     */
    void plot (std::string const & filename, double threshold = 0.) const
    {
#ifdef WITH_PNG
        // create output file
        std::ofstream out (filename, std::ios_base::out | std::ios_base::binary);

        // create PNG data generator
        PngGenerator png (this, threshold);

        // write PNG file
        png.write(out);
#else
        HexException("The program was not built with PNG++ library (use -DWITH_PNG).");
#endif // WITH_PNG
    }

    /**
     * @brief Link to a disk file.
     * 
     * Set a default I/O file that will be used if any of the functions
     * @ref hdfload or @ref hdfsave will be used without an explicit filename.
     */
    void hdflink (std::string name) { name_ = name; }
    void unlink () { name_.clear(); }

    /**
     * Get linked HDF name.
     */
    std::string hdfname () const { return name_; }

    /**
     * Save matrix to HDF file.
     */
    //@{
    bool hdfsave () const { return hdfsave(name_); }
    bool hdfsave (std::string name) const
    {
        HDFFile hdf(name, HDFFile::overwrite);

        // write dimensions
        hdf.write("m", &m_, 1);
        hdf.write("n", &n_, 1);

        // write indices
        if (not p_.empty())
            hdf.write("p", &(p_[0]), p_.size());
        if (not i_.empty())
            hdf.write("i", &(i_[0]), i_.size());

        // write data
        if (not x_.empty())
        {
            hdf.write
            (
                "x",
                reinterpret_cast<typename typeinfo<DataT>::cmpttype const*>(&(x_[0])),
                x_.size() * typeinfo<DataT>::ncmpt
            );
        }

        return true;
    }
    //@}

    /**
     * Load matrix from HDF file.
     */
    //@{
    bool hdfload () { return hdfload(name_); }
    bool hdfload (std::string name)
    {
        HDFFile hdf(name, HDFFile::readonly);
        if (not hdf.valid()) return false;

        // read dimensions
        if (not hdf.read("m", &m_, 1)) return false;
        if (not hdf.read("n", &n_, 1)) return false;

        // read indices
        if (p_.resize(hdf.size("p")) and not hdf.read("p", &(p_[0]), p_.size()))
            return false;
        if (i_.resize(hdf.size("i")) and not hdf.read("i", &(i_[0]), i_.size()))
            return false;

        // read data
        if
        (
            x_.resize(hdf.size("x") / 2)
            and not hdf.read
            (
                "x",
                reinterpret_cast<typename typeinfo<DataT>::cmpttype*>(&(x_[0])),
                x_.size() * typeinfo<DataT>::ncmpt
            )
        )
            return false;

        return true;
    }
    //@}

    /**
     * Element-wise access (const).
     * If a non-existing element is referenced, zero is returned.
     */
    DataT operator() (IdxT i, IdxT j) const
    {
        // get all column indices, which have nonzero element in row "i"
        IdxT idx1 = p_[i];
        IdxT idx2 = p_[i + 1];

        // find the correct column ("j")
        auto it = std::lower_bound(i_.begin() + idx1, i_.begin() + idx2, j);

        if (it == i_.end() or *it != j)
            return 0.;

        // return the value
        return x_[it - i_.begin()];
    }

    /**
     * Multiplication by a number.
     */
    CsrMatrix<IdxT,DataT> & operator *= (DataT r)
    {
        std::size_t N = i_.size();
        for (std::size_t i = 0; i < N; i++)
            x_[i] *= r;
        return *this;
    }

    /**
     * Addition of another CSR matrix.
     * The matrices MUST HAVE THE SAME SPARSE STRUCTURE, as no indices are
     * checked.
     */
    CsrMatrix<IdxT,DataT> & operator &= (CsrMatrix<IdxT,DataT> const & B)
    {
        std::size_t N = i_.size();

        // check at least dimensions and non-zero element count
        assert(m_ == B.m_);
        assert(n_ == B.n_);
        assert(N == B.i_.size());

        for (std::size_t i = 0; i < N; i++)
            x_[i] += B.x_[i];

        return *this;
    }

    /**
     * Subtraction of another CSR matrix.
     * The matrices MUST HAVE THE SAME SPARSE STRUCTURE, as no indices are
     * checked.
     */
    CsrMatrix<IdxT,DataT> & operator ^= (CsrMatrix<IdxT,DataT> const & B)
    {
        std::size_t N = i_.size();

        assert(m_ == B.m_);
        assert(n_ == B.n_);
        assert(N == B.i_.size());

        for (std::size_t i = 0; i < N; i++)
            x_[i] -= B.x_[i];

        return *this;
    }

    /**
     * Sets fill-in elements so that the storage structure of this matrix
     * will be identical to that of the other CSR matrix. We can the use
     * “fast” arithmetic operators & and ^.
     * @return self
     */
    CsrMatrix<IdxT,DataT> sparse_like (CsrMatrix<IdxT,DataT> const & B) const
    {
        // check dimensions
        assert(m_ == B.m_);
        assert(n_ == B.n_);

        // prepare zero matrix with the same storage pattern the matrix B has
        CsrMatrix<IdxT,DataT> A = B;
        std::memset(A.x_.data(), 0, A.x_.size() * sizeof(DataT));

        // copy all nonzero elements of "this" matrix
        for (IdxT row = 0; row < m_; row++)
        {
            IdxT idx1 = A.p_[row];
            IdxT idx2 = A.p_[row+1];

            for (std::size_t idx = idx1; idx < idx2; idx++)
            {
                IdxT col = A.i_[idx];
                A.x_[idx] = (*this)(row,col);
            }
        }

        // return temporary matrix
        return A;
    }

    /**
     * Return dense array with diagonal elements of the matrix.
     */
    NumberArray<DataT> diag () const
    {
        NumberArray<DataT> D (std::min(m_, n_));

        for (int irow = 0; irow < m_; irow++)
        for (int idx = p_[irow]; idx < p_[irow+1]; idx++)
        if (i_[idx] == irow)
            D[irow] = x_[idx];

        return D;
    }

    /**
     * Convert to COO format.
     */
    CooMatrix<IdxT,DataT> tocoo () const
    {
        if (p().empty() or i().empty() or x().empty())
            return CooMatrix<IdxT,DataT>();

        IdxT nz = p().back();

        NumberArray<IdxT>  I; I.reserve(nz);
        NumberArray<IdxT>  J; J.reserve(nz);
        NumberArray<DataT> V; V.reserve(nz);

        for (IdxT irow = 0; irow < (IdxT)rows(); irow++)
        for (IdxT idx = p()[irow]; idx < p()[irow + 1]; idx++)
        {
            I.push_back(irow);
            J.push_back(i()[idx]);
            V.push_back(x()[idx]);
        }

        return CooMatrix<IdxT,DataT>
        (
            rows(),
            cols(),
            std::move(I),
            std::move(J),
            std::move(V)
        );
    }

    /**
     * Convert to dense matrix.
     */
    RowMatrix<DataT> torow () const
    {
        RowMatrix<DataT> M (rows(),cols());
        for (IdxT irow = 0; irow < rows(); irow++)
        {
            IdxT rptr_begin = p_[irow];
            IdxT rptr_end = p_[irow + 1];

            for (IdxT idx = rptr_begin; idx < rptr_end; idx++)
            {
                IdxT icol = i_[idx];
                Complex x = x_[idx];

                M(irow, icol) = x;
            }
        }

        return M;
    }

    /**
     * @brief Upper triangular solve.
     * 
     * Solves upper triangular system of equations using backsubstitution.
     * The matrix needs not be triangular, but elements under the diagonal
     * won't be used or changed.
     * 
     * @param b Right-hand side.
     * @param omega Over-relaxation factor.
     */
    NumberArray<DataT> upperSolve (ArrayView<DataT> const & b, DataT omega = 1) const
    {
        // check size
        IdxT N = b.size();
        assert(m_ == N);
        assert(n_ == N);

        // create output array
        NumberArray<DataT> x (N);

        // loop over rows
        for (IdxT i = 0; i < N; i++)
        {
            IdxT row = N - 1 - i;
            DataT accum = 0.;

            // get relevant columns of the sparse matrix
            IdxT idx1 = p_[row];
            IdxT idx2 = p_[row + 1];

            // diagonal element of the matrix
            DataT a = 0.;

            // loop over the columns
            for (IdxT idx = idx1; idx < idx2; idx++)
            {
                // which column is this?
                IdxT col = i_[idx];

                // diagonal element will be useful in a moment, store it
                if (col == row)
                    a = x_[idx] / omega;

                // backsubstitute
                else if (col > row)
                    accum += x_[idx] * x[col];
            }

            // triangular matrix, in order to be regular, needs nonzero diagonal elements
            assert(a != 0.);

            // compute and store the new root
            x[row] = (b[row] - accum) / a;
        }

        return x;
    }

    /**
     * @brief Lower triangular solve.
     * 
     * Solves lower triangular system of equations using backsubstitution.
     * The matrix needs not be triangular, but elements above the diagonal
     * won't be used or changed.
     * 
     * @param b Right-hand side.
     * @param omega Over-relaxation factor.
     */
    NumberArray<DataT> lowerSolve (ArrayView<DataT> const & b, DataT omega = 1) const
    {
        // check size
        IdxT N = b.size();
        assert(m_ == N);
        assert(n_ == N);

        // create output array
        NumberArray<DataT> x (N);

        // loop over rows
        for (IdxT row = 0; row < N; row++)
        {
            DataT accum = 0.;

            // get relevant columns of the sparse matrix
            IdxT idx1 = p_[row];
            IdxT idx2 = p_[row + 1];

            // diagonal element of the matrix
            DataT a = 0.;

            // loop over the columns
            for (IdxT idx = idx1; idx < idx2; idx++)
            {
                // which column is this?
                IdxT col = i_[idx];

                // diagonal element will be useful in a moment, store it
                if (col == row)
                    a = x_[idx] / omega;

                // backsubstitute
                else if (col < row)
                    accum += x_[idx] * x[col];
            }

            // triangular matrix, in order to be regular, needs nonzero diagonal elements
            assert(a != 0.);

            // compute and store the new root
            x[row] = (b[row] - accum) / a;
        }

        return x;
    }

    /**
     * Applies a user transformation on <b>nonzero</b> matrix elements.
     * @param f A functor compatible with following declaration:
     * @code
     * Complex (*f) (size_t i, size_t j, Complex z);
     * @endcode
     */
    template <class Functor> CsrMatrix nzTransform (Functor f) const
    {
        // create output matrix
        CsrMatrix A = *this;

        // loop over rows
        for (size_t row = 0; row < rows(); row++)
        {
            // get relevant columns
            size_t idx1 = p_[row];
            size_t idx2 = p_[row + 1];

            // loop over columns
            for (size_t idx = idx1; idx < idx2; idx++)
            {
                // which column is this?
                size_t col = i_[idx];

                // transform the output matrix
                A.x_[idx] = f(row, col, x_[idx]);
            }
        }

        return A;
    }

private:

    // dimensions
    IdxT m_;
    IdxT n_;

    // representation
    NumberArray<IdxT> p_;
    NumberArray<IdxT> i_;
    NumberArray<DataT> x_;

    // linked HDF file
    std::string name_;
};

/**
 * Computes a sum of two csr-matrices OF THE SAME SPARSE STRUCTURE.
 */
template <class IdxT, class DataT>
CsrMatrix<IdxT,DataT> operator & (CsrMatrix<IdxT,DataT> const & A, CsrMatrix<IdxT,DataT> const & B)
{
    CsrMatrix<IdxT,DataT> C = A;
    return C &= B;
}

/**
 * Computes a difference of two csr-matrices OF THE SAME SPARSE STRUCTURE.
 */
template <class IdxT, class DataT>
CsrMatrix<IdxT,DataT> operator ^ (CsrMatrix<IdxT,DataT> const & A, CsrMatrix<IdxT,DataT> const & B)
{
    CsrMatrix<IdxT,DataT> C = A;
    return C ^= B;
}

/**
 * Multiplication of csr-matrix by a number.
 */
template <class IdxT, class DataT>
CsrMatrix<IdxT,DataT> operator * (DataT z, CsrMatrix<IdxT,DataT> const &  B)
{
    CsrMatrix<IdxT,DataT> C = B;
    return C *= z;
}

/**
 * Multiplication of csr-matrix by a number.
 */
template <class IdxT, class DataT>
CsrMatrix<IdxT,DataT> operator * (CsrMatrix<IdxT,DataT> const & A, double r)
{
    CsrMatrix<IdxT,DataT> C = A;
    return C *= r;
}

/**
 * @brief Kronecker product.
 * 
 * Calculates the Kronecker product @f$ w = (A \otimes B) \cdot u @f$, or in
 * components
 * @f[
 *     w_{i r_B + j} = A_{ik} B_{jl} u_{k c_B + l}
 * @f]
 * where @f$ u @f$ is a column vector and the matrices @f$ A @f$ and @f$ B @f$
 * are of the CSR type. The dimensions of the matrices are @f$ r_A \times c_A @f$
 * and @f$ r_B \times c_B @f$, respectively. The result is a column vector, again.
 */
template <class IdxT, class DataT>
NumberArray<DataT> kron_dot (CsrMatrix<IdxT,DataT> const & A, CsrMatrix<IdxT,DataT> const & B, NumberArray<DataT> const & u)
{
    // dimensions
    std::size_t rA = A.rows(), rB = B.rows();
    std::size_t cA = A.cols(), cB = B.cols();

    // check input size
    assert(u.size() == cA * cB);

    // V[col] = B U[col]
    NumberArray<DataT> V (rB * cA);
    for (std::size_t i = 0; i < cA; i++)
        ArrayView<DataT>(V, i * rB, rB) = B.dot(ArrayView<DataT>(u, i * cB, cB));
    transpose(V, rB, cA);

    // W = A V
    NumberArray<DataT> W (rA * rB);
    for (std::size_t i = 0; i < rB; i++)
        ArrayView<DataT>(W, i * rA, rA) = A.dot(ArrayView<DataT>(u, i * cA, cA));
    transpose(W, rA, rB);

    return W;
}

// --------------------------------------------------------------------------------- //

#endif // HEX_CSRMATRIX_H
