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

#ifndef HEX_DENSEMATRIX_H
#define HEX_DENSEMATRIX_H

#ifdef WITH_PNG
#include <png++/png.hpp>
#endif

#include "hex-arrays.h"
#include "hex-blas.h"
#include "hex-matrix.h"
#include "hex-misc.h"

/**
 * @brief DenseMatrix.
 * 
 * Base class for row- and column- oriented matrices. The data
 * are stored in a one-dimensional array. The class has a basic
 * data access interface; however, the main purpose is that it
 * is the base class of the storage-order-specialized classes
 * RowMatrix and ColMatrix. It is recommended to use the latter
 * two when working with dense matrices.
 */
template <class T> class DenseMatrixView
{
    public:
        
        // default constructor
        DenseMatrixView () : rows_(0), cols_(0), data_()
        {
            // do nothing
        }
        
        // size constructor
        DenseMatrixView (std::size_t rows, std::size_t cols) : rows_(rows), cols_(cols), data_()
        {
            // do nothing
        }
        
        // view constructor
        DenseMatrixView (std::size_t rows, std::size_t cols, const ArrayView<T> data) : rows_(rows), cols_(cols), data_(data)
        {
            // check that the sizes match
            assert(data.size() == rows * cols);
        }
        
        // get elements array (read-write access)
        ArrayView<T> data () { return data_; }
        
        // get elements array (read-only access)
        const ArrayView<T> data () const { return data_; }
        
        // throw away all data
        void drop () { rows_ = cols_ = 0; }
        
        // use new data
        virtual void update (std::size_t rows, std::size_t cols, const ArrayView<T> elems)
        {
            rows_ = rows;
            cols_ = cols;
            data_.reset(rows * cols, elems.data());
        }
        
        // get element
        virtual T & operator () (std::size_t i, std::size_t j) = 0;
        virtual T operator () (std::size_t i, std::size_t j) const = 0;
        
        // get number of elements
        std::size_t size () const { return rows_ * cols_; }
        
        // get number of columns
        int cols () const { return cols_; }
        
        // get number of rows
        int rows () const { return rows_; }
        
        // get pointer to the beginning of the elements array (read-write access)
        T * begin () { return data_.begin(); }
        
        // get pointer to the beginning of the elements array (read-only access)
        T const * begin () const { return data_.begin(); }
        
        // layout (row- or column- major)
        virtual char layout () const
        {
            HexException("General DenseMatrix has undefined layout type.");
            return 0;
        }
        
        // leading dimension
        virtual std::size_t ld () const
        {
            HexException("General DenseMatrix has undefined leading dimension.");
            return 0;
        }
        
    protected:
        
        /// Row count.
        int rows_;
        
        /// Column count.
        int cols_;
        
        /// Matrix elements, consecutive rows joined to one array.
        ArrayView<T> data_;
};

/**
 * @brief DenseMatrix.
 * 
 * Base class for row- and column- oriented matrices. The data
 * are stored in a one-dimensional array. The class has a basic
 * data access interface; however, the main purpose is that it
 * is the base class of the storage-order-specialized classes
 * RowMatrix and ColMatrix. It is recommended to use the latter
 * two when working with dense matrices.
 */
template <class T> class DenseMatrix : public DenseMatrixView<T>
{
    public:
        
        enum Layouts
        {
            RowMajorLayout = 'T',
            ColumnMajorLayout = 'N'
        };
        
        // default constructor
        DenseMatrix () : DenseMatrixView<T>(), storage_()
        {
            // does nothing
        }
        
        // size constructor
        DenseMatrix (std::size_t rows, std::size_t cols) : DenseMatrixView<T>(rows, cols), storage_(rows * cols)
        {
            // let the view point to the storage
            DenseMatrixView<T>::data_.reset(rows * cols, storage_.data());
        }
        
        DenseMatrix (std::size_t rows, std::size_t cols, const ArrayView<T> data) : DenseMatrixView<T>(rows, cols), storage_(data)
        {
            // check that the sizes match
            assert(data.size() == rows * cols);
            
            // let the view point to the storage
            DenseMatrixView<T>::data_.reset(rows * cols, storage_.data());
        }
        
        // get matrix elements array (read-write access)
        ArrayView<T> data ()
        {
            return DenseMatrixView<T>::data();
        }
        
        // get matrix elements array (read-only access)
        const ArrayView<T> data () const
        {
            return DenseMatrixView<T>::data();
        }
        
        // use new data
        virtual void update (std::size_t rows, std::size_t cols, const ArrayView<T> elems)
        {
            assert(rows * cols == elems.size());
            storage_ = elems;
            DenseMatrixView<T>::update(rows, cols, storage_);
        }
        
        // use new data
        virtual void update (std::size_t rows, std::size_t cols, NumberArray<T> && elems)
        {
            assert(rows * cols == elems.size());
            storage_ = std::move(elems);
            DenseMatrixView<T>::update(rows, cols, storage_);
        }
        
        // throw away all data
        void drop ()
        {
            DenseMatrixView<T>::drop();
            storage_.drop();
        }
        
        // retrieve size of the matrix (number of elements)
        std::size_t size () const
        {
            return DenseMatrixView<T>::rows() * DenseMatrixView<T>::cols();
        }
        
        // get number of columns
        int cols () const { return DenseMatrixView<T>::cols_; }
        
        // get number of rows
        int rows () const { return DenseMatrixView<T>::rows_; }
        
        // get pointer to the beginning of the elements array (read-write access)
        T * begin () { return DenseMatrixView<T>::data_.begin(); }
        
        // get pointer to the beginning of the elements array (read-only access)
        T const * begin () const { return DenseMatrixView<T>::data_.begin(); }
        
        // layout (row- or column- major)
        virtual char layout () const
        {
            HexException("General DenseMatrix has undefined layout type.");
            return 0;
        }
        
        // leading dimension
        virtual std::size_t ld () const
        {
            HexException("General DenseMatrix has undefined leading dimension.");
            return 0;
        }
        
    protected:
        
        /// Matrix elements, consecutive rows joined to one array.
        NumberArray<T> storage_;
};

/**
 * @brief Dense (column-oriented) matrix.
 * 
 * This class represents a generally non-symmetric dense matrix. It is
 * an encapsulation of the NumberArray class with some operators
 * bundled to it.
 */
template <class Type, class Base> class ColMatrix : public Base
{
    public:
        
        // constructor
        ColMatrix ()
            : Base(), ld_(0) {}
        ColMatrix (int size)
            : Base(size, size), ld_(size) {}
        ColMatrix (int rows, int cols)
            : Base(rows, cols), ld_(rows) {}
        ColMatrix (int rows, int cols, const ArrayView<Type> data, int ld = -1)
            : Base(rows, cols, data), ld_(ld < 0 ? rows : ld) {}
        ColMatrix (ColMatrix<Type> const & m)
            : Base(m.rows(), m.cols(), m.data()), ld_(m.ld()) { }
        explicit ColMatrix (RowMatrix<Type> const & m)
            : Base(m.rows(), m.cols(), m.data()), ld_(m.rows()) { reorder_(); }
        
        /**
         * @brief Assignment operator.
         */
        ColMatrix<Type> & operator = (ColMatrix<Type> const & A)
        {
            this->update(A.rows_, A.cols_, A.data_);
            ld_   = A.ld_;
            return *this;
        }
        
        /**
         * @brief Move assignment.
         */
        ColMatrix<Type> & operator = (ColMatrix<Type> const && A)
        {
            this->update(A.rows_, A.cols_, std::move(A.data_));
            ld_   = A.ld_;
            return *this;
        }
        
        /**
         * @brief Populator.
         * 
         * This method accepts a function object that is called
         * for every matrix element.
         */
        template <class Function> void populate (Function f)
        {
            // pointer to data
            Type * p = this->begin();
            
            // fill the matrix
            for (int icol = 0; icol < this->cols(); icol++)
            for (int irow = 0; irow < this->rows(); irow++)
                *(p + icol * ld_ + irow) = f(irow,icol);
        }
        
        /**
         * @brief Matrix transpose.
         * 
         * Return transposed matrix. The result is of the type RowMatrix,
         * so that no reordering of the entries is necessary.
         */
        RowMatrix<Type> T () const
        {
            return RowMatrix<Type> (this->cols(), this->rows(), this->data(), ld_);
        }
        
        /**
         * @brief Matrix hermitian conjugate.
         * 
         * Return hermitian conjugated matrix. The result is of the type RowMatrix,
         * so that no reordering of the entries is necessary.
         */
        RowMatrix<Type> H () const
        {
            return RowMatrix<Type> (this->cols(), this->rows(), NumberArray<Type>(this->data()).conj(), ld_);
        }
        
        /**
         * @brief Matrix column.
         *
         * Return shallow copy of the chosen matrix column.
         */
        //@{
        ArrayView<Type> col (int i)
        {
            return ArrayView<Type> (this->data(), i * ld_, this->rows());
        }
        const ArrayView<Type> col (int i) const
        {
            return ArrayView<Type> (this->data(), i * ld_, this->rows());
        }
        //@}
        
        /// Element access.
        //@{
        virtual Type operator() (std::size_t i, std::size_t j) const { return col(j)[i]; }
        virtual Type & operator() (std::size_t i, std::size_t j) { return col(j)[i]; }
        //@}
        
        /// Inversion.
        void invert (ColMatrix<Type> & inv) const;
        
        /// Diagonalization.
        void diagonalize
        (
            NumberArray<Type> & eigval,
            ColMatrix<Type> * eigvecL = nullptr,
            ColMatrix<Type> * eigvecR = nullptr
        ) const;
        
        /// Storage layout.
        virtual char layout () const { return DenseMatrix<Type>::ColumnMajorLayout; }
        
        /// Leading dimension.
        virtual std::size_t ld () const { return ld_; }
        
        void write (std::ostream & out, std::string const & pre = "", std::string const & pos = "") const
        {
            // data pointer
            Type const * ptr = this->data().begin();
            
            for (int irow = 0; irow < this->rows(); irow++)
            {
                out << pre;
                for (int icol = 0; icol < this->cols(); icol++)
                {
                    Type x = *(ptr + icol * ld_ + irow);
                    
                    if (x == std::abs(x))
                    {
                        // positive entry
                        out << " " << std::abs(x) << " ";
                    }
                    else
                    {
                        // negative entry
                        out << x << " ";
                    }
                }
                
                out << pos << "\n";
            }
        }
    
    private:
        
        /// Change "data" from row-oriented to column-oriented.
        void reorder_ ()
        {
            NumberArray<Type> new_data (this->rows() * this->cols());
            
            for (int irow = 0; irow < this->rows(); irow++)
            for (int icol = 0; icol < this->cols(); icol++)
                new_data[icol * this->rows() + irow] = this->data_[irow * this->cols() + icol];
            
            this->data_ = new_data;
        }
        
        /// Leading dimension.
        std::size_t ld_;
};

/**
 * @brief Dense (row-oriented) matrix.
 * 
 * This class represents a generally non-symmetric dense matrix. It is
 * an encapsulation of the cArray class with some operators bundled to it.
 */
template <class Type, class Base> class RowMatrix : public Base
{
    public:
        
        // constructor
        RowMatrix ()
            : Base(), ld_(0) {}
        RowMatrix (int size)
            : Base(size, size), ld_(size) {}
        RowMatrix (int rows, int cols)
            : Base(rows, cols), ld_(cols) {}
        RowMatrix (int rows, int cols, const ArrayView<Type> data, int ld = -1)
            : Base(rows, cols, data), ld_(ld < 0 ? cols : ld) {}
        RowMatrix (RowMatrix<Type> const & m)
            : Base(m.rows(), m.cols(), m.data()), ld_(m.ld()) { }
        explicit RowMatrix (ColMatrix<Type> const & m)
            : Base(m.rows(), m.cols(), m.data()), ld_(m.cols()) { reorder_(); }
        
        /**
         * @brief Assignment operator.
         */
        RowMatrix<Type> & operator = (RowMatrix<Type> const & A)
        {
            this->update(A.rows_, A.cols_, A.data_);
            ld_  = A.ld_;
            return *this;
        }
        
        /**
         * @brief Move assignment.
         */
        RowMatrix<Type> & operator = (RowMatrix<Type> const && A)
        {
            this->update(A.rows_, A.cols_, std::move(A.data_));
            ld_  = A.ld_;
            return *this;
        }
        
        /**
         * @brief Populator.
         * 
         * This method accepts a function object that is called
         * for every matrix element. Use functions with signature
           @code
               double (*) (int i, int j);
           @endcode
         */
        template <class Function> void populate (Function f)
        {
            // pointer to data
            Type * p = this->begin();
            
            // fill the matrix
            for (int irow = 0; irow < this->rows(); irow++)
            for (int icol = 0; icol < this->cols(); icol++)
                *(p + irow * ld_ + icol) = f(irow,icol);
        }
        
        /**
         * @brief Transformator.
         * 
         * This method accepts a function object that is called
         * for every matrix element. Use functions with signature
           @code
               void (*) (int i, int j, T& x);
           @endcode
         */
        template <class Function> void transform (Function f)
        {
            // pointer to data
            Type * p = this->begin();
            
            // transform the matrix
            for (int irow = 0; irow < this->rows(); irow++)
            for (int icol = 0; icol < this->cols(); icol++)
                f(irow,icol,*(p + irow * ld_ + icol));
        }
        
        /**
         * @brief Matrix transpose.
         * 
         * Return transposed matrix. The result is of the type ColMatrix,
         * so that no reordering of the entries is necessary.
         */
        ColMatrix<Type> T () const
        {
            return ColMatrix<Type> (this->cols(), this->rows(), this->data(), ld_);
        }
        
        /**
         * @brief Matrix row.
         * 
         * Return shallow copy of a matrix row.
         */
        //@{
        ArrayView<Type> row (int i)
        {
            return ArrayView<Type> (this->data(), i * ld_, this->cols());
        }
        const ArrayView<Type> row (int i) const
        {
            return ArrayView<Type> (this->data(), i * ld_, this->cols());
        }
        //@}
        
        /// Element access.
        //@{
        virtual Type operator() (std::size_t i, std::size_t j) const { return row(i)[j]; }
        virtual Type & operator() (std::size_t i, std::size_t j) { return row(i)[j]; }
        //@}
        
        /**
         * @brief Bilinear form.
         * 
         * This function computes the scalar product
         * @f[
         *     u^T \mathsf{A} v \ .
         * @f]
         */
        Type operator() (const ArrayView<Type> u, const ArrayView<Type> v) const
        {
            // return value
            Type sum = 0;
            
            // pointer to matrix elements
            Type const * restrict pA = this->data().data();
            
            // pointers to vector elements
            Type const * const restrict pu = u.data();
            Type const * const restrict pv = v.data();
            
            // for all elements
            for (int i = 0; i < this->rows(); i++)
            for (int j = 0; j < this->cols(); j++)
                sum += pu[i] * (*(pA + i * ld_ + j)) * pv[j];
            
            return sum;
        }
        
        /// Storage layout.
        virtual char layout () const { return DenseMatrix<Type>::RowMajorLayout; }
        
        /// Leading dimension.
        virtual std::size_t ld () const { return ld_; }
        
        /// Identity matrix.
        static RowMatrix<Type> Eye (int size)
        {
            RowMatrix<Type> M (size, size);
            
            for (int i = 0; i < 0; i++)
                M.data()[i * size + i] = Type(1);
            
            return M;
        }
        
        // arithmetic operators
        RowMatrix<Type> const & operator += (DenseMatrix<Type> const & A)
        {
            assert(this->rows() == A.rows());
            assert(this->cols() == A.cols());
            
            this->data_ += A.data();
            
            return *this;
        }
        RowMatrix<Type> const & operator -= (DenseMatrix<Type> const & A)
        {
            assert(this->rows() == A.rows());
            assert(this->cols() == A.cols());
            
            this->data_ -= A.data();
            
            return *this;
        }
        RowMatrix<Type> const & operator *= (Type x)
        {
            this->data_ *= x;
            return *this;
        }
        
        /**
         * @brief Output to file.
         * 
         * Write matrix to a text file. The format is influenced by
         * the "pre" and "pos" characters supplied to the routine.
         * The "pre" string will be present on the beginning of every
         * row of the matrix, the "pos" string on the end. For example,
         * the call
         * @code
         *     RowMatrix (
         *         2, 2, rArray({1., 2., -2., 1.})
         *     ).write (std::cout, "[", "]")
         * @endcode
         * will result in the following output to the STDOUT:
         * @verbatim
         *     [ 1 2 ]
         *     [-2 1 ]
         * @endverbatim
         */
        void write (std::ostream & out, std::string const & pre = "", std::string const & pos = "") const
        {
            // data pointer
            Type const * ptr = this->data().begin();
            
            for (int irow = 0; irow < this->rows(); irow++)
            {
                out << pre;
                for (int icol = 0; icol < this->cols(); icol++)
                {
                    Type x = *(ptr + irow * ld_ + icol);
                    
                    if (x == std::abs(x))
                    {
                        // positive entry
                        out << " " << std::abs(x) << " ";
                    }
                    else
                    {
                        // negative entry
                        out << x << " ";
                    }
                }
                
                out << pos << "\n";
            }
        }
        
#ifdef WITH_PNG
        /**
         * @brief Plot to PNG file.
         * 
         * This will write PNG file data to a supplied stream. The data
         * written are the gray-scale representation of the absolute
         * values of the dense matrix entries.
         * 
         * @note This function uses PNG++. Use compile flag -DWITH_PNG to
         * enable it.
         */
        void plot_abs (std::ofstream & out) const
        {
            // create empty gray-scale 1-bit image
            png::image<png::gray_pixel_16> image (this->cols_, this->rows_);
            
            // skip empty matrix
            if (this->data_.size() > 0)
            {
                // find minimal and maximal value
                double min = std::abs(this->data_[0]), max = std::abs(this->data_[0]);
                for (Type const & x : this->data_)
                {
                    double absx = std::abs(x);
                    if (absx < min) min = absx;
                    if (absx > max) max = absx;
                }
                
                // for all elements
                for (int y = 0; y < this->rows(); y++)
                {
                    for (int x = 0; x < this->cols(); x++)
                    {
                        double fraction = 1. - (std::abs(this->data_[y * ld_ + x]) - min) / (max - min);
                        image.set_pixel(x, y, std::ceil(65535 * fraction));
                    }
                }
            }
            
            // save image
            image.write_stream(out);
        }
#endif
        
    private:
        
        /// Leading dimension.
        std::size_t ld_;
        
        /// Change "data" from column-oriented to row-oriented.
        void reorder_ ()
        {
            NumberArray<Type> new_data (this->rows_ * this->cols_);
            
            for (int irow = 0; irow < this->rows_; irow++)
            for (int icol = 0; icol < this->cols_; icol++)
                new_data[irow * this->cols_ + icol] = this->data_[icol * this->rows() + irow];
            
            this->data_ = new_data;
        }
};

/// Sum two row-major matrices of the same type.
template <class Type, class Base1, class Base2>
RowMatrix<Type> operator + (RowMatrix<Type,Base1> const & A, RowMatrix<Type,Base2> const & B)
{
    RowMatrix<Type> C (A);
    C.data() += B.data();
    return C;
}

/// Subtract two row-major matrices of the same type.
template <class Type, class Base1, class Base2>
RowMatrix<Type> operator - (RowMatrix<Type,Base2> const & A, RowMatrix<Type,Base2> const & B)
{
    RowMatrix<Type> C (A);
    C.data() -= B.data();
    return C;
}

/// Subtract two column-major matrices of the same type.
template <class Type, class Base1, class Base2>
ColMatrix<Type> operator - (ColMatrix<Type,Base1> const & A, ColMatrix<Type,Base2> const & B)
{
    ColMatrix<Type> C (A);
    C.data() -= B.data();
    return C;
}

/// Multiply row-major matrix (from left) by a scalar of the same type.
template <class Type, class Base>
RowMatrix<Type> operator * (Type x, RowMatrix<Type,Base> const & A)
{
    RowMatrix<Type> B (A);
    B.data() *= x;
    return B;
}

/// Multiply row-major matrix (from right) by a scalar of the same type.
template <class Type, class Base>
RowMatrix<Type> operator * (RowMatrix<Type> const & A, Type x)
{
    RowMatrix<Type> B (A);
    B.data() *= x;
    return B;
}

/// Divide row-major maxtrix by a scalar of the same type.
template <class Type, class Base>
RowMatrix<Type> operator / (RowMatrix<Type> const & A, Type x)
{
    RowMatrix<Type> B(A);
    B.data() *= (1./x);
    return B;
}

/**
 * @brief Dot product of kronecker product and a vector
 * 
 * This function will compute the following expression, given two matrices and a vector,
 * @f[
 *     \mathbf{w} = (\mathsf{A} \otimes \mathsf{B}) \cdot \mathbf{v} \,,
 * @f]
 * witnout the need of evaluating (and storing) the Kronecker product.
 */
template <class Base1, class Base2>
cArray kron_dot (RowMatrix<Complex,Base1> const & A, RowMatrix<Complex,Base2> const & B, cArrayView const v)
{
    // allocate output array and a temporary array
    cArray w (A.rows() * B.rows());
    cArray z (A.rows() * B.cols());
    
    // reshape vectors
    RowMatrixView<Complex> V (A.cols(), B.cols(), v);
    RowMatrixView<Complex> W (A.rows(), B.rows(), w);
    RowMatrixView<Complex> Z (A.rows(), B.cols(), z);
    
    // calculate the contraction
    blas::gemm(1., A, V, 0., Z);
    blas::gemm(1., Z, B.T(), 0., W);
    
    // return result
    return w;
}

/**
 * @brief Dot product of kronecker product and a vector
 * 
 * This function will compute the following expression, given two matrices and a vector,
 * @f[
 *     \mathbf{w} = (\mathsf{A} \otimes \mathsf{B}) \cdot \mathbf{v} \,,
 * @f]
 * witnout the need of evaluating (and storing) the Kronecker product.
 */
template <class Base1, class Base2>
cArray kron_dot (ColMatrix<Complex,Base1> const & A, ColMatrix<Complex,Base2> const & B, cArrayView const v)
{
    // allocate output array and a temporary array
    cArray w (A.rows() * B.rows());
    cArray z (A.rows() * B.cols());
    
    // reshape vectors
    RowMatrixView<Complex> V (A.cols(), B.cols(), v);
    RowMatrixView<Complex> W (A.rows(), B.rows(), w);
    RowMatrixView<Complex> Z (A.rows(), B.cols(), z);
    
    // calculate the contraction
    blas::gemm(1., A, V, 0., Z);
    blas::gemm(1., Z, B.T(), 0., W);
    
    // return result
    return w;
}

template <class Type>
NumberArray<Type> operator | (DenseMatrixView<Type> const & A, ArrayView<Type> v)
{
    NumberArray<Type> w (A.rows());
    blas::gemv(1., A, v, 0., w);
    return w;
}

template <class Type>
NumberArray<Type> operator * (DenseMatrixView<Type> const & A, ArrayView<Type> v)
{
    NumberArray<Type> w (A.rows());
    blas::gemv(1., A, v, 0., w);
    return w;
}

template <class Type>
NumberArray<Type> operator | (ArrayView<Type> v, RowMatrixView<Type> const & A)
{
    NumberArray<Type> w (A.cols());
    blas::gemv(1., A.T(), v, 0., w);
    return w;
}

template <class Type>
NumberArray<Type> operator * (ArrayView<Complex> v, RowMatrixView<Complex> const & A)
{
    NumberArray<Type> w (A.cols());
    blas::gemv(1., A.T(), v, 0., w);
    return w;
}

/**
 * @brief General dense matrix transposition.
 * 
 * Transposes a dense matrix with given dimensions.
 * 
 * @param A Dense matrix elements array to transpose of length @c ldA0 x @c ldA.
 * @param ldA0 Original leading dimension.
 * @param ldA New leading dimensions.
 */
template <class Type>
void transpose (ArrayView<Type> A, std::size_t ldA0, std::size_t ldA)
{
    // the matrix size must be divisible by the leading dimension
    assert(A.size() == ldA * ldA0);
    
    // backup the original matrix
    NumberArray<Type> A0 = A;
    
    // fill transposed elements
    for (std::size_t i = 0; i < ldA0; i++)
    for (std::size_t j = 0; j < ldA; j++)
        A[i * ldA + j] = A0[j * ldA0 + i];
}

#endif // HEX_DENSEMATRIX_H
