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

#ifndef NO_PNG
#include <png++/png.hpp>
#endif

#include "arrays.h"
#include "misc.h"

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
        
        // constructors
        DenseMatrixView()
            : rows_(0), cols_(0), data_() {}
        DenseMatrixView(std::size_t rows, std::size_t cols)
            : rows_(rows), cols_(cols), data_(rows * cols) {}
        DenseMatrixView(std::size_t rows, std::size_t cols, const ArrayView<T> data)
            : rows_(rows), cols_(cols), data_(data) { assert(data.size() == rows * cols); }
        
        // explicit conversion to cArrayView
        ArrayView<T> data () { return data_; }
        const ArrayView<T> data () const { return data_; }
        
        // getters
        void drop () { rows_ = cols_ = 0; data_.drop(); }
        std::size_t size () const { return rows_ * cols_; }
        int cols () const { return cols_; }
        int rows () const { return rows_; }
        T * begin () { return data_.begin(); }
        T const * begin () const { return data_.begin(); }
        
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
template <class T> class DenseMatrix
{
    public:
        
        // constructors
        DenseMatrix()
            : rows_(0), cols_(0), data_() {}
        DenseMatrix(std::size_t rows, std::size_t cols)
            : rows_(rows), cols_(cols), data_(rows * cols) {}
        DenseMatrix(std::size_t rows, std::size_t cols, const ArrayView<T> data)
            : rows_(rows), cols_(cols), data_(data) { assert(data.size() == rows * cols); }
        
        // explicit conversion to cArrayView
        ArrayView<T> data () { return data_; }
        const ArrayView<T> data () const { return data_; }
        
        // getters
        void drop () { rows_ = cols_ = 0; data_.drop(); }
        std::size_t size () const { return rows_ * cols_; }
        int cols () const { return cols_; }
        int rows () const { return rows_; }
        T * begin () { return data_.begin(); }
        T const * begin () const { return data_.begin(); }
        
    protected:
        
        /// Row count.
        int rows_;
        
        /// Column count.
        int cols_;
        
        /// Matrix elements, consecutive rows joined to one array.
        NumberArray<T> data_;
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
            : Base() {}
        ColMatrix (int size)
            : Base(size, size) {}
        ColMatrix (int rows, int cols)
            : Base(rows, cols) {}
        ColMatrix (int rows, int cols, const ArrayView<Type> data)
            : Base(rows, cols, data) {}
        ColMatrix (ColMatrix<Type> const & m)
            : Base(m.rows(), m.cols(), m.data()) { }
        ColMatrix (RowMatrix<Type> const & m)
            : Base(m.rows(), m.cols(), m.data()) { reorder_(); }
        
        /**
         * @brief Assignment operator.
         */
        ColMatrix<Type> & operator = (ColMatrix<Type> const & A)
        {
            this->rows_ = A.rows_;
            this->cols_ = A.cols_;
            this->data_ = A.data_;
            return *this;
        }
        
        /**
         * @brief Move assignment.
         */
        ColMatrix<Type> & operator = (ColMatrix<Type> const && A)
        {
            this->rows_ = std::move(A.rows_);
            this->cols_ = std::move(A.cols_);
            this->data_ = std::move(A.data_);
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
            Type * ptr = this->begin();
            
            // fill the matrix
            for (int icol = 0; icol < this->cols(); icol++)
            for (int irow = 0; irow < this->rows(); irow++)
                *(ptr++) = f(irow,icol);
        }
        
        /**
         * @brief Matrix transpose.
         * 
         * Return transposed matrix. The result is of the type RowMatrix,
         * so that no reordering of the entries is necessary.
         */
        RowMatrix<Type> T () const
        {
            return RowMatrix<Type>
            (
                this->cols(),
                this->rows(),
                this->data()
            );
        }
        
        /**
         * @brief Matrix hermitian conjugate.
         * 
         * Return hermitian conjugated matrix. The result is of the type RowMatrix,
         * so that no reordering of the entries is necessary.
         */
        RowMatrix<Type> H () const
        {
            return RowMatrix<Type>
            (
                this->cols(),
                this->rows(),
                NumberArray<Type>(this->data()).conj()
            );
        }
        
        /**
         * @brief Matrix column.
         *
         * Return shallow copy of the chosen matrix column.
         */
        //@{
        ArrayView<Type> col (int i)
        {
            return ArrayView<Type>
            (
                this->data(),
                i * this->rows(),
                this->rows()
            );
        }
        const ArrayView<Type> col (int i) const
        {
            return ArrayView<Type>
            (
                this->data(),
                i * this->rows(),
                this->rows()
            );
        }
        //@}
        
        /// Element access.
        //@{
        Type operator() (int i, int j) const { return col(j)[i]; }
        Type & operator() (int i, int j) { return col(j)[i]; }
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
    
    private:
        
        /// Change "data" from row-oriented to column-oriented.
        void reorder_ ()
        {
            NumberArray<Type> new_data (this->data().size());
            
            for (int irow = 0; irow < this->rows(); irow++)
            for (int icol = 0; icol < this->cols(); icol++)
                new_data[icol * this->rows() + irow] = this->data()[irow * this->cols() + icol];
            
            this->data() = new_data;
        }
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
            : Base() {}
        RowMatrix (int size)
            : Base(size, size) {}
        RowMatrix (int rows, int cols)
            : Base(rows, cols) {}
        RowMatrix (int rows, int cols, const ArrayView<Type> data)
            : Base(rows, cols, data) {}
        RowMatrix (RowMatrix<Type> const & m)
            : Base(m.rows(), m.cols(), m.data()) { }
        RowMatrix (ColMatrix<Type> const & m)
            : Base(m.rows(), m.cols(), m.data()) { reorder_(); }
        
        /**
         * @brief Assignment operator.
         */
        RowMatrix<Type> & operator = (RowMatrix<Type> const & A)
        {
            this->rows_ = A.rows_;
            this->cols_ = A.cols_;
            this->data_ = A.data_;
            return *this;
        }
        
        /**
         * @brief Move assignment.
         */
        RowMatrix<Type> & operator = (RowMatrix<Type> const && A)
        {
            this->rows_ = std::move(A.rows_);
            this->cols_ = std::move(A.cols_);
            this->data_ = std::move(A.data_);
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
            Type * ptr = this->begin();
            
            // fill the matrix
            for (int irow = 0; irow < this->rows(); irow++)
            for (int icol = 0; icol < this->cols(); icol++)
                *(ptr++) = f(irow,icol);
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
            Type * ptr = this->begin();
            
            // transform the matrix
            for (int irow = 0; irow < this->rows(); irow++)
            for (int icol = 0; icol < this->cols(); icol++)
                f(irow,icol,*(ptr++));
        }
        
        /**
         * @brief Matrix transpose.
         * 
         * Return transposed matrix. The result is of the type ColMatrix,
         * so that no reordering of the entries is necessary.
         */
        ColMatrix<Type> T() const
        {
            return ColMatrix<Type>
            (
                this->cols(),
                this->rows(),
                this->data()
            );
        }
        
        /**
         * @brief Matrix row.
         * 
         * Return shallow copy of a matrix row.
         */
        //@{
        ArrayView<Type> row (int i)
        {
            return ArrayView<Type>
            (
                this->data(),
                i * this->cols(),
                this->cols()
            );
        }
        const ArrayView<Type> row (int i) const
        {
            return ArrayView<Type>
            (
                this->data(),
                i * this->cols(),
                this->cols()
            );
        }
        //@}
        
        /// Element access.
        //@{
        Type operator() (int i, int j) const { return row(i)[j]; }
        Type & operator() (int i, int j) { return row(i)[j]; }
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
                sum += pu[i] * (*(pA++)) * pv[j];
            
            return sum;
        }
        
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
            assert (this->rows() == A.rows());
            assert (this->cols() == A.cols());
            
            this->data() += A.data();
            
            return *this;
        }
        RowMatrix<Type> const & operator -= (DenseMatrix<Type> const & A)
        {
            assert (this->rows() == A.rows());
            assert (this->cols() == A.cols());
            
            this->data() -= A.data();
            
            return *this;
        }
        RowMatrix<Type> const & operator *= (Type x)
        {
            this->data() *= x;
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
                    if (*ptr == std::abs(*ptr))
                    {
                        // positive entry
                        out << " " << std::abs(*ptr++) << " ";
                    }
                    else
                    {
                        // negative entry
                        out << *ptr++ << " ";
                    }
                }
                
                out << pos << "\n";
            }
        }
        
#ifndef NO_PNG
        /**
         * @brief Plot to PNG file.
         * 
         * This will write PNG file data to a supplied stream. The data
         * written are the gray-scale representation of the absolute
         * values of the dense matrix entries.
         * 
         * @note This function uses PNG++ and can be disables if the macro
         * NO_PNG is defined.
         */
        void plot_abs (std::ofstream & out) const
        {
            // create empty gray-scale 1-bit image
            png::image<png::gray_pixel_16> image (this->cols(), this->rows());
            
            // skip empty matrix
            if (this->data().size() > 0)
            {
                // find minimal and maximal value
                double min = std::abs(this->data()[0]), max = std::abs(this->data()[0]);
                for (Type const & x : this->data())
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
                        double fraction = 1. - (std::abs(this->data()[y * this->cols() + x]) - min) / (max - min);
                        image.set_pixel(x, y, std::ceil(65535 * fraction));
                    }
                }
            }
            
            // save image
            image.write_stream(out);
        }
#endif
        
    private:
        
        /// Change "data" from column-oriented to row-oriented.
        void reorder_ ()
        {
            NumberArray<Type> new_data (this->data().size());
            
            for (int irow = 0; irow < this->rows(); irow++)
            for (int icol = 0; icol < this->cols(); icol++)
                new_data[irow * this->cols() + icol] = this->data()[icol * this->rows() + irow];
            
            this->data() = new_data;
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
 * 
 * @note This routine avoid reallocating of its intermediate work matrix
 * by making it static and only reallocating if the dimensions change. If you
 * need to free the memory occupied by the work matrix, call this function
 * with zero/nullptr arguments.
 */
void dense_kron_dot
(
    int A_rows, int A_cols, Complex const * A_data,
    int B_rows, int B_cols, Complex const * B_data,
    Complex const * v_data,
    Complex       * w_data
);

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
    // allocate output array
    cArray w (v.size());
    
    // calculate the kronecker contraction
    dense_kron_dot
    (
        A.rows(), A.cols(), A.data().data(),
        B.rows(), B.cols(), B.data().data(),
        v.data(), w.data()
    );
    
    // return result
    return w;
}

/**
 * @brief Matrix-vector multiplication.
 * 
 * Multiplies row-major matrix times colum vector, producing a column vector.
 * Internally, the xGEMV BLAS function is used.
 * 
 * @note This is a general template, which fails. However, it is overloaded with supported type combinations.
 */
template <class Type, class Base>
NumberArray<Type> operator * (RowMatrix<Type,Base> const & A, const ArrayView<Type> v)
{
    HexException("Don't know how to multipy matrix times vector of type %s.", typeid(Type).name);
}

/**
 * @brief Matrix-vector multiplication.
 * 
 * Multiplies complex row-major matrix times complex colum vector, producing a complex column vector.
 * Internally, the ZGEMV BLAS function is used.
 */
template <class Base>
rArray operator * (RowMatrix<double,Base> const & A, const rArrayView v)
{
    // output array
    rArray w (v.size());
    
    // auxiliary variables
    char trans = 'T';
    int m = A.rows(), n = A.cols(), inc = 1;
    double alpha = 1, beta = 0;
    
    // calculate the product
    dgemv_
    (
        &trans, &m, &n, &alpha, const_cast<double*>(A.data().data()), &m,
        const_cast<double*>(v.data()), &inc, &beta, w.data(), &inc
    );
    
    // return the result
    return w;
}

/**
 * @brief Matrix-vector multiplication.
 * 
 * Multiplies complex row-major matrix times complex colum vector, producing a complex column vector.
 * Internally, the ZGEMV BLAS function is used.
 */
template <class Base>
cArray operator * (RowMatrix<Complex,Base> const & A, const cArrayView v)
{
    // output array
    cArray w (v.size());
    
    // auxiliary variables
    char trans = 'T';
    int m = A.rows(), n = A.cols(), inc = 1;
    Complex alpha = 1, beta = 0;
    
    // calculate the product
    zgemv_
    (
        &trans, &m, &n, &alpha, const_cast<Complex*>(A.data().data()), &m,
        const_cast<Complex*>(v.data()), &inc, &beta, w.data(), &inc
    );
    
    // return the result
    return w;
}

/**
 * @brief Matrix-matrix multiplication.
 * 
 * Multiplies row-major matrix times colum-major matrix, producing a row-major matrix.
 * Internally, the xGEMM BLAS function is used.
 * 
 * @note This is a general template, which fails. However, it is overloaded with supported type combinations.
 */
template <class Type, class Base1, class Base2>
RowMatrix<Type,DenseMatrix<Type>> operator * (RowMatrix<Type,Base1> const & A, ColMatrix<Type,Base2> const & B)
{
    HexException("Don't know how to multipy matrices of type %s.", typeid(Type).name);
}

/**
 * @brief Matrix-matrix multiplication.
 * 
 * Multiplies real row-major matrix times real colum-major matrix, producing a row-major matrix.
 * Internally, the DGEMM BLAS function is used.
 */
template <class Base1, class Base2>
RowMatrix<double> operator * (RowMatrix<double,Base1> const & A, ColMatrix<double,Base2> const & B)
{
    // C(row) = A(row) B(col) = C'(col) = A'(col) B(col)
    // =>
    // C(col) = (A'(col) B(col))' = B'(col) A(col)
    
    if (A.cols() != B.rows())
        HexException("Matrix multiplication requires A.cols() == B.rows(), but %d != %d.", A.cols(), B.rows());
    
    // dimensions
    int m = B.cols(), n = A.cols(), k = A.rows();
    
    // create output matrix
    RowMatrix<double> C (A.rows(), B.cols());
    
    // use the BLAS-3 routine
    char transA = 'T', transB = 'N';
    double alpha = 1, beta = 0;
    dgemm_
    (
        &transA, &transB, &m, &n, &k,
        &alpha, const_cast<double*>(B.data().data()), &k,
                const_cast<double*>(A.data().data()), &k,
        &beta, C.data().data(), &m
    );

    // return result
    return C;
}

/**
 * @brief Matrix-matrix multiplication.
 * 
 * Multiplies complex row-major matrix times complex colum-major matrix, producing a row-major matrix.
 * Internally, the ZGEMM BLAS function is used.
 */
template <class Base1, class Base2>
RowMatrix<Complex> operator * (RowMatrix<Complex,Base1> const & A, ColMatrix<Complex,Base2> const & B)
{
    // C(row) = A(row) B(col) = C'(col) = A'(col) B(col)
    // =>
    // C(col) = (A'(col) B(col))' = B'(col) A(col)
    
    if (A.cols() != B.rows())
        HexException("Matrix multiplication requires A.cols() == B.rows(), but %d != %d.", A.cols(), B.rows());
    
    // dimensions
    int m = B.cols(), n = A.cols(), k = B.rows();
    
    // output matrix
    RowMatrix<Complex> C (A.rows(), B.cols());
    
    // use the BLAS-3 routine
    char transA = 'T', transB = 'N';
    Complex alpha = 1, beta = 0;
    zgemm_
    (
        &transA, &transB, &m, &n, &k,
        &alpha, const_cast<Complex*>(B.data().data()), &k,
                const_cast<Complex*>(A.data().data()), &k,
        &beta, C.data().data(), &m
    );
    
    // return resulting matrix
    return C;
}

/**
 * @brief Matrix-matrix multiplication.
 * 
 * Multiplies column-major matrix times colum-major matrix, producing a column-major matrix.
 * Internally, the xGEMM BLAS function is used.
 * 
 * @note This is a general template, which fails. However, it is overloaded with supported type combinations.
 */
template <class Type, class Base1, class Base2>
ColMatrix<Type> operator * (ColMatrix<Type> const & A, ColMatrix<Type> const & B)
{
    HexException("Don't know how to multipy matrices of type %s.", typeid(Type).name);
}

/**
 * @brief Matrix-matrix multiplication.
 * 
 * Multiplies real column-major matrix times real colum-major matrix, producing a real column-major matrix.
 * Internally, the ZGEMM BLAS function is used.
 */
template <class Base1, class Base2>
ColMatrix<double> operator * (ColMatrix<double,Base1> const & A, ColMatrix<double,Base2> const & B)
{
    // C(col) = A(col) B(col)
    
    if (A.cols() != B.rows())
        HexException("Matrix multiplication requires A.cols() == B.rows(), but %d != %d.", A.cols(), B.rows());
    
    // dimensions
    int m = A.rows(), n = B.cols(), k = A.cols();
    
    // output matrix
    ColMatrix<double> C (A.rows(), B.cols());
    
    // use the BLAS-3 routine
    char transA = 'N', transB = 'N';
    double alpha = 1, beta = 0;
    dgemm_
    (
        &transA, &transB, &m, &n, &k,
        &alpha, const_cast<double*>(A.data().data()), &k, // FIXME: Not sure about the indices in
                const_cast<double*>(B.data().data()), &k, //        the non-square matrix case.
        &beta, C.data().data(), &m
    );
    
    // return resulting matrix
    return C;
}

/**
 * @brief Matrix-matrix multiplication.
 * 
 * Multiplies complex column-major matrix times complex colum-major matrix, producing a column-major matrix.
 * Internally, the ZGEMM BLAS function is used.
 */
template <class Base1, class Base2>
ColMatrix<Complex> operator * (ColMatrix<Complex,Base1> const & A, ColMatrix<Complex,Base2> const & B)
{
    // C(col) = A(col) B(col)
    
    if (A.cols() != B.rows())
        HexException("Matrix multiplication requires A.cols() == B.rows(), but %d != %d.", A.cols(), B.rows());
    
    // dimensions
    int m = A.rows(), n = B.cols(), k = A.cols();
    
    // output matrix
    ColMatrix<Complex> C (A.rows(), B.cols());
    
    // use the BLAS-3 routine
    char transA = 'N', transB = 'N';
    Complex alpha = 1, beta = 0;
    zgemm_
    (
        &transA, &transB, &m, &n, &k,
        &alpha, const_cast<Complex*>(A.data().data()), &k,  // FIXME: Not sure about the indices in
                const_cast<Complex*>(B.data().data()), &k,  //        the non-square matrix case.
        &beta, C.data().data(), &m
    );
    
    // return resulting matrix
    return C;
}

#endif // HEX_DENSEMATRIX_H
