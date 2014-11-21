//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2014, Jakub Benda, Charles University in Prague                    //
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

#ifndef HEX_MATRIX
#define HEX_MATRIX

#include <cassert>
#include <complex>
#include <cstring>
#include <iostream>
#include <vector>

#ifndef NO_PNG
#include <png++/png.hpp>
#endif

#include <umfpack.h>

#include "arrays.h"
#include "complex.h"

// forward declaration of classes in order to enable them as return types
// of other classes before proper definition
template <class T> class RowMatrix;
template <class T> class ColMatrix;
class CooMatrix;
class CscMatrix;
class CsrMatrix;
class SymDiaMatrix;

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
        DenseMatrix(int rows, int cols)
            : rows_(rows), cols_(cols), data_(rows * cols) {}
        DenseMatrix(int rows, int cols, const ArrayView<T> data)
            : rows_(rows), cols_(cols), data_(data) { assert(data.size() == (size_t)(rows * cols)); }
        
        // explicit conversion to cArrayView
        ArrayView<T> data () { return data_; }
        const ArrayView<T> data () const { return data_; }
        
        // getters
        size_t size () const { return rows_ * cols_; }
        int cols () const { return cols_; }
        int rows () const { return rows_; }
        T * begin () { return data_.begin(); }
        T const * begin () const { return data_.begin(); }
        
    private:
        
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
template <class Type> class ColMatrix : public DenseMatrix<Type>
{
    public:
        
        // constructor
        ColMatrix ()
            : DenseMatrix<Type>() {}
        ColMatrix (int size)
            : DenseMatrix<Type>(size, size) {}
        ColMatrix (int rows, int cols)
            : DenseMatrix<Type>(rows, cols) {}
        ColMatrix (int rows, int cols, const ArrayView<Type> data)
            : DenseMatrix<Type>(rows, cols, data) {}
        ColMatrix (RowMatrix<Type> const & m)
            : DenseMatrix<Type>(m.rows(), m.cols(), m.data()) { reorder_(); }
        
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
        ColMatrix<Type> invert () const;
        
        /// Diagonalization.
        std::tuple<NumberArray<Type>,ColMatrix<Type>,ColMatrix<Type>> diagonalize () const;
    
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
template <class Type> class RowMatrix : public DenseMatrix<Type>
{
    public:
        
        // constructor
        RowMatrix ()
            : DenseMatrix<Type>() {}
        RowMatrix (int size)
            : DenseMatrix<Type>(size, size) {}
        RowMatrix (int rows, int cols)
            : DenseMatrix<Type>(rows, cols) {}
        RowMatrix (int rows, int cols, const ArrayView<Type> data)
            : DenseMatrix<Type>(rows, cols, data) {}
        RowMatrix (ColMatrix<Type> const & m)
            : DenseMatrix<Type>(m.rows(), m.cols(), m.data()) { reorder_(); }
        
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

template <class T> RowMatrix<T> operator + (RowMatrix<T> const & A, RowMatrix<T> const & B) { RowMatrix<T> C(A); C.data() += B.data(); return C; }
template <class T> RowMatrix<T> operator - (RowMatrix<T> const & A, RowMatrix<T> const & B) { RowMatrix<T> C(A); C.data() -= B.data(); return C; }
template <class T> RowMatrix<T> operator * (T x, RowMatrix<T> const & A) { RowMatrix<T> B(A); B.data() *= x; return B; }
template <class T> RowMatrix<T> operator * (RowMatrix<T> const & A, T x) { RowMatrix<T> B(A); B.data() *= x; return B; }
template <class T> RowMatrix<T> operator / (RowMatrix<T> const & A, T x) { RowMatrix<T> B(A); B.data() /= x; return B; }

/**
 * @brief Dot product of kronecker product and a vector
 * 
 * This function will compute the following expression, given two matrices and a vector:
 * @f[
 *     (A \otimes B) \cdot v
 * @f]
 * witnout the need of evaluating (and storing) the Kronecker product.
 */
cArray kron_dot (RowMatrix<Complex> const & A, RowMatrix<Complex> const & B, cArrayView const v);

/**
 * @brief Dense matrix multiplication.
 */
template <class Type> RowMatrix<Type> operator * (RowMatrix<Type> const & A, ColMatrix<Type> const & B)
{
    assert(A.cols() == B.rows());
    
    // sizes
    int rows = A.rows();
    int cols = B.cols();
    int comm = A.cols();
    
    // output matrix, initialized to zero
    RowMatrix<Type> C(A.rows(), B.cols());
    
    // data pointers
    Type const * restrict pA = A.begin();
    Type const * restrict pB = B.begin();
    Type       * restrict pC = C.begin();
    
    // for all rows of A
    for (int irow = 0; irow < rows; irow++)
    {
        // for all columns of B
        for (int icol = 0; icol < cols; icol++)
        {
            // compute the scalar product "row * column"
            for (int k = 0; k < comm; k++)
                (*pC) += (*(pA + k)) * (*pB++);
            
            // move to next element of C
            pC++;
        }
        
        // move to the next row of A
        pA += A.cols();
        
        // reset B's data pointer
        pB = B.begin();
    }
    
    // return result
    return C;
}

/**
 * @brief Matrix parts.
 * 
 * This enumeration contains some matrix part identificators
 * that can be used e.g. together with matrix-vector multiplication
 * routines (if implemented). The multiplication is done, then, as if
 * the other parts were zero. For example
   @code
       y = A.dot(x, strict_upper | strict_lower)
   @endcode
 * will multiply vector @c x as if the diagonal of @c A were zero.
 */
typedef enum {
    none         = 0,  // b000,
    strict_lower = 1,  // b001, strict_lower
    strict_upper = 2,  // b010,                           strict_upper
    strict_both  = 3,  // b011, strict_lower |            strict_upper
    diagonal     = 4,  // b100,                diagonal
    lower        = 5,  // b101, strict_lower | diagonal
    upper        = 6,  // b110,                diagonal | strict_upper
    both         = 7   // b111, strict_lower | diagonal | strict_upper
} MatrixTriangle;

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
    
    CscMatrix()
        : m_(0), n_(0) {}
    CscMatrix(size_t m, size_t n)
        : m_(m), n_(n) {}
    CscMatrix(CscMatrix const & A)
        : m_(A.m_), n_(A.n_), p_(A.p_), i_(A.i_), x_(A.x_) {}
    CscMatrix(size_t m, size_t n, const lArrayView p, const lArrayView i, const cArrayView x)
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
    long m_;
    long n_;
    
    // representation
    lArray p_;
    lArray i_;
    cArray x_;
    
};

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
class CsrMatrix {
    
public:
    
    // Constructors
    
    CsrMatrix()
        : m_(0), n_(0), name_() {}
    CsrMatrix(size_t m, size_t n)
        : m_(m), n_(n), name_() {}
    CsrMatrix(CsrMatrix const & A)
        : m_(A.m_), n_(A.n_), p_(A.p_), i_(A.i_), x_(A.x_), name_() {}
    CsrMatrix(size_t m, size_t n, lArrayView const & p, lArrayView const & i, cArrayView const & x)
        : m_(m), n_(n), p_(p), i_(i), x_(x), name_() {}
    
    // Destructor
    
    ~CsrMatrix() {}
    
    /**
     * @brief Clear data.
     */
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
    cArray dot (cArrayView const & b) const;
    
    // Getters
    
    std::size_t size () const { return i_.size(); }
    std::size_t rows () const { return m_; }
    std::size_t cols () const { return n_; }
    lArray const & p () const { return p_; }
    lArray const & i () const { return i_; }
    cArray const & x () const { return x_; }
    
    // return absolute value of the (in absolute value) largest element
    double norm() const;
    
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
    void write(const char* filename) const;
    
#ifndef NO_PNG
    /**
     * PNG row data generator for use in \ref plot function.
     * 
     * @note This class requires PNG++. It can be disabled if the macro
     * NO_PNG is defined.
     */
    class PngGenerator : public png::generator<png::gray_pixel_1,PngGenerator>
    {
        typedef png::generator<png::gray_pixel_1,PngGenerator> base_t;
        typedef png::packed_pixel_row<png::gray_pixel_1> row;
        typedef png::row_traits<row> row_traits;
        
        public:
            
            PngGenerator (CsrMatrix const * mat, double threshold);
            ~PngGenerator ();
            
            png::byte* get_next_row (size_t pos);
            
        private:
            
            CsrMatrix const * M_;
            row buff_;
            double threshold_;
    };

    /**
     * Save matrix structure as a black-and-white image.
     * @param filename File name.
     * @param threshold Largest absolute value represented by white colour.
     */
    void plot (const char* filename, double threshold = 0.) const;
#endif

#ifndef NO_UMFPACK
    /**
     * @brief LU factorization object
     * 
     * This class is returned by the function CsrMatrix::factorize() and
     * it provides some functions that can be used when solving equations with
     * that LU factorization. Also, it is possible to store the decomposition
     * to disk (link,save), destro the data (drop) and load later when needed
     * (load). The most important function is "solve".
     */
    class LUft
    {
    public:
        
        /// Default constructor.
        LUft ()
            : numeric_(nullptr), matrix_(nullptr), filename_(), info_(UMFPACK_INFO) {}
        
        /// Copy constructor.
        LUft (LUft const & lu)
            : numeric_(lu.numeric_), matrix_(lu.matrix_), filename_(lu.filename_), info_(lu.info_) {}
        
        /// Move constructor.
        LUft (LUft && lu)
            : numeric_(lu.numeric_), matrix_(lu.matrix_), filename_(lu.filename_), info_(lu.info_) { lu.numeric_ = nullptr; }
        
        /// Initialize the structure using the matrix and its numeric decomposition.
        LUft (const CsrMatrix* matrix, void* numeric)
            : numeric_(numeric), matrix_(matrix), filename_(), info_(UMFPACK_INFO) {}
        
        /// Destructor.
        ~LUft ()
        {
            drop ();
        }
        
        /**
         * @brief Transfer data from another factorization object.
         * 
         * Move contents of another LUft object to this one. If there are already some
         * data in this object, delete them.
         */
        void transfer (LUft && lu)
        {
            // clean this object
            drop();
            
            // transfer data
            numeric_  = lu.numeric_;
            matrix_   = lu.matrix_;
            filename_ = lu.filename_;
            info_     = lu.info_;
            
            // clean source
            lu.numeric_  = nullptr;
        }
        
        /**
         * @brief Size of the numerical data.
         * 
         * Return the number of bytes occupied by the stored elements
         * of the LU-factorization. This doesn't contain any other structural data.
         */
        std::size_t size () const
        {
            if (numeric_ == nullptr)
                return 0;
            
            long lnz, unz, m, n, nz_udiag;
            long status = umfpack_zl_get_lunz
            (
                &lnz, &unz, &m, &n, &nz_udiag, numeric_
            );
            return status == 0 ? (lnz + unz) * 16 : 0; // Byte count
        }
        
        /**
         * @brief Solve equations.
         * 
         * The parameter "b" is assumed to contain several right hand
         * side vectors (their count is supplied as the optional parameter
         * "eqs"). The results are stored in "x", which has the same size
         * as "b".
         */
        //@{
        cArray solve (const cArrayView b, unsigned eqs = 1) const
        {
            // reserve space for the solution
            cArray x(b.size());
            
            // solve
            solve(b, x, eqs);
            
            // return the result
            return x;
        }
        void solve (const cArrayView b, cArrayView x, int eqs = 1) const
        {
            // check sizes
            assert (eqs * matrix_->n_ == (int)x.size());
            assert (eqs * matrix_->n_ == (int)b.size());
            
            // solve for all RHSs
            for (int eq = 0; eq < eqs; eq++)
            {
                // solve for current RHS
                long status = umfpack_zl_solve
                (
                    UMFPACK_Aat,
                    matrix_->p_.data(), matrix_->i_.data(),
                    reinterpret_cast<const double*>(matrix_->x_.data()), nullptr,
                    reinterpret_cast<double*>(&x[0] + eq * matrix_->n_), nullptr,
                    reinterpret_cast<const double*>(&b[0] + eq * matrix_->n_), nullptr,
                    numeric_, nullptr, &info_[0]
                );
                
                // check output
                if (status != UMFPACK_OK)
                {
                    std::cerr << "\n[CsrMatrix::LUft::solve] Exit status " << status << std::endl;
                    umfpack_zl_report_status(0, status);
                }
            }
        }
        //@}
        
        /**
         * @brief Get info array.
         * 
         * Get UMFPACK "info" array.
         */
        rArray const & info () const { return info_; }
        
        /**
         * @brief Link to a disk file.
         * 
         * This function will set a filename that will be used if
         * any of the functions @ref save or @ref load is used without
         * a specific filename.
         */
        void link (std::string name) { filename_ = name; }
        void unlink () { filename_.clear(); }
        
        /**
         * @brief Name of the linked disk file.
         */
        std::string name () const { return filename_; }
        
        /**
         * @brief Save Numeric object to a disk file.
         * 
         * Stores the LU-factorization data in the native UMFPACK format
         * to a disk file.
         */
        void save (std::string name) const
        {
            long err = umfpack_zl_save_numeric (numeric_, const_cast<char*>(filename_.c_str()));
            
            if (err == UMFPACK_ERROR_invalid_Numeric_object)
                throw exception ("[LUft::save] Invalid numeric object.");
            
            if (err == UMFPACK_ERROR_file_IO)
                throw exception ("[LUft::save] Failed to save LU object \"%s\" (size = %ld).", name.c_str(), size());
        }
        void save () const
        {
            save (filename_);
        }
        
        /**
         * @brief Load Numeric object from a disk file.
         * 
         * The expected format is the format of umfpack_zl_save_numeric.
         */
        //@{
        void load (std::string name, bool throw_on_io_failure = true)
        {
            long err = umfpack_zl_load_numeric (&numeric_, const_cast<char*>(filename_.c_str()));
            
            if (err == UMFPACK_ERROR_out_of_memory)
                throw exception ("[LUft::load] Out of memory.");
            
            if (err == UMFPACK_ERROR_file_IO and throw_on_io_failure)
                throw exception ("[LUft::save] Failed to load LU object \"%s\".", name.c_str());
        }
        void load ()
        {
            load (filename_, true);
        }
        void silent_load ()
        {
            load (filename_, false);
        }
        //@}
        
        /**
         * @brief Free memory.
         * 
         * Release memory occupied by the LU-factorization numeric object.
         */
        void drop ()
        {
            if (numeric_ != nullptr)
            {
                umfpack_zl_free_numeric (&numeric_);
                numeric_ = nullptr;
            }
        }
        
    private:
        
        /// Numeric decomposition as produced by UMFPACK.
        void* numeric_;
        
        /// Pointer to the matrix that has been factorized. Necessary for validity of @ref numeric_.
        const CsrMatrix* matrix_;
        
        /// Linked HDF file name.
        std::string filename_;
        
        /// Set of status flags produced by UMFPACK.
        mutable rArray info_;
        
        // Disable bitwise copy
        LUft const & operator= (LUft const &);
    };

    /**
     * @brief Compute (incomplete) LU factorization.
     * 
     * This function computes the LU factorization of the matrix. It uses
     * the free UMFPACK library.
     * 
     * @param droptol Drop tolerance. If an element of the factorization matrices
     *                should be in absolute value smaller than "droptol" the
     *                library will omit it completely effectively making the result
     *                more sparse.
     * @return Special structure of the LUft type, holding opaque information about
     *         the factorization.
     */
    LUft factorize (double droptol = 0) const;
#endif
    
    /**
     * @brief Solve the Ax = b problem, where "b" can be a matrix.
     * 
     * Uses UMFPACK through the LUft::solve function.
     * 
     * @param b Complex vector containing column-major ordered data; it may be
     *          a flattened matrix.
     * @param eqs Number of columns.
     * @return Array of roots in the same shape as "b".
     */
    cArray solve (const cArrayView b, size_t eqs = 1) const;
    
    /**
     * @brief Link to a disk file.
     * 
     * Set a default I/O file that will be used if any of the functions
     * @ref hdfload or @ref hdfsave will be used without an explicit filename.
     */
    void hdflink (std::string name);
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
    bool hdfsave (std::string name) const;
    //@}
    
    /**
     * Load matrix from HDF file.
     */
    //@{
    bool hdfload () { return hdfload(name_); }
    bool hdfload (std::string name);
    //@}
    
    /**
     * Element-wise access (const).
     * If a non-existing element is referenced, zero is returned.
     */
    Complex operator() (unsigned i, unsigned j) const;
    
    /**
     * Multiplication by a number.
     */
    CsrMatrix & operator *= (Complex r);
    
    /**
     * Addition of another CSR matrix.
     * The matrices MUST HAVE THE SAME SPARSE STRUCTURE, as no indices are
     * checked.
     */
    CsrMatrix & operator &= (CsrMatrix const & B);
    
    /**
     * Subtraction of another CSR matrix.
     * The matrices MUST HAVE THE SAME SPARSE STRUCTURE, as no indices are
     * checked.
     */
    CsrMatrix& operator ^= (CsrMatrix const & B);
    
    /**
     * Sets fill-in elements so that the storage structure of this matrix
     * will be identical to that of the other CSR matrix. We can the use
     * “fast” arithmetic operators & and ^.
     * @return self
     */
    CsrMatrix sparse_like (CsrMatrix const & B) const;
    
    /**
     * Return dense array with diagonal elements of the matrix.
     */
    cArray diag() const;
    
    /**
     * Convert to COO format.
     */
    CooMatrix tocoo() const;
    
    /**
     * Convert to dense matrix.
     */
    RowMatrix<Complex> torow() const;
    
    /**
     * Solves upper triangular system of equations using backsubstitution.
     * The matrix needs not be triangular, but elements under the diagonal
     * won't be used or changed.
     * @param b Right-hand side.
     */
    cArray upperSolve (cArrayView const & b) const;
    
    /**
     * Solves lower triangular system of equations using backsubstitution.
     * The matrix needs not be triangular, but elements above the diagonal
     * won't be used or changed.
     * @param b Right-hand side.
     */
    cArray lowerSolve (cArrayView const & b) const;
    
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
    long m_;
    long n_;
    
    // representation
    lArray p_;
    lArray i_;
    cArray x_;
    
    // linked HDF file
    std::string name_;
};


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
class CooMatrix {
    
public:
    
    // Empty constructors
    
    CooMatrix()
        : m_(0), n_(0), sorted_(true) {}
    CooMatrix(size_t m, size_t n)
        : m_(m), n_(n), sorted_(true) {}
    CooMatrix(CooMatrix const & A)
        : m_(A.m_), n_(A.n_), i_(A.i_), j_(A.j_), x_(A.x_), sorted_(false) {}
    CooMatrix(size_t m, size_t n, NumberArray<long> const & i, NumberArray<long> const & j, NumberArray<Complex> const & x)
        : m_(m), n_(n), i_(i), j_(j), x_(x), sorted_(false) {}
    
    /**
     * Copy constructor initialized from dense array.
     * @param m Row counz of new matrix.
     * @param n Column count of new matrix.
     * @param a Column-major ordered dense array with matrix elements.
     *          Only nonzero elements are copied into internal storage.
     */
    template <class T> CooMatrix (size_t m, size_t n, T a) : m_(m), n_(n), sorted_(false)
    {
        // initialize from column-major formatted input
        size_t i = 0;
        for (size_t col = 0; col < n; col++)
        {
            for (size_t row = 0; row < m; row++)
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
    ~CooMatrix() {}
    
    /// Convert 1×1 matrix to a complex number.
    operator Complex () const
    {
        if (m_ == 1 and n_ == 1)
        {
            // get number of nonzero matrix elements
            size_t elems = this->shake().x_.size();
            if (elems > 1)
                throw "[CooMatrix::operator Complex] more elements than nominal volume!\n";
            else if (elems == 1)
                return this->shake().x_[0];
            else
                return 0.;
        }
        else
            throw "[CooMatrix::operator Complex] matrix is not 1×1!\n";
    }
    
    // Getters
    
    size_t rows() const { return m_; }
    size_t cols() const { return n_; }
    size_t size() const { return i_.size(); }
    lArray const & i() const { return i_; }
    lArray const & j() const { return j_; }
    cArray const & v() const { return x_; }
    
    /// Index operator. Returns the existing value or zero.
    Complex operator() (size_t ix, size_t iy) const
    {
        for (size_t n = 0; n < i_.size(); n++)
             if ((size_t)i_[n] == ix and (size_t)j_[n] == iy)
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
    template <class Functor> CooMatrix& symm_populate_band (size_t d, Functor f)
    {
        Complex val;
        
        for (size_t row = 0; row < m_; row++)
        {
            for (size_t col = row; col < n_ and col - row <= d; col++)
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
        
        for (size_t row = 0; row < m_; row++)
        {
            for (size_t col = 0; col < n_; col++)
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
    void add (long i, long j, Complex v)
    {
        i_.push_back(i);
        j_.push_back(j);
        x_.push_back(v);
        
        sorted_ = false;
    }
    
    /// Transposition, implemented as an interchange of "i" and "j" data.
    CooMatrix transpose() const
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
        for (size_t i = prev_size; i < x_.size(); i++)
            x_[i] = -x_[i];
        
        sorted_ = false;
        
        return *this;
    }
    
    /// Element-wise multiplication by complex number.
    CooMatrix& operator *= (Complex c)
    {
        long nz = i_.size();
        for (long i = 0; i < nz; i++)
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
    void resize (size_t m, size_t n)
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
    CooMatrix reshape (size_t m, size_t n) const;
    
    /// Convert matrix to dense column-major ordered 1D-array.
    cArray todense () const;
    
    /// Sort indices (by i_, then by j_)
    void sort();
    bool sorted() const { return sorted_; }
    
    /// Convert to CSC matrix.
    CscMatrix tocsc() const;
    
    /// Convert to CSR matrix.
    CsrMatrix tocsr() const;
    
    /// Convert to dense matrix of a given underlying type.
    template <typename DenseMatrixType> DenseMatrixType todense() const
    {
        DenseMatrixType M (rows(), cols());
        for (unsigned idx = 0; idx < x_.size(); idx++)
            M (i_[idx], j_[idx]) = x_[idx];
        return M;
    }
    
    /// Convert to dense matrix (row-ordered).
    RowMatrix<Complex> torow() const
    {
        return todense<RowMatrix<Complex>>();
    }
    
    /// Convert to dense matrix (column-ordered).
    ColMatrix<Complex> tocol() const
    {
        return todense<ColMatrix<Complex>>();
    }
    
    /**
     * @brief Convert to symmetric DIA format.
     * 
     * Use "lower" for conversion of lower triangle,
     * "upper" for conversion of upper triangle and "both" for
     * conversion of the main diagonal only.
     */
    SymDiaMatrix todia(MatrixTriangle triangle = lower) const;
    
    /**
     * @brief Solve matrix equation.
     * 
     * Solve the Ax = b problem, where "b" can be a matrix.
     * @param b Complex vector containing column-major ordered data; it may be
     *          a flattened matrix.
     * @param eqs Number of columns.
     * @return Array of roots in the same shape as "b".
     */
    cArray solve (const cArrayView b, size_t eqs = 1) const
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
 * @brief Symmetric diagonal matrix.
 * 
 * A multi-diagonal matrix is a sparse matrix of the structure
 * @f[
 *     A = \pmatrix {
 *            \ast   & \dots  & \ast   &        &        &        \cr
 *            \vdots & \ddots &        & \ddots &        &        \cr
 *            \ast   &        & \ddots &        & \ddots &        \cr
 *                   & \ddots &        & \ddots &        & \ast   \cr
 *                   &        & \ddots &        & \ddots & \vdots \cr
 *                   &        &        & \ast   & \dots  & \ast   \cr
 *     } \ ,
 * @f]
 * i.e. it is banded, with nonzero elements only near to the diagonal.
 * Some of the diagonals may be identically zero. This class holds all
 * nonzero main and upper diagonals (lower diagonals are not necessary
 * in the symmetric case). The diagonal storage has several advantages:
 * 
 * - Matrix-vector multiplication can be vectorized, in contrast to the
 *   CSR-matrix-vector multiplication which requires a strongly irregular
 *   memory access.
 * - Matrix-matrix multiplication conserves the DIA structure, with just
 *   a simple increase of bandwidth. It holds that the bandwidth of the
 *   matrix C = AB is equal to the sum of the bandwidths of the factors
 *   decreased by one.
 */
class SymDiaMatrix
{

public:
    //
    // Constructors
    //

    /// Empty constructor.
    SymDiaMatrix ();
    
    /// Size constructor.
    SymDiaMatrix (int n);
    
    /**
     * @brief Data constructor.
     * 
     * @param n Size of the matrix.
     * @param id Identifyiers of the diagonals (positive integers expected).
     */
    SymDiaMatrix (int n, const iArrayView id);
    
    /**
     * @brief Data constructor.
     * 
     * @param n Size of the matrix.
     * @param id Identifyiers of the diagonals (positive integers expected).
     * @param v Stacked (and padded if necessary) diagonals.
     */
    SymDiaMatrix (int n, const iArrayView id, const cArrayView v);
    
    /// Copy constructor.
    SymDiaMatrix (SymDiaMatrix const & A);

    /// Move constructor.
    SymDiaMatrix (SymDiaMatrix && A);
    
    /// Constructor - HDF loader.
    SymDiaMatrix (std::string filename);
    
    /**
     * @brief Plain symmetrical populator.
     *
     * Given a functor of the signature
       @code
           Complex (*) (int, int);
       @endcode
     * the function will call the functor with row and column number of every
     * element that is to be set.
     * 
     * @param d How many upper diagonals to populate. The main diagonal will
     *          be populated always.
     * @param f The functor that will compute the matrix elements.
     */
    template <class Functor> SymDiaMatrix & populate (unsigned d, Functor f)
    {
        // throw away old data
        elems_.resize(0);
        
        // for all diagonals
        for (size_t id = 0; id <= d; id++)
        {
            // add this diagonal
            idiag_.push_back(id);
            
            // for all elements of the diagonal
            for (int icol = id; icol < n_; icol++)
            {
                // get the row index, too
                int irow = icol - id;
                
                // evaluate the element
                elems_.push_back(f(irow,icol));
            }
        }
        
        // update diagonal pointers
        setup_dptrs_();
        
        return *this;
    }

    //
    // Destructor
    //

    ~SymDiaMatrix () {}
    
    /// Free all fields, set dimensions to zero.
    void drop ()
    {
        n_ = 0;
        elems_.drop();
        idiag_.drop();
        dptrs_.clear();
    }
    
    //
    // Getters
    //
    
    /**
     * @brief Diagonal indices.
     * 
     * Return array of indices of the stored diagonals. These are always
     * only main and upper diagonals, so all numbers are non-negative.
     * The array is sorted and begins with zero. Its length is always
     * larger than zero.
     */
    //@{
    iArray const & diag () const { return idiag_; }
    iArray & diag () { return idiag_; }
    //@}
    
    /**
     * @brief Diagonal index.
     * 
     * Return diagonal index for of i-th stored diagonal. The zero-th stored
     * diagonal is always the main diagonal (= 0), but it doesn't have to hold
     * for next diagonals.
     */
    int diag (int i) const { return idiag_[i]; }
    
    /**
     * @brief Main diagonal.
     * 
     * Return direct-access view of the main diagonal.
     */
    //@{
    cArrayView main_diagonal () const { return cArrayView(elems_, 0, n_); }
    cArrayView main_diagonal () { return cArrayView(elems_, 0, n_); }
    //@}
    
    /**
     * @brief Data pointer.
     * 
     * Return direct-access data pointer.
     */
    //@{
    cArray const & data () const { return elems_; }
    cArray & data () { return elems_; }
    //@}
    
    /**
     * @brief Pointer to diagonal data.
     * 
     * @param i Index of the diagonal in the "idiag_" array.
     *          The maximal value is thus less then the number stored
     *          diagonals.
     */
    //@{
    Complex const * dptr (int i) const { return dptrs_[i]; }
    Complex * dptr (int i) { return dptrs_[i]; }
    //@}
    
    /**
     * @brief Matrix dimension.
     * 
     * Return row/column count. The matrix is symmetric and so both
     * counts are equal.
     */
    std::size_t size () const { return n_; }
    
    /**
     * @brief Bandwidth.
     * 
     * Return the bandwidth of the matrix, i.e. number of all (upper, main an lower)
     * diagonals that would have to be stored in a full banded-matrix format.
     */
    int bandwidth () const { return 1 + 2 * idiag_.back(); }
    
    /**
     * @brief Check compatibility of matrices.
     * 
     * Check that the matrix B has the same dimensions as *this matrix and
     * also that they keep the same diagonals. Such matrices can be very effectively
     * summed and subtracted -- just by doing the operation on the stored element arrays.
     */
    bool is_compatible (SymDiaMatrix const & B) const;
    
    //
    // Arithmetic and other operators
    //
    
    SymDiaMatrix const & operator = (SymDiaMatrix && A);
    SymDiaMatrix const & operator = (SymDiaMatrix const & A);
    
    SymDiaMatrix const & operator += (SymDiaMatrix const & B);
    SymDiaMatrix const & operator -= (SymDiaMatrix const & B);
    
    /**
     * @brief Dot product.
     *
     * This is a key member of the structure, defining e.g. the speed of conjugate
     * gradients and evaluation of the scattering amplitudes.
     * 
     * @param B Dense matrix. It is supposed to be stored by columns and to have
     *          dimensions n times k, where n is the column count of (*this) matrix.
     *          Also, though only a view of the array is required, it is assumed
     *          that B is actually NumberArray, i.e. that it is aligned with the alignment
     *          2*sizeof(T).
     * @param triangle Whether to use only the upper or only the lower or both triangles
     *                 of the othwerwise symmetric matrix.
     * @param parallelize Whether to use OpenMP to parallelize the SpMV operation.
     */
    cArray dot (const cArrayView B, MatrixTriangle triangle = both, bool parallelize = false) const;
    
    /**
     * @brief Back-substitution (lower).
     * 
     * Assume the matrix is normalized lower-triangular (i.e. has unit main diagonal
     * and zero upper triangle) and do the triangular solve.
     * 
     * @param b Right hand side of the triangular system.
     */
    cArray lowerSolve (const cArrayView b) const;
    
    /**
     * @brief Back-substitution (upper).
     * 
     * Assume the matrix is normalized upper-triangular (i.e. has unit main diagonal
     * and zero lower triangle) and do the triangular solve.
     * 
     * @param b Right hand side of the triangular system.
     */
    cArray upperSolve (const cArrayView b) const;
    
    /**
     * @brief Kronecker product.
     * 
     * Compute Kronecker product with other matrix.
     */
    SymDiaMatrix kron (SymDiaMatrix const & B) const;

    //
    // HDF interface
    //
    
    /// Link matrix to a disk file.
    void hdflink (std::string name);
    
    /// Return the name of the linked disk file.
    std::string hdfname () const { return name_; }
    
    /// Return content of the 'name' file as a new SymDiaMatrix object.
    SymDiaMatrix hdfget () const { return SymDiaMatrix(name_); }
    
    /**
     * @brief Load from file.
     * 
     * Load the matrix from a HDF5 file created by the routine @ref hdfsave.
     * 
     * @return True on successful read, false otherwise (mostly when doesn't exist).
     */
    //@{
    bool hdfload () { return hdfload (name_); }
    bool hdfload (std::string name);
    //@}
    
    /**
     * @brief Save data to file.
     * 
     * Save the matrix to the HDF5 file as a set of four datasets:
     * - idiag - identifiers of (positive) diagonals that contain nonzeros,
     * - n - single number containing size of the matrix (rows or columns)
     * - x - concatenated diagonal elements (interleaved complex values)
     * - zero_blocks - information on the compression of the "x" array, see @ref NumberArray::compress.
     * 
     * @return True on successful write, false otherwise.
     */
    //@{
    bool hdfsave () const { return hdfsave (name_); }
    bool hdfsave (std::string name, bool docompress = false, int consec = 10) const;
    //@}

    //
    // Conversions to other formats
    //
    
    /**
     * @brief Zero-pad rows.
     * 
     * Pad rows with zeros as below:
     * @f[
     *      \left(
     *           \matrix {
     *               \ast & \ast &      &      &      \cr
     *               \ast & \ast & \ast &      &      \cr
     *               \ast & \ast & \ast & \ast &      \cr
     *               \ast & \ast & \ast & \ast & \ast \cr
     *                    & \ast & \ast & \ast & \ast \cr
     *                    &      & \ast & \ast & \ast \cr
     *                    &      &      & \ast & \ast \cr
     *           }
     *      \right)
     *      \longrightarrow
     *      \matrix {
     *           0    & 0    & 0 \cr
     *                & 0    & 0 \cr
     *                &      & 0 \cr
     *                &      & . \cr
     *                &      & . \cr
     *                &      & . \cr
     *                &      & . \cr
     *      }
     *      \left(
     *           \matrix {
     *               \ast & \ast &      &      &      \cr
     *               \ast & \ast & \ast &      &      \cr
     *               \ast & \ast & \ast & \ast &      \cr
     *               \ast & \ast & \ast & \ast & \ast \cr
     *                    & \ast & \ast & \ast & \ast \cr
     *                    &      & \ast & \ast & \ast \cr
     *                    &      &      & \ast & \ast \cr
     *           }
     *      \right)
     *      \matrix {
     *           .    &      &      \cr
     *           .    &      &      \cr
     *           .    &      &      \cr
     *           .    &      &      \cr
     *           0    &      &      \cr
     *           0    & 0    &      \cr
     *           0    & 0    & 0    \cr
     *      }
     *      \longrightarrow
     *      \matrix {
     *          0    & 0    & 0    & \ast & \ast \cr
     *          0    & 0    & \ast & \ast & \ast \cr
     *          0    & \ast & \ast & \ast & \ast \cr
     *          \ast & \ast & \ast & \ast & \ast \cr
     *          \ast & \ast & \ast & \ast & 0    \cr
     *          \ast & \ast & \ast & 0    & 0    \cr
     *          \ast & \ast & 0    & 0    & 0    \cr
     *      }
     * @f]
     * Then concatenate rows and return as a single array.
     */
    cArray toPaddedRows () const;
    
    /**
     * @brief Zero-pad columns.
     * 
     * Pad rows with zeros as in @ref toPaddedRows, then concatenate
     * columns and return as a single array.
     */
    cArray toPaddedCols () const;
    
    /// Convert matrix part to CooMatrix.
    CooMatrix tocoo (MatrixTriangle triangle = both) const;
    
    /// Convert matrix part to RowMatrix.
    RowMatrix<Complex> torow (MatrixTriangle triangle = both) const;
    
    /// Output to a text stream.
    friend std::ostream & operator << (std::ostream & out, SymDiaMatrix const & A);
    
private:

    // dimension (only square matrices allowed)
    int n_;

    // diagonals: concatenated diagonals starting from the longest
    // to the shortest (i.e. with rising right index)
    cArray elems_;
    
    // upper diagonal indices starting from zero
    iArray idiag_;
    
    // name of linked HDF file
    std::string name_;
    
    // diagonal data pointers
    std::vector<Complex*> dptrs_;
    
    /** 
     * @brief Setup diagonal data pointers
     * 
     * Make the pointer list dptrs_ point to the beginnings of the diagonal
     * data. If there are 2d+1 diagonals in the martix, the pointers are
     * assigned in the following way:
     * @verbatim
     * dptrs_[0] ... main diagonal
     * dptrs_[1] ... first upper diagonal (identical to 1st lower diagonal)
     * ...
     * dptrs_[d] ... d-th upper diagonal (identical to d-th lower diagonal)
     * @endverbatim
     */
    void setup_dptrs_();
};

// --------------------------------------------------------------------------//
// ---- Binary arithmetic operators ---------------------------------------- //
// --------------------------------------------------------------------------//

/**
 * Computes a sum of two csr-matrices OF THE SAME SPARSE STRUCTURE.
 */
inline CsrMatrix operator & (const CsrMatrix& A, const CsrMatrix& B)
{
    CsrMatrix C = A;
    return C &= B;
}

/**
 * Computes a difference of two csr-matrices OF THE SAME SPARSE STRUCTURE.
 */
inline CsrMatrix operator ^ (const CsrMatrix& A, const CsrMatrix& B)
{
    CsrMatrix C = A;
    return C ^= B;
}

/**
 * Multiplication of csr-matrix by a number.
 */
inline CsrMatrix operator * (Complex z, const CsrMatrix&  B)
{
    CsrMatrix C = B;
    return C *= z;
}

/**
 * Multiplication of csr-matrix by a number.
 */
inline CsrMatrix operator * (const CsrMatrix& A, double r)
{
    CsrMatrix C = A;
    return C *= r;
}

/**
 * Computes a sum of two csc-matrices OF THE SAME SPARSE STRUCTURE.
 */
inline CscMatrix operator & (const CscMatrix& A, const CscMatrix& B)
{
    CscMatrix C = A;
    return C &= B;
}

/**
 * Computes a difference of two csc-matrices OF THE SAME SPARSE STRUCTURE.
 */
inline CscMatrix operator ^ (const CscMatrix& A, const CscMatrix& B)
{
    CscMatrix C = A;
    return C ^= B;
}

/**
 * Multiplication of csc-matrix by a number.
 */
inline CscMatrix operator * (double r, const CscMatrix& B)
{
    CscMatrix C = B;
    return C *= r;
}

/**
 * Multiplication of csc-matrix by a number.
 */
inline CscMatrix operator * (const CscMatrix& A, double r)
{
    CscMatrix C = A;
    return C *= r;
}

/**
 * Computes a sum of two coo-matrices.
 * @param A COO-matrix.
 * @param B COO-matrix.
 */
inline CooMatrix operator + (const CooMatrix& A, const CooMatrix& B)
{
    CooMatrix C = A;
    return C += B;
}

/**
 * Computes a difference of two coo-matrices.
 * @param A COO-matrix.
 * @param B COO-matrix.
 */
inline CooMatrix operator - (const CooMatrix& A, const CooMatrix& B)
{
    CooMatrix C = A;
    return C -= B;
}

/**
 * Computes product of a matrix and a complex number.
 * @param z Number.
 * @param B Matrix.
 */
inline CooMatrix operator * (const Complex& z, const CooMatrix& B)
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
inline CooMatrix operator * (const CooMatrix& A, const Complex& z)
{
    CooMatrix C = A;
    return C *= z;
}

SymDiaMatrix operator + (SymDiaMatrix const & A, SymDiaMatrix const & B);
SymDiaMatrix operator - (SymDiaMatrix const & A, SymDiaMatrix const & B);
SymDiaMatrix operator * (SymDiaMatrix const & A, SymDiaMatrix const & B);
SymDiaMatrix operator * (double z, SymDiaMatrix const & A);
SymDiaMatrix operator * (Complex z, SymDiaMatrix const & A);

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
SymDiaMatrix kron (SymDiaMatrix const & A, SymDiaMatrix const & B);

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

#endif // HEX_MATRIX
