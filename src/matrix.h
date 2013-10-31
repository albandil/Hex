/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2013                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HEX_SPMATRIX
#define HEX_SPMATRIX

#include <cassert>
#include <complex>
#include <cstring>
#include <iostream>
#include <vector>

#ifndef NO_PNG
#include <png++/png.hpp>
#endif

#include <unistd.h>
#include <umfpack.h>

#include "arrays.h"
#include "complex.h"

// declaration of classes in order to enable them as return types
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
 * Base class for row- and column- oriented matrices.
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
 * an encapsulation of the cArray class with some operators bundled to it.
 */
template <class Type> class ColMatrix : public DenseMatrix<Type>
{
    public:
        
        // constructor
        ColMatrix ()
            : DenseMatrix<Type>() {}
        ColMatrix (int rows, int cols)
            : DenseMatrix<Type>(rows, cols) {}
        ColMatrix (int rows, int cols, const ArrayView<Type> data)
            : DenseMatrix<Type>(rows, cols, data) {}
        
        // transpose
        RowMatrix<Type> T () const
        {
            return RowMatrix<Type>
            (
                this->cols(),
                this->rows(),
                this->data()
            );
        }
        
        // column views
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
        
        // element access
        Type operator() (int i, int j) const { return col(j)[i]; }
        Type & operator() (int i, int j) { return col(j)[i]; }
        
    private:
        
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
        RowMatrix (int rows, int cols)
            : DenseMatrix<Type>(rows, cols) {}
        RowMatrix (int rows, int cols, const ArrayView<Type> data)
            : DenseMatrix<Type>(rows, cols, data) {}
        RowMatrix (ColMatrix<Type> const & m)
            : DenseMatrix<Type>(m.rows(), m.cols(), m.data()) { reorder_(); }
        
        // transpose
        ColMatrix<Type> T() const
        {
            return ColMatrix<Type>
            (
                this->cols(),
                this->rows(),
                this->data()
            );
        }
        
        // row views
        ArrayView<Type> row(int i)
        {
            return ArrayView<Type>
            (
                this->data(),
                i * this->cols(),
                this->cols()
            );
        }
        const ArrayView<Type> row(int i) const
        {
            return ArrayView<Type>
            (
                this->data(),
                i * this->cols(),
                this->cols()
            );
        }
        
        // element access
        Type operator() (int i, int j) const { return row(i)[j]; }
        Type & operator() (int i, int j) { return row(i)[j]; }
        
        // write to file
        void write (std::ofstream & out) const
        {
            // data pointer
            Type const * ptr = this->data().begin();
            
            for (int irow = 0; irow < this->rows(); irow++)
            {
                for (int icol = 0; icol < this->cols(); icol++)
                    out << *ptr++ << " ";
                
                out << "\n";
            }
        }
        
        void plot_abs (std::ofstream & out) const
        {
            // create empty gray-scale 1-bit image
            png::image<png::gray_pixel_1> image (this->cols(), this->rows());
            
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
        
    private:
        
        /// Change "data" from column-oriented to row-oriented.
        void reorder_ ()
        {
            cArray new_data (this->data().size());
            
            for (int irow = 0; irow < this->rows(); irow++)
            for (int icol = 0; icol < this->cols(); icol++)
                new_data[irow * this->cols() + icol] = this->data()[icol * this->rows() + irow];
            
            this->data() = new_data;
        }
};

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
    Complex const * restrict pA = A.begin();
    Complex const * restrict pB = B.begin();
    Complex       * restrict pC = C.begin();
    
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
 * @code
 *     y = A.dot(x, strict_upper | strict_lower)
 * @endcode
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
    
    ~CscMatrix() {}
    
    /**
     * Convert to COO-matrix.
     */
    CooMatrix tocoo() const;
    
    /**
     * Matrix-vector product using a transposed matrix, \f$ A^T \cdot b \f$.
     * For ordinary matrix vector product convert this matrix first to CSR
     * format.
     * @param b Vector to multiply with.
     */
    cArray dotT(const cArrayView b) const;
    
    // Getters
    
    size_t size() const { return i_.size(); }
    size_t rows() const { return m_; }
    size_t cols() const { return n_; }
    lArray const & p() const { return p_; }
    lArray const & i() const { return i_; }
    cArray const & x() const { return x_; }
    
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
    bool hdfsave(const char* name) const;
    
    /**
     * Load matrix from HDF file.
     * @param name Filename.
     */
    bool hdfload(const char* name);
    
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
        : m_(0), n_(0) {}
    CsrMatrix(size_t m, size_t n)
        : m_(m), n_(n) {}
    CsrMatrix(CsrMatrix const & A)
        : m_(A.m_), n_(A.n_), p_(A.p_), i_(A.i_), x_(A.x_) {}
    CsrMatrix(size_t m, size_t n, lArrayView const & p, lArrayView const & i, cArrayView const & x)
        : m_(m), n_(n), p_(p), i_(i), x_(x) {}
    
    // Destructor
    
    ~CsrMatrix() {}
    
    /**
     * Ordinary matrix-vector product, \f$ A\cdot b \f$.
     * @param b Vector to multiply with.
     */
    cArray dot(cArrayView const & b) const;
    
    // Getters
    
    size_t size() const { return i_.size(); }
    size_t rows() const { return m_; }
    size_t cols() const { return n_; }
    lArray const & p() const { return p_; }
    lArray const & i() const { return i_; }
    cArray const & x() const { return x_; }
    
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
     * LU factorization
     */
    class LUft
    {
    public:
        
        /// Default constructor.
        LUft () :
            numeric_(0), matrix_(0), filename_(0), info_(UMFPACK_INFO) {}
        
        /// Copy constructor.
        LUft (LUft const & lu) :
            numeric_(lu.numeric_), matrix_(lu.matrix_), filename_(lu.filename_), info_(lu.info_) {}
        
        /// Initialize the structure using the matrix and its numeric decomposition.
        LUft (const CsrMatrix* matrix, void* numeric) : 
            numeric_(numeric), matrix_(matrix), filename_(0), info_(UMFPACK_INFO) {}
        
        void free ()
        { 
            if (numeric_ != 0)
            {
                umfpack_zl_free_numeric(&numeric_);
                numeric_ = 0;
            }
        }
        
        size_t size () const
        {
            long lnz, unz, m, n, nz_udiag;
            long status = umfpack_zl_get_lunz (
                &lnz, &unz, &m, &n, &nz_udiag, numeric_
            );
            return status == 0 ? (lnz + unz) * 16 : 0; // Byte count
        }
        
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
                long status = umfpack_zl_solve (
                    UMFPACK_Aat,
                    matrix_->p_.data(), matrix_->i_.data(),
                    reinterpret_cast<const double*>(matrix_->x_.data()), 0,
                    reinterpret_cast<double*>(&x[0] + eq * matrix_->n_), 0,
                    reinterpret_cast<const double*>(&b[0] + eq * matrix_->n_), 0,
                    numeric_, nullptr, &info_[0]
                );
                
                // check output
                if (status != 0)
                {
                    std::cerr << "\n[CsrMatrix::LUft::solve] Exit status " << status << "\n";
                    umfpack_zl_report_status(0, status);
                }
            }
        }
        
        rArray const & info() const { return info_; }
        
    private:
        
        /// Numeric decomposition as produced by UMFPACK.
        void* numeric_;
        
        /// Pointer to the matrix that has been factorized. Necessary for validity of @ref numeric_.
        const CsrMatrix* matrix_;
        
        /// (not used at the moment)
        char* filename_;
        
        /// Set of status flags produced by UMFPACK.
        mutable rArray info_;
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
     * Solve the Ax = b problem, where "b" can be a matrix. Uses UMFPACK.
     * @param b Complex vector containing column-major ordered data; it may be
     *          a flattened matrix.
     * @param eqs Number of columns.
     * @return Array of roots in the same shape as "b".
     */
    cArray solve(const cArrayView b, size_t eqs = 1) const;
    
    /**
     * Save matrix to HDF file.
     * @param name Filename.
     */
    bool hdfsave(const char* name) const;
    
    /**
     * Load matrix from HDF file.
     * @param name Filename.
     */
    bool hdfload(const char* name);
    
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
    CsrMatrix sparse_like(CsrMatrix const & B) const;
    
    /**
     * Return dense array with diagonal elements of the matrix.
     */
    cArray diag() const;
    
    /**
     * Convert to COO format.
     */
    CooMatrix tocoo() const;
    
    /**
     * Solves upper triangular system of equations using backsubstitution.
     * The matrix needs not be triangular, but elements under the diagonal
     * won't be used or changed.
     * @param b Right-hand side.
     */
    cArray upperSolve(cArrayView const & b) const;
    
    /**
     * Solves lower triangular system of equations using backsubstitution.
     * The matrix needs not be triangular, but elements above the diagonal
     * won't be used or changed.
     * @param b Right-hand side.
     */
    cArray lowerSolve(cArrayView const & b) const;
    
    /**
     * Applies a user transformation on <b>nonzero</b> matrix elements.
     * @param f A functor compatible with following declaration:
     * @code
     * Complex (*f) (size_t i, size_t j, Complex z);
     * @endcode
     */
    template <class Functor> CsrMatrix nzTransform(Functor f) const
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
    template <class T> CooMatrix(size_t m, size_t n, T a) : m_(m), n_(n), sorted_(false)
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
    
    // Conversions
    
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
    
    Complex operator() (size_t ix, size_t iy) const
    {
        for (size_t n = 0; n < i_.size(); n++)
             if ((size_t)i_[n] == ix and (size_t)j_[n] == iy)
                return x_[n];
        return 0.;
    }
    
    /**
     * Symmetrical band Populator
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
    template <class Functor> CooMatrix& symm_populate_band(size_t d, Functor f)
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
     * Full populator
     * 
     * Sets all values. The matrix can become easily dense!
     * @param f Object compatible with signature:
     * @code
     *  Complex (*) (long, long)
     * @endcode
     */
    template <class Functor> CooMatrix& populate(Functor f)
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
     * Addition of an element to matrix. If an existing coordinated are used,
     * the numbers will be summed.
     */
    void add(long i, long j, Complex v)
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
    
    // Addition.
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
    
    // Subtraction.
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
    
    // Element-wise multiplication by complex number.
    CooMatrix& operator *= (Complex c)
    {
        long nz = i_.size();
        for (long i = 0; i < nz; i++)
            x_[i] *= c;
        
        return *this;
    }
    
    // multiplication
    CooMatrix dot (const cArrayView B) const;
    
    /**
     * Double inner matrix-matrix product, \f$ A : B \f$.
     * \note Works only on sorted data.
     * @param B Other matrix.
     */
    Complex ddot(CooMatrix const & B) const;
    
    /**
     * Matrix multiplication
     * @param B Dense column-major ordered 1D-array. It's row count
     *          is assumed to be equal to the column count of *this
     *          and the column count is computed from the size of the
     *          array, which must be integer multiple of *this's column
     *          count.
     */
    CooMatrix& operator *= (const cArrayView B);
    
    // Element-wise divide by a complex number
    CooMatrix& operator /= (Complex c)
    {
        return *this *= 1./c;
    }
    
    /**
     * Change dimension of the matrix.
     * \warning No row/column index range checking.
     */
    void resize(size_t m, size_t n)
    {
        m_ = m;
        n_ = n;
    }
    
    /**
     * Change matrix shape. Conserves volume, i.e. it holds
     * @f[
     *                m \cdot n = m_0 \cdot n_0
     * @f]
     * @param m New row count.
     * @param n New column coount.
     */
    CooMatrix reshape(size_t m, size_t n) const;
    
    // Convert matrix to dense column-major ordered 1D-array.
    cArray todense() const;
    
    // sort indices (by i_, then by j_)
    void sort();
    bool sorted() const { return sorted_; }
    
    // Convert to CSC matrix.
    CscMatrix tocsc() const;
    
    // Convert to CSR matrix.
    CsrMatrix tocsr() const;
    
    /**
     * @brief Convert to symmetric DIA format.
     * 
     * Use "lower" for conversion of lower triangle,
     * "upper" for conversion of upper triangle and "both" for
     * conversion of the main diagonal only.
     */
    SymDiaMatrix todia(MatrixTriangle triangle = lower) const;
    
    /**
     * Solve the Ax = b problem, where "b" can be a matrix.
     * @param b Complex vector containing column-major ordered data; it may be
     *          a flattened matrix.
     * @param eqs Number of columns.
     * @return Array of roots in the same shape as "b".
     */
    cArray solve(const cArrayView b, size_t eqs = 1) const
    {
        // COO format is not optimal for solving -> covert to CSC
        return tocsr().solve(b, eqs);
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
    void write(const char* filename) const;
    
    // Shake the content, i.e. sum same element entries.
    CooMatrix shake() const;
    
    /**
     * Save matrix to HDF file.
     * @param name Filename.
     */
    bool hdfsave(const char* name) const;
    
    /**
     * Load matrix from HDF file.
     * @param name Filename.
     */
    bool hdfload(const char* name);

private:

    // dimensions
    size_t m_, n_;

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
class SymDiaMatrix {

public:
    //
    // Constructors
    //

    /// Empty constructor.
    SymDiaMatrix();
    
    /// Size constructor.
    SymDiaMatrix(int n);
    
    /**
     * @brief Data constructor.
     * 
     * @param n Size of the matrix.
     * @param id Identifyiers of the diagonals (positive integers expected).
     * @param v Stacked (and padded if ncessary) diagonals.
     */
    SymDiaMatrix(int n, const iArrayView id, const cArrayView v);
    
    /// Copy constructor.
    SymDiaMatrix(SymDiaMatrix const & A);

    /// Move constructor.
    SymDiaMatrix(SymDiaMatrix && A);

    /**
     * @brief Plain symmetrical populator.
     *
     * Given a functor of the signature
     * @code
     *     Complex (*) (int, int);
     * @endcode
     * the function will call the functor with row and column number of every
     * element that is to be set.
     * 
     * @param d How many upper diagonals to populate. The main diagonal will
     *          be populated always.
     * @param f The functor that will compute the matrix elements.
     */
    template <class Functor> SymDiaMatrix & populate(unsigned d, Functor f)
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
        
        return *this;
    }

    //
    // Destructor
    //

    ~SymDiaMatrix() {}
    
    // free all fields, set dimensions to zero
    void clear () { *this = std::move(SymDiaMatrix()); }
    
    //
    // Getters
    //
    
    iArray const & diag() const { return idiag_; }
    int diag(int i) const { return idiag_[i]; }
    cArray const & data() const { return elems_; }
    cArrayView main_diagonal() const { return cArrayView(elems_, 0, n_); }
    Complex const * dptr(int i) const { return dptrs_[i]; }
    
    iArray & diag() { return idiag_; }
    cArray & data() { return elems_; }
    cArrayView main_diagonal() { return cArrayView(elems_, 0, n_); }
    Complex * dptr(int i) { return dptrs_[i]; }
    
    size_t size() const { return n_; }
    int bandwidth() const { return 1 + 2 * idiag_.back(); }
    
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
     * @param triangle Whether to use only the upper or only the lower or both triangles
     *                 of the othwerwise symmetric matrix.
     */
    cArray dot (const cArrayView B, MatrixTriangle triangle = both) const;
    
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
     */
    SymDiaMatrix kron (SymDiaMatrix const & B) const;

    //
    // HDF interface
    //
    
    /**
     * @brief Load from file.
     * 
     * Load the matrix from a HDF5 file created by the routine @ref hdfsave.
     * 
     * @return True on successful read, false otherwise (mostly when doesn't exist).
     */
    bool hdfload(const char* name);
    
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
    bool hdfsave(const char* name, bool docompress = false, int consec = 10) const;

    //
    // Conversions to other formats
    //
    
    CooMatrix tocoo (MatrixTriangle triangle = both) const;
    RowMatrix<Complex> torow (MatrixTriangle triangle = both) const;
    
    //
    // Output
    //
    
    friend std::ostream & operator << (std::ostream & out, SymDiaMatrix const & A);
    
private:

    // dimension (only square matrices allowed)
    int n_;

    // diagonals: concatenated diagonals starting from the longest
    // to the shortest (i.e. with rising right index)
    cArray elems_;
    
    // upper diagonal indices starting from zero
    iArray idiag_;
    
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
CooMatrix kron(const CooMatrix& A, const CooMatrix& B);

/**
 * Identity matrix.
 * @param N Dimension.
 */
CooMatrix eye(size_t N);

/**
 * Diagonal matrix with diagonal consisting of 0, 1, 2, 3, ..., N - 1.
 * @param N Dimension.
 */
CooMatrix stairs(size_t N);

#endif
