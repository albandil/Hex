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

#ifdef WITH_PNGPP
#include <png++/png.hpp>
#endif

#include <unistd.h>
#include <umfpack.h>

#include "arrays.h"
#include "complex.h"

// declaration of classes in order to enable them as return types
// of other classes before proper definition
class CooMatrix;
class CscMatrix;
class CsrMatrix;
class SymDiaMatrix;

typedef enum {
    none  = 0,  // b00
    lower = 1,  // b01
    upper = 2,  // b10
    both  = 3   // b11
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
    
    // Friends
    
    friend class CooMatrix;
    
    // Constructors
    
    CscMatrix()
        : _m_(0), _n_(0) {}
    CscMatrix(size_t m, size_t n)
        : _m_(m), _n_(n) {}
    CscMatrix(const CscMatrix& A)
        : _m_(A._m_), _n_(A._n_), _p_(A._p_), _i_(A._i_), _x_(A._x_) {}
    CscMatrix(size_t m, size_t n, lArrayView const & p, lArrayView const & i, cArrayView const & x)
        : _m_(m), _n_(n), _p_(p), _i_(i), _x_(x) {}
    
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
    cArray dotT(cArrayView const &  b) const;
    
    // Getters
    
    size_t size() const { return _i_.size(); }
    size_t rows() const { return _m_; }
    size_t cols() const { return _n_; }
    
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
    long _m_;
    long _n_;
    
    // representation
    lArray _p_;
    lArray _i_;
    cArray _x_;
    
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
    
    // Friends
    
    friend class CooMatrix;
    
    // Constructors
    
    CsrMatrix()
        : _m_(0), _n_(0) {}
    CsrMatrix(size_t m, size_t n)
        : _m_(m), _n_(n) {}
    CsrMatrix(CsrMatrix const & A)
        : _m_(A._m_), _n_(A._n_), _p_(A._p_), _i_(A._i_), _x_(A._x_) {}
    CsrMatrix(size_t m, size_t n, lArrayView const & p, lArrayView const & i, cArrayView const & x)
        : _m_(m), _n_(n), _p_(p), _i_(i), _x_(x) {}
    
    // Destructor
    
    ~CsrMatrix() {}
    
    /**
     * Ordinary matrix-vector product, \f$ A\cdot b \f$.
     * @param b Vector to multiply with.
     */
    cArray dot(cArrayView const & b) const;
    
    // Getters
    
    size_t size() const { return _i_.size(); }
    size_t rows() const { return _m_; }
    size_t cols() const { return _n_; }
    lArray const & p() const { return _p_; }
    lArray const & i() const { return _i_; }
    cArray const & x() const { return _x_; }
    
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
    
#ifdef WITH_PNGPP
    
    /**
     * PNG row data generator for use in \ref plot function.
     */
    class PngGenerator : public png::generator<png::gray_pixel_1,PngGenerator>
    {
        typedef png::generator<png::gray_pixel_1,PngGenerator> base_t;
        typedef png::packed_pixel_row<png::gray_pixel_1> row;
        typedef png::row_traits<row> row_traits;
        
        public:
            
            PngGenerator(CsrMatrix const * mat, double threshold);
            ~PngGenerator();
            
            png::byte* get_next_row(size_t pos);
            
        private:
            
            CsrMatrix const * Mat;
            row buffer;
            double Threshold;
    };
    
    /**
     * Save matrix structure as a black-and-white image.
     * @param filename File name.
     * @param threshold Largest absolute value represented by white colour.
     */
    void plot(const char* filename, double threshold = 0.) const;
    
#endif
    
    /**
     * LU factorization
     */
    class LUft
    {
    public:
        
        LUft() :
            __Numeric__(0), __matrix__(0), filename(0) {}
        
        LUft(const LUft& lu) :
            __Numeric__(lu.__Numeric__), __matrix__(lu.__matrix__), filename(lu.filename) {}
        
        LUft(const CsrMatrix* matrix, void* Numeric) : 
            __Numeric__(Numeric), __matrix__(matrix) {}
        
        void free()
        { 
            if (__Numeric__ != 0)
            {
                umfpack_zl_free_numeric(&__Numeric__);
                __Numeric__ = 0;
            }
        }
        
        size_t size() const
        {
            long lnz, unz, m, n, nz_udiag;
            long status = umfpack_zl_get_lunz (
                &lnz, &unz, &m, &n, &nz_udiag, __Numeric__
            );
            return status == 0 ? (lnz + unz) * 16 : 0; // Byte count
        }
        
        cArray solve(cArrayView const & b, unsigned eqs = 1)
        {
            // reserve space for the solution
            cArray x(b.size());
            
            // solve
            solve(b, x, eqs);
            
            // return the result
            return x;
        }
        
        void solve(cArrayView const & b, cArrayView & x, unsigned eqs = 1)
        {
            // solve for all RHSs
            for (size_t eq = 0; eq < eqs; eq++)
            {
                // is there enough RHSs?
                if ((eq + 1) * __matrix__->_n_ > b.size())
                    break;
                
                // solve for current RHS
                long status = umfpack_zl_solve (
                    UMFPACK_Aat,
                    __matrix__->_p_.data(), __matrix__->_i_.data(),
                    reinterpret_cast<const double*>(__matrix__->_x_.data()), 0,
                    reinterpret_cast<double*>(&x[0] + eq * __matrix__->_n_), 0,
                    reinterpret_cast<const double*>(&b[0] + eq * __matrix__->_n_), 0,
                    __Numeric__, 0, 0
                );
                
                // check output
                if (status != 0)
                {
                    std::cerr << "\n[CsrMatrix::LUft::solve] Exit status " << status << "\n";
                    umfpack_zl_report_status(0, status);
                }
            }
        }
        
    private:
        void* __Numeric__;
        const CsrMatrix* __matrix__;
        char* filename;
    };
    
    /**
     * Compute LU factorization. Uses UMFPACK.
     */
    LUft factorize() const;
    
    /**
     * Solve the Ax = b problem, where "b" can be a matrix. Uses UMFPACK.
     * @param b Complex vector containing column-major ordered data; it may be
     *          a flattened matrix.
     * @param eqs Number of columns.
     * @return Array of roots in the same shape as "b".
     */
    cArray solve(cArray const & b, size_t eqs = 1) const;
    
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
            size_t idx1 = _p_[row];
            size_t idx2 = _p_[row + 1];
            
            // loop over columns
            for (size_t idx = idx1; idx < idx2; idx++)
            {
                // which column is this?
                size_t col = _i_[idx];
                
                // transform the output matrix
                A._x_[idx] = f(row, col, _x_[idx]);
            }
        }
        
        return A;
    }
    
private:
    
    // dimensions
    long _m_;
    long _n_;
    
    // representation
    lArray _p_;
    lArray _i_;
    cArray _x_;
    
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
    
    // Friends
    
    friend CooMatrix kron(const CooMatrix& A, const CooMatrix& B);
    
    // Empty constructors
    
    CooMatrix()
        : _m_(0), _n_(0), sorted_(true) {}
    CooMatrix(size_t m, size_t n)
        : _m_(m), _n_(n), sorted_(true) {}
    CooMatrix(CooMatrix const & A)
        : _m_(A._m_), _n_(A._n_), _i_(A._i_), _j_(A._j_), _x_(A._x_), sorted_(false) {}
    CooMatrix(size_t m, size_t n, NumberArray<long> const & i, NumberArray<long> const & j, NumberArray<Complex> const & x)
        : _m_(m), _n_(n), _i_(i), _j_(j), _x_(x), sorted_(false) {}
    
    /**
     * Copy constructor initialized from dense array.
     * @param m Row counz of new matrix.
     * @param n Column count of new matrix.
     * @param a Column-major ordered dense array with matrix elements.
     *          Only nonzero elements are copied into internal storage.
     */
    template <class T> CooMatrix(size_t m, size_t n, T a) : _m_(m), _n_(n), sorted_(false)
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
                    _i_.push_back(row);
                    _j_.push_back(col);
                    _x_.push_back(val);
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
        if (_m_ == 1 and _n_ == 1)
        {
            // get number of nonzero matrix elements
            size_t elems = this->shake()._x_.size();
            if (elems > 1)
                throw "[CooMatrix::operator Complex] more elements than nominal volume!\n";
            else if (elems == 1)
                return this->shake()._x_[0];
            else
                return 0.;
        }
        else
            throw "[CooMatrix::operator Complex] matrix is not 1×1!\n";
    }
    
    // Getters
    
    size_t rows() const { return _m_; }
    size_t cols() const { return _n_; }
    size_t size() const { return _i_.size(); }
    lArray const & i() const { return _i_; }
    lArray const & j() const { return _j_; }
    cArray const & v() const { return _x_; }
    
    Complex operator() (size_t ix, size_t iy) const
    {
        for (size_t n = 0; n < _i_.size(); n++)
             if ((size_t)_i_[n] == ix and (size_t)_j_[n] == iy)
                return _x_[n];
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
        
        for (size_t row = 0; row < _m_; row++)
        {
            for (size_t col = row; col < _n_ and col - row <= d; col++)
            {
                val = f(row,col);
                
                if (val != 0.)
                {
                    _i_.push_back(row);
                    _j_.push_back(col);
                    _x_.push_back(val);
                    
                    if (row != col)
                    {
                        _i_.push_back(col);
                        _j_.push_back(row);
                        _x_.push_back(val);
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
        
        for (size_t row = 0; row < _m_; row++)
        {
            for (size_t col = 0; col < _n_; col++)
            {
                val = f(row,col);
                
                if (val != 0.)
                {
                    _i_.push_back(row);
                    _j_.push_back(col);
                    _x_.push_back(val);
                }
            }
        }
        
        sorted_ = false;
        
        return *this;
    }
    
    // Assignment
    CooMatrix & operator = (CooMatrix const & A)
    {
        _m_ = A._m_;
        _n_ = A._n_;
        _i_ = A._i_;
        _j_ = A._j_;
        _x_ = A._x_;
        
        sorted_ = A.sorted_;
        
        return *this;
    }
    
    /**
     * Addition of an element to matrix. If an existing coordinated are used,
     * the numbers will be summed.
     */
    void add(long i, long j, Complex v)
    {
        _i_.push_back(i);
        _j_.push_back(j);
        _x_.push_back(v);
        
        sorted_ = false;
    }
    
    /// Transposition, implemented as an interchange of "i" and "j" data.
    CooMatrix transpose() const
    {
        CooMatrix tr;
        
        tr._m_ = _n_;
        tr._n_ = _m_;
        tr._i_ = _j_;
        tr._j_ = _i_;
        tr._x_ = _x_;
        
        tr.sorted_ = false;
        
        return tr;
    }
    
    // Addition.
    CooMatrix& operator += (const CooMatrix& A)
    {
        assert(_m_ == A._m_);
        assert(_n_ == A._n_);
        
        _i_.append(A._i_.begin(), A._i_.end());
        _j_.append(A._j_.begin(), A._j_.end());
        _x_.append(A._x_.begin(), A._x_.end());
        
        sorted_ = false;
        
        return *this;
    }
    
    // Subtraction.
    CooMatrix& operator -= (const CooMatrix& A)
    {
        assert(_m_ == A._m_);
        assert(_n_ == A._n_);
        
        size_t prev_size = _x_.size();
        
        _i_.append(A._i_.begin(), A._i_.end());
        _j_.append(A._j_.begin(), A._j_.end());
        _x_.append(A._x_.begin(), A._x_.end());
        
        // negate the newly added elements
        for (size_t i = prev_size; i < _x_.size(); i++)
            _x_[i] = -_x_[i];
        
        sorted_ = false;
        
        return *this;
    }
    
    // Element-wise multiplication by complex number.
    CooMatrix& operator *= (Complex c)
    {
        long nz = _i_.size();
        for (long i = 0; i < nz; i++)
            _x_[i] *= c;
        
        return *this;
    }
    
    // multiplication
    CooMatrix dot (cArrayView const &  B) const;
    
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
    CooMatrix& operator *= (cArray const &  B);
    
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
        _m_ = m;
        _n_ = n;
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
    
    // sort indices (by _i_, then by _j_)
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
    cArray solve(cArray const & b, size_t eqs = 1) const
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
    size_t _m_, _n_;

    // ijv-representation
    lArray _i_, _j_;
    cArray _x_;
    
    bool sorted_;
};

/**
 * @brief Symmetric diagonal matrix.
 * 
 * Only the main and upper diagonals are stored.
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
    SymDiaMatrix(int n, iArrayView const & id, cArrayView const & v);
    
    /// Copy constructor.
    SymDiaMatrix(SymDiaMatrix const & A);

    /// Move constructor.
    SymDiaMatrix(SymDiaMatrix && A);

    /**
     * @brief Plain symmetrical populator.
     *
     * @param d How many upper diagonals to populate. The main diagonal will
     *          be populated always.
     */
    template <class Functor> SymDiaMatrix & populate(unsigned d, Functor f)
    {
        // empty
        elems_.clear();
        
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
    
    //
    // Getters
    //
    
    iArray const & diag() const { return idiag_; }
    cArray const & data() const { return elems_; }
    
    size_t size() const { return n_; }
    
    cArrayView main_diagonal() const { return cArrayView(elems_, 0, n_); }
    
    bool is_compatible (SymDiaMatrix const & B) const;
    
    //
    // Arithmetic and other operators
    //
    
    SymDiaMatrix const & operator = (SymDiaMatrix && A);
    SymDiaMatrix const & operator = (SymDiaMatrix const & A);
    
    friend SymDiaMatrix operator + (SymDiaMatrix const & A, SymDiaMatrix const & B);
    friend SymDiaMatrix operator - (SymDiaMatrix const & A, SymDiaMatrix const & B);
    friend SymDiaMatrix operator * (Complex z, SymDiaMatrix const & A);

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
    cArray dot (cArrayView const & B, MatrixTriangle triangle = both) const;

    /**
     * @brief Kronecker product.
     *
     */
    SymDiaMatrix kron (SymDiaMatrix const & B) const;

    //
    // HDF interface
    //

    bool hdfload(const char* name);

    bool hdfsave(const char* name, bool docompress = false, int consec = 10) const;

    //
    // Conversions to other formats
    //
    
    CooMatrix tocoo(MatrixTriangle triangle = both) const;
    
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
    
    // diagonal right indices starting from zero
    iArray idiag_;
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
inline CooMatrix operator * (const CooMatrix& A, cArray const & B)
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
 * Identity
 * @param N Dimension.
 */
CooMatrix eye(size_t N);

/**
 * Diagonal matrix with diagonal consisting of 0, 1, 2, 3, ..., N - 1.
 * @param N Dimension.
 */
CooMatrix stairs(size_t N);

/**
 * @brief Sparse incomplete Cholesky decomposition.
 * 
 * This routine computes the LDL-decomposition of a symmetric matrix,
 * @f[
 *     A = L D L^T \ ,
 * @f]
 * where @f$ L @f$ is a lower triangular matrix normalized so that it has
 * units on the diagonal and @f$ D @f$ is a diagonal matrix.
 * 
 * @param A Matrix elements in the form of a consecutive array @f$ \left\{a_i\right\}_{i=1}^N @f$
 *          as in
 *          @f[
 *                   \pmatrix
 *                   {
 *                       a_1 &     &     &        &       \cr
 *                       a_2 & a_3 &     &        &       \cr
 *                       a_4 & a_5 & a_6 &        &       \cr
 *                       a_7 & a_8 & a_9 & a_{10} &       \cr
 *                       \vdots &   &     &        & \ddots \cr
 *                   }
 *          @f]
 *          Whenever @f$ a_k @f$ is equal to zero, it is (or can be) omitted from the input array.
 * @param I %Array of column indices (one for every element of A).
 * @param P %Array of row pointers, i.e. starting positions of rows of A. For dense matrix
 *          it would be 0, 1, 3, 6, 10, ... The last element must be equal to the length
 *          of both A and I.
 * 
 * @return The elements of @f$ L @f$ (below diagonal) and @f$ D @f$ (at diagonal) with the
 *         exact sparse pattern as the input array A, i.e. specified by I and P arrays.
 */
cArray iChol(cArrayView const & A, lArrayView const & I, lArrayView const & P);

#endif
