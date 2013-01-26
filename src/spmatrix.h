/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2012                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HEX_SPMATRIX
#define HEX_SPMATRIX

#include <assert.h>
#include <complex>
#include <cstring>
#include <vector>
#include <png++/png.hpp>

#ifdef WITH_MPI
	#include <mpi.h>
#endif

#include <unistd.h>

/*
 * Sparse matrixc operations are done using UMFPACK, of which following
 * functions are used:
 * - umfpack_zl_triplet_to_col
 * - umfpack_zl_col_to_triplet
 * - umfpack_zl_symbolic
 * - umfpack_zl_numeric
 * - umfpack_zl_solve
 * - umfpack_zl_free_symbolic
 * - umfpack_zl_free_numeric
 * - umfpack_zl_report_status
 */
#include <suitesparse/umfpack.h>

#include "arrays.h"
#include "complex.h"

// declaration of classes in order to enable them as return types
// of other classes before proper definition
class CooMatrix;
class CscMatrix;
class CsrMatrix;

/**
 * \brief Complex CSC matrix.
 * 
 * Complex compressed-storage-column-major-ordered sparse matrix.
 * The data are stored in three arrays,
 * \code
 * p, i, x
 * \endcode
 * The format is that used in UMFPACK with interleaved real and comples parts;
 * for explanation see UMFPACK Uses's Guide.
 */
class CscMatrix {
	
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
	CscMatrix(size_t m, size_t n, const std::vector<long>& p, const std::vector<long>& i, const std::vector<Complex>& x)
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
	 * \param b Vector to multiply with.
	 */
	cArray dotT(cArray const &  b) const;
	
	// Getters
	
	size_t size() const { return _i_.size(); }
	size_t rows() const { return _m_; }
	size_t cols() const { return _n_; }
	
	/**
	 * Multiplication by a real number.
	 */
	CscMatrix& operator *= (double r)
	{
		size_t N = _i_.size();
		
		for (size_t i = 0; i < N; i++)
			_x_[i] *= r;
		
		return *this;
	}
	
	/**
	 * Addition of another CSC matrix.
	 * The matrices MUST HAVE THE SAME SPARSE STRUCTURE, as no indices are
	 * checked.
	 */
	CscMatrix& operator &= (const CscMatrix&  B)
	{
		size_t N = _i_.size();
		
		assert(_m_ == B._m_);
		assert(_n_ == B._n_);
		assert(N == B._i_.size());
		
		for (size_t i = 0; i < N; i++)
		{
			assert(_i_[i] == B._i_[i]);
			
			_x_[i] += B._x_[i];
		}
		
		return *this;
	}
	
	/**
	 * Subtraction of another CSC matrix.
	 * The matrices MUST HAVE THE SAME SPARSE STRUCTURE, as no indices are
	 * checked.
	 */
	CscMatrix& operator ^= (const CscMatrix&  B)
	{
		size_t N = _i_.size();
		
		assert(_m_ == B._m_);
		assert(_n_ == B._n_);
		assert(N == B._i_.size());
		
		for (size_t i = 0; i < N; i++)
		{
			assert(_i_[i] == B._i_[i]);
			
			_x_[i] -= B._x_[i];
		}
		
		return *this;
	}
	
	/**
	 * Save matrix to HDF file.
	 * \param name Filename.
	 */
	bool hdfsave(const char* name) const;
	
	/**
	 * Load matrix from HDF file.
	 * \param name Filename.
	 */
	bool hdfload(const char* name);
	
private:
	
	// dimensions
	long _m_;
	long _n_;
	
	// representation
	std::vector<long> _p_;
	std::vector<long> _i_;
	std::vector<Complex> _x_;
	
};

/**
 * \brief Complex CSR matrix.
 * 
 * Complex compressed-storage-row-major-ordered sparse matrix.
 * The data are stored in four arrays,
 * \code
 * p, i, x, z
 * \endcode
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
	CsrMatrix(const CsrMatrix& A)
		: _m_(A._m_), _n_(A._n_), _p_(A._p_), _i_(A._i_), _x_(A._x_) {}
	CsrMatrix(size_t m, size_t n, const std::vector<long>& p, const std::vector<long>& i, const std::vector<Complex>& x)
		: _m_(m), _n_(n), _p_(p), _i_(i), _x_(x) {}
	
	// Destructor
	
	~CsrMatrix() {}
	
	/**
	 * Ordinary matrix-vector product, \f$ A\cdot b \f$.
	 * \param b Vector to multiply with.
	 */
	cArray dot(cArray const &  b) const;
	cArray sdot(cArray const &  b) const;
	
	// Getters
	
	size_t size() const { return _i_.size(); }
	size_t rows() const { return _m_; }
	size_t cols() const { return _n_; }
	const long* p() const { return _p_.data(); }
	const long* i() const { return _i_.data(); }
	const Complex* x() const { return _x_.data(); }
	
	// return absolute value of the (in absolute value) largest element
	double norm() const;
	
	/**
	 * Write the matrix data to a file.
	 * 
	 * \param filename Filename.
	 * 
	 * Format of the fields i,j,x,z is:
	 * \code
	 * %d\t%d\t%g\t%g
	 * \endcode
	 */
	void write(const char* filename) const;
	
	/**
	 * Returns a submatrix.
	 * \param a Starting row.
	 * \param b Ending row plus one.
	 * \param c Starting column.
	 * \param d Ending column plus one.
	 */
	CsrMatrix submatrix(unsigned a, unsigned b, unsigned c, unsigned d) const;
	
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
	 * \param filename File name.
	 * \param threshold Largest absolute value represented by white colour.
	 */
	void plot(const char* filename, double threshold = 0.) const;
	
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
		
		cArray solve(cArray const & b, unsigned eqs = 1)
		{
			// reserve space for the solution
			cArray x(b.size());
			
			// solve for all RHSs
			for (size_t eq = 0; eq < eqs; eq++)
			{
				// is there enough RHSs?
				if ((eq + 1) * __matrix__->_n_ > b.size())
					break;
				
				// solve for current RHS
				long status = umfpack_zl_solve(
					UMFPACK_Aat,
					__matrix__->_p_.data(), __matrix__->_i_.data(),
					reinterpret_cast<const double*>(__matrix__->_x_.data()), 0,
					reinterpret_cast<double*>(x.data() + eq * __matrix__->_n_), 0,
					reinterpret_cast<const double*>(&b[0] + eq * __matrix__->_n_), 0,
					__Numeric__, 0, 0
				);
				
				// check output
				if (status != 0)
				{
					fprintf(stderr, "\n[CsrMatrix::LUft::solve] %ld ", status);
					umfpack_zl_report_status(0, status);
				}
			}
			
			return x;
		}
		
// 	private:
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
	 * \param b Complex vector containing column-major ordered data; it may be
	 *          a flattened matrix.
	 * \param eqs Number of columns.
	 * \return Array of roots in the same shape as "b".
	 */
	cArray solve(cArray const & b, size_t eqs = 1) const;
	
	/**
	 * Structure with information defining a preconditioner.
	 */
	typedef struct {
		
		/// preconditioner identification
		unsigned preconditioner;
		
		/// parameter of SSOR iterations
		double omega;
		
		/// number of equal blocks on diagonal (for block inversion preconditioner)
		unsigned Nblock;
		
	} PreconditionerInfo;
	
	/**
	 * Solve a system of equations with this matrix using conjugate
	 * gradient method.
	 * \param b Vector of right hand sides.
	 * \param x On entry initial guess, on return solution. Use zeros if no
	 *          idea.
	 * \param eps Relative stop precision.
	 * \param min_iterations Minimal number of iterations.
	 * \param max_iterations Maximal number of iterations.
	 * \param pi Preconditioner parameters.
	 * \return Number of iterations.
	 */
	unsigned cg(
		cArray const & b,
		cArray& x,
		double eps,
		unsigned min_iterations,
		unsigned max_iterations,
		const PreconditionerInfo& pi
	) const;
	
	/**
	 * Save matrix to HDF file.
	 * \param name Filename.
	 */
	bool hdfsave(const char* name) const;
	
	/**
	 * Load matrix from HDF file.
	 * \param name Filename.
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
	CsrMatrix& operator *= (Complex r)
	{
		size_t N = _i_.size();
		
		for (size_t i = 0; i < N; i++)
			_x_[i] *= r;
		
		return *this;
	}
	
	/**
	 * Addition of another CSR matrix.
	 * The matrices MUST HAVE THE SAME SPARSE STRUCTURE, as no indices are
	 * checked.
	 */
	CsrMatrix& operator &= (const CsrMatrix&  B)
	{
		size_t N = _i_.size();
		
		// check at least dimensions and non-zero element count
		assert(_m_ == B._m_);
		assert(_n_ == B._n_);
		assert(N == B._i_.size());
		
		for (size_t i = 0; i < N; i++)
			_x_[i] += B._x_[i];
		
		return *this;
	}
	
	/**
	 * Subtraction of another CSR matrix.
	 * The matrices MUST HAVE THE SAME SPARSE STRUCTURE, as no indices are
	 * checked.
	 */
	CsrMatrix& operator ^= (const CsrMatrix&  B)
	{
		size_t N = _i_.size();
		
		assert(_m_ == B._m_);
		assert(_n_ == B._n_);
		assert(N == B._i_.size());
		
		for (size_t i = 0; i < N; i++)
			_x_[i] -= B._x_[i];
		
		return *this;
	}
	
	/**
	 * Sets fill-in elements so that the storage structure of this matrix
	 * will be identical to that of the other CSR matrix. We can the use
	 * “fast” arithmetic operators & and ^.
	 * \return self
	 */
	CsrMatrix sparse_like(const CsrMatrix& B) const;
	
	/**
	 * Return dense array with diagonal elements of the matrix.
	 */
	cArray diag() const;
	
	/**
	 * Solves upper triangular system of equations using backsubstitution.
	 * The matrix needs not be triangular, but elements under the diagonal
	 * won't be used or changed.
	 * \param b Right-hand side.
	 */
	cArray upperSolve(cArray const & b) const;
	
	/**
	 * Solves lower triangular system of equations using backsubstitution.
	 * The matrix needs not be triangular, but elements above the diagonal
	 * won't be used or changed.
	 * \param b Right-hand side.
	 */
	cArray lowerSolve(cArray const & b) const;
	
	/**
	 * Applies a user transformation on <b>nonzero</b> matrix elements.
	 * \param f A functor compatible with following declaration:
	 * \code
	 * Complex (*f) (size_t i, size_t j, Complex z);
	 * \endcode
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
	std::vector<long> _p_;
	std::vector<long> _i_;
	std::vector<Complex> _x_;
	
};


/**
 * \brief Complex COO matrix.
 * 
 * Coordinate (ijv) storage format sparse matrix.
 * The data are stored in four arrays,
 * \code
 * i, j, x, z
 * \endcode
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
		: _m_(0), _n_(0) {}
	CooMatrix(size_t m, size_t n)
		: _m_(m), _n_(n) {}
	CooMatrix(const CooMatrix& A)
		: _m_(A._m_), _n_(A._n_), _i_(A._i_), _j_(A._j_), _x_(A._x_) {}
	CooMatrix(size_t m, size_t n, const std::vector<long>& i, const std::vector<long>& j, const std::vector<Complex>& x)
		: _m_(m), _n_(n), _i_(i), _j_(j), _x_(x) {}
	
	/**
	 * Copy constructor initialized from dense array.
	 * \param m Row counz of new matrix.
	 * \param n Column count of new matrix.
	 * \param a Column-major ordered dense array with matrix elements.
	 *          Only nonzero elements are copied into internal storage.
	 */
	template <class T> CooMatrix(size_t m, size_t n, T a) : _m_(m), _n_(n)
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
	 * \param d Band halfwidth (full width is 2 * d + 1).
	 * \param f Object compatible with signature:
	 * \code
	 *  Complex (*) (long, long)
	 * \endcode
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
		
		return *this;
	}
	
	/**
	 * Full populator
	 * 
	 * Sets all values. The matrix can become easily dense!
	 * \param f Object compatible with signature:
	 * \code
	 *  Complex (*) (long, long)
	 * \endcode
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
		
		return *this;
	}
	
	// Assignment
	CooMatrix& operator = (const CooMatrix& A)
	{
		_m_ = A._m_;
		_n_ = A._n_;
		_i_ = A._i_;
		_j_ = A._j_;
		_x_ = A._x_;
		
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
		
		return tr;
	}
	
	// Addition.
	CooMatrix& operator += (const CooMatrix& A)
	{
		assert(_m_ == A._m_);
		assert(_n_ == A._n_);
		
		_i_.insert(_i_.end(), A._i_.begin(), A._i_.end());
		_j_.insert(_j_.end(), A._j_.begin(), A._j_.end());
		_x_.insert(_x_.end(), A._x_.begin(), A._x_.end());
		
		return *this;
	}
	
	// Subtraction.
	CooMatrix& operator -= (const CooMatrix& A)
	{
		assert(_m_ == A._m_);
		assert(_n_ == A._n_);
		
		size_t prev_size = _x_.size();
		
		_i_.insert(_i_.end(), A._i_.begin(), A._i_.end());
		_j_.insert(_j_.end(), A._j_.begin(), A._j_.end());
		_x_.insert(_x_.end(), A._x_.begin(), A._x_.end());
		
		// negate the newly added elements
		for (size_t i = prev_size; i < _x_.size(); i++)
			_x_[i] = -_x_[i];
		
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
	CooMatrix dot (cArray const &  B);
	
	/**
	 * Matrix multiplication
	 * \param B Dense column-major ordered 1D-array. It's row count
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
	 */
	void resize(size_t m, size_t n)
	{
		_m_ = m;
		_n_ = n;
	}
	
	/**
	 * Change matrix shape. Conserves volume, i.e. it holds
	 * \f[
	 * m \cdot n = m_0 \cdot n_0
	 * \f]
	 * \param m New row count.
	 * \param n New column coount.
	 */
	CooMatrix reshape(size_t m, size_t n) const;
	
	// Convert matrix to dense column-major ordered 1D-array.
	cArray todense() const;
	
	// sort indices (by _i_, then by _j_)
	void sort();
	
	// Convert to CSC matrix.
	CscMatrix tocsc() const;
	
	// Convert to CSR matrix.
	CsrMatrix tocsr() const;
	
	/**
	 * Solve the Ax = b problem, where "b" can be a matrix.
	 * \param b Complex vector containing column-major ordered data; it may be
	 *          a flattened matrix.
	 * \param eqs Number of columns.
	 * \return Array of roots in the same shape as "b".
	 */
	cArray solve(cArray const & b, size_t eqs = 1) const
	{
		// COO format is not optimal for solving -> covert to CSC
		return tocsr().solve(b, eqs);
	}
	
	/**
	 * Write the matrix data to a file.
	 * 
	 * \param filename Filename.
	 * 
	 * Format of the fields i,j,x,z is:
	 * \code
	 * %d\t%d\t%g\t%g
	 * \endcode
	 */
	void write(const char* filename) const;
	
	// Shake the content, i.e. sum same element entries.
	CooMatrix shake() const;
	
	/**
	 * Save matrix to HDF file.
	 * \param name Filename.
	 */
	bool hdfsave(const char* name) const;
	
	/**
	 * Load matrix from HDF file.
	 * \param name Filename.
	 */
	bool hdfload(const char* name);

private:

	// dimensions
	size_t _m_, _n_;

	// ijv-representation
	std::vector<long> _i_, _j_;
	std::vector<Complex> _x_;
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
 * \param A COO-matrix.
 * \param B COO-matrix.
 */
inline CooMatrix operator + (const CooMatrix& A, const CooMatrix& B)
{
	CooMatrix C = A;
	return C += B;
}

/**
 * Computes a difference of two coo-matrices.
 * \param A COO-matrix.
 * \param B COO-matrix.
 */
inline CooMatrix operator - (const CooMatrix& A, const CooMatrix& B)
{
	CooMatrix C = A;
	return C -= B;
}

/**
 * Computes product of a matrix and a complex number.
 * \param z Number.
 * \param B Matrix.
 */
inline CooMatrix operator * (const Complex& z, const CooMatrix& B)
{
	CooMatrix C = B;
	return C *= z;
}

/**
 * Computes product of a sparse matrix and a dense matrix.
 * \param A COO-Matrix.
 * \param B Dense matrix (as a column-major ordered array).
 */
inline CooMatrix operator * (const CooMatrix& A, cArray const & B)
{
	CooMatrix C = A;
	return C *= B;
}

/**
 * Computes product of a matrix and a complex number.
 * \param A COO-Matrix.
 * \param z Number.
 */
inline CooMatrix operator * (const CooMatrix& A, const Complex& z)
{
	CooMatrix C = A;
	return C *= z;
}

/**
 * Kronecker product (“flattened tensor product”).
 * \param A First matrix
 * \param B Second matrix.
 */
CooMatrix kron(const CooMatrix& A, const CooMatrix& B);

/**
 * Identity
 * \param N Dimension.
 */
CooMatrix eye(size_t N);

/**
 * Diagonal matrix with diagonal consisting of 0, 1, 2, 3, ..., N - 1.
 * \param N Dimension.
 */
CooMatrix stairs(size_t N);

/**
 * Callback-based Conjugate gradients.
 * \param b Right hand side.
 * \param x Output vector. An initial guess may be present at beginning.
 * \param eps Relative tolerance.
 * \param min_iterations Minimal iterations count.
 * \param max_iterations Maximal iterations count.
 * \param apply_preconditioner Functor compatible with void(*)(const Array&, Array&) prototype.
 *                             Should apply a custom preconditioner to a given vector.
 * \param matrix_multiply Functor compatible with void(*)(const Array&, Array&) prototype.
 *                        Should multiply the given vector by the matrix of the equation set.
 * \return Iteration count.
 */
template <typename TFunctor1, typename TFunctor2>
unsigned cg_callbacks(
	cArray const & b, cArray & x,
	double eps,
	unsigned min_iterations, unsigned max_iterations,
	TFunctor1 apply_preconditioner,
	TFunctor2 matrix_multiply
) {
	// some arrays (search directions etc.)
	size_t N = b.size();
	cArray p(N), q(N), z(N);
	
	// residual; initialized to starting residual using the initial guess
	cArray r(N);
#ifdef WITH_MPI
	fprintf(stdout, "\t[cg:%d] Computing residual...\n", MPI::COMM_WORLD.Get_rank()); fflush(stdout);
#endif
	matrix_multiply(x, r);
	r = b - r;

	// some other scalar variables
	Complex rho_new;		// contains inner product r_i^T · r_i
	Complex rho_old;		// contains inner product r_{i-1}^T · r_{i-1}
	Complex alpha, beta;	// contains projection ratios
	
	// Iterate
	
	unsigned k;
	for/*ever*/ (k = 0; ; k++)
	{
#ifdef WITH_MPI
		fprintf(stdout, "\t[cg:%d] Residual relative magnitude: %g\n", MPI::COMM_WORLD.Get_rank(), r.norm() / b.norm()); fflush(stdout);
#else
		fprintf(stdout, "\t[cg] Residual relative magnitude after %d iterations: %g\n", k, r.norm() / b.norm()); fflush(stdout);
#endif
		
		// apply desired preconditioner
#ifdef WITH_MPI
		fprintf(stdout, "\t[cg:%d] Applying preconditioner\n", MPI::COMM_WORLD.Get_rank()); fflush(stdout);
#endif
		apply_preconditioner(r, z);
		
		// compute projection ρ = r·z
		rho_new = (r|z);
		
		// setup search direction p
		if (k == 0)
		{
			p = z;
		}
		else
		{
			beta = rho_new / rho_old;
			p = z + beta * p;
		}
		
		// move to next Krylov subspace by multiplying A·p
#ifdef WITH_MPI
		fprintf(stdout, "\t[cg:%d] Krylov space update...\n", MPI::COMM_WORLD.Get_rank()); fflush(stdout);
#endif
		matrix_multiply(p, q);
		
		// compute projection ratio α
		alpha = rho_new / (p|q);
		
		// update the solution and the residual
		x += alpha * p;
		r -= alpha * q;
		
		// once in a while check convergence but do at least "min_iterations" iterations
		if (k >= min_iterations and k % 4 == 0 and r.norm() / b.norm() < eps)
			break;
		
		// check iteration limit (stop at "max_iterations" iterations)
		if (k >= max_iterations)
		{
#ifdef WITH_MPI
			fprintf(stdout, "\t[cg:%d] Iteration limit %d reached.\n", MPI::COMM_WORLD.Get_rank(), max_iterations); fflush(stdout);
#endif
			break;
		}
		
		// move to the next iteration: store previous projection
		rho_old = rho_new;
	}
	
	return k;
}

/**
 * Callback-based BiCG-STAB.
 * \param b Right hand side.
 * \param x Output vector. An initial guess may be present at beginning.
 * \param eps Relative tolerance.
 * \param min_iterations Minimal iterations count.
 * \param max_iterations Maximal iterations count.
 * \param apply_preconditioner Functor compatible with void(*)(const Array&, Array&) prototype.
 *                             Should apply a custom preconditioner to a given vector.
 * \param matrix_multiply Functor compatible with void(*)(const Array&, Array&) prototype.
 *                        Should multiply the given vector by the matrix of the equation set.
 * \return Iteration count.
 */
template <typename TFunctor1, typename TFunctor2>
unsigned bicgstab_callbacks(
	cArray const & b, cArray & x,
	double eps,
	unsigned min_iterations, unsigned max_iterations,
	TFunctor1 apply_preconditioner,
	TFunctor2 matrix_multiply
) {
	int N = b.size();
	
	cArray r_im1(N), rt(N), p_i(N), p_im1(N), v_i(N), v_im1(N), phat(N), s(N), shat(N), t(N), x_im1(N), x_i(N), r_i(N);
	Complex rho_im1, rho_im2, alpha_i, alpha_im1, omega_i, omega_im1;
	
	x_im1 = x;
	matrix_multiply(x_im1,r_im1);
	rt = r_im1 = b - r_im1;
	
	for/*ever*/ (int i = 1; ; i++)
	{
		rho_im1 = (rt | r_im1);
		if (rho_im1 == 0.)
		{
			fprintf(stderr, "BiCG failed, ρ = 0\n");
			return i-1;
		}
		
		if (i == 1)
		{
			p_i = r_im1;
		}
		else
		{
			Complex beta_im1 = (rho_im1 / rho_im2) * (alpha_im1 / omega_im1);
			p_i = r_im1 + beta_im1 * (p_im1 - omega_im1 * v_im1);
		}
		
		apply_preconditioner(p_i, phat);
		matrix_multiply(phat, v_i);
		alpha_i = rho_im1 / (rt|v_i);
		s = r_im1 - alpha_i * v_i;
		
		fprintf(stdout, "\t[bcg] s-residual relative magnitude after %d iterations: %g\n", i, s.norm() / b.norm()); fflush(stdout);
		
		if (s.norm() < eps * b.norm())
		{
			x = x_im1 + alpha_i * phat;
			return i-1;
		}
		
		apply_preconditioner(s, shat);
		matrix_multiply(shat, t);
		omega_i = (t|s) / (t|t);
		x_i = x_im1 + alpha_i * phat + omega_i * shat;
		r_i = s - omega_i * t;
		
		fprintf(stdout, "\t[bcg] r-residual relative magnitude after %d iterations: %g\n", i, r_i.norm() / b.norm()); fflush(stdout);
		
		if (r_i.norm() < eps * b.norm())
		{
			x = x_i;
			return i;
		}
		
		if (omega_i == 0.)
		{
			fprintf(stderr, "BiCG failed, ω = 0\n");
			return i;
		}
		
		// update
		x_im1 = x_i;
		r_im1 = r_i;
		rho_im2 = rho_im1;
		p_im1 = p_i;
		v_im1 = v_i;
		alpha_im1 = alpha_i;
		omega_im1 = omega_i;
	}
}

#endif