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

#include <algorithm>
#include <complex>
#include <cstdio>
#include <cstring>
#include <set>
#include <vector>

#ifdef WITH_PNGPP
	#include <png++/png.hpp>
#endif

#include "hdffile.h"

// preconditioners
#define P_NONE			0
#define P_JACOBI		1
#define P_SSOR			2
#define P_BLOCK_INV		3

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
#include <umfpack.h>

#include "arrays.h"
#include "spmatrix.h"

CooMatrix kron(const CooMatrix& A, const CooMatrix& B)
{
	// shorthands
	size_t Csize = A._i_.size() * B._i_.size();
	size_t Brows = B._m_;
	size_t Bcols = B._n_;
	
	// create temporary matrix to hold the Kronecker product
	CooMatrix C;
	
	// set correct dimensions, pre-allocate space
	C._m_ = A._m_ * B._m_;
	C._n_ = A._n_ * B._n_;
	C._i_ = std::vector<long>(Csize);
	C._j_ = std::vector<long>(Csize);
	C._x_ = std::vector<Complex>(Csize);
	
	// get iterators
	size_t ic = 0;
	
	// loop over A data
	size_t Asize = A._i_.size();
	for (size_t ia = 0; ia < Asize; ia++)
	{
		// loop over B data
		size_t Bsize = B._i_.size();
		for (size_t ib = 0; ib < Bsize; ib++)
		{
			// compute new row index
			C._i_[ic] = A._i_[ia] * Brows + B._i_[ib];
			
			// compute new column index
			C._j_[ic] = A._j_[ia] * Bcols + B._j_[ib];
			
			// compute product of the two elements
			C._x_[ic] = A._x_[ia] * B._x_[ib];
			
			// move to next value of C
			ic++;
		}
	}
	
	return C;
}

CooMatrix eye(size_t N)
{
	return CooMatrix(N,N).symm_populate_band(
		0,
		[](size_t i, size_t j) -> Complex {
			return 1.;
		}
	);
}

CooMatrix stairs(size_t N)
{
	return CooMatrix(N,N).symm_populate_band(
		0,
		[](size_t i, size_t j) -> Complex {
			return (i == j) ? i : 0.;
		}
	);
}



// ------------------------------------------------------------------------- //

// -- CSC matrix methods --------------------------------------------------- //

// ------------------------------------------------------------------------- //

CscMatrix & CscMatrix::operator *= (double r)
{
	size_t N = _i_.size();
	
	for (size_t i = 0; i < N; i++)
		_x_[i] *= r;
	
	return *this;
}

CscMatrix & CscMatrix::operator &= (const CscMatrix&  B)
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

CscMatrix & CscMatrix::operator ^= (const CscMatrix&  B)
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

cArray CscMatrix::dotT(const cArrayView&  b) const
{
	// create output array
	cArray c (_n_);
	
	// the matrix "*this" is actually transposed
	for (unsigned icol = 0; icol < _n_; icol++)
	{
		size_t idx1 = _p_[icol];
		size_t idx2 = _p_[icol+1];
		
		// for all nonzero elements in this column
		for (size_t idx = idx1; idx < idx2; idx++)
		{
			// get row number
			unsigned irow = _i_[idx];
			
			// store product
			c[icol] += _x_[idx] * b[irow];
		}
	}
	
	return c;
}

CooMatrix CscMatrix::tocoo() const
{
	// reserve space for the auxiliary (__j) and the output (Ti,Tj,Tx,Tz) arrays
	size_t N = _x_.size();
	std::vector<long> Ti(N), Tj(N), __j(N);
	std::vector<Complex> Tx(N);
	
	// do we have any elements at all?
	if (N != 0)
	{
	
		// do the conversion
		long status = umfpack_zl_col_to_triplet(_n_, _p_.data(), __j.data());
		
		// check success
		if (status != 0)
		{
			fprintf(stderr, "\n[CscMatrix::tocoo] %ld ", status);
			umfpack_zl_report_status(0, status);
		}
		
		// copy only non-zero entries to output arrays
		size_t nz = 0;
		for (size_t i = 0; i < N; i++)
		{
			if (_x_[i] != 0.)
			{
				Ti[nz] = _i_[i];
				Tj[nz] = __j[i];
				Tx[nz] = _x_[i];
				nz++;
			}
		}
		
		// crop the output arrays
		Ti.resize(nz);
		Tj.resize(nz);
		Tx.resize(nz);
		
	}
	
	// return new CooMatrix
	return CooMatrix (_m_, _n_, Ti, Tj, Tx);
}

bool CscMatrix::hdfsave(const char* name) const
{
#ifndef NO_HDF
	try
	{
		HDFFile file(name, HDFFile::overwrite);
		
		// write dimensions
		file.write("m", &_m_, 1);
		file.write("n", &_n_, 1);
		
		// write indices
		if (not _p_.empty())
			file.write("p", &(_p_[0]), _p_.size());
		if (not _i_.empty())
			file.write("i", &(_i_[0]), _i_.size());
		
		// write complex data as a "double" array
		if (not _x_.empty())
		{
			file.write (
				"x",
				reinterpret_cast<double const*>( &(_x_[0]) ),
				_x_.size() * 2
			);
		}
		
		return true;
	}
	catch (...)
	{
		return false;
	}
#else /* NOHDF */
	return false;
#endif
}

bool CscMatrix::hdfload(const char* name)
{
#ifndef NO_HDF
	try
	{
		HDFFile hdf(name, HDFFile::readonly);
		
		// read dimensions
		hdf.read("m", &_m_, 1);
		hdf.read("n", &_n_, 1);
		
		// read indices
		if (_p_.resize(hdf.size("p")))
			hdf.read("p", &(_p_[0]), _p_.size());
		if (_i_.resize(hdf.size("i")))
			hdf.read("i", &(_i_[0]), _i_.size());
		
		// read data
		if (_x_.resize(hdf.size("x") / 2))
		{
			hdf.read (
				"x",
				reinterpret_cast<double*>(&(_x_[0])),
				_x_.size() * 2
			);
		}
		
		return true;
	}
	catch (...)
	{
		return false;
	}
#else /* NO_HDF */
	return false;
#endif
}


// ------------------------------------------------------------------------- //

// -- CSR matrix methods --------------------------------------------------- //

// ------------------------------------------------------------------------- //

CsrMatrix & CsrMatrix::operator *= (Complex r)
{
	size_t N = _i_.size();
	
	for (size_t i = 0; i < N; i++)
		_x_[i] *= r;
	
	return *this;
}

CsrMatrix & CsrMatrix::operator &= (CsrMatrix const &  B)
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

CsrMatrix & CsrMatrix::operator ^= (CsrMatrix const &  B)
{
	size_t N = _i_.size();
	
	assert(_m_ == B._m_);
	assert(_n_ == B._n_);
	assert(N == B._i_.size());
	
	for (size_t i = 0; i < N; i++)
		_x_[i] -= B._x_[i];
	
	return *this;
}

cArray CsrMatrix::dot(const cArrayView& b) const
{
	// create output array
	cArray c(_m_);
	
	for (unsigned irow = 0; irow < _m_; irow++)
	{
		size_t idx1 = _p_[irow];
		size_t idx2 = _p_[irow+1];
		
		// for all nonzero elements in this row
		for (size_t idx = idx1; idx < idx2; idx++)
		{
			// get column number
			unsigned icol = _i_[idx];
			
			// store product
			c[irow] += _x_[idx] * b[icol];
		}
	}
	
	return c;
}

#ifdef WITHPNG
CsrMatrix::PngGenerator::PngGenerator(const CsrMatrix* mat, double threshold)
	: base_t(mat->cols(), mat->rows()), Mat(mat), buffer(mat->cols()), Threshold(threshold)
{
}

CsrMatrix::PngGenerator::~PngGenerator()
{
}

png::byte* CsrMatrix::PngGenerator::get_next_row(size_t irow)
{
	// get column indices
	int idx_min = Mat->_p_[irow];
	int idx_max = Mat->_p_[irow + 1];
	
	// clear memory
	for (int icol = 0; icol < (int)Mat->cols(); icol++)
		buffer[icol] = 1;
	
	// for all nonzero columns
	for (int idx = idx_min; idx < idx_max; idx++)
		if (abs(Mat->_x_[idx]) > Threshold)
			buffer[Mat->_i_[idx]] = 0;
		
	// pass the buffer
	return reinterpret_cast<png::byte*>(row_traits::get_data(buffer));
}


void CsrMatrix::plot(const char* filename, double threshold) const
{
	// create output file
	std::ofstream out(filename, std::ios_base::out | std::ios_base::binary);
	
	// create PNG data generator
	PngGenerator png(this, threshold);
	
	// write PNG file
	png.write(out);
}
#endif

CsrMatrix::LUft CsrMatrix::factorize() const
{
	// Use standard UMFPACK sequence
	void *Symbolic, *Numeric;
	long status;
	
	// analyze the sparse structure
	status = umfpack_zl_symbolic(
		_m_, _n_,					// matrix dimensions
		_p_.data(), _i_.data(),		// column and row indices
		reinterpret_cast<const double*>(_x_.data()), 0,	// matrix data
		&Symbolic, 0, 0				// UMFPACK internals
	);
	if (status != 0)
	{
		fprintf(stderr, "\n[CscMatrix::factorize] %ld ", status);
		umfpack_zl_report_status(0, status);
		abort();
	}
	
	// do some factorizations
	status = umfpack_zl_numeric(
		_p_.data(), _i_.data(),	// column and row indices
		reinterpret_cast<const double*>(_x_.data()), 0,	// matrix data
		Symbolic, &Numeric, 0, 0	// UMFPACK internals
	);
	if (status != 0)
	{
		fprintf(stderr, "\n[CscMatrix::factorize] %ld ", status);
		umfpack_zl_report_status(0, status);
		abort();
	}
	
	// release unused data
	umfpack_zl_free_symbolic(&Symbolic);
	return LUft(this, Numeric);
}

cArray CsrMatrix::solve(const cArray&  b, size_t eqs) const
{
	// only square matrices are allowed
	assert(_m_ == _n_);
	
	// compute the LU factorization
	LUft luft = factorize();
	  
	// solve the equations
	cArray solution = luft.solve(b, eqs);
	  
	// cleanup and return
	luft.free();
	return solution;
}

void CsrMatrix::write(const char* filename) const
{
	FILE *f = fopen(filename, "w");
	fprintf(f, "# Matrix %ld × %ld with %ld nonzero elements:\n\n", _m_, _n_, _x_.size());
	for (unsigned irow = 0; irow < _m_; irow++)
	{
		size_t idx1 = _p_[irow];
		size_t idx2 = _p_[irow + 1];
		
		for (size_t idx = idx1; idx < idx2; idx++)
			fprintf(f, "%d\t%ld\t%g\t%g\n", irow, _i_[idx], _x_[idx].real(), _x_[idx].imag());
	}
	fclose(f);
}

bool CsrMatrix::hdfsave(const char* name) const
{
#ifndef NO_HDF
	try
	{
		HDFFile hdf(name, HDFFile::overwrite);
		
		// write dimensions
		hdf.write("m", &_m_, 1);
		hdf.write("n", &_n_, 1);
		
		// write indices
		if (not _p_.empty())
			hdf.write("p", &(_p_[0]), _p_.size());
		if (not _i_.empty())
			hdf.write("i", &(_i_[0]), _i_.size());
		
		// write data
		if (not _x_.empty())
		{
			hdf.write (
				"x",
				reinterpret_cast<double const*>(&(_x_[0])),
				_x_.size() * 2
			);
		}
		
		return true;
	}
	catch (...)
	{
		return false;
	}
#else /* NO_HDF */
	return false;
#endif
}

bool CsrMatrix::hdfload(const char* name)
{
#ifndef NO_HDF
	try
	{
		HDFFile hdf(name, HDFFile::readonly);
		
		// read dimensions
		hdf.read("m", &_m_, 1);
		hdf.read("n", &_n_, 1);
		
		// read indices
		if (_p_.resize(hdf.size("p")))
			hdf.read("p", &(_p_[0]), _p_.size());
		if (_i_.resize(hdf.size("i")))
			hdf.read("i", &(_i_[0]), _i_.size());
		
		// read data
		if (_x_.resize(hdf.size("x") / 2))
		{
			hdf.read (
				"x",
				reinterpret_cast<double*>(&(_x_[0])),
				_x_.size() * 2
			);
		}
		
		return true;
	}
	catch (...)
	{
		return false;
	}
#else /* NO_HDF */
	return false;
#endif
}

double CsrMatrix::norm() const
{
	size_t N = _i_.size();
	double res = 0.;
	
	// return the abs(largest element)
	for (size_t i = 0; i < N; i++)
	{
		// compute the absolute value
		double val = abs(_x_[i]);
		
		// update the winner
		if (val > res)
			res = val;
	}
	
	return res;
}

cArray CsrMatrix::upperSolve(cArrayView const &  b) const
{
	// check size
	size_t N = b.size();
	assert((size_t)_m_ == N);
	assert((size_t)_n_ == N);

	// create output array
	cArray x(N);
	
	// loop over rows
	for (size_t i = 0; i < N; i++)
	{
		size_t row = N - 1 - i;
		Complex accum = 0.;
		
		// get relevant columns of the sparse matrix
		size_t idx1 = _p_[row];
		size_t idx2 = _p_[row + 1];
		
		// diagonal element of the matrix
		Complex a = 0.;
		
		// loop over the columns
		for (size_t idx = idx1; idx < idx2; idx++)
		{
			// which column is this?
			size_t col = _i_[idx];
			
			// diagonal element will be useful in a moment, store it
			if (col == row)
				a = _x_[idx];
			
			// backsubstitute
			else if (col > row)
				accum += _x_[idx] * x[col];
		}
		
		// triangular matrix, in order to be regular, needs nonzero diagonal elements
		assert(a != 0.);
		
		// compute and store the new root
		x[row] = (b[row] - accum) / a;
	}

	return x;
}

cArray CsrMatrix::lowerSolve(cArrayView const & b) const
{
	// check size
	size_t N = b.size();
	assert((size_t)_m_ == N);
	assert((size_t)_n_ == N);
	
	// create output array
	cArray x(N);
	
	// loop over rows
	for (size_t row = 0; row < N; row++)
	{
		Complex accum = 0.;
		
		// get relevant columns of the sparse matrix
		size_t idx1 = _p_[row];
		size_t idx2 = _p_[row + 1];
		
		// diagonal element of the matrix
		Complex a = 0.;
		
		// loop over the columns
		for (size_t idx = idx1; idx < idx2; idx++)
		{
			// which column is this?
			size_t col = _i_[idx];
			
			// diagonal element will be useful in a moment, store it
			if (col == row)
				a = _x_[idx];
			
			// backsubstitute
			else if (col < row)
				accum += _x_[idx] * x[col];
		}
		
		// triangular matrix, in order to be regular, needs nonzero diagonal elements
		assert(a != 0.);
		
		// compute and store the new root
		x[row] = (b[row] - accum) / a;
	}
	
	return x;
}

cArray CsrMatrix::diag() const
{
	cArray D ( std::min(_m_, _n_) );
	
	for (size_t irow = 0; irow < (size_t)_m_; irow++)
		for (size_t idx = _p_[irow]; idx < (size_t)_p_[irow+1]; idx++)
			if ((size_t)_i_[idx] == irow)
				D[irow] = _x_[idx];
	
	return D;
}

CsrMatrix CsrMatrix::sparse_like(const CsrMatrix& B) const
{
	// check dimensions
	assert(_m_ == B._m_);
	assert(_n_ == B._n_);
	
	// prepare zero matrix with the same storage pattern the matrix B has
	CsrMatrix A = B;
	memset(A._x_.data(), 0, A._x_.size() * sizeof(Complex));
	
	// copy all nonzero elements of "this" matrix
	for (unsigned row = 0; row < _m_; row++)
	{
		size_t idx1 = A._p_[row];
		size_t idx2 = A._p_[row+1];
		
		for (size_t idx = idx1; idx < idx2; idx++)
		{
			unsigned col = A._i_[idx];
			A._x_[idx] = (*this)(row,col);
		}
	}
	
	// return temporary matrix
	return A;
}

Complex CsrMatrix::operator() (unsigned i, unsigned j) const
{
	// get all column indices, which have nonzero element in row "i"
	size_t idx1 = _p_[i];
	size_t idx2 = _p_[i + 1];
	
	// find the correct column ("j")
	auto it = std::lower_bound(_i_.begin() + idx1, _i_.begin() + idx2, j);

	if (it == _i_.end() or *it != j)
		return 0.;
	
	// return the value
	return _x_[it - _i_.begin()];
}


// ------------------------------------------------------------------------- //

// -- COO matrix methods --------------------------------------------------- //

// ------------------------------------------------------------------------- //



CscMatrix CooMatrix::tocsc() const
{
	size_t nz = _x_.size();
	
	// CSC matrix data
	std::vector<long> Ap(_n_ + 1), Ai(nz);
	std::vector<Complex> Ax(nz);
	
	// do we have any elements at all?
	if (nz != 0)
	{
	
		long status = umfpack_zl_triplet_to_col (
			_m_,			// rows
			_n_,			// cols
			nz,				// data length
			_i_.data(),		// row indices
			_j_.data(), 	// column indices
			reinterpret_cast<const double *>(_x_.data()), 0, 	// interleaved data
			Ap.data(),		// column pointers
			Ai.data(),		// row pointers
			reinterpret_cast<double *>(Ax.data()), 0,	// interleaved data
			0
		);
		
		// check success
		if (status != 0)
		{
			fprintf(stderr, "\n[CooMatrix::tocsc] %ld ", status);
			umfpack_zl_report_status(0, status);
		}
		
		// crop storage
		size_t N = Ap[_n_];
		Ai.resize(N);
		Ax.resize(N);
	}
	
	return CscMatrix(_m_, _n_, Ap, Ai, Ax);
}


CsrMatrix CooMatrix::tocsr() const
{
	size_t nz = _x_.size();
	
	// CSC matrix data
	std::vector<long> Ap(_n_ + 1), Ai(nz);
	std::vector<Complex> Ax(nz);
	
	// do we have any elements at all?
	if (nz != 0)
	{
	
		long status = umfpack_zl_triplet_to_col (
			_n_,			// cols (rows of transposed matrix)
			_m_,			// rows (cols of transposed matrix)
			nz,				// data length
			_j_.data(), 	// column indices (rows of transposed matrix)
			_i_.data(),		// row indices (cols of transposed matrix)
			reinterpret_cast<const double *>(_x_.data()), 0,	// interleaved data
			Ap.data(),		// row pointers
			Ai.data(),		// column indices
			reinterpret_cast<double *>(Ax.data()), 0,	// interleaved data
			0
		);
		
		// check success
		if (status != 0)
		{
			fprintf(stderr, "\n[CooMatrix::tocsr] %ld ", status);
			umfpack_zl_report_status(0, status);
		}
		
		// crop storage
		size_t N = Ap[_m_];
		Ai.resize(N);
		Ax.resize(N);
	}
	
	return CsrMatrix(_m_, _n_, Ap, Ai, Ax);
}

SymDiaMatrix CooMatrix::todia() const
{
	// diagonal indices
	std::set<int> diags;
	
	// elements per diagonal
	cArrays elems(_m_);
	
	// for all nonzero elements
	for (size_t i = 0; i < _x_.size(); i++)
	{
		// get row and column index
		int irow = _i_[i];
		int icol = _j_[i];
		
		// get diagonal
		int diagonal = icol - irow;
		
		// we only want to store the upper triangle
		if (diagonal < 0)
			continue;
		
		// add this diagonal (if not already added)
		diags.insert(diagonal);
		
		// reserve space for this diagonal if not already done
		if (elems[diagonal].size() == 0)
			elems[diagonal].resize(_m_ - diagonal);
		
		// add element
		elems[diagonal][irow] = _x_[i];
	}
	
	// concatenate the diagonals, construct matrix object
	return SymDiaMatrix (
		_m_,
		Array<int>(diags.cbegin(), diags.cend()),
		join(elems)
	);
}

void CooMatrix::write(const char* filename) const
{
	std::ofstream f(filename);
	
	f << "# Matrix " << _m_ << " × " << _n_ << " with " << _x_.size() << " nonzero elements:\n\n";
	
	for (size_t i = 0; i < _i_.size(); i++)
		f << _i_[i] << "\t" << _j_[i] << "\t" << _x_[i].real() << "\t" << _x_[i].imag() << "\n";
}


CooMatrix CooMatrix::reshape(size_t m, size_t n) const
{
	CooMatrix C = *this;
	
	// conserved dimensions
	size_t N = C._i_.size();
	size_t H = C._m_;
	
	// check dimensions
	assert(m * n == C._m_ * C._n_);
	
	// reshape
	for (size_t i = 0; i < N; i++)
	{
		// conserved position in column-ordered array
		size_t idx = C._i_[i] + C._j_[i] * H;
		
		// new coordinates
		size_t row = idx % m;
		size_t col = idx / m;
		
		// update values
		C._i_[i] = row;
		C._j_[i] = col;
	}
	
	C._m_ = m;
	C._n_ = n;
	return C;
}

cArray CooMatrix::todense() const
{
	// return column-major dense representations
	cArray v (_m_ * _n_);
	
	size_t N = _i_.size();
	for (size_t i = 0; i < N; i++)
		v[_i_[i] + _j_[i] * _m_] += _x_[i];
	
	return v;
}

CooMatrix& CooMatrix::operator *= (cArray const &  B)
{
	return *this = this->dot(B);
}

CooMatrix CooMatrix::dot(cArrayView const & B) const
{
	// FIXME: This is a memory INEFFICIENT method.
	// NOTE: Row-major storage assumed for B.
	
	// volumes
	size_t A_vol = _x_.size();
	size_t B_vol = B.size();
	
	// check B shape
	assert(B_vol % _n_ == 0);
	
	// create output matrix
	unsigned C_rows = _m_;
	unsigned C_cols = B_vol / _n_;
	CooMatrix C(C_rows, C_cols);
	
	// for all elements of A
	for (size_t i = 0; i < A_vol; i++)
	{
		unsigned row = _i_[i];
		unsigned col = _j_[i];
		
		// for all columns of B
		for (unsigned icol = 0; icol < C_cols; icol++)
		{
			C.add (
				row, icol,
				_x_[i] * B[col*C_cols + icol]
			);
		}
	}
	
	// summation is done by shaking
	return C.shake();
}

Complex CooMatrix::ddot(CooMatrix const & B) const
{
	assert(_m_ == B._m_);
	assert(_n_ == B._n_);
	
	// sort by _i_ and _j_
	if (not sorted() or not B.sorted())
		throw exception("[CooMatrix] Sort matrices before ddot!");
		
	Complex result = 0;
	
	auto Ai = _i_.begin();
	auto Aj = _j_.begin();
	auto Av = _x_.begin();
	
	auto Bi = B._i_.begin();
	auto Bj = B._j_.begin();
	auto Bv = _x_.begin();
	
	while (Av != _x_.end() and Bv != B._x_.end())
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

CooMatrix CooMatrix::shake() const
{
	// ugly and memory inefficient method... FIXME
	return tocsc().tocoo();
}

bool CooMatrix::hdfsave(const char* name) const
{
#ifndef NO_HDF
	try
	{
		HDFFile hdf(name, HDFFile::overwrite);
		
		// write dimensions
		hdf.write("m", &_m_, 1);
		hdf.write("n", &_n_, 1);
		
		// write indices
		if (not _i_.empty())
			hdf.write("i", &(_i_[0]), _i_.size());
		if (not _j_.empty())
			hdf.write("j", &(_j_[0]), _j_.size());
		
		// write data
		if (not _x_.empty())
		{
			hdf.write (
				"x",
				reinterpret_cast<double const*>(&_x_),
				_x_.size() * 2
			);
		}
		
		return true;
	}
	catch (...)
	{
		return false;
	}
#else /* NO_HDF */
	return false;
#endif
}

bool CooMatrix::hdfload(const char* name)
{
	sorted_ = false;
	
#ifndef NO_HDF
	try
	{
		HDFFile hdf(name, HDFFile::readonly);
		
		// read dimensions
		hdf.read("m", &_m_, 1);
		hdf.read("n", &_n_, 1);
		
		// read indices
		if (_i_.resize(hdf.size("i")))
			hdf.read("i", &(_i_[0]), _i_.size());
		if (_j_.resize(hdf.size("j")))
			hdf.read("j", &(_j_[0]), _j_.size());
		
		// read data
		if (_x_.resize(hdf.size("x") / 2))
		{
			hdf.read (
				"x",
				reinterpret_cast<double*>(&(_x_[0])),
				_x_.size() * 2
			);
		}
		
		return true;
	}
	catch (...)
	{
		return false;
	}
#else /* NO_HDF */
	return false;
#endif
}

SymDiaMatrix::SymDiaMatrix() : n_(0), elems_(0), idiag_(0) {}

SymDiaMatrix::SymDiaMatrix(int n) : n_(n), elems_(0), idiag_(0) {}

SymDiaMatrix::SymDiaMatrix(int n, ArrayView<int> const & id, ArrayView<Complex> const & v)
	: n_(n), elems_(v), idiag_(id)
{
// 	std::cout << "Constructing from arrays.\n";
}

SymDiaMatrix::SymDiaMatrix(SymDiaMatrix const & A)
	: n_(A.n_), elems_(A.elems_), idiag_(A.idiag_)
{
// 	std::cout << "Constructing from l-value const reference.\n";
}

SymDiaMatrix::SymDiaMatrix(SymDiaMatrix&& A)
{
// 	std::cout << "Constructing from r-value reference.\n";
	
	n_ = std::move(A.n_);
	elems_ = std::move(A.elems_);
	idiag_ = std::move(A.idiag_);
}

SymDiaMatrix const & SymDiaMatrix::operator += (SymDiaMatrix const & B)
{
	if (is_compatible(B))
	{
		for (size_t i = 0; i < elems_.size(); i++)
			elems_[i] += B.elems_[i];
	}
	return *this;
}

SymDiaMatrix const & SymDiaMatrix::operator -= (SymDiaMatrix const & B)
{
	if (is_compatible(B))
	{
		for (size_t i = 0; i < elems_.size(); i++)
			elems_[i] -= B.elems_[i];
	}
	return *this;
}

SymDiaMatrix const & SymDiaMatrix::operator = (SymDiaMatrix&& A)
{
// 	std::cout << "[SymDiaMatrix::operator=] (from r-value reference)\n";
	
	n_ = std::move(A.n_);
	elems_ = std::move(A.elems_);
	idiag_ = std::move(A.idiag_);
	return *this;
}

SymDiaMatrix const & SymDiaMatrix::operator = (SymDiaMatrix const & A)
{
// 	std::cout << "[SymDiaMatrix::operator=] (from l-value const reference)\n";
	
	n_ = A.n_;
	elems_ = A.elems_;
	idiag_ = A.idiag_;
	return *this;
}

SymDiaMatrix operator + (SymDiaMatrix const & A, SymDiaMatrix const & B)
{
	A.is_compatible(B);
	return SymDiaMatrix(A.n_, A.idiag_, A.elems_ + B.elems_);
}

SymDiaMatrix operator - (SymDiaMatrix const & A, SymDiaMatrix const & B)
{
	A.is_compatible(B);
	return SymDiaMatrix(A.n_, A.idiag_, A.elems_ - B.elems_);
}

SymDiaMatrix operator * (Complex z, SymDiaMatrix const & A)
{
	return SymDiaMatrix(A.n_, A.idiag_, z * A.elems_);
}

bool SymDiaMatrix::is_compatible(const SymDiaMatrix& B) const
{
	if (n_ != B.n_)
		throw exception ("[SymDiaMatrix::operator+=] Unequal ranks.");
	if (idiag_.size() != B.idiag_.size())
		throw exception ("[SymDiaMatrix::operator+=] Unequal number of diagonals (%d != %d).", idiag_.size(), B.idiag_.size());
	for (size_t i = 0; i < idiag_.size(); i++)
		if (idiag_[i] != B.idiag_[i])
			throw exception ("[SymDiaMatrix::operator+=] Unequal distribution of diagonals.");
	return true;
}

bool SymDiaMatrix::hdfload(const char* name)
{
#ifndef NO_HDF
	try
	{
		HDFFile hdf(name, HDFFile::readonly);
		if (not hdf.valid())
			return false;
		
		// read dimension
		if (not hdf.read("n", &n_, 1))
			return false;
		
		// read non-zero diagonal identificators
		if (idiag_.resize(hdf.size("idiag")))
			if (not hdf.read("idiag", &(idiag_[0]), idiag_.size()))
				return false;
		
		// comperssed array info
		NumberArray<int> zero_blocks;
		NumberArray<Complex> elements;
		
		if (zero_blocks.resize(hdf.size("zero_blocks")))
			if (not hdf.read("zero_blocks", &(zero_blocks[0]), zero_blocks.size()))
				return false;
		
		// load compressed elements
		if (elements.resize(hdf.size("x") / 2))
			if (not hdf.read("x", &(elements[0]), elements.size()))
				return false;
		
		// decompress
		elems_ = elements.decompress(zero_blocks);
		return true;
	}
	catch (H5::FileIException exc)
	{
		return false;
	}
#else /* NO_HDF */
	return false;
#endif
}

bool SymDiaMatrix::hdfsave(const char* name, bool docompress, int consec) const
{
#ifndef NO_HDF
	HDFFile hdf(name, HDFFile::overwrite);
		
	if (not hdf.valid())
	{
		return false;
// 		throw exception ("Unable to save HDF file \"%s\".", name);
	}
		
	// write dimension and diagonal info
	if (not hdf.write("n", &n_, 1))
		return false;
	if (not hdf.write("idiag", &(idiag_[0]), idiag_.size()))
		return false;
	
	// compress elements array
	NumberArray<int> zero_blocks;
	NumberArray<Complex> elements;
	
	if (docompress)
		std::tie(zero_blocks, elements) = elems_.compress(consec);
	else
		elements = elems_;
	
	// write compressed elements array
	if (not zero_blocks.empty())
	{
		if (not hdf.write (
			"zero_blocks",
			&(zero_blocks[0]),
			zero_blocks.size()
		)) return false;
	}
	if (not elements.empty())
	{
		if (not hdf.write (
			"x",
			&(elements[0]),
			elements.size()
		)) return false;
	}
	
	return true;
#else /* NO_HDF */
	return false;
#endif
}

CooMatrix SymDiaMatrix::tocoo() const
{
	std::vector<long> i, j;
	std::vector<Complex> v;
	
	auto el = elems_.begin();
	
	// for all diagonals
	for (auto id : idiag_)
	{
		// for all elements in this diagonal
		for (int iel = 0; iel < n_ - id; iel++)
		{
			// skip zero elements
			if (*el == 0.)
			{
				el++;
				continue;
			}
			
			// add this element to COO
			i.push_back(iel);
			j.push_back(iel+id);
			v.push_back(*el);
			
			// main diagonal shall be added only once
			if (id == 0)
			{
				el++;
				continue;
			}
			
			// and also its symmetric counterpart
			i.push_back(iel+id);
			j.push_back(iel);
			v.push_back(*el);
			
			// move on to the next element
			el++;
		}
	}
	
	return CooMatrix(n_, n_, i, j, v);
}

cArray SymDiaMatrix::dot(cArrayView const & B) const
{
	// check dimensions
	if ((int)B.size() != n_)
		throw exception ("[SymDiaMatrix::dot] Incompatible dimensions.");
	
	// the result
	cArray res(n_);
	
	// data pointers
	// - "restricted" and "aligned" for maximization of the cache usage
	// - "aligned" to convince the auto-vectorizer that vectorization is worth
	// NOTE: cArray (= NumberArray<Complex>) is aligned on sizeof(Complex) boundary
	// NOTE: GCC needs -ffast-math (included in -Ofast) to auto-vectorize both the ielem-loops below
	Complex       *       __restrict rp_res    = (Complex*)__builtin_assume_aligned(&res[0],    sizeof(Complex));
	Complex const * const __restrict rp_elems_ = (Complex*)__builtin_assume_aligned(&elems_[0], sizeof(Complex));
	Complex const * const __restrict rp_B      = (Complex*)__builtin_assume_aligned(&B[0],      sizeof(Complex));
	
	// for all elements in the main diagonal
	for (int ielem = 0; ielem < n_; ielem++)
		rp_res[ielem] = rp_elems_[ielem] * rp_B[ielem];
	
	// beginning of the current diagonal
	size_t beg = n_;
	
	// for all other diagonals
	for (unsigned id = 1; id < idiag_.size(); id++)
	{
		// index of this diagonal
		int idiag = idiag_[id];
		
		// number of elements in the current diagonal
		int Nelem = n_ - idiag;
		
		// for all elements of the current diagonal
		for (int ielem = 0; ielem < Nelem; ielem++)
		{
			rp_res[ielem]         += rp_elems_[beg + ielem] * rp_B[ielem + idiag];
			rp_res[ielem + idiag] += rp_elems_[beg + ielem] * rp_B[ielem];
		}
		
		// move to the beginning of the next diagonal
		beg += Nelem;
	}
	
	return res;
}

SymDiaMatrix SymDiaMatrix::kron (SymDiaMatrix const & B) const
{
	// FIXME
	return ::kron (this->tocoo(), B.tocoo()).todia();
}

std::ostream & operator << (std::ostream & out, SymDiaMatrix const & A)
{
	Array<Complex>::const_iterator iter = A.elems_.begin();
	
	for (auto id : A.idiag_)
	{
		out << "[" << id << "]: ";
		for (int i = 0; i < A.n_ - id; i++)
			out << *(iter++) << " ";
		
		std::cout << "\n";
	}
	
	return out;
}
