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
#include <vector>

#ifdef WITH_PNGPP
	#include <png++/png.hpp>
#endif

#include <H5Cpp.h>

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
#include <suitesparse/umfpack.h>

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
		H5::Exception::dontPrint();
		H5::H5File h5file(name, H5F_ACC_TRUNC);
		
		// all arrays are rank-1
		int rank = 1;
		
		// length of array
		hsize_t length;
		
		// save row count
		length = 1;
		H5::DataSpace dspc_m(rank, &length);
		H5::IntType dtype_m(H5::PredType::NATIVE_ULONG);
		H5::DataSet dset_m = h5file.createDataSet("m", dtype_m, dspc_m);
		dset_m.write(&(_m_), H5::PredType::NATIVE_ULONG);
		
		// save column count
		length = 1;
		H5::DataSpace dspc_n(rank, &length);
		H5::IntType dtype_n(H5::PredType::NATIVE_ULONG);
		H5::DataSet dset_n = h5file.createDataSet("n", dtype_n, dspc_n);
		dset_n.write(&(_n_), H5::PredType::NATIVE_ULONG);
		
		// save column pointers
		length = _p_.size();
		H5::DataSpace dspc_p(rank, &length);
		H5::IntType dtype_p(H5::PredType::NATIVE_LONG);
		H5::DataSet dset_p = h5file.createDataSet("p", dtype_p, dspc_p);
		dset_p.write(_p_.data(), H5::PredType::NATIVE_LONG);
		
		// save row indices
		length = _i_.size();
		H5::DataSpace dspc_i(rank, &length);
		H5::IntType dtype_i(H5::PredType::NATIVE_LONG);
		H5::DataSet dset_i = h5file.createDataSet("i", dtype_i, dspc_i);
		dset_i.write(_i_.data(), H5::PredType::NATIVE_LONG);
		
		// save interleaved data
		length = 2 * _x_.size();
		H5::DataSpace dspc_x(rank, &length);
		H5::FloatType dtype_x(H5::PredType::NATIVE_DOUBLE);
		H5::DataSet dset_x = h5file.createDataSet("x", dtype_x, dspc_x);
		dset_x.write(_x_.data(), H5::PredType::NATIVE_DOUBLE);
		
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
		H5::Exception::dontPrint();
		H5::H5File h5file(name, H5F_ACC_RDONLY);
		
		// read row count
		H5::DataSet dset_m = h5file.openDataSet("m");
		H5::DataSpace dspc_m = dset_m.getSpace();
		dset_m.read(&(_m_), H5::PredType::NATIVE_ULONG, dspc_m, dspc_m);
		
		// read column count
		H5::DataSet dset_n = h5file.openDataSet("n");
		H5::DataSpace dspc_n = dset_n.getSpace();
		dset_m.read(&(_n_), H5::PredType::NATIVE_ULONG, dspc_n, dspc_n);
		
		// read column pointers
		H5::DataSet dset_p = h5file.openDataSet("p");
		H5::DataSpace dspc_p = dset_p.getSpace();
		_p_.resize( dspc_p.getSimpleExtentNpoints() );
		dset_p.read(&(_p_[0]), H5::PredType::NATIVE_LONG, dspc_p, dspc_p);
		
		// read row indices
		H5::DataSet dset_i = h5file.openDataSet("i");
		H5::DataSpace dspc_i = dset_i.getSpace();
		_i_.resize( dspc_i.getSimpleExtentNpoints() );
		dset_i.read(&(_i_[0]), H5::PredType::NATIVE_LONG, dspc_i, dspc_i);
		
		// read interleaved data
		H5::DataSet dset_x = h5file.openDataSet("x");
		H5::DataSpace dspc_x = dset_x.getSpace();
		_x_.resize( dspc_x.getSimpleExtentNpoints() / 2);
		dset_x.read(&(_x_[0]), H5::PredType::NATIVE_DOUBLE, dspc_x, dspc_x);
		
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

CsrMatrix CsrMatrix::submatrix(unsigned a, unsigned b, unsigned c, unsigned d) const
{
	// if empty, return empty
	if (b <= a or d <= c)
		return CsrMatrix();
	
	CsrMatrix subm(b-a, d-c);
	subm._p_.push_back(0);
	
	// for all relevant rows
	for (unsigned irow = a; irow < b; irow++)
	{
		// get all columns
		unsigned idx1 = _p_[irow];
		unsigned idx2 = _p_[irow + 1];
		
		// get beginning of relevant row fragment
		std::vector<long>::const_iterator i_iter_begin = std::lower_bound(
			_i_.begin() + idx1,
			_i_.begin() + idx2,
			c
		);
		
		// get end of relevant row fragment
		std::vector<long>::const_iterator i_iter_end = std::lower_bound(
			_i_.begin() + idx1,
			_i_.begin() + idx2,
			d
		);
		
		// data length
		size_t n = i_iter_end - i_iter_begin;
		subm._p_.push_back(subm._p_.back() + n);
		
		// copy column indices (shifted by "c")
		subm._i_.resize(subm._i_.size() + n);
		std::transform(
			i_iter_begin,
			i_iter_end,
			subm._i_.end() - n,
			[ c ](long col) -> long { return col - c; }
		);
		
		// copy data
		subm._x_.insert(
			subm._x_.end(),
			_x_.data() + (i_iter_begin - _i_.begin()),
			_x_.data() + (i_iter_end   - _i_.begin())
		);
	}
	
	return subm;
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

// unsigned CsrMatrix::cg(
// 	const cArray&  b, cArray&  x,
// 	double eps,
// 	unsigned min_iterations, unsigned max_iterations,
// 	const PreconditionerInfo& pi
// ) const
// {
// 	// check dimensions
// 	size_t N = b.size();
// 	assert(rows() == N);
// 	assert(cols() == N);
// 	assert(x.size() == N);
// 
// 	//
// 	// 1) Set up desired preconditioner
// 	//
// 	
// 	// if Jacobi preconditioner is to be used, extract inverse diagonal
// 	cArray Dm1;	// A's inverse diagonal as 1D-array
// 	if (pi.preconditioner == P_JACOBI)
// 	{
// 		Dm1 = this->diag().transform(
// 			[](Complex z) -> Complex {
// 				return 1. / z;
// 			}
// 		);
// 	}
// 	
// 	// if SSOR preconditioner is to be used, scale diagonals by omega
// 	CsrMatrix A;		// this matrix with ω-scaled diagonal
// 	cArray factor_D;		// the scaled diagonal itself as a 1D-array
// 	Complex factor = (2. - pi.omega) / pi.omega;	// some other scaling factor
// 	if (pi.preconditioner == P_SSOR)
// 	{
// 		A = this->nzTransform(
// 			[ pi ](size_t i, size_t j, Complex z) -> Complex {
// 				return (i == j) ? z / pi.omega : z;
// 			}
// 		);
// 		factor_D = factor * A.diag();
// 	}
// 
// 	// if block inversion preconditioner is to be used:
// 	std::vector<CsrMatrix> blocks(pi.Nblock);
// 	std::vector<CsrMatrix::LUft> lufts(pi.Nblock);
// 	if (pi.preconditioner == P_BLOCK_INV)
// 	{
// 		// the dimension of a diagonal block
// 		unsigned blocksize = N / pi.Nblock;
// 		
// 		// for all diagonal blocks
// 		for (unsigned iblock = 0; iblock < pi.Nblock; iblock++)
// 		{
// 			// get a copy of iblock-th block
// 			blocks[iblock] = submatrix(
// 					iblock * blocksize, (iblock + 1) * blocksize,
// 					iblock * blocksize, (iblock + 1) * blocksize
// 			);
// 			
// 			// store the block's factorization
//  			lufts[iblock] = blocks[iblock].factorize();
// 		}
// 	}
// 	
// 	//
// 	// 2) Declare/initialize used variables
// 	//
// 	
// 	// some arrays (search directions etc.)
// 	cArray p(N), q(N), z(N);
// 	
// 	// residual; initialized to starting residual using the initial guess
// 	cArray r = b - this->dot(x);
// 
// 	// some other scalar variables
// 	Complex rho_new;		// contains inner product r_i^T · r_i
// 	Complex rho_old;		// contains inner product r_{i-1}^T · r_{i-1}
// 	Complex alpha, beta;	// contains projection ratios
// 	
// 	//
// 	// 3) Iterate
// 	//
// 	
// 	unsigned k;
// 	for/*ever*/ (k = 0; ; k++)
// 	{
// 		// apply desired preconditioner
// 		if (pi.preconditioner == P_NONE)
// 		{
// 			z = r;
// 		}
// 		else if (pi.preconditioner == P_JACOBI)
// 		{
// 			z = Dm1 * r;
// 		}
// 		else if (pi.preconditioner == P_SSOR)
// 		{
// 			z = A.upperSolve( factor_D * A.lowerSolve(r) );
// 		}
// 		else if (pi.preconditioner == P_BLOCK_INV)
// 		{
// 			// solve the preconditioner equation set Mz = r
// 			# pragma omp parallel for
// 			for (unsigned iblock = 0; iblock < pi.Nblock; iblock++)
// 			{
// 				size_t chunksize = cols() / pi.Nblock;
// 				
// 				// create a copy of a RHS segment
// 				cArray r_block(chunksize);
// 				memcpy(
// 					&r_block[0],				// dest
// 					&r[0] + iblock * chunksize,	// src
// 					chunksize * sizeof(Complex)	// n
// 				);
// 				
// 				// multiply by an inverted block
// 				cArray z_block = lufts[iblock].solve(r_block, 1);
// 				
// 				// copy output segment to the whole array
// 				memcpy(
// 					&z[0] + iblock * chunksize,	// dest
// 					&z_block[0],				// src
// 					chunksize * sizeof(Complex)	// n
// 				);
// 			}
// 		}
// 		
// 		// compute projection ρ = r·z
// 		rho_new = (r|z);
// 		
// 		// setup search direction p
// 		if (k == 0)
// 		{
// 			p = z;
// 		}
// 		else
// 		{
// 			beta = rho_new / rho_old;
// 			p = z + beta * p;
// 		}
// 		
// 		// move to next Krylov subspace by multiplying A·p
// 		q = this->dot(p);
// 		
// 		// compute projection ratio α
// 		alpha = rho_new / (p|q);
// 		
// 		// update the solution and the residual
// 		x += alpha * p;
// 		r -= alpha * q;
// 		
// 		// once in a while check convergence but do at least "min_iterations" iterations
// 		if (k >= min_iterations and k % 4 == 0 and r.norm() / b.norm() < eps)
// 			break;
// 		
// 		// check iteration limit (stop at "max_iterations" iterations)
// 		if (k >= max_iterations)
// 		{
// 			printf("[CsrMatrix::cg] Iteration limit %d reached.\n", max_iterations);
// 			break;
// 		}
// 		
// 		// move to the next iteration: store previous projection
// 		rho_old = rho_new;
// 	}
// 
// 	// cleanup
// 	
// 	if (pi.preconditioner == P_BLOCK_INV)
// 		for (unsigned iblock = 0; iblock < pi.Nblock; iblock++)
// 			lufts[iblock].free();
// 	
// 	return k;
// }

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
		H5::Exception::dontPrint();
		H5::H5File h5file(name, H5F_ACC_TRUNC);
		
		// all arrays are rank-1
		int rank = 1;
		
		// length of array
		hsize_t length;
		
		// save row count
		length = 1;
		H5::DataSpace dspc_m(rank, &length);
		H5::IntType dtype_m(H5::PredType::NATIVE_ULONG);
		H5::DataSet dset_m = h5file.createDataSet("m", dtype_m, dspc_m);
		dset_m.write(&_m_, H5::PredType::NATIVE_ULONG);
		
		// save column count
		length = 1;
		H5::DataSpace dspc_n(rank, &length);
		H5::IntType dtype_n(H5::PredType::NATIVE_ULONG);
		H5::DataSet dset_n = h5file.createDataSet("n", dtype_n, dspc_n);
		dset_n.write(&_n_, H5::PredType::NATIVE_ULONG);
		
		// save row pointers
		length = _p_.size();
		H5::DataSpace dspc_p(rank, &length);
		H5::IntType dtype_p(H5::PredType::NATIVE_LONG);
		H5::DataSet dset_p = h5file.createDataSet("p", dtype_p, dspc_p);
		dset_p.write(_p_.data(), H5::PredType::NATIVE_LONG);
		
		// save column indices
		length = _i_.size();
		H5::DataSpace dspc_i(rank, &length);
		H5::IntType dtype_i(H5::PredType::NATIVE_LONG);
		H5::DataSet dset_i = h5file.createDataSet("i", dtype_i, dspc_i);
		dset_i.write(_i_.data(), H5::PredType::NATIVE_LONG);
		
		// save interleaved data
		length = 2 * _x_.size();
		H5::DataSpace dspc_x(rank, &length);
		H5::FloatType dtype_x(H5::PredType::NATIVE_DOUBLE);
		H5::DataSet dset_x = h5file.createDataSet("x", dtype_x, dspc_x);
		dset_x.write(_x_.data(), H5::PredType::NATIVE_DOUBLE);
		
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
		H5::Exception::dontPrint();
		H5::H5File h5file(name, H5F_ACC_RDONLY);
		
		// read row count
		H5::DataSet dset_m = h5file.openDataSet("m");
		H5::DataSpace dspc_m = dset_m.getSpace();
		dset_m.read(&_m_, H5::PredType::NATIVE_ULONG, dspc_m, dspc_m);
		
		// read column count
		H5::DataSet dset_n = h5file.openDataSet("n");
		H5::DataSpace dspc_n = dset_n.getSpace();
		dset_m.read(&_n_, H5::PredType::NATIVE_ULONG, dspc_n, dspc_n);
		
		// read row pointers
		H5::DataSet dset_p = h5file.openDataSet("p");
		H5::DataSpace dspc_p = dset_p.getSpace();
		_p_.resize( dspc_p.getSimpleExtentNpoints() );
		dset_p.read(&_p_[0], H5::PredType::NATIVE_LONG, dspc_p, dspc_p);
		
		// read colunm indices
		H5::DataSet dset_i = h5file.openDataSet("i");
		H5::DataSpace dspc_i = dset_i.getSpace();
		_i_.resize( dspc_i.getSimpleExtentNpoints() );
		dset_i.read(&_i_[0], H5::PredType::NATIVE_LONG, dspc_i, dspc_i);
		
		// read interleaved data
		H5::DataSet dset_x = h5file.openDataSet("x");
		H5::DataSpace dspc_x = dset_x.getSpace();
		_x_.resize( dspc_x.getSimpleExtentNpoints() / 2);
		dset_x.read(&_x_[0], H5::PredType::NATIVE_DOUBLE, dspc_x, dspc_x);
		
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

cArray CsrMatrix::upperSolve(const cArray&  b) const
{
	// check size
	size_t N = b.size();
	assert((size_t)_m_ == N);
	assert((size_t)_n_ == N);

	// create output array
	cArray x(N);
	
// #ifdef __INTEL_MKL__
// 	
// 	//
// 	// test environment so that we can use 'long' == '_INTEGER_t'
// 	//
// #ifdef MKL_ILP64
// 	// MKL uses 64bit integers, i.e. "long long". We have "long". Is it interchangeable?
// 	if (sizeof(long) !=  sizeof(long long))
// 	{
// 		// no, notify the user
// 		printf("[Intel DSS] Error. Hex uses 'long' as internal representation of "
// 		       "matrix indices. Your ILP64 version of MKL uses 'long long', which is "
// 			   "larger than 'long' on your specific system. Please recompile Hex WITHOUT the "
// 			   "-DMKL_ILP64 flag.\n");
// 		abort();
// 	}
// #else // MKL_ILP64
// 	// MKL uses 32bit integers, i.e. "int". We have "long". Is it interchangeable?
// 	if (sizeof(long) !=  sizeof(int))
// 	{
// 		// no, notify the user
// 		printf("[Intel DSS] Error. Hex uses 'long' as internal representation of "
// 		       "matrix indices. Your LP64 version of MKL uses 'int', which is "
// 			   "smaller than 'long' on your specific system. Please recompile Hex WITH the "
// 			   "-DMKL_ILP64 flag.\n");
// 		abort();
// 	}
// #endif // MKL_ILP64
// 	
// 	long __m = this->rows();
// 	mkl_cspblas_zcsrtrsv(
// 		"U",	// use upper triangle
// 		"N",	// do not transpose
// 		"N",	// this matrix is not unit triangular
// 		&__m,	// row count
// 		const_cast<Complex*>(&_x_[0]), 	// matrix elements
// 		const_cast<long*>(&_p_[0]),		// row pointers
// 		const_cast<long*>(&_i_[0]),		// column indices
// 		const_cast<Complex*>(&b[0]),	// rhs
// 		&x[0]	// solution
// 	);
// 	
// #else // __INTEL_MKL__
	
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
	
// #endif // __INTEL_MKL__

	return x;
}

cArray CsrMatrix::lowerSolve(const cArray& b) const
{
	// check size
	size_t N = b.size();
	assert((size_t)_m_ == N);
	assert((size_t)_n_ == N);
	
	// create output array
	cArray x(N);
	
// #ifdef __INTEL_MKL__
// 
// 	//
// 	// test environment so that we can use 'long' == '_INTEGER_t'
// 	//
// #ifdef MKL_ILP64
// 	// MKL uses 64bit integers, i.e. "long long". We have "long". Is it interchangeable?
// 	if (sizeof(long) !=  sizeof(long long))
// 	{
// 		// no, notify the user
// 		printf("[Intel DSS] Error. Hex uses 'long' as internal representation of "
// 		       "matrix indices. Your ILP64 version of MKL uses 'long long', which is "
// 			   "larger than 'long' on your specific system. Please recompile Hex WITHOUT the "
// 			   "-DMKL_ILP64 flag.\n");
// 		abort();
// 	}
// #else // MKL_ILP64
// 	// MKL uses 32bit integers, i.e. "int". We have "long". Is it interchangeable?
// 	if (sizeof(long) !=  sizeof(int))
// 	{
// 		// no, notify the user
// 		printf("[Intel DSS] Error. Hex uses 'long' as internal representation of "
// 		       "matrix indices. Your LP64 version of MKL uses 'int', which is "
// 			   "smaller than 'long' on your specific system. Please recompile Hex WITH the "
// 			   "-DMKL_ILP64 flag.\n");
// 		abort();
// 	}
// #endif // MKL_ILP64
// 	
// 	long __m = this->rows();
// 	mkl_cspblas_zcsrtrsv(
// 		"L",	// use lower triangle
// 		"N",	// do not transpose
// 		"N",	// this matrix is not unit triangular
// 		&__m,	// row count
// 		const_cast<Complex*>(&_x_[0]), 	// matrix elements
// 		const_cast<long*>(&_p_[0]),		// row pointers
// 		const_cast<long*>(&_i_[0]),		// column indices
// 		const_cast<Complex*>(&b[0]),	// rhs
// 		&x[0]	// solution
// 	);
// 	
// #else // __INTEL_MKL__
	
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
	
// #endif // __INTEL_MKL__
	
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
	
		long status = umfpack_zl_triplet_to_col(
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


void CooMatrix::write(const char* filename) const
{
	FILE *f = fopen(filename, "w");
	fprintf(f, "# Matrix %ld × %ld with %ld nonzero elements:\n\n", _m_, _n_, _x_.size());
	for (size_t i = 0; i < _i_.size(); i++)
		fprintf(f, "%ld\t%ld\t%g\t%g\n", _i_[i], _j_[i], _x_[i].real(), _x_[i].imag());
	fclose(f);
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
		H5::Exception::dontPrint();
		H5::H5File h5file(name, H5F_ACC_TRUNC);
		
		// all arrays are rank-1
		int rank = 1;
		
		// length of array
		hsize_t length;
		
		// save row count
		length = 1;
		H5::DataSpace dspc_m(rank, &length);
		H5::IntType dtype_m(H5::PredType::NATIVE_ULONG);
		H5::DataSet dset_m = h5file.createDataSet("m", dtype_m, dspc_m);
		dset_m.write(&_m_, H5::PredType::NATIVE_ULONG);
		
		// save column count
		length = 1;
		H5::DataSpace dspc_n(rank, &length);
		H5::IntType dtype_n(H5::PredType::NATIVE_ULONG);
		H5::DataSet dset_n = h5file.createDataSet("n", dtype_n, dspc_n);
		dset_n.write(&_n_, H5::PredType::NATIVE_ULONG);
		
		// save row indices
		length = _i_.size();
		H5::DataSpace dspc_i(rank, &length);
		H5::IntType dtype_i(H5::PredType::NATIVE_LONG);
		H5::DataSet dset_i = h5file.createDataSet("i", dtype_i, dspc_i);
		dset_i.write(_i_.data(), H5::PredType::NATIVE_LONG);
		
		// save column indices
		length = _j_.size();
		H5::DataSpace dspc_j(rank, &length);
		H5::IntType dtype_j(H5::PredType::NATIVE_LONG);
		H5::DataSet dset_j = h5file.createDataSet("j", dtype_j, dspc_j);
		dset_j.write(_j_.data(), H5::PredType::NATIVE_LONG);
		
		// save interleaved data
		length = 2 * _x_.size();
		H5::DataSpace dspc_x(rank, &length);
		H5::FloatType dtype_x(H5::PredType::NATIVE_DOUBLE);
		H5::DataSet dset_x = h5file.createDataSet("x", dtype_x, dspc_x);
		dset_x.write(_x_.data(), H5::PredType::NATIVE_DOUBLE);
		
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
		H5::Exception::dontPrint();
		H5::H5File h5file(name, H5F_ACC_RDONLY);
		
		// read row count
		H5::DataSet dset_m = h5file.openDataSet("m");
		H5::DataSpace dspc_m = dset_m.getSpace();
		dset_m.read(&_m_, H5::PredType::NATIVE_ULONG, dspc_m, dspc_m);
		
		// read column count
		H5::DataSet dset_n = h5file.openDataSet("n");
		H5::DataSpace dspc_n = dset_n.getSpace();
		dset_m.read(&_n_, H5::PredType::NATIVE_ULONG, dspc_n, dspc_n);
		
		// read row indices
		H5::DataSet dset_i = h5file.openDataSet("i");
		H5::DataSpace dspc_i = dset_i.getSpace();
		_i_.resize( dspc_i.getSimpleExtentNpoints() );
		dset_i.read(&_i_[0], H5::PredType::NATIVE_LONG, dspc_i, dspc_i);
		
		// read colunm indices
		H5::DataSet dset_j = h5file.openDataSet("j");
		H5::DataSpace dspc_j = dset_j.getSpace();
		_j_.resize( dspc_j.getSimpleExtentNpoints() );
		dset_j.read(&_j_[0], H5::PredType::NATIVE_LONG, dspc_j, dspc_j);
		
		// read interleaved data
		H5::DataSet dset_x = h5file.openDataSet("x");
		H5::DataSpace dspc_x = dset_x.getSpace();
		_x_.resize( dspc_x.getSimpleExtentNpoints() / 2);
		dset_x.read(&_x_[0], H5::PredType::NATIVE_DOUBLE, dspc_x, dspc_x);
		
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
