/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2014                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <algorithm>
#include <complex>
#include <cstdio>
#include <cstring>
#include <set>
#include <vector>

#ifndef NO_PNG
#include <png++/png.hpp>
#endif

#ifndef NO_UMFPACK
#include <umfpack.h>
#endif

#ifndef NO_HDF
#include "hdffile.h"
#endif

#include "arrays.h"
#include "matrix.h"
#include "misc.h"

//
// Dense matrix routines
//

#ifndef NO_BLAS
extern "C" void zgemm_
(
    char * TRANSA, 
    char * TRANSB, 
    int * M, 
    int * N, 
    int * K, 
    Complex * ALPHA, 
    Complex * A, 
    int * LDA, 
    Complex * B, 
    int * LDB, 
    Complex * BETA, 
    Complex * C, 
    int * LDC 
);
#endif

#ifndef NO_LAPACK
// Calculates LU decomposition of a complex dense matrix.
extern "C" void zgetrf_
(
    int * M,
    int * N,
    Complex * A,
    int * LDA,
    int * IPIV,
    int * INFO
);
// Inverts a complex dense matrix using precomputed LU decomposition.
extern "C" void zgetri_
(
    int * N,
    Complex * A,
    int * LDA,
    int * IPIV,
    Complex * WORK,
    int * LWORK,
    int * INFO
);
// Eigen-diagonalization of a general complex matrix.
extern "C" void zgeev_
(
    char * JOBVL,
    char * JOBVR,
    int * N,
    Complex * A,
    int * LDA,
    Complex * W,
    Complex * VL,
    int * LDVL,
    Complex * VR,
    int * LDVR,
    Complex * WORK,
    int * LWORK,
    double * RWORK,
    int * INFO
);
#endif

// specialization of multiplication of real dense row-major ("C-like") matrices
template<> RowMatrix<double> operator * (RowMatrix<double> const & A, ColMatrix<double> const & B)
{
    // check sanity
    if (A.cols() != B.rows());
        throw exception ("Matrix multiplication requires A.cols() == B.rows(), but %d != %d.", A.cols(), B.rows());
    
    // create output matrix
    RowMatrix<double> C (A.rows(), B.cols());
    
    // get sizes
    int m = A.rows(), n = B.cols(), r = A.cols();
    
    // get restricted and aligned pointers for fast access
#if defined(__GNUC__) || defined(__INTEL_COMPILER)
    double const * const restrict pA = (double*) __builtin_assume_aligned (A.data().data(), 32u);
    double const * const restrict pB = (double*) __builtin_assume_aligned (B.data().data(), 32u);
    double       * const restrict pC = (double*) __builtin_assume_aligned (C.data().data(), 32u);
#else
    double const * const restrict pA = A.data().data();
    double const * const restrict pB = B.data().data();
    double       * const restrict pC = C.data().data();
#endif /* defined(__GNUC__) || defined(__INTEL_COMPILER) */
    
    // multiply
    for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
    for (int k = 0; k < r; k++)
        pC[i * n + j] += pA[i * r + k] * pB[j * r + k];
    
    // return result
    return C;
}

template<> ColMatrix<Complex> ColMatrix<Complex>::invert () const
{
#ifndef NO_LAPACK
    if (rows() != cols())
        throw exception ("Only square matrices can be diagonalized.");
    int N = rows();
    
    ColMatrix<Complex> inv (*this);
    iArray IPIV(N);
    int LWORK = -1;
    cArray WORK(1);
    int INFO;
    
    // compute the LU factorization
    zgetrf_
    (
        &N,             // number of rows in the matrix
        &N,             // number of columns of the matrix
        inv.begin(),    // matrix data (on return the LU factors)
        &N,             // leading dimension of A
        &IPIV[0],       // pivot indices
        &INFO           // diagnostic information
    );
    
    // query for optimal work size for inversion
    zgetri_
    (
        &N,             // order of the matrix
        inv.begin(),    // LU-factors (on return the inverse matrix)
        &N,             // leading dimension of A
        &IPIV[0],       // pivot indices
        &WORK[0],       // work array
        &LWORK,         // work array length
        &INFO           // diagnostic information
    );
    
    // resize work array
    LWORK = WORK.front().real();
    WORK.resize(LWORK);
    
    // compute the inverse
    zgetri_
    (
        &N,             // order of the matrix
        inv.begin(),    // LU-factors (on return the inverse matrix)
        &N,             // leading dimension of A
        &IPIV[0],       // pivot indices
        &WORK[0],       // work array
        &LWORK,         // work array length
        &INFO           // diagnostic information
    );
    
    return inv;
#else
    throw exception ("Cannot invert dense matrix without LAPACK support.");
#endif
}

// specialization of diagonalization of complex dense column-major ("Fortran-like") matrices
template<> std::tuple<cArray,ColMatrix<Complex>,ColMatrix<Complex>> ColMatrix<Complex>::diagonalize () const
{
#ifndef NO_LAPACK
    if (rows() != cols())
        throw exception ("Only square matrices can be diagonalized.");
    int N = rows();
    
    ColMatrix<Complex> VL(N,N), VR(N,N);
    cArray W(N), WORK(1);
    rArray RWORK(2*N);
    int LWORK = -1, INFO;
    char JOB = 'V';
    
    // copy matrix data (it will be overwritten)
    cArray A = data();
    
    // get work size
    zgeev_
    (
        &JOB,       // compute left eigenvectors
        &JOB,       // compute right eigenvectors
        &N,         // order of the matrix
        &A[0],      // matrix A elements
        &N,         // leading dimension of the matrix A
        &W[0],      // eigenvalues
        VL.begin(), // left eigenvectors elements
        &N,         // leading dimension of the matrix VL
        VR.begin(), // right eigenvectors elements
        &N,         // leading dimension of the matrix VR
        &WORK[0],   // work array
        &LWORK,     // length of the work array
        &RWORK[0],  // work array
        &INFO       // diagnostic information
    );
    
    // resize work array
    LWORK = WORK[0].real();
    WORK.resize(LWORK);
    
    // run the diagonalization
    zgeev_
    (
        &JOB,       // compute left eigenvectors
        &JOB,       // compute right eigenvectors
        &N,         // order of the matrix
        &A[0],      // matrix A elements
        &N,         // leading dimension of the matrix A
        &W[0],      // eigenvalues
        VL.begin(), // left eigenvectors elements
        &N,         // leading dimension of the matrix VL
        VR.begin(), // right eigenvectors elements
        &N,         // leading dimension of the matrix VR
        &WORK[0],   // work array
        &LWORK,     // length of the work array
        &RWORK[0],  // work array
        &INFO       // diagnostic information
    );
    
    if (INFO != 0)
    {
        std::cerr << "Diagonalization error: " << INFO << std::endl;
    }
    
    return std::make_tuple(W,VL,VR);
#else
    throw exception ("Cannot diagonalize dense matrix without LAPACK support.");
#endif
}

cArray kron_dot (RowMatrix<Complex> const & A, RowMatrix<Complex> const & B, cArrayView const v)
{
    assert(A.cols() * B.cols() == v.size());
    
    // return vector
    NumberArray<Complex> w(A.rows() * B.rows());
    
    // auxiliary matrix
    ColMatrix<Complex> C(B.rows(),A.cols());
    Complex alpha = 1, beta = 0;
    char norm = 'N', trans = 'T';
    
    // C = V * A^T
    {
        int m = B.rows(), k = B.cols(), n = A.rows();
        zgemm_(&norm, &norm, &m, &n, &k, &alpha, const_cast<Complex*>(v.data()), &n, const_cast<Complex*>(A.data().data()), &k, &beta, C.data().begin(), &m);
    }
    
    // W = B * C
    {
        int m = A.rows(), k = A.cols(), n = C.cols();
        zgemm_(&trans, &norm, &m, &n, &k, &alpha, const_cast<Complex*>(B.data().begin()), &n, C.data().begin(), &n, &beta, w.data(), &m);
    }
    
    return w;
}

//
// Sparse matrix routines.
//

CooMatrix kron (const CooMatrix& A, const CooMatrix& B)
{
    // shorthands
    size_t Asize = A.v().size();
    size_t Bsize = B.v().size();
    size_t Csize = Asize * Bsize;
    size_t Brows = B.rows();
    size_t Bcols = B.cols();
    
    // set correct dimensions, pre-allocate space
    int m = A.rows() * B.rows();
    int n = A.cols() * B.cols();
    lArray C_i (Csize), C_j (Csize);
    cArray C_x (Csize);
    
    // get pointers
    long const * restrict pA_i = A.i().data();
    long const * restrict pB_i = B.i().data();
    long       * restrict pC_i = C_i.data();
    long const * restrict pA_j = A.j().data();
    long const * restrict pB_j = B.j().data();
    long       * restrict pC_j = C_j.data();
    Complex const * restrict pA_x = A.v().data();
    Complex const * restrict pB_x = B.v().data();
    Complex       * restrict pC_x = C_x.data();
    
    // loop over A data
    for (size_t ia = 0; ia < Asize; ia++)
    {
        // loop over B data
        for (size_t ib = 0; ib < Bsize; ib++)
        {
            // compute new row index
            *pC_i = pA_i[ia] * Brows + pB_i[ib];
            
            // compute new column index
            *pC_j = pA_j[ia] * Bcols + pB_j[ib];
            
            // compute product of the two elements
            *pC_x = pA_x[ia] * pB_x[ib];
            
            // move to next value of C
            pC_i++; pC_j++; pC_x++;
        }
    }
    
    // return the Kronecker product
    return CooMatrix(m, n, C_i, C_j, C_x);
}

CooMatrix eye (size_t N)
{
    return CooMatrix(N,N).symm_populate_band
    (
        0,
        [](size_t i, size_t j) -> Complex { return 1.; }
    );
}

CooMatrix stairs (size_t N)
{
    return CooMatrix(N,N).symm_populate_band
    (
        0,
        [](size_t i, size_t j) -> Complex { return (i == j) ? i : 0.; }
    );
}

// 
// CSC matrix routines
// 

CscMatrix & CscMatrix::operator *= (double r)
{
    size_t N = i_.size();
    
    for (size_t i = 0; i < N; i++)
        x_[i] *= r;
    
    return *this;
}

CscMatrix & CscMatrix::operator &= (const CscMatrix&  B)
{
    size_t N = i_.size();
    
    assert(m_ == B.m_);
    assert(n_ == B.n_);
    assert(N == B.i_.size());
    
    for (size_t i = 0; i < N; i++)
    {
        assert(i_[i] == B.i_[i]);
        
        x_[i] += B.x_[i];
    }
    
    return *this;
}

CscMatrix & CscMatrix::operator ^= (const CscMatrix&  B)
{
    size_t N = i_.size();
    
    assert(m_ == B.m_);
    assert(n_ == B.n_);
    assert(N == B.i_.size());
    
    for (size_t i = 0; i < N; i++)
    {
        assert(i_[i] == B.i_[i]);
        
        x_[i] -= B.x_[i];
    }
    
    return *this;
}

cArray CscMatrix::dotT (const cArrayView b) const
{
    // create output array
    cArray c (n_);
    
    // the matrix "*this" is actually transposed
    for (unsigned icol = 0; icol < n_; icol++)
    {
        size_t idx1 = p_[icol];
        size_t idx2 = p_[icol+1];
        
        // for all nonzero elements in this column
        for (size_t idx = idx1; idx < idx2; idx++)
        {
            // get row number
            unsigned irow = i_[idx];
            
            // store product
            c[icol] += x_[idx] * b[irow];
        }
    }
    
    return c;
}

#ifndef NO_UMFPACK
CooMatrix CscMatrix::tocoo () const
{
    // reserve space for the auxiliary (J) and the output (Ti,Tj,Tx,Tz) arrays
    size_t N = x_.size();
    std::vector<long> Ti(N), Tj(N), J(N);
    std::vector<Complex> Tx(N);
    
    // do we have any elements at all?
    if (N != 0)
    {
    
        // do the conversion
        long status = umfpack_zl_col_to_triplet(n_, p_.data(), J.data());
        
        // check success
        if (status != 0)
        {
            std::cerr << "\n[CscMatrix::tocoo] Exit code " << status << "\n";
            umfpack_zl_report_status(0, status);
        }
        
        // copy only non-zero entries to output arrays
        size_t nz = 0;
        for (size_t i = 0; i < N; i++)
        {
            if (x_[i] != 0.)
            {
                Ti[nz] = i_[i];
                Tj[nz] = J[i];
                Tx[nz] = x_[i];
                nz++;
            }
        }
        
        // crop the output arrays
        Ti.resize(nz);
        Tj.resize(nz);
        Tx.resize(nz);
        
    }
    
    // return new CooMatrix
    return CooMatrix (m_, n_, Ti, Tj, Tx);
}
#endif

bool CscMatrix::hdfsave (const char* name) const
{
#ifndef NO_HDF
    try
    {
        HDFFile file(name, HDFFile::overwrite);
        
        // write dimensions
        file.write("m", &m_, 1);
        file.write("n", &n_, 1);
        
        // write indices
        if (not p_.empty())
            file.write("p", &(p_[0]), p_.size());
        if (not i_.empty())
            file.write("i", &(i_[0]), i_.size());
        
        // write complex data as a "double" array
        if (not x_.empty())
        {
            file.write
            (
                "x",
                reinterpret_cast<double const*>( &(x_[0]) ),
                x_.size() * 2
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

bool CscMatrix::hdfload (const char* name)
{
#ifndef NO_HDF
    try
    {
        HDFFile hdf(name, HDFFile::readonly);
        
        // read dimensions
        hdf.read("m", &m_, 1);
        hdf.read("n", &n_, 1);
        
        // read indices
        if (p_.resize(hdf.size("p")))
            hdf.read("p", &(p_[0]), p_.size());
        if (i_.resize(hdf.size("i")))
            hdf.read("i", &(i_[0]), i_.size());
        
        // read data
        if (x_.resize(hdf.size("x") / 2))
        {
            hdf.read
            (
                "x",
                reinterpret_cast<double*>(&(x_[0])),
                x_.size() * 2
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

//
// CSR matrix routines
//

CsrMatrix & CsrMatrix::operator *= (Complex r)
{
    size_t N = i_.size();
    
    for (size_t i = 0; i < N; i++)
        x_[i] *= r;
    
    return *this;
}

CsrMatrix & CsrMatrix::operator &= (CsrMatrix const &  B)
{
    size_t N = i_.size();
    
    // check at least dimensions and non-zero element count
    assert(m_ == B.m_);
    assert(n_ == B.n_);
    assert(N == B.i_.size());
    
    for (size_t i = 0; i < N; i++)
        x_[i] += B.x_[i];
    
    return *this;
}

CsrMatrix & CsrMatrix::operator ^= (CsrMatrix const &  B)
{
    size_t N = i_.size();
    
    assert(m_ == B.m_);
    assert(n_ == B.n_);
    assert(N == B.i_.size());
    
    for (size_t i = 0; i < N; i++)
        x_[i] -= B.x_[i];
    
    return *this;
}

cArray CsrMatrix::dot (cArrayView const & b) const
{
    // create output array
    cArray c(m_);
    
    for (unsigned irow = 0; irow < m_; irow++)
    {
        size_t idx1 = p_[irow];
        size_t idx2 = p_[irow+1];
        
        // for all nonzero elements in this row
        for (size_t idx = idx1; idx < idx2; idx++)
        {
            // get column number
            unsigned icol = i_[idx];
            
            // store product
            c[irow] += x_[idx] * b[icol];
        }
    }
    
    return c;
}

#ifndef NO_PNG
CsrMatrix::PngGenerator::PngGenerator (CsrMatrix const * mat, double threshold)
    : base_t(mat->cols(), mat->rows()), M_(mat), buff_(mat->cols()), threshold_(threshold)
{
}

CsrMatrix::PngGenerator::~PngGenerator ()
{
}

png::byte* CsrMatrix::PngGenerator::get_next_row (size_t irow)
{
    // get column indices
    int idx_min = M_->p_[irow];
    int idx_max = M_->p_[irow + 1];
    
    // clear memory
    for (int icol = 0; icol < (int)M_->cols(); icol++)
        buff_[icol] = 1;
    
    // for all nonzero columns
    for (int idx = idx_min; idx < idx_max; idx++)
        if (std::abs(M_->x_[idx]) > threshold_)
            buff_[M_->i_[idx]] = 0;
        
    // pass the buffer
    return reinterpret_cast<png::byte*>(row_traits::get_data(buff_));
}


void CsrMatrix::plot (const char* filename, double threshold) const
{
    // create output file
    std::ofstream out(filename, std::ios_base::out | std::ios_base::binary);
    
    // create PNG data generator
    PngGenerator png(this, threshold);
    
    // write PNG file
    png.write(out);
}
#endif

#ifndef NO_UMFPACK
CsrMatrix::LUft CsrMatrix::factorize (double droptol) const
{
    // Use standard UMFPACK sequence
    void *Symbolic, *Numeric;
    long status;
    
    // get default setting
    double Control[UMFPACK_CONTROL];
    umfpack_zl_defaults(Control);
    
    // modify the drop tolerance
    Control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_SYMMETRIC;
    Control[UMFPACK_DROPTOL] = droptol;
    
    // analyze the sparse structure
    status = umfpack_zl_symbolic
    (
        m_, n_,                    // matrix dimensions
        p_.data(), i_.data(),        // column and row indices
        reinterpret_cast<const double*>(x_.data()), 0,    // matrix data
        &Symbolic, Control, nullptr                // UMFPACK internals
    );
    if (status != 0)
    {
        std::cerr << "\n[CscMatrix::factorize] Exit status " << status << "\n";
        umfpack_zl_report_status(0, status);
        abort();
    }
    
    // do some factorizations
    status = umfpack_zl_numeric
    (
        p_.data(), i_.data(),    // column and row indices
        reinterpret_cast<const double*>(x_.data()), 0,    // matrix data
        Symbolic, &Numeric, Control, nullptr    // UMFPACK internals
    );
    if (status != 0)
    {
        std::cerr << "\n[CscMatrix::factorize] Exit status " << status << "\n";
        umfpack_zl_report_status(0, status);
        abort();
    }
    
    // release unused data
    umfpack_zl_free_symbolic(&Symbolic);
    return LUft (this, Numeric);
}
#endif

#ifndef NO_UMFPACK
cArray CsrMatrix::solve (const cArrayView b, size_t eqs) const
{
    // only square matrices are allowed
    assert(m_ == n_);
    
    // compute the LU factorization
    LUft luft ( std::move(factorize()) );
    
    // solve the equations
    cArray solution = luft.solve(b, eqs);
    
    // return
    return solution;
}
#endif

void CsrMatrix::write (const char* filename) const
{
    std::ofstream out (filename);
    out << "# Matrix " << m_ << " × " << n_ << " with " << x_.size() << " nonzero elements:\n\n";
    for (unsigned irow = 0; irow < m_; irow++)
    {
        size_t idx1 = p_[irow];
        size_t idx2 = p_[irow + 1];
        
        for (size_t idx = idx1; idx < idx2; idx++)
            out << irow << "\t" << i_[idx] << "\t" << x_[idx].real() << "\t" << x_[idx].imag() << "\n";
    }
    out.close();
}

void CsrMatrix::link (std::string name)
{
    name_ = name;
}

bool CsrMatrix::hdfsave (std::string name) const
{
#ifndef NO_HDF
    try
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
                reinterpret_cast<double const*>(&(x_[0])),
                x_.size() * 2
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

bool CsrMatrix::hdfload (std::string name)
{
#ifndef NO_HDF
    try
    {
        HDFFile hdf(name, HDFFile::readonly);
        
        // read dimensions
        hdf.read("m", &m_, 1);
        hdf.read("n", &n_, 1);
        
        // read indices
        if (p_.resize(hdf.size("p")))
            hdf.read("p", &(p_[0]), p_.size());
        if (i_.resize(hdf.size("i")))
            hdf.read("i", &(i_[0]), i_.size());
        
        // read data
        if (x_.resize(hdf.size("x") / 2))
        {
            hdf.read
            (
                "x",
                reinterpret_cast<double*>(&(x_[0])),
                x_.size() * 2
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

double CsrMatrix::norm () const
{
    size_t N = i_.size();
    double res = 0.;
    
    // return the abs(largest element)
    for (size_t i = 0; i < N; i++)
    {
        // compute the absolute value
        double val = abs(x_[i]);
        
        // update the winner
        if (val > res)
            res = val;
    }
    
    return res;
}

cArray CsrMatrix::upperSolve (cArrayView const &  b) const
{
    // check size
    size_t N = b.size();
    assert((size_t)m_ == N);
    assert((size_t)n_ == N);

    // create output array
    cArray x(N);
    
    // loop over rows
    for (size_t i = 0; i < N; i++)
    {
        size_t row = N - 1 - i;
        Complex accum = 0.;
        
        // get relevant columns of the sparse matrix
        size_t idx1 = p_[row];
        size_t idx2 = p_[row + 1];
        
        // diagonal element of the matrix
        Complex a = 0.;
        
        // loop over the columns
        for (size_t idx = idx1; idx < idx2; idx++)
        {
            // which column is this?
            size_t col = i_[idx];
            
            // diagonal element will be useful in a moment, store it
            if (col == row)
                a = x_[idx];
            
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

cArray CsrMatrix::lowerSolve (cArrayView const & b) const
{
    // check size
    size_t N = b.size();
    assert((size_t)m_ == N);
    assert((size_t)n_ == N);
    
    // create output array
    cArray x(N);
    
    // loop over rows
    for (size_t row = 0; row < N; row++)
    {
        Complex accum = 0.;
        
        // get relevant columns of the sparse matrix
        size_t idx1 = p_[row];
        size_t idx2 = p_[row + 1];
        
        // diagonal element of the matrix
        Complex a = 0.;
        
        // loop over the columns
        for (size_t idx = idx1; idx < idx2; idx++)
        {
            // which column is this?
            size_t col = i_[idx];
            
            // diagonal element will be useful in a moment, store it
            if (col == row)
                a = x_[idx];
            
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

cArray CsrMatrix::diag () const
{
    cArray D ( std::min(m_, n_) );
    
    for (size_t irow = 0; irow < (size_t)m_; irow++)
        for (size_t idx = p_[irow]; idx < (size_t)p_[irow+1]; idx++)
            if ((size_t)i_[idx] == irow)
                D[irow] = x_[idx];
    
    return D;
}

#ifndef NO_UMFPACK
CooMatrix CsrMatrix::tocoo () const
{
    // reserve space for the auxiliary (__j) and the output (Ti,Tj,Tx,Tz) arrays
    size_t N = x_.size();
    lArray Ti(N), Tj(N), __j(N);
    cArray Tx(N);
    
    // do we have any elements at all?
    if (N != 0)
    {
    
        // do the conversion
        long status = umfpack_zl_col_to_triplet(n_, p_.data(), __j.data());
        
        // check success
        if (status != 0)
        {
            std::cerr << "\n[CscMatrix::tocoo] Exit code " << status << "\n";
            umfpack_zl_report_status(0, status);
        }
        
        // copy only non-zero entries to output arrays
        size_t nz = 0;
        for (size_t i = 0; i < N; i++)
        {
            if (x_[i] != 0.)
            {
                Ti[nz] = i_[i];
                Tj[nz] = __j[i];
                Tx[nz] = x_[i];
                nz++;
            }
        }
        
        // crop the output arrays
        Ti.resize(nz);
        Tj.resize(nz);
        Tx.resize(nz);
    }
    
    // return new CooMatrix
    return CooMatrix (m_, n_, Tj, Ti, Tx);
}
#endif

RowMatrix<Complex> CsrMatrix::torow () const
{
    RowMatrix<Complex> M (rows(),cols());
    for (unsigned irow = 0; irow < rows(); irow++)
    {
        unsigned rptr_begin = p_[irow];
        unsigned rptr_end = p_[irow + 1];
        
        for (unsigned idx = rptr_begin; idx < rptr_end; idx++)
        {
            unsigned icol = i_[idx];
            Complex x = x_[idx];
            
            M(irow, icol) = x;
        }
    }
    
    return M;
}

CsrMatrix CsrMatrix::sparse_like (const CsrMatrix& B) const
{
    // check dimensions
    assert(m_ == B.m_);
    assert(n_ == B.n_);
    
    // prepare zero matrix with the same storage pattern the matrix B has
    CsrMatrix A = B;
    memset(A.x_.data(), 0, A.x_.size() * sizeof(Complex));
    
    // copy all nonzero elements of "this" matrix
    for (unsigned row = 0; row < m_; row++)
    {
        size_t idx1 = A.p_[row];
        size_t idx2 = A.p_[row+1];
        
        for (size_t idx = idx1; idx < idx2; idx++)
        {
            unsigned col = A.i_[idx];
            A.x_[idx] = (*this)(row,col);
        }
    }
    
    // return temporary matrix
    return A;
}

Complex CsrMatrix::operator () (unsigned i, unsigned j) const
{
    // get all column indices, which have nonzero element in row "i"
    size_t idx1 = p_[i];
    size_t idx2 = p_[i + 1];
    
    // find the correct column ("j")
    auto it = std::lower_bound(i_.begin() + idx1, i_.begin() + idx2, j);

    if (it == i_.end() or *it != j)
        return 0.;
    
    // return the value
    return x_[it - i_.begin()];
}

Complex sparse_row_scalar_product
(
  int n1, Complex const * const restrict x1, long const * const restrict i1,
  int n2, Complex const * const restrict x2, long const * const restrict i2
)
{
    Complex ssp = 0.;
    if (n1 < 1 or n2 < 1)
        return ssp;
    
    for (int pos1 = 0, pos2 = 0; pos1 < n1 and pos2 < n2;)
    {
        if (i1[pos1] < i2[pos2])
        {
            pos1++;
        }
        else if (i1[pos1] > i2[pos2])
        {
            pos2++;
        }
        else
        {
            ssp += x1[pos1] * x2[pos2];
            pos1++;
            pos2++;
        }
    }
    
    return ssp;
}

CsrMatrix operator * (CsrMatrix const & A, CscMatrix const & B)
{
    // check compatibility
    if (A.cols() != B.rows())
        throw exception ("Wrong sizes for matrix multiplication: %d %d, %d %d.", A.rows(), A.cols(), B.rows(), B.cols());
    
    // result arrays
    lArray Cp = { 0 }, Ci;
    cArray Cx;
    
    // for all rows of the resulting matrix
    for (int irow = 0; irow < (int)A.rows(); irow++)
    {
        // for all columns of the resulting matrix 
        for (int icol = 0; icol < (int)B.cols(); icol++)
        {
            // compute scalar product of the rows
            Complex ssp = sparse_row_scalar_product
            (
                A.p()[irow+1] - A.p()[irow],   // nnz count
                A.x().data() + A.p()[irow],    // elements
                A.i().data() + A.p()[irow],    // column indices
                B.p()[icol+1] - B.p()[icol],   // nnz count
                B.x().data() + B.p()[icol],    // elements
                B.i().data() + B.p()[icol]     // row indices
            );
            
            // continue if zero
            if (ssp == 0.)
                continue;
            
            // otherwise store the value (and column index)
            Cx.push_back(ssp);
            Ci.push_back(icol);
        }
        
        // move to the next row
        Cp.push_back(Cx.size());
    }
    
    return CsrMatrix(A.rows(), B.cols(), Cp, Ci, Cx);
}


//
// COO matrix routines
//

#ifndef NO_UMFPACK
CscMatrix CooMatrix::tocsc () const
{
    size_t nz = x_.size();
    
    // CSC matrix data
    lArray Ap(n_ + 1), Ai(nz);
    cArray Ax(nz);
    
    // do we have any elements at all?
    if (nz != 0)
    {
    
        long status = umfpack_zl_triplet_to_col
        (
            m_,            // rows
            n_,            // cols
            nz,                // data length
            i_.data(),        // row indices
            j_.data(),     // column indices
            reinterpret_cast<const double *>(x_.data()), 0,     // interleaved data
            Ap.data(),        // column pointers
            Ai.data(),        // row pointers
            reinterpret_cast<double *>(Ax.data()), 0,    // interleaved data
            0
        );
        
        // check success
        if (status != 0)
        {
            std::cerr << "\n[CooMatrix::tocsc] Exit code " << status << "\n";
            umfpack_zl_report_status(0, status);
        }
        
        // crop storage
        size_t N = Ap[n_];
        Ai.resize(N);
        Ax.resize(N);
    }
    
    return CscMatrix(m_, n_, Ap, Ai, Ax);
}
#endif

#ifndef NO_UMFPACK
CsrMatrix CooMatrix::tocsr () const
{
    size_t nz = x_.size();
    
    // CSC matrix data
    lArray Ap(n_ + 1), Ai(nz);
    cArray Ax(nz);
    
    // do we have any elements at all?
    if (nz != 0)
    {
    
        long status = umfpack_zl_triplet_to_col
        (
            n_,            // cols (rows of transposed matrix)
            m_,            // rows (cols of transposed matrix)
            nz,             // data length
            j_.data(),     // column indices (rows of transposed matrix)
            i_.data(),     // row indices (cols of transposed matrix)
            
            // interleaved data
            reinterpret_cast<const double *>(x_.data()),
            nullptr,
            
            Ap.data(),      // row pointers
            Ai.data(),      // column indices
            
            // interleaved data
            reinterpret_cast<double *>(Ax.data()),
            nullptr,
         
            // map
            nullptr
        );
        
        // check success
        if (status != 0)
        {
            std::cerr << "\n[CooMatrix::tocsr] Exit code " << status << "\n";
            umfpack_zl_report_status(0, status);
        }
        
        // crop storage
        size_t N = Ap[m_];
        Ai.resize(N);
        Ax.resize(N);
    }
    
    return CsrMatrix(m_, n_, Ap, Ai, Ax);
}
#endif

SymDiaMatrix CooMatrix::todia (MatrixTriangle triangle) const
{
    // diagonal indices
    std::set<int> diags;
    
    // elements per diagonal
    cArrays elems(m_);
    
    // for all nonzero elements
    for (size_t i = 0; i < x_.size(); i++)
    {
        // get row and column index
        int irow = i_[i];
        int icol = j_[i];
        
        // get diagonal
        int d = icol - irow;
        
        // do we want the lower triangle?
        if ((d < 0) and not (triangle & strict_lower))
            continue;
        
        // do we want the upper triangle?
        if ((d > 0) and not (triangle & strict_upper))
            continue;
        
        // do we want the diagonal?
        if ((d == 0) and not (triangle & diagonal))
            continue;
        
        // separate sign and value of the diagonal index
        MatrixTriangle dtri = (d > 0 ? upper : lower);
        d = std::abs(d);
        
        // add this diagonal (if not already added)
        diags.insert(d);
        
        // reserve space for this diagonal if not already done
        if (elems[d].size() == 0)
            elems[d].resize(m_ - d);
        
        // add element numbered by the correct index (row/column)
        if (dtri == upper)
            elems[d][irow] = x_[i];
        else // dtri == lower
            elems[d][icol] = x_[i];
    }
    
    // concatenate the diagonals, construct matrix object
    return SymDiaMatrix
    (
        m_,
        Array<int>(diags.cbegin(), diags.cend()),
        join(elems)
    );
}

void CooMatrix::write (const char* filename) const
{
    std::ofstream f(filename);
    
    f << "# Matrix " << m_ << " × " << n_ << " with " << x_.size() << " nonzero elements:\n\n";
    
    for (size_t i = 0; i < i_.size(); i++)
        f << i_[i] << "\t" << j_[i] << "\t" << x_[i].real() << "\t" << x_[i].imag() << "\n";
}


CooMatrix CooMatrix::reshape (size_t m, size_t n) const
{
    CooMatrix C = *this;
    
    // conserved dimensions
    size_t N = C.i_.size();
    size_t H = C.m_;
    
    // check dimensions
    assert(m * n == C.m_ * C.n_);
    
    // reshape
    for (size_t i = 0; i < N; i++)
    {
        // conserved position in column-ordered array
        size_t idx = C.i_[i] + C.j_[i] * H;
        
        // new coordinates
        size_t row = idx % m;
        size_t col = idx / m;
        
        // update values
        C.i_[i] = row;
        C.j_[i] = col;
    }
    
    C.m_ = m;
    C.n_ = n;
    return C;
}

cArray CooMatrix::todense () const
{
    // return column-major dense representations
    cArray v (m_ * n_);
    
    size_t N = i_.size();
    for (size_t i = 0; i < N; i++)
        v[i_[i] + j_[i] * m_] += x_[i];
    
    return v;
}

#ifndef NO_UMFPACK
CooMatrix& CooMatrix::operator *= (const cArrayView B)
{
    return *this = this->dot(B);
}
#endif

#ifndef NO_UMFPACK
CooMatrix CooMatrix::dot (const cArrayView B) const
{
    // FIXME: This is a memory INEFFICIENT method.
    // NOTE: Row-major storage assumed for B.
    
    // volumes
    size_t A_vol = x_.size();
    size_t B_vol = B.size();
    
    // check B shape
    assert(B_vol % n_ == 0);
    
    // create output matrix
    unsigned C_rows = m_;
    unsigned C_cols = B_vol / n_;
    CooMatrix C(C_rows, C_cols);
    
    // for all elements of A
    for (size_t i = 0; i < A_vol; i++)
    {
        unsigned row = i_[i];
        unsigned col = j_[i];
        
        // for all columns of B
        for (unsigned icol = 0; icol < C_cols; icol++)
        {
            C.add
            (
                row, icol,
                x_[i] * B[col*C_cols + icol]
            );
        }
    }
    
    // summation is done by shaking
    return C.shake();
}
#endif

Complex CooMatrix::ddot (CooMatrix const & B) const
{
    assert(m_ == B.m_);
    assert(n_ == B.n_);
    
    // sort by i_ and j_
    if (not sorted() or not B.sorted())
        throw exception("[CooMatrix] Sort matrices before ddot!");
        
    Complex result = 0;
    
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

#ifndef NO_UMFPACK
CooMatrix CooMatrix::shake () const
{
    // ugly and memory inefficient method... FIXME
    return tocsc().tocoo();
}
#endif

bool CooMatrix::hdfsave (const char* name) const
{
#ifndef NO_HDF
    try
    {
        HDFFile hdf(name, HDFFile::overwrite);
        
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
                reinterpret_cast<double const*>(&x_),
                x_.size() * 2
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

bool CooMatrix::hdfload (const char* name)
{
    sorted_ = false;
    
#ifndef NO_HDF
    try
    {
        HDFFile hdf(name, HDFFile::readonly);
        
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
                reinterpret_cast<double*>(&(x_[0])),
                x_.size() * 2
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

SymDiaMatrix::SymDiaMatrix ()
    : n_(0), elems_(0), idiag_(0), name_(), dptrs_(0)
{
    
}

SymDiaMatrix::SymDiaMatrix (int n)
    : n_(n), elems_(n), idiag_(0), name_(), dptrs_(0)
{
    
}

SymDiaMatrix::SymDiaMatrix (int n, const iArrayView id)
    : n_(n), idiag_(id), name_()
{
    // compute needed number of elements
    std::size_t vol = 0;
    for (int d : idiag_)
        vol += n*n - d;
    
    // resize storage
    elems_.resize(vol);
    
    // setup diagonal pointers
    setup_dptrs_();
}

SymDiaMatrix::SymDiaMatrix (int n, const iArrayView id, const cArrayView v)
    : n_(n), elems_(v), idiag_(id), name_()
{
    // setup diagonal pointers
    setup_dptrs_();
}

SymDiaMatrix::SymDiaMatrix(SymDiaMatrix const & A)
    : n_(A.n_), elems_(A.elems_), idiag_(A.idiag_), name_()
{
    // setup diagonal pointers
    setup_dptrs_();
}

SymDiaMatrix::SymDiaMatrix(SymDiaMatrix&& A)
    : n_(std::move(A.n_)), elems_(std::move(A.elems_)), idiag_(std::move(A.idiag_)), name_(A.name_)
{
    // setup diagonal pointers
    setup_dptrs_();
}

void SymDiaMatrix::setup_dptrs_()
{
    // set data pointers for all (upper and lower) diagonals
    int d = idiag_.size() - 1;
    dptrs_.resize(d + 1);
    Complex * ptr = elems_.data();
    for (int i = 0; i <= d; i++)
    {
        dptrs_[i] = ptr;
        ptr += n_ - idiag_[i];
    }
}

SymDiaMatrix const & SymDiaMatrix::operator += (SymDiaMatrix const & B)
{
    // check sizes
    if (size() != B.size())
        throw exception ("[SymDiaMatrix::operator+=] Unequal sizes!");
    
    // check diagonals
    if (diag().size() == B.diag().size() and all(diag() == B.diag()))
    {
        // fast version
        for (size_t i = 0; i < elems_.size(); i++)
            elems_[i] += B.elems_[i];
    }
    else
    {
        // general version
        cArrays diags (std::max(diag().size(), B.diag().size()));
        std::set<int> idiags;
        
        // add all A's diagonals
        int dataptr = 0;
        for (auto id : idiag_)
        {
            idiags.insert(id);
            diags[id] = cArrayView(elems_, dataptr, n_ - id);
            dataptr += n_ - id;
        }
        
        // add all B's diagonals
        dataptr = 0;
        for (auto id : B.idiag_)
        {
            idiags.insert(id);
            
            if (diags[id].size() == 0)
                diags[id] = cArrayView(elems_, dataptr, n_ - id);
            else
                diags[id] += cArrayView(elems_, dataptr, n_ - id);
            
            dataptr += n_ - id;
        }
        
        // use new diagonals
        idiag_ = iArray(idiags.begin(), idiags.end());
        
        // use new data
        elems_ = join(diags);
        setup_dptrs_();
    }
    
    return *this;
}

SymDiaMatrix const & SymDiaMatrix::operator -= (SymDiaMatrix const & B)
{
    // check sizes
    if (size() != B.size())
        throw exception ("[SymDiaMatrix::operator-=] Unequal sizes!");
    
    // check diagonals
    if (diag().size() == B.diag().size() and all(diag() == B.diag()))
    {
        // fast version
        for (size_t i = 0; i < elems_.size(); i++)
            elems_[i] -= B.elems_[i];
    }
    else
    {
        // general version
        cArrays diags(std::max(diag().size(), B.diag().size()));
        std::set<int> idiags;
        
        // add all A's diagonals
        int dataptr = 0;
        for (auto id : idiag_)
        {
            idiags.insert(id);
            
            diags[id] = cArray(elems_.data() + dataptr, elems_.begin() + dataptr + n_ - id);
            
            dataptr += n_ - id;
        }
        
        // add all B's diagonals
        dataptr = 0;
        for (auto id : B.idiag_)
        {
            idiags.insert(id);
            
            if (diags[id].size() == 0)
                diags[id] = -cArray(elems_.data() + dataptr, elems_.begin() + dataptr + n_ - id);
            else
                diags[id] -= cArray(elems_.data() + dataptr, elems_.begin() + dataptr + n_ - id);
            
            dataptr += n_ - id;
        }
        
        // use new diagonals
        idiag_ = iArray(idiags.begin(), idiags.end());
        
        // use new data
        elems_ = join(diags);
        setup_dptrs_();
    }
    
    return *this;
}

SymDiaMatrix const & SymDiaMatrix::operator = (SymDiaMatrix&& A)
{
    n_ = std::move(A.n_);
    elems_ = std::move(A.elems_);
    idiag_ = std::move(A.idiag_);
    setup_dptrs_();
    return *this;
}

SymDiaMatrix const & SymDiaMatrix::operator = (SymDiaMatrix const & A)
{
    n_ = A.n_;
    elems_ = A.elems_;
    idiag_ = A.idiag_;
    setup_dptrs_();
    return *this;
}

SymDiaMatrix operator + (SymDiaMatrix const & A, SymDiaMatrix const & B)
{
    A.is_compatible(B);
    return SymDiaMatrix(A.size(), A.diag(), A.data() + B.data());
}

SymDiaMatrix operator - (SymDiaMatrix const & A, SymDiaMatrix const & B)
{
    A.is_compatible(B);
    return SymDiaMatrix(A.size(), A.diag(), A.data() - B.data());
}

SymDiaMatrix operator * (double z, SymDiaMatrix const & A)
{
    return SymDiaMatrix(A.size(), A.diag(), z * A.data());
}

SymDiaMatrix operator * (Complex z, SymDiaMatrix const & A)
{
    return SymDiaMatrix(A.size(), A.diag(), z * A.data());
}

SymDiaMatrix operator * (SymDiaMatrix const & A, SymDiaMatrix const & B)
{
    // FIXME : write an optimized routine
//     return (A.tocoo().tocsr() * B.tocoo().tocsc()).tocoo().todia();
    
    // check dimensions
    if (A.size() != B.size())
        throw exception ("Cannot multiply matrices: A's column count != B's row count.");
    
    // compute bandwidth, half-bandwidth and other parameters of the matrices
    int C_size = A.size();
    int C_bw = std::min(A.bandwidth() + B.bandwidth() - 1, 2 * C_size - 1);
    int C_hbw = C_bw / 2;
    int C_nnz = (C_hbw + 1) * (C_size + C_size - C_hbw) / 2;
    int A_Ndiag = A.diag().size();
    
    // create the output matrix
    SymDiaMatrix C (
        A.size(),
        linspace<int>(0, C_hbw, C_hbw + 1),
        cArray(C_nnz)
    );
    
    // for all (main or upper) output diagonals
    for (int d = 0; d <= C_hbw; d++)
    {
        // for all A's diagonals
        for (int idA = -(A_Ndiag - 1); idA <= A_Ndiag - 1; idA++)
        {
            // get diagonal label
            int dA = A.diag(std::abs(idA)) * signum(idA);
            
            // determine corresponding B's diagonal (shifted from A's by 'd'), and its label
            int dB = dA - d;
            int const * bptr = std::find(B.diag().begin(), B.diag().end(), std::abs(dB));
            if (bptr == B.diag().end())
                continue;
            int idB = signum(dB) * (bptr - B.diag().begin());
            
            // get starting/ending columns of the diagonals
            int start_column = std::max
            (
                (dA > 0 ? dA : 0),
                (dB > 0 ? dB : 0)
            );
            int end_column = std::min
            (
                (dA > 0 ? A.size() - 1 : A.size() - 1 + dA),
                (dB > 0 ? B.size() - 1 : B.size() - 1 + dB)
            );
            int Nelems = end_column - start_column + 1;
            if (Nelems <= 0)
                continue;
            
            // get pointers
            Complex const * restrict pdA = A.dptr(std::abs(idA)) + (dA > 0 ? start_column - dA : start_column);
            Complex const * restrict pdB = B.dptr(std::abs(idB)) + (dB > 0 ? start_column - dB : start_column);
            Complex * restrict pdC = C.dptr(d) + start_column - dA;
            
            // add multiple
            for (int i = 0; i < Nelems; i++)
                pdC[i] += pdA[i] * pdB[i];
        }
    }
    
    return C;
}

bool SymDiaMatrix::is_compatible (SymDiaMatrix const & B) const
{
    if (n_ != B.n_)
        throw exception ("[SymDiaMatrix] Unequal ranks.");
    if (idiag_.size() != B.idiag_.size())
        throw exception ("[SymDiaMatrix] Unequal number of diagonals (%d != %d).", idiag_.size(), B.idiag_.size());
    for (size_t i = 0; i < idiag_.size(); i++)
        if (idiag_[i] != B.idiag_[i])
            throw exception ("[SymDiaMatrix] Unequal distribution of diagonals.");
    return true;
}

bool SymDiaMatrix::hdfload (std::string name)
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
        
        // compressed array info
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
        setup_dptrs_();
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

bool SymDiaMatrix::hdfsave (std::string name, bool docompress, int consec) const
{
#ifndef NO_HDF
    HDFFile hdf(name, HDFFile::overwrite);
        
    if (not hdf.valid())
        return false;
        
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
    
    // write compressed elements array (if non-empty); check result
    if (not zero_blocks.empty() and
        not hdf.write("zero_blocks", &(zero_blocks[0]), zero_blocks.size()))
        return false;
    
    if (not elements.empty() and
        not hdf.write("x", &(elements[0]), elements.size()))
        return false;
    
    return true;
#else /* NO_HDF */
    return false;
#endif
}

CooMatrix SymDiaMatrix::tocoo (MatrixTriangle triangle) const
{
    lArray i, j;
    cArray v;
    
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
            
            // add this element to COO (upper triangle)
            if ((id != 0) and (triangle & strict_upper))
            {
                i.push_back(iel);
                j.push_back(iel+id);
                v.push_back(*el);
            }
            
            // add this element to COO (upper triangle)
            if ((id != 0) and (triangle & strict_lower))
            {
                i.push_back(iel+id);
                j.push_back(iel);
                v.push_back(*el);
            }
            
            // main diagonal
            if ((id == 0) and (triangle & diagonal))
            {
                i.push_back(iel);
                j.push_back(iel);
                v.push_back(*el);
            }
            
            // move on to the next element
            el++;
        }
    }
    
    return CooMatrix(n_, n_, i, j, v);
}

cArray SymDiaMatrix::dot (const cArrayView B, MatrixTriangle triangle, bool parallelize) const
{
    // check dimensions
    if ((int)B.size() != n_)
        throw exception ("[SymDiaMatrix::dot] Incompatible dimensions.");
    
    // size of the matrix
    int Nrows = n_;
    int Ndiag = idiag_.size();
    
    // the result
    cArray res(Nrows);
    
    // data pointers
    // - "restricted" for maximization of the cache usage
    // - "aligned" to convince the auto-vectorizer of the worth of the vectorization using SIMD
    // Note that:
    //    - cArray (= NumberArray<Complex>) is aligned on 2*sizeof(Complex) boundary
    //    - GCC needs -fcx-limited-range to auto-vectorize 'complex' operations
    // The option -fcx-limited-range will inhibit some run-time range checking, so, technically,
    // some 'complex' operations may overflow. However, e.g. GNU Fortran compiler never(!) checks
    // for overflows, so why should we here, when we do not expect any.
#if defined(__GNUC__) || defined(__INTEL_COMPILER)
    Complex       *       restrict rp_res    = (Complex*)__builtin_assume_aligned(res.data(),    32u);
    Complex const *       restrict rp_elems_ = (Complex*)__builtin_assume_aligned(elems_.data(), 32u);
    Complex const * const restrict rp_B      = (Complex*)__builtin_assume_aligned(B.data(),      32u);
#else
    Complex       *       restrict rp_res    = &res[0];
    Complex const *       restrict rp_elems_ = &elems_[0];
    Complex const * const restrict rp_B      = &B[0];
#endif
    
    // for all elements in the main diagonal
    if (triangle & diagonal)
    {
       # pragma omp parallel for default (none) schedule (static) firstprivate (Nrows,rp_res,rp_elems_,rp_B) if (parallelize)
        for (int ielem = 0; ielem < Nrows; ielem++)
            rp_res[ielem] = rp_elems_[ielem] * rp_B[ielem];
    }
    
    // if only diagonal multiplication has been requested, return the result
    if (not (triangle & strict_upper) and not (triangle & strict_lower))
        return res;

    // beginning of the first non-main diagonal
    rp_elems_ += Nrows;
    
    // for all other diagonals
    for (int id = 1; id < Ndiag; id++)
    {
        // index of this diagonal
        int idiag = idiag_[id];
        
        // number of elements in the current diagonal
        int Nelem = Nrows - idiag;
        
        // for all elements of the current diagonal
        if (triangle & strict_upper)
        {
            # pragma omp parallel for default (none) schedule (static) firstprivate (Nelem,rp_res,rp_elems_,rp_B,idiag) if (parallelize)
            for (int ielem = 0; ielem < Nelem; ielem++)
                rp_res[ielem]         += rp_elems_[ielem] * rp_B[ielem + idiag];
        }
        if (triangle & strict_lower)
        {
            # pragma omp parallel for default (none) schedule (static) firstprivate (Nelem,rp_res,rp_elems_,rp_B,idiag) if (parallelize)
            for (int ielem = 0; ielem < Nelem; ielem++)
                rp_res[ielem + idiag] += rp_elems_[ielem] * rp_B[ielem];
        }
        
        // move to the beginning of the next diagonal
        rp_elems_ += Nelem;
    }
    
    return res;
}

SymDiaMatrix kron (SymDiaMatrix const & A, SymDiaMatrix const & B)
{
    iArray Cdiags;
    
    // compose diagonals of the result
    for (int i = 0; i < (int)A.diag().size(); i++)
    for (int j = -(int)B.diag().size() + 1; j < (int)B.diag().size(); j++)
    {
        // A's diagonal label
        int dA = A.diag(i);
        
        // B's diagonal label
        int dB = (j < 0 ? -B.diag(std::abs(j)) : B.diag(j));
        
        // C's diagonal label
        int dC = dA * B.size() + dB;
        
        // add the new diagonal
        if (dC >= 0)
            Cdiags.push_back(dC);
    }
    
    // create result matrix
    SymDiaMatrix C (A.size() * B.size(), Cdiags);
    
    // for all A's and B's diagonals
    # pragma omp parallel for collapse (2)
    for (int i = 0; i < (int)A.diag().size(); i++)
    for (int j = -(int)B.diag().size() + 1; j < (int)B.diag().size(); j++)
    {
        // A's diagonal label and pointer
        int dA = A.diag(i);
        Complex const * restrict pA = A.dptr(i);
        
        // B's diagonal label and pointer
        int dB = (j < 0 ? -B.diag(std::abs(j)) : B.diag(j));
        Complex const * restrict pB = B.dptr(std::abs(j));
        
        // C's diagonal label
        int dC = dA * B.size() + dB;
        
        // compute elements on C's current diagonal dC
        if (dC >= 0)
        {
            // C's diagonal pointer
            Complex * const restrict pC = C.dptr(dC);
            
            // for all elements on A's diagonal dA
            for (std::size_t ia = 0; ia < A.size() - dA; ia++)
            {
                // for all elements on B's diagonal dB
                for (std::size_t ib = 0; ib < B.size() - dB; ib++)
                {
                    // get position on the C's diagonal
                    std::size_t ic = ia * B.size() + ib;
                    
                    // compute element
                    pC[ic] = pA[ia] * pB[ib];
                }
            }
        }
    }
    
    // return final result
    return C;
}

SymDiaMatrix SymDiaMatrix::kron (SymDiaMatrix const & B) const
{
    // FIXME this is ugly and inefficient
//     return ::kron (this->tocoo(), B.tocoo()).todia();
    
    return ::kron (*this, B);
}

cArray SymDiaMatrix::lowerSolve (const cArrayView b) const
{
    assert(size() == b.size());
    
    // the solution; handle the unit diagonal by using "b" right away
    cArray x = b;
    
    // pointer to the solution
    Complex * const px = x.data();
    
    // diagonal count and row count
    register int Nrows = size();
    register int Ndiag = diag().size();
    
    // pointer to the matrix diagonal labels
    int const * const pAd = idiag_.data();
    
    // for all matrix rows (or solution elements) starting from the top
    for (register int irow = 0; irow < Nrows; irow++)
    {
        // pointer to the matrix elements
        Complex const * pA = elems_.data() + Nrows;
        
        // for all diagonals (except the main diagonal, which has been already taken care of)
        for (register int id = 1; id < Ndiag; id++)
        {
            // skip diagonals that are not relevant for this matrix row
            if (pAd[id] > irow)
                break;
            
            // pick the correct element on this diagonal corresponding to row "irow"
            px[irow] -= px[irow - pAd[id]] * pA[irow - pAd[id]];
            
            // shift the pointer to the beginning of the next diagonal
            pA += Nrows - pAd[id];
        }
    }
    
    return x;
}

cArray SymDiaMatrix::upperSolve (const cArrayView b) const
{
    assert(size() == b.size());
    
    // the solution; handle the unit diagonal by using "b" right away
    cArray x = b;
    
    // pointer to the solution
    Complex * const px = x.data();
    
    // diagonal count and row count
    register int Nrows = size();
    register int Ndiag = diag().size();
    
    // pointer to the matrix diagonal labels
    int const * const pAd = idiag_.data();
    
    // for all matrix rows (or solution elements) starting from the bottom
    for (register int irow = 0; irow < Nrows; irow++) // "irow" numbered from bottom!
    {
        // pointer to the matrix elements
        Complex const * pA = elems_.data() + Nrows;
        
        // for all diagonals (except the main diagonal, which has been already taken care of)
        for (register int id = 1; id < Ndiag; id++)
        {
            // skip diagonals that are not relevant for this matrix row
            if (pAd[id] > irow)
                break;
            
            // pick the correct element on this diagonal corresponding to row "irow"
            px[Nrows - 1 - irow] -= px[Nrows - 1 - irow + pAd[id]] * pA[Nrows - 1 - irow];
            
            // shift the pointer to the beginning of the next diagonal
            pA += Nrows - pAd[id];
        }
    }
    
    return x;
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

RowMatrix<Complex> SymDiaMatrix::torow (MatrixTriangle triangle) const
{
    RowMatrix<Complex> M(size(), size());
    
    // for all diagonals
    for (int d : diag())
    {
        // main diagonal
        if ((d == 0) and (triangle & diagonal))
        {
            for (unsigned i = 0; i < size(); i++)
                M(i,i) = main_diagonal()[i];
        }
        
        // upper triangle
        if ((d != 0) and (triangle & strict_upper))
        {
            for (unsigned i = 0; i < size() - d; i++)
                M(i,i+d) = dptr(d)[i];
        }
        
        // lower triangle
        if ((d != 0) and (triangle & strict_lower))
        {
            for (unsigned i = 0; i < size() - d; i++)
                M(i+d,i) = dptr(d)[i];
        }
    }
    
    return M;
}

cArray SymDiaMatrix::toPaddedRows () const
{
    // stored diagonals count
    int ndiag = diag().size();
    
    // all diagonals count
    int Ndiag = 2 * ndiag - 1;
    
    // row count
    int Nrows = n_;
    
    // create zero array of the right length
    cArray padr (Nrows * Ndiag);
    
    // add main diagonal
    for (int i = 0; i < Nrows; i++)
        padr[ndiag - 1 + i * Ndiag] = data()[i];
    
    // add non-main diagonals
    for (int idiag = 1; idiag < ndiag; idiag++)
    {
        // get data pointer for this diagonal
        Complex const * const pa = dptr(idiag);
        
        // get diagonal label
        int ldiag = diag(idiag);
        
        // store the upper diagonal elements
        for (int irow = 0; irow + ldiag < Nrows; irow++)
            padr[ndiag - 1 + idiag + irow * Ndiag] = pa[irow];
        
        // store the lower diagonal elements
        for (int irow = ldiag; irow < Nrows; irow++)
            padr[ndiag - 1 - idiag + irow * Ndiag] = pa[irow - ldiag];
    }
    
    return padr;
}

cArray SymDiaMatrix::toPaddedCols () const
{
    // stored diagonals count
    int ndiag = diag().size();
    
    // all diagonals count
    int Ndiag = 2 * ndiag - 1;
    
    // row count
    int Nrows = n_;
    
    // create zero array of the right length
    cArray padr (Nrows * Ndiag);
    
    // add main diagonal
    for (int i = 0; i < Nrows; i++)
        padr[(ndiag - 1) * Nrows + i] = data()[i];
    
    // add non-main diagonals
    for (int idiag = 1; idiag < ndiag; idiag++)
    {
        // get data pointer for this diagonal
        Complex const * const pa = dptr(idiag);
        
        // get diagonal label
        int ldiag = diag(idiag);
        
        // store the upper diagonal elements
        for (int irow = 0; irow + ldiag < Nrows; irow++)
            padr[(ndiag - 1 + idiag) * Nrows + irow] = pa[irow];
        
        // store the lower diagonal elements
        for (int irow = ldiag; irow < Nrows; irow++)
            padr[(ndiag - 1 - idiag) * Nrows + irow] = pa[irow - ldiag];
    }
    
    return padr;
}

void SymDiaMatrix::link (std::string name)
{
    name_ = name;
}
