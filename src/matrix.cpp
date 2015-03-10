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

#ifdef _OPENMP
    #include <omp.h>
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

template<> void ColMatrix<Complex>::invert (ColMatrix<Complex> & inv) const
{
#ifndef NO_LAPACK
    if (rows() != cols())
        HexException("Only square matrices can be diagonalized.");
    
    int N = rows();
    iArray IPIV(N);
    int LWORK = -1;
    cArray WORK(1);
    int INFO;
    
    // create copy of current matrix (will be overwritten)
    inv = *this;
    
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
    
#else
    Exception("Cannot invert dense matrix without LAPACK support.");
#endif
}

// specialization of diagonalization of complex dense column-major ("Fortran-like") matrices
template<> void ColMatrix<Complex>::diagonalize
(
    NumberArray<Complex> & W,
    ColMatrix<Complex> * VL,
    ColMatrix<Complex> * VR
) const
{
#ifndef NO_LAPACK
    if (rows() != cols())
        HexException("Only square matrices can be diagonalized.");
    
    int N = rows();
    cArray WORK(1);
    rArray RWORK(2*N);
    int LWORK = -1, INFO;
    char JOBL = 'N', JOBR = 'N';
    Complex *pVL = nullptr, *pVR = nullptr;
    
    W.resize(N);
    
    if (VL != nullptr)
    {
        // compute left eigenvectors
        JOBL = 'V';
        
        // resize output data
        if (VL->rows() != N or VL->cols() != N)
            *VL = std::move(ColMatrix<Complex>(N,N));
        pVL = VL->begin();
    }
    
    if (VR != nullptr)
    {
        // compute right eigenvectors
        JOBR = 'V';
        
        // resize output data
        if (VR->rows() != N or VR->cols() != N)
            *VR = std::move(ColMatrix<Complex>(N,N));
        pVR = VR->begin();
    }
    
    // copy matrix data (it will be overwritten)
    cArray A = data();
    
    // get work size
    zgeev_
    (
        &JOBL,      // compute left eigenvectors
        &JOBR,      // compute right eigenvectors
        &N,         // order of the matrix
        &A[0],      // matrix A elements
        &N,         // leading dimension of the matrix A
        &W[0],      // eigenvalues
        pVL,        // left eigenvectors elements
        &N,         // leading dimension of the matrix VL
        pVR,        // right eigenvectors elements
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
        &JOBL,      // compute left eigenvectors
        &JOBR,      // compute right eigenvectors
        &N,         // order of the matrix
        &A[0],      // matrix A elements
        &N,         // leading dimension of the matrix A
        &W[0],      // eigenvalues
        pVL,        // left eigenvectors elements
        &N,         // leading dimension of the matrix VL
        pVR,        // right eigenvectors elements
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
#else
    Exception("Cannot diagonalize dense matrix without LAPACK support.");
#endif
}

void dense_kron_dot
(
    int A_rows, int A_cols, Complex const * A_data,
    int B_rows, int B_cols, Complex const * B_data,
    Complex const * v_data,
    Complex       * w_data
)
{
    // work matrix
    static int C_rows = 0, C_cols = 0;
    static Complex * C_data = nullptr;
    
    // realloc the work matrix if necessary
    if (B_rows * A_cols != C_rows * C_cols)
    {
        delete [] C_data;
        C_data = new Complex [B_rows * A_cols];
    }
    
    // update sizes of the intermediate matrix
    C_rows = B_rows; C_cols = A_cols;
    
    // skip invalid inputs (used to dealloc the work matrix)
    if (A_rows == 0 or A_cols == 0 or A_data == nullptr or
        B_rows == 0 or B_cols == 0 or B_data == nullptr or
        v_data == nullptr or w_data == nullptr)
        return;
    
    // auxiliary variables
    Complex alpha = 1, beta = 0;
    char norm = 'N', trans = 'T';
    
    // C = V * A^T
    {
        int m = B_rows, k = B_cols, n = A_rows;
        zgemm_
        (
            &norm, &norm, &m, &n, &k,
            &alpha, const_cast<Complex*>(v_data), &n,
            const_cast<Complex*>(A_data), &k,
            &beta, C_data, &m
        );
    }
    
    // W = B * C
    {
        int m = A_rows, k = A_cols, n = C_cols;
        zgemm_
        (
            &trans, &norm, &m, &n, &k,
            &alpha, const_cast<Complex*>(B_data), &n,
            C_data, &n,
            &beta, w_data, &m
        );
    }
}

cArray kron_dot (SymBandMatrix const & A, SymBandMatrix const & B, const cArrayView v)
{
    cArray w (v.size());
    
    std::size_t A_size = A.size();
    std::size_t B_size = B.size();
    
    # pragma omp parallel for collapse (2)
    for (std::size_t i = 0; i < A_size; i++)
    for (std::size_t j = 0; j < B_size; j++)
    {
        // iteration bounds
        std::size_t kmin = (i >= A.halfbw() ? i - A.halfbw() : 0);
        std::size_t lmin = (j >= B.halfbw() ? j - B.halfbw() : 0);
        std::size_t kmax = std::min(A.size() - 1, i + A.halfbw() - 1);
        std::size_t lmax = std::min(B.size() - 1, j + B.halfbw() - 1);
        
        // calculate the scalar product
        Complex res = 0;
        for (std::size_t k = kmin; k <= kmax; k++)
        for (std::size_t l = lmin; l <= lmax; l++)
            res += A(i,k) * B(j,l) * v[k * B.size() + l];
        
        // save result
        w[i * B.size() + j] = res;
    }
    
    return w;
}

//
// Sparse matrix routines.
//

CooMatrix kron (const CooMatrix& A, const CooMatrix& B)
{
    // shorthands
    std::size_t Asize = A.v().size();
    std::size_t Bsize = B.v().size();
    std::size_t Csize = Asize * Bsize;
    std::size_t Brows = B.rows();
    std::size_t Bcols = B.cols();
    
    // set correct dimensions, pre-allocate space
    int m = A.rows() * B.rows();
    int n = A.cols() * B.cols();
    lArray C_i (Csize), C_j (Csize);
    cArray C_x (Csize);
    
    // get pointers
    std::int64_t const * restrict pA_i = A.i().data();
    std::int64_t const * restrict pB_i = B.i().data();
    std::int64_t       * restrict pC_i = C_i.data();
    std::int64_t const * restrict pA_j = A.j().data();
    std::int64_t const * restrict pB_j = B.j().data();
    std::int64_t       * restrict pC_j = C_j.data();
    Complex const * restrict pA_x = A.v().data();
    Complex const * restrict pB_x = B.v().data();
    Complex       * restrict pC_x = C_x.data();
    
    // loop over A data
    for (std::size_t ia = 0; ia < Asize; ia++)
    {
        // loop over B data
        for (std::size_t ib = 0; ib < Bsize; ib++)
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

CooMatrix eye (std::size_t N)
{
    return CooMatrix(N,N).symm_populate_band
    (
        0,
        [](unsigned i, unsigned j) -> Complex { return 1.; }
    );
}

CooMatrix stairs (std::size_t N)
{
    return CooMatrix(N,N).symm_populate_band
    (
        0,
        [](unsigned i, unsigned j) -> Complex { return (i == j) ? i : 0.; }
    );
}

// 
// CSC matrix routines
// 

CscMatrix & CscMatrix::operator *= (double r)
{
    std::size_t N = i_.size();
    
    for (std::size_t i = 0; i < N; i++)
        x_[i] *= r;
    
    return *this;
}

CscMatrix & CscMatrix::operator &= (const CscMatrix&  B)
{
    std::size_t N = i_.size();
    
    assert(m_ == B.m_);
    assert(n_ == B.n_);
    assert(N == B.i_.size());
    
    for (std::size_t i = 0; i < N; i++)
    {
        assert(i_[i] == B.i_[i]);
        
        x_[i] += B.x_[i];
    }
    
    return *this;
}

CscMatrix & CscMatrix::operator ^= (const CscMatrix&  B)
{
    std::size_t N = i_.size();
    
    assert(m_ == B.m_);
    assert(n_ == B.n_);
    assert(N == B.i_.size());
    
    for (std::size_t i = 0; i < N; i++)
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
        std::size_t idx1 = p_[icol];
        std::size_t idx2 = p_[icol+1];
        
        // for all nonzero elements in this column
        for (std::size_t idx = idx1; idx < idx2; idx++)
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
    std::size_t N = x_.size();
    std::vector<std::int64_t> Ti(N), Tj(N), J(N);
    std::vector<Complex> Tx(N);
    
    // do we have any elements at all?
    if (N != 0)
    {
    
        // do the conversion
        std::int64_t status = umfpack_zl_col_to_triplet(n_, p_.data(), J.data());
        
        // check success
        if (status != 0)
        {
            std::cerr << "\n[CscMatrix::tocoo] Exit code " << status << "\n";
            umfpack_zl_report_status(0, status);
        }
        
        // copy only non-zero entries to output arrays
        std::size_t nz = 0;
        for (std::size_t i = 0; i < N; i++)
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
    std::size_t N = i_.size();
    
    for (std::size_t i = 0; i < N; i++)
        x_[i] *= r;
    
    return *this;
}

CsrMatrix & CsrMatrix::operator &= (CsrMatrix const &  B)
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

CsrMatrix & CsrMatrix::operator ^= (CsrMatrix const &  B)
{
    std::size_t N = i_.size();
    
    assert(m_ == B.m_);
    assert(n_ == B.n_);
    assert(N == B.i_.size());
    
    for (std::size_t i = 0; i < N; i++)
        x_[i] -= B.x_[i];
    
    return *this;
}

cArray CsrMatrix::dot (cArrayView const & b) const
{
    // create output array
    cArray c(m_);
    
    for (unsigned irow = 0; irow < m_; irow++)
    {
        std::size_t idx1 = p_[irow];
        std::size_t idx2 = p_[irow+1];
        
        // for all nonzero elements in this row
        for (std::size_t idx = idx1; idx < idx2; idx++)
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
    std::int64_t status;
    
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
        std::cerr << "\n[CsrMatrix::factorize] Exit status " << status << "\n";
        umfpack_zl_report_status(0, status);
        std::exit(EXIT_FAILURE);
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
        std::exit(EXIT_FAILURE);
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
        std::size_t idx1 = p_[irow];
        std::size_t idx2 = p_[irow + 1];
        
        for (std::size_t idx = idx1; idx < idx2; idx++)
            out << irow << "\t" << i_[idx] << "\t" << x_[idx].real() << "\t" << x_[idx].imag() << "\n";
    }
    out.close();
}

void CsrMatrix::hdflink (std::string name)
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
    std::size_t N = i_.size();
    double res = 0.;
    
    // return the abs(largest element)
    for (std::size_t i = 0; i < N; i++)
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
    std::size_t N = b.size();
    assert((std::size_t)m_ == N);
    assert((std::size_t)n_ == N);

    // create output array
    cArray x(N);
    
    // loop over rows
    for (std::size_t i = 0; i < N; i++)
    {
        std::size_t row = N - 1 - i;
        Complex accum = 0.;
        
        // get relevant columns of the sparse matrix
        std::size_t idx1 = p_[row];
        std::size_t idx2 = p_[row + 1];
        
        // diagonal element of the matrix
        Complex a = 0.;
        
        // loop over the columns
        for (std::size_t idx = idx1; idx < idx2; idx++)
        {
            // which column is this?
            std::size_t col = i_[idx];
            
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
    std::size_t N = b.size();
    assert((std::size_t)m_ == N);
    assert((std::size_t)n_ == N);
    
    // create output array
    cArray x(N);
    
    // loop over rows
    for (std::size_t row = 0; row < N; row++)
    {
        Complex accum = 0.;
        
        // get relevant columns of the sparse matrix
        std::size_t idx1 = p_[row];
        std::size_t idx2 = p_[row + 1];
        
        // diagonal element of the matrix
        Complex a = 0.;
        
        // loop over the columns
        for (std::size_t idx = idx1; idx < idx2; idx++)
        {
            // which column is this?
            std::size_t col = i_[idx];
            
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
    cArray D (std::min(m_, n_));
    
    for (int irow = 0; irow < m_; irow++)
    for (int idx = p_[irow]; idx < p_[irow+1]; idx++)
    if (i_[idx] == irow)
        D[irow] = x_[idx];
    
    return D;
}

#ifndef NO_UMFPACK
CooMatrix CsrMatrix::tocoo () const
{
    // reserve space for the auxiliary (__j) and the output (Ti,Tj,Tx,Tz) arrays
    std::size_t N = x_.size();
    lArray Ti(N), Tj(N), __j(N);
    cArray Tx(N);
    
    // do we have any elements at all?
    if (N != 0)
    {
    
        // do the conversion
        std::int64_t status = umfpack_zl_col_to_triplet(n_, p_.data(), __j.data());
        
        // check success
        if (status != 0)
        {
            std::cerr << "\n[CscMatrix::tocoo] Exit code " << status << "\n";
            umfpack_zl_report_status(0, status);
        }
        
        // copy only non-zero entries to output arrays
        std::size_t nz = 0;
        for (std::size_t i = 0; i < N; i++)
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
    std::memset(A.x_.data(), 0, A.x_.size() * sizeof(Complex));
    
    // copy all nonzero elements of "this" matrix
    for (unsigned row = 0; row < m_; row++)
    {
        std::size_t idx1 = A.p_[row];
        std::size_t idx2 = A.p_[row+1];
        
        for (std::size_t idx = idx1; idx < idx2; idx++)
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
    std::size_t idx1 = p_[i];
    std::size_t idx2 = p_[i + 1];
    
    // find the correct column ("j")
    auto it = std::lower_bound(i_.begin() + idx1, i_.begin() + idx2, j);

    if (it == i_.end() or *it != j)
        return 0.;
    
    // return the value
    return x_[it - i_.begin()];
}

Complex sparse_row_scalar_product
(
  int n1, Complex const * const restrict x1, std::int64_t const * const restrict i1,
  int n2, Complex const * const restrict x2, std::int64_t const * const restrict i2
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
        HexException("Wrong sizes for matrix multiplication: %d %d, %d %d.", A.rows(), A.cols(), B.rows(), B.cols());
    
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
    std::size_t nz = x_.size();
    
    // CSC matrix data
    lArray Ap(n_ + 1), Ai(nz);
    cArray Ax(nz);
    
    // do we have any elements at all?
    if (nz != 0)
    {
    
        std::int64_t status = umfpack_zl_triplet_to_col
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
        std::size_t N = Ap[n_];
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
    
        std::int64_t status = umfpack_zl_triplet_to_col
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
        std::size_t N = Ap[m_];
        Ai.resize(N);
        Ax.resize(N);
    }
    
    return CsrMatrix(m_, n_, Ap, Ai, Ax);
}
#endif

/*SymBandMatrix CooMatrix::tosymband (MatrixTriangle triangle) const
{
    // diagonal indices
    std::set<int> diags;
    
    // elements per diagonal
    cArrays elems(m_);
    
    // for all nonzero elements
    for (std::size_t i = 0; i < x_.size(); i++)
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
    return SymBandMatrix
    (
        m_,
        Array<int>(diags.cbegin(), diags.cend()),
        join(elems)
    );
}*/

void CooMatrix::write (const char* filename) const
{
    std::ofstream f(filename);
    
    f << "# Matrix " << m_ << " × " << n_ << " with " << x_.size() << " nonzero elements:\n\n";
    
    for (size_t i = 0; i < i_.size(); i++)
        f << i_[i] << "\t" << j_[i] << "\t" << x_[i].real() << "\t" << x_[i].imag() << "\n";
}


CooMatrix CooMatrix::reshape (std::size_t m, std::size_t n) const
{
    CooMatrix C = *this;
    
    // conserved dimensions
    std::size_t N = C.i_.size();
    std::size_t H = C.m_;
    
    // check dimensions
    assert(m * n == C.m_ * C.n_);
    
    // reshape
    for (std::size_t i = 0; i < N; i++)
    {
        // conserved position in column-ordered array
        std::size_t idx = C.i_[i] + C.j_[i] * H;
        
        // new coordinates
        std::size_t row = idx % m;
        std::size_t col = idx / m;
        
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
    
    std::size_t N = i_.size();
    for (std::size_t i = 0; i < N; i++)
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
    std::size_t A_vol = x_.size();
    std::size_t B_vol = B.size();
    
    // check B shape
    assert(B_vol % n_ == 0);
    
    // create output matrix
    unsigned C_rows = m_;
    unsigned C_cols = B_vol / n_;
    CooMatrix C(C_rows, C_cols);
    
    // for all elements of A
    for (std::size_t i = 0; i < A_vol; i++)
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
        HexException("[CooMatrix] Sort matrices before ddot!");
        
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

SymBandMatrix::SymBandMatrix ()
    : n_(0), d_(0), elems_(0), name_()
{}

SymBandMatrix::SymBandMatrix (std::size_t n, std::size_t d)
    : n_(n), d_(d), elems_(n * d), name_()
{}

SymBandMatrix::SymBandMatrix (std::size_t n, std::size_t d, const cArrayView v)
    : n_(n), d_(d), elems_(v), name_()
{}

SymBandMatrix::SymBandMatrix (SymBandMatrix const & A)
    : n_(A.n_), d_(A.d_), elems_(A.elems_), name_()
{}

SymBandMatrix::SymBandMatrix (SymBandMatrix&& A)
    : n_(std::move(A.n_)), d_(std::move(A.d_)), elems_(std::move(A.elems_)), name_(A.name_)
{
}

SymBandMatrix::SymBandMatrix (std::string filename)
    : name_(filename)
{
    # pragma omp critical
    if (not hdfload())
        HexException("Unable to load from %s.", filename.c_str());
}

SymBandMatrix const & SymBandMatrix::operator += (SymBandMatrix const & B)
{
    is_compatible(B);
    elems_ += B.elems_;
    return *this;
}

SymBandMatrix const & SymBandMatrix::operator -= (SymBandMatrix const & B)
{
    is_compatible(B);
    elems_ -= B.elems_;
    return *this;
}

SymBandMatrix const & SymBandMatrix::operator = (SymBandMatrix&& A)
{
    n_ = std::move(A.n_);
    d_ = std::move(A.d_);
    elems_ = std::move(A.elems_);
    return *this;
}

SymBandMatrix const & SymBandMatrix::operator = (SymBandMatrix const & A)
{
    n_ = A.n_;
    d_ = A.d_;
    elems_ = A.elems_;
    return *this;
}

SymBandMatrix operator + (SymBandMatrix const & A, SymBandMatrix const & B)
{
    A.is_compatible(B);
    return SymBandMatrix(A.size(), A.halfbw(), A.data() + B.data());
}

SymBandMatrix operator - (SymBandMatrix const & A, SymBandMatrix const & B)
{
    A.is_compatible(B);
    return SymBandMatrix(A.size(), A.halfbw(), A.data() - B.data());
}

SymBandMatrix operator * (double z, SymBandMatrix const & A)
{
    return SymBandMatrix(A.size(), A.halfbw(), z * A.data());
}

SymBandMatrix operator * (Complex z, SymBandMatrix const & A)
{
    return SymBandMatrix(A.size(), A.halfbw(), z * A.data());
}

bool SymBandMatrix::is_compatible (SymBandMatrix const & B) const
{
    if (n_ != B.n_)
        HexException("Unequal ranks (%d != %d).", n_, B.n_);
    
    if (d_ != B.d_)
        HexException("Unequal half-bandwidths (%d != %d).", d_, B.d_);
    
    return true;
}

#ifndef NO_HDF
bool SymBandMatrix::hdfload (HDFFile & hdf, std::string prefix)
{
    // check the HDF file
    if (not hdf.valid())
        return false;
    
    // set prefix
    if (not prefix.empty())
        hdf.prefix = prefix;
    
    // read dimension
    if (not hdf.read("n", &n_, 1))
        return false;
    
    // read dimension
    if (not hdf.read("d", &d_, 1))
        return false;
    
    // compressed array info
    iArray zero_blocks_re, zero_blocks_im;
    rArray elements_re, elements_im;
    
    // read data from file
    if ((zero_blocks_re.resize(hdf.size("zero_blocks_re")) and not hdf.read("zero_blocks_re", &(zero_blocks_re[0]), zero_blocks_re.size())) or
        (zero_blocks_im.resize(hdf.size("zero_blocks_im")) and not hdf.read("zero_blocks_im", &(zero_blocks_im[0]), zero_blocks_im.size())) or
        (elements_re.resize(hdf.size("re")) and not hdf.read("re", &(elements_re[0]), elements_re.size())) or
        (elements_im.resize(hdf.size("im")) and not hdf.read("im", &(elements_im[0]), elements_im.size())))
        return false;
    
    // decompress
    elems_ = std::move
    (
        interleave
        (
            elements_re.decompress(zero_blocks_re),
            elements_im.decompress(zero_blocks_im)
        )
    );
    
    return true;
}
#endif

bool SymBandMatrix::hdfload (std::string name)
{
#ifndef NO_HDF
    // open the HDF file
    HDFFile hdf(name, HDFFile::readonly);
    if (not hdf.valid())
        return false;
    
    return hdfload(hdf);
#else /* NO_HDF */
    return false;
#endif
}

bool SymBandMatrix::hdfsave (std::string name, HDFFile::FileAccess flags, bool docompress, std::size_t consec) const
{
#ifndef NO_HDF
    HDFFile hdf(name, flags);
    
    if (not hdf.valid())
        return false;
    
    // write dimension and diagonal info
    if (not hdf.write("n", &n_, 1) or
        not hdf.write("d", &d_, 1))
        return false;
    
    // compress elements array
    iArray zero_blocks_re, zero_blocks_im;
    rArray elements_re, elements_im;
    
    if (docompress)
    {
        std::tie(zero_blocks_re, elements_re) = realpart(elems_).compress(consec);
        std::tie(zero_blocks_im, elements_im) = imagpart(elems_).compress(consec);
    }
    else
    {
        elements_re = std::move(realpart(elems_));
        elements_im = std::move(imagpart(elems_));
    }
    
    // write compressed elements array (if non-empty); check result
    if ((not zero_blocks_re.empty() and
         not hdf.write("zero_blocks_re", &(zero_blocks_re[0]), zero_blocks_re.size())) or
       ((not zero_blocks_im.empty() and
         not hdf.write("zero_blocks_im", &(zero_blocks_im[0]), zero_blocks_im.size()))) or
       ((not elements_re.empty() and
         not hdf.write("re", &(elements_re[0]), elements_re.size()))) or
       ((not elements_im.empty() and
        not hdf.write("im", &(elements_im[0]), elements_im.size()))))
        return false;
    
    return true;
#else /* NO_HDF */
    return false;
#endif
}

CooMatrix SymBandMatrix::tocoo (MatrixTriangle triangle) const
{
    lArray I, J;
    cArray V;
    
    Complex const * el = elems_.data();
    
    // for all elements
    for (std::size_t i = 0; i < n_; i++)
    for (std::size_t d = 0; d < d_; d++)
    {
        if (i + d < n_)
        {
            // skip zero elements
            if (*el == 0.)
            {
                el++;
                continue;
            }
            
            // add this element to COO (upper triangle)
            if ((d != 0) and (triangle & strict_upper))
            {
                I.push_back(i);
                J.push_back(i+d);
                V.push_back(*el);
            }
            
            // add this element to COO (upper triangle)
            if ((d != 0) and (triangle & strict_lower))
            {
                I.push_back(i+d);
                J.push_back(i);
                V.push_back(*el);
            }
            
            // main diagonal
            if ((d == 0) and (triangle & diagonal))
            {
                I.push_back(i);
                J.push_back(i);
                V.push_back(*el);
            }
        }
        
        // move on to the next element
        el++;
    }
    
    return CooMatrix(n_, n_, I, J, V);
}

cArray SymBandMatrix::sym_band_dot (int n, int d, const cArrayView M, const cArrayView X)
{
    // check dimensions
    if (X.size() % n != 0)
        HexException("Incompatible dimensions: %d (mat) × %ld (vec). You are probably mixing radial data for different grids.", n, X.size());
    
    // get number of source vectors
    int Nvec = X.size() / n;
    
    // allocate the same number of destination vectors
    cArray Y (X.size());
    
    // skip if there are no data
    if (n == 0 or d == 0)
        return Y;
    
    // auxiliary variables
    char uplo = 'L';
    int k = d - 1, inc = 1;
    Complex alpha = 1, beta = 1;
    
    // for all source vectors (columns)
    for (int j = 0; j < Nvec; j++)
    {
        zsbmv_
        (
            &uplo, &n, &k, &alpha, const_cast<Complex*>(M.data()), &d,
            const_cast<Complex*>(X.data()) + j * n, &inc, &beta,
            Y.data() + j * n, &inc
        );
    }
    
    return Y;
}

cArray SymBandMatrix::dot (const cArrayView X) const
{
    return sym_band_dot(n_, d_, elems_, X);
}

RowMatrix<Complex> SymBandMatrix::torow (MatrixTriangle triangle) const
{
    RowMatrix<Complex> M(n_, n_);
    
    Complex const * el = elems_.data();
    
    // for all elements
    for (std::size_t i = 0; i < n_; i++)
    for (std::size_t d = 0; d < d_; d++)
    {
        if (i + d < n_)
        {
            // main diagonal
            if ((d == 0) and (triangle & diagonal))
                M(i,i) = *el;
            
            // upper triangle
            if ((d != 0) and (triangle & strict_upper))
                M(i,i+d) = *el;
            
            // lower triangle
            if ((d != 0) and (triangle & strict_lower))
                M(i+d,i) = *el;
        }
        // move on to the next elements
        el++;
    }
    
    return M;
}

void SymBandMatrix::hdflink (std::string name)
{
    name_ = name;
}

bool BlockSymBandMatrix::hdfcheck () const
{
    return HDFFile(diskfile_, HDFFile::readonly).valid();
}

cArray BlockSymBandMatrix::dot (cArrayView v, bool parallelize) const
{
    // check vector size
    if (v.size() != size_ * size_)
        HexException("[BlockSymBandMatrix::dot] Different size of matrix and vector in multiplication routine: %ld (mat) != %ld (vec).", size_ * size_, v.size());
    
    // create output vector
    cArray w(v.size());
    
    // check if matrix is empty
    if (halfbw_ == 0)
        return w;
    
#ifndef NO_HDF
    // open data file for reading
    HDFFile * hdf = nullptr;
    if (not inmemory_)
        hdf = new HDFFile (diskfile_, HDFFile::readonly);
#endif
    
    // data
    cArray diskdata;
    
    // write locks
    omp_lock_t lock[size_];
    for (std::size_t i = 0; i < size_; i++)
        omp_init_lock(&lock[i]);
    
    // for all blocks
    for (std::size_t d = 0; d < halfbw_; d++)
    {
        // paralle processing of blocks on this diagonal (only if requested, and they are present in memory)
        # pragma omp parallel for schedule (dynamic,1) if (parallelize && inmemory_)
        for (std::size_t i = 0; i < size_; i++)
        if (i + d < size_)
        {
            // block volume and offset
            std::size_t vol = size_ * halfbw_;
            std::size_t offset = (i * halfbw_ + d) * vol;
            
            // data view of this block
            cArrayView view (data_, offset, vol);
            
#ifndef NO_HDF
            if (not inmemory_)
            {
                // read data from the disk
                diskdata.resize(vol);
                if (not hdf->read("data", &diskdata[0], vol, offset))
                    HexException("Failed to read HDF file \"%s\".\nHDF error stack:\n%s", diskfile_.c_str(), hdf->error().c_str());
                
                // reset view to the new data
                view.reset(vol, diskdata.data());
            }
#endif
            
            cArray product;
            
            // multiply by diagonal block
            if (d == 0)
            {
                product = SymBandMatrix::sym_band_dot
                (
                    size_, halfbw_, view,
                    cArrayView(v, i * size_, size_)
                );
                
                omp_set_lock(&lock[i]);
                cArrayView(w, i * size_, size_) += product;
                omp_unset_lock(&lock[i]);
            }
            
            // multiply by the other diagonals (both symmetries)
            if (d != 0)
            {
                product = SymBandMatrix::sym_band_dot
                (
                    size_, halfbw_, view,
                    cArrayView(v, (i + d) * size_, size_)
                );
                
                omp_set_lock(&lock[i]);
                cArrayView(w, i * size_, size_) += product;
                omp_unset_lock(&lock[i]);
                
                product = SymBandMatrix::sym_band_dot
                (
                    size_, halfbw_, view,
                    cArrayView(v, i * size_, size_)
                );
                
                omp_set_lock(&lock[i + d]);
                cArrayView(w, (i + d) * size_, size_) += product;
                omp_unset_lock(&lock[i + d]);
            }
        }
    }
    
#ifndef NO_HDF
    if (not inmemory_)
        delete hdf;
#endif
    
    for (std::size_t i = 0; i < size_; i++)
        omp_destroy_lock(&lock[i]);
    
    return w;
}

CooMatrix BlockSymBandMatrix::tocoo () const
{
    // number of structurally non-zero blocks (both upper and lower)
    std::size_t nblocks = size_ * (2 * halfbw_ - 1);
    
    // number of structurally non-zero elements (the structure is recursive)
    std::size_t nelem = nblocks * nblocks;
    
    // allocate the ijv arrays
    lArray I; I.reserve(nelem);
    lArray J; J.reserve(nelem); 
    cArray V; V.reserve(nelem);
    
#ifndef NO_HDF
    // open data file for reading
    HDFFile * hdf = nullptr;
    if (not inmemory_)
        hdf = new HDFFile (diskfile_, HDFFile::readonly);
#endif
    
    // for all blocks
    for (std::size_t i = 0; i < size_; i++)
    for (std::size_t d = 0; d < halfbw_; d++)
    if (i + d < size_)
    {
        // data view of this block diagonal
        cArrayView view (data_, (i * halfbw_ + d) * size_ * halfbw_, size_ * halfbw_);
        
#ifndef NO_HDF
        // it may be necessary to load the data from disk
        cArray diskdata;
        if (not inmemory_)
        {
            // read data from the disk
            diskdata.resize(size_ * halfbw_);
            if (not hdf->read("data", &diskdata[0], size_ * halfbw_, (i * halfbw_ + d) * size_ * halfbw_))
                HexException("Failed to read HDF file \"%s\".\nHDF error stack:\n%s", diskfile_.c_str(), hdf->error().c_str());
            
            // reset view to the new data
            view.reset(size_ * halfbw_, diskdata.data());
        }
#endif
        
        // convert the block to COO format
        CooMatrix coo = SymBandMatrix(size_, halfbw_, view).tocoo();
        
        // copy all elements to whole-matrix arrays
        for (std::size_t k = 0; k < coo.v().size(); k++)
        {
            I.push_back( coo.i()[k] +  i      * size_ );
            J.push_back( coo.j()[k] + (i + d) * size_ );
            V.push_back( coo.v()[k]                   );
        }
        if (d != 0)
        for (std::size_t k = 0; k < coo.v().size(); k++)
        {
            I.push_back( coo.i()[k] + (i + d) * size_ );
            J.push_back( coo.j()[k] +  i      * size_ );
            V.push_back( coo.v()[k]                   );
        }
    }
    
#ifndef NO_HDF
    if (not inmemory_)
        delete hdf;
#endif
    
    // compose the final matrix
    return CooMatrix
    (
        size_ * size_,  // number of rows
        size_ * size_,  // number of columns
        std::move(I),   // row indices
        std::move(J),   // column indices
        std::move(V)    // structurally nonzero matrix entries
    );
}

bool BlockSymBandMatrix::hdfinit () const
{
    // create a new file
    HDFFile hdf (diskfile_, HDFFile::overwrite);
    if (not hdf.valid())
        HexException("Cannot open file \"%s\" for writing.\nHDF error stack:\n%s", diskfile_.c_str(), hdf.error().c_str());
    
    // initialize the dataset to its full length
    if (not hdf.write("data", (Complex*)nullptr, size_ * halfbw_ * size_ * halfbw_))
        HexException("Failed to initialize file \"%s\" (size %.1f).\n%s", diskfile_.c_str(), double(size_ * halfbw_)/std::pow(2,30), hdf.error().c_str());
    
    return true;
}

bool BlockSymBandMatrix::hdfsave () const
{
    // open disk file
    HDFFile hdf (diskfile_, HDFFile::overwrite);
    
    // check success
    if (not hdf.valid())
        return false;
    
    // save data
    if (not hdf.write("data", &data_[0], data_.size()))
        return false;
    
    return true;
}


bool BlockSymBandMatrix::hdfload ()
{
    // open disk file
    HDFFile hdf (diskfile_, HDFFile::readonly);
    
    // check success
    if (not hdf.valid())
        return false;
    
    // allocate memory
    data_.resize(size_ * size_ * halfbw_ * halfbw_);
    
    // read data
    if (not hdf.read("data", &data_[0], data_.size()))
        return false;
    
    // turn on the "in memory" flag
    inmemory_ = true;
    
    return true;
}

void BlockSymBandMatrix::drop ()
{
    // release memory
    data_.drop();
    
    // turn off the "in memory" flag
    inmemory_ = false;
}

cArray BlockSymBandMatrix::getBlock (int i) const
{
    if (inmemory_)
    {
        return cArrayView (data_, i * size_ * halfbw_, size_ * halfbw_);
    }
    else
    {
        // open disk file
        HDFFile hdf (diskfile_, HDFFile::readonly);
        
        // check success
        if (not hdf.valid())
            HexException("Cannot open file \"%s\".\nHDF error stack:\n%s", diskfile_.c_str(), hdf.error().c_str());
        
        // create output array
        cArray data (size_ * halfbw_);
        
        // read data
        if (not hdf.read("data", &data[0], size_ * halfbw_, i * size_ * halfbw_))
            HexException("Cannot access block %d in file \"%s\".\nHDF error stack:\n%s", i, diskfile_.c_str(), hdf.error().c_str());
        
        return data;
    }
}

void BlockSymBandMatrix::setBlock (int i, const cArrayView data)
{
    if (data.size() != size_ * halfbw_)
        HexException("Wrong dimensions: %ld != %ld.", data.size(), size_ * halfbw_);
    
    if (inmemory_)
    {
        cArrayView(data_, i * size_ * halfbw_, size_ * halfbw_) = data;
    }
    else
    {
        // check that the file exists
        HDFFile * phdf = new HDFFile (diskfile_, HDFFile::readwrite);
        if (not phdf->valid())
        {
            // create a new file
            phdf = new HDFFile (diskfile_, HDFFile::overwrite);
            if (not phdf->valid())
                HexException("Cannot open file \"%s\" for writing.\nHDF error stack:\n%s", diskfile_.c_str(), phdf->error().c_str());
            
            // initialize the dataset to its full length
            if (not phdf->write("data", (Complex*)nullptr, size_ * halfbw_ * size_ * halfbw_))
                HexException("Failed to initialize file \"%s\" (size %.1f).\n%s", diskfile_.c_str(), double(size_ * halfbw_)/std::pow(2,30), phdf->error().c_str());
        }
        
        // write data
        if (not phdf->write("data", &data[0], size_ * halfbw_, i * size_ * halfbw_))
            HexException("Cannot access block %d in file \"%s\".\nHDF error stack:\n%s", i, diskfile_.c_str(), phdf->error().c_str());
        
        delete phdf;
    }
}
