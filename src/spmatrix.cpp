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

#ifndef NO_PNG
    #include <png++/png.hpp>
#endif

#ifndef NO_UMFPACK
    #include <umfpack.h>
#endif

#include "arrays.h"
#include "hdffile.h"
#include "spmatrix.h"

CooMatrix kron(const CooMatrix& A, const CooMatrix& B)
{
    // shorthands
    size_t Csize = A.i_.size() * B.i_.size();
    size_t Brows = B.m_;
    size_t Bcols = B.n_;
    
    // create temporary matrix to hold the Kronecker product
    CooMatrix C;
    
    // set correct dimensions, pre-allocate space
    C.m_ = A.m_ * B.m_;
    C.n_ = A.n_ * B.n_;
    C.i_ = lArray(Csize);
    C.j_ = lArray(Csize);
    C.x_ = cArray(Csize);
    
    // get iterators
    size_t ic = 0;
    
    // loop over A data
    size_t Asize = A.i_.size();
    for (size_t ia = 0; ia < Asize; ia++)
    {
        // loop over B data
        size_t Bsize = B.i_.size();
        for (size_t ib = 0; ib < Bsize; ib++)
        {
            // compute new row index
            C.i_[ic] = A.i_[ia] * Brows + B.i_[ib];
            
            // compute new column index
            C.j_[ic] = A.j_[ia] * Bcols + B.j_[ib];
            
            // compute product of the two elements
            C.x_[ic] = A.x_[ia] * B.x_[ib];
            
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

cArray CscMatrix::dotT(const cArrayView&  b) const
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
CooMatrix CscMatrix::tocoo() const
{
    // reserve space for the auxiliary (__j) and the output (Ti,Tj,Tx,Tz) arrays
    size_t N = x_.size();
    std::vector<long> Ti(N), Tj(N), __j(N);
    std::vector<Complex> Tx(N);
    
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
    return CooMatrix (m_, n_, Ti, Tj, Tx);
}
#endif

bool CscMatrix::hdfsave(const char* name) const
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
            file.write (
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

bool CscMatrix::hdfload(const char* name)
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
            hdf.read (
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


// ------------------------------------------------------------------------- //

// -- CSR matrix methods --------------------------------------------------- //

// ------------------------------------------------------------------------- //

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

cArray CsrMatrix::dot(const cArrayView& b) const
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
    int idx_min = Mat->p_[irow];
    int idx_max = Mat->p_[irow + 1];
    
    // clear memory
    for (int icol = 0; icol < (int)Mat->cols(); icol++)
        buffer[icol] = 1;
    
    // for all nonzero columns
    for (int idx = idx_min; idx < idx_max; idx++)
        if (std::abs(Mat->x_[idx]) > Threshold)
            buffer[Mat->i_[idx]] = 0;
        
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

#ifndef NO_UMFPACK
CsrMatrix::LUft CsrMatrix::factorize(double droptol) const
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
    status = umfpack_zl_symbolic (
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
    status = umfpack_zl_numeric (
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
    return LUft(this, Numeric);
}
#endif

#ifndef NO_UMFPACK
cArray CsrMatrix::solve(const cArray&  b, size_t eqs) const
{
    // only square matrices are allowed
    assert(m_ == n_);
    
    // compute the LU factorization
    LUft luft = factorize();
    
    // solve the equations
    cArray solution = luft.solve(b, eqs);
    
    // cleanup and return
    luft.free();
    return solution;
}
#endif

void CsrMatrix::write(const char* filename) const
{
    FILE *f = fopen(filename, "w");
    fprintf(f, "# Matrix %ld × %ld with %ld nonzero elements:\n\n", m_, n_, x_.size());
    for (unsigned irow = 0; irow < m_; irow++)
    {
        size_t idx1 = p_[irow];
        size_t idx2 = p_[irow + 1];
        
        for (size_t idx = idx1; idx < idx2; idx++)
            fprintf(f, "%d\t%ld\t%g\t%g\n", irow, i_[idx], x_[idx].real(), x_[idx].imag());
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
            hdf.write (
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

bool CsrMatrix::hdfload(const char* name)
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
            hdf.read (
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

double CsrMatrix::norm() const
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

cArray CsrMatrix::upperSolve(cArrayView const &  b) const
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

cArray CsrMatrix::lowerSolve(cArrayView const & b) const
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

cArray CsrMatrix::diag() const
{
    cArray D ( std::min(m_, n_) );
    
    for (size_t irow = 0; irow < (size_t)m_; irow++)
        for (size_t idx = p_[irow]; idx < (size_t)p_[irow+1]; idx++)
            if ((size_t)i_[idx] == irow)
                D[irow] = x_[idx];
    
    return D;
}

#ifndef NO_UMFPACK
CooMatrix CsrMatrix::tocoo() const
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

CsrMatrix CsrMatrix::sparse_like(const CsrMatrix& B) const
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

Complex CsrMatrix::operator() (unsigned i, unsigned j) const
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


// ------------------------------------------------------------------------- //

// -- COO matrix methods --------------------------------------------------- //

// ------------------------------------------------------------------------- //


#ifndef NO_UMFPACK
CscMatrix CooMatrix::tocsc() const
{
    size_t nz = x_.size();
    
    // CSC matrix data
    lArray Ap(n_ + 1), Ai(nz);
    cArray Ax(nz);
    
    // do we have any elements at all?
    if (nz != 0)
    {
    
        long status = umfpack_zl_triplet_to_col (
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
CsrMatrix CooMatrix::tocsr() const
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

SymDiaMatrix CooMatrix::todia(MatrixTriangle triangle) const
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
    return SymDiaMatrix (
        m_,
        Array<int>(diags.cbegin(), diags.cend()),
        join(elems)
    );
}

void CooMatrix::write(const char* filename) const
{
    std::ofstream f(filename);
    
    f << "# Matrix " << m_ << " × " << n_ << " with " << x_.size() << " nonzero elements:\n\n";
    
    for (size_t i = 0; i < i_.size(); i++)
        f << i_[i] << "\t" << j_[i] << "\t" << x_[i].real() << "\t" << x_[i].imag() << "\n";
}


CooMatrix CooMatrix::reshape(size_t m, size_t n) const
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

cArray CooMatrix::todense() const
{
    // return column-major dense representations
    cArray v (m_ * n_);
    
    size_t N = i_.size();
    for (size_t i = 0; i < N; i++)
        v[i_[i] + j_[i] * m_] += x_[i];
    
    return v;
}

#ifndef NO_UMFPACK
CooMatrix& CooMatrix::operator *= (cArray const &  B)
{
    return *this = this->dot(B);
}
#endif

#ifndef NO_UMFPACK
CooMatrix CooMatrix::dot(cArrayView const & B) const
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
            C.add (
                row, icol,
                x_[i] * B[col*C_cols + icol]
            );
        }
    }
    
    // summation is done by shaking
    return C.shake();
}
#endif

Complex CooMatrix::ddot(CooMatrix const & B) const
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
CooMatrix CooMatrix::shake() const
{
    // ugly and memory inefficient method... FIXME
    return tocsc().tocoo();
}
#endif

bool CooMatrix::hdfsave(const char* name) const
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
            hdf.write (
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

bool CooMatrix::hdfload(const char* name)
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
            hdf.read (
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

SymDiaMatrix::SymDiaMatrix() : n_(0), elems_(0), idiag_(0) {}

SymDiaMatrix::SymDiaMatrix(int n) : n_(n), elems_(0), idiag_(0) {}

SymDiaMatrix::SymDiaMatrix(int n, ArrayView<int> const & id, ArrayView<Complex> const & v)
    : n_(n), elems_(v), idiag_(id) {}

SymDiaMatrix::SymDiaMatrix(SymDiaMatrix const & A)
    : n_(A.n_), elems_(A.elems_), idiag_(A.idiag_) {}

SymDiaMatrix::SymDiaMatrix(SymDiaMatrix&& A)
    : n_(std::move(A.n_)), elems_(std::move(A.elems_)), idiag_(std::move(A.idiag_)) {}

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
    n_ = std::move(A.n_);
    elems_ = std::move(A.elems_);
    idiag_ = std::move(A.idiag_);
    return *this;
}

SymDiaMatrix const & SymDiaMatrix::operator = (SymDiaMatrix const & A)
{
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
            if ((id > 0) and (triangle & strict_upper))
            {
                i.push_back(iel);
                j.push_back(iel+id);
                v.push_back(*el);
            }
            
            // add this element to COO (upper triangle)
            if ((id > 0) and (triangle & strict_lower))
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

cArray SymDiaMatrix::dot(cArrayView const & B, MatrixTriangle triangle) const
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
    if (triangle & diagonal)
    {
        for (int ielem = 0; ielem < n_; ielem++)
            rp_res[ielem] = rp_elems_[ielem] * rp_B[ielem];
    }
    
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
            if (triangle & strict_upper)
            {
                rp_res[ielem]         += rp_elems_[beg + ielem] * rp_B[ielem + idiag];
            }
            if (triangle & strict_lower)
            {
                rp_res[ielem + idiag] += rp_elems_[beg + ielem] * rp_B[ielem];
            }
        }
        
        // move to the beginning of the next diagonal
        beg += Nelem;
    }
    
    return res;
}

SymDiaMatrix SymDiaMatrix::kron (SymDiaMatrix const & B) const
{
    // FIXME this is ugly and inefficient
    return ::kron (this->tocoo(), B.tocoo()).todia();
}

cArray SymDiaMatrix::lowerSolve(cArrayView const & b) const
{
    assert(size() == b.size());
    
    // the solution; handle the unit diagonal by using "b" right away
    cArray x = b;
    
    // for all matrix rows (or solution elements) starting from the top
    for (size_t irow = 0; irow < b.size(); irow++)
    {
        // pointer to the beginning of the diagonal data
        size_t dptr = size();
        
        // for all diagonals (except the main diagonal, which has been already taken care of)
        for (size_t id = 1; id < idiag_.size(); id++)
        {
            // skip diagonals that are not relevant for this matrix row
            if (idiag_[id] > (int)irow)
                break;
            
            // pick the correct element on this diagonal corresponding to row "irow"
            x[irow] -= x[irow - idiag_[id]] * elems_[dptr + irow - idiag_[id]];
            
            // shift the pointer to the beginning of the next diagonal
            dptr += size() - idiag_[id];
        }
    }
    
    return x;
}

cArray SymDiaMatrix::upperSolve(cArrayView const & b) const
{
    assert(size() == b.size());
    
    // the solution; handle the unit diagonal by using "b" right away
    cArray x = b;
    
    // for all matrix rows (or solution elements) starting from the bottom
    for (size_t irow = 0; irow < b.size(); irow++) // "irow" numbered from bottom!
    {
        // pointer to the beginning of the diagonal data
        size_t dptr = size();
        
        // for all diagonals (except the main diagonal, which has been already taken care of)
        for (size_t id = 1; id < idiag_.size(); id++)
        {
            // skip diagonals that are not relevant for this matrix row
            if (idiag_[id] > (int)irow)
                break;
            
            // pick the correct element on this diagonal corresponding to row "irow"
            x[size() - 1 - irow] -= x[size() - 1 - irow + idiag_[id]] * elems_[dptr + size() - 1 - irow];
            
            // shift the pointer to the beginning of the next diagonal
            dptr += size() - idiag_[id];
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

cArray iChol(cArrayView const & A, lArrayView const & I, lArrayView const & P)
{
    // this will be returned
    cArray LD(A.size());
    
    // check lengths
    assert(A.size() == (size_t)P.back());
    assert(I.size() == (size_t)P.back());
    
    // current row
    int irow = 0;
    
    // for all elements of the output array
    for (int pos = 0; pos < (int)LD.size(); pos++)
    {
        // get column index of this element
        int icol = I[pos];
        
        // is this a diagonal?
        if (icol == irow)
        {
            // Compute an element of D
            // - start by copying corresponding coefficient from A
            LD[pos] = A[pos];
            
            // - continue by subtracting all existing contributions from the current row
            //   (loop over ELEMENTS)
            for (int ielem = P[irow]; ielem < pos; ielem++)
                LD[pos] -= LD[ielem] * LD[ielem] * LD[P[I[ielem]+1]-1];
            
            irow++;
        }
        else
        {
            // Compute an element of L
            // - start by copying corresponding coefficient from A
            LD[pos] = A[pos];
            
            // - continue by subtracting all existing contributions
            //   (loop over COLUMNS)
            int pos1 = P[irow], pos2 = P[icol];
            while (pos1 < P[irow+1] and I[pos1] < icol and pos2 < P[icol+1] and I[pos2] < icol)
            {
                if (I[pos1] < I[pos2])
                {
                    pos1++;
                }
                else if (I[pos1] > I[pos2])
                {
                    pos2++;
                }
                else
                {
                    LD[pos] -= LD[pos1] * LD[pos2] * LD[P[I[pos1]+1]-1];
                    pos1++; pos2++;
                }
            }
            
            // - finish by dividing by the diagonal element
            LD[pos] /= LD[P[icol+1]-1];
        }
    }
    
    return LD;
}
