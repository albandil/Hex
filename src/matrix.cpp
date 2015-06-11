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
#include <memory>
#include <set>
#include <vector>

#include "arrays.h"
#include "hdffile.h"
#include "matrix.h"
#include "misc.h"

// -------------------------------------------------------------------------------------
// Dense matrix routines
//

template<>
void ColMatrix<Complex>::invert (ColMatrix<Complex> & inv) const
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
    HexException("Cannot invert dense matrix without LAPACK support.");
#endif
}

// specialization of diagonalization of complex dense column-major ("Fortran-like") matrices
template<>
void ColMatrix<Complex>::diagonalize
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
    HexException("Cannot diagonalize dense matrix without LAPACK support.");
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

// -------------------------------------------------------------------------------------
// Sparse matrix routines.
//

template<>
cArray SymBandMatrix<Complex>::sym_band_dot (int n, int d, const cArrayView M, const cArrayView X)
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

// -------------------------------------------------------------------------------------
// UMFPACK-dependent functions.
//

#ifdef WITH_UMFPACK
#include <umfpack.h>

template<>
std::shared_ptr<LUft<LU_int_t,Complex>> CsrMatrix<LU_int_t,Complex>::factorize_umfpack (double droptol, void * data) const
{
    // Use standard UMFPACK sequence
    void *Symbolic, *Numeric;
    LU_int_t status;
    
    // get default setting
    double Control[UMFPACK_CONTROL];
    UMFPACK_DEFAULTS_F(Control);
    
    // modify the drop tolerance
    Control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_SYMMETRIC;
    Control[UMFPACK_DROPTOL] = droptol;
    
    // analyze the sparse structure
    status = UMFPACK_SYMBOLIC_F
    (
        m_, n_,                    // matrix dimensions
        p_.data(), i_.data(),        // column and row indices
        reinterpret_cast<const double*>(x_.data()), 0,    // matrix data
        &Symbolic, Control, nullptr                // UMFPACK internals
    );
    if (status != 0)
    {
        std::cerr << "\n[CsrMatrix::factorize] Exit status " << status << std::endl;
        UMFPACK_REPORT_STATUS_F(0, status);
        std::exit(EXIT_FAILURE);
    }
    
    // do some factorizations
    status = UMFPACK_NUMERIC_F
    (
        p_.data(), i_.data(),    // column and row indices
        reinterpret_cast<const double*>(x_.data()), 0,    // matrix data
        Symbolic, &Numeric, Control, nullptr    // UMFPACK internals
    );
    if (status != 0)
    {
        std::cerr << "\n[CsrMatrix::factorize] Exit status " << status << std::endl;
        UMFPACK_REPORT_STATUS_F(0, status);
        std::exit(EXIT_FAILURE);
    }
    
    // release unused data
    UMFPACK_FREE_SYMBOLIC_F(&Symbolic);
    
    // create a new LU factorization container
    LUft<LU_int_t,Complex> * lu_ptr = new LUft_UMFPACK<LU_int_t,Complex>(this, Numeric);
    
    // wrap the pointer into smart pointer
    return std::shared_ptr<LUft<LU_int_t,Complex>>(lu_ptr);
}

template<>
CooMatrix<LU_int_t,Complex> CsrMatrix<LU_int_t,Complex>::tocoo () const
{
    // reserve space for the auxiliary (jj) and the output (Ti,Tj,Tx,Tz) arrays
    std::size_t N = x_.size();
    NumberArray<LU_int_t> Ti(N), Tj(N), jj(N);
    cArray Tx(N);
    
    // do we have any elements at all?
    if (N != 0)
    {
    
        // do the conversion
        LU_int_t status = UMFPACK_COL_TO_TRIPLET_F(n_, p_.data(), jj.data());
        
        // check success
        if (status != 0)
        {
            std::cerr << "\n[CscMatrix::tocoo] Exit code " << status << std::endl;
            UMFPACK_REPORT_STATUS_F(0, status);
        }
        
        // copy only non-zero entries to output arrays
        std::size_t nz = 0;
        for (std::size_t i = 0; i < N; i++)
        {
            if (x_[i] != 0.)
            {
                Ti[nz] = i_[i];
                Tj[nz] = jj[i];
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
    return CooMatrix<LU_int_t,Complex> (m_, n_, Tj, Ti, Tx);
}

template<>
CsrMatrix<LU_int_t,Complex> CooMatrix<LU_int_t,Complex>::tocsr () const
{
    std::size_t nz = x_.size();
    
    // CSC matrix data
    NumberArray<LU_int_t> Ap(n_ + 1), Ai(nz);
    cArray Ax(nz);
    
    // do we have any elements at all?
    if (nz != 0)
    {
    
        LU_int_t status = UMFPACK_TRIPLET_TO_COL_F
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
            UMFPACK_REPORT_STATUS_F(0, status);
        }
        
        // crop storage
        std::size_t N = Ap[m_];
        Ai.resize(N);
        Ax.resize(N);
    }
    
    return CsrMatrix<LU_int_t,Complex> (m_, n_, Ap, Ai, Ax);
}

#endif // WITH_UMFPACK


// -------------------------------------------------------------------------------------
// SUPERLU-dependent functions.
//

#ifdef WITH_SUPERLU
#include <slu_zdefs.h>

template<>
std::shared_ptr<LUft<int,Complex>> CsrMatrix<int,Complex>::factorize_superlu (double droptol, void * data) const
{
    //
    // Create matrix of the system.
    //
    
        NRformat AStore;
        AStore.nnz    = this->x().size();   // number of non-zero elements
        AStore.nzval  = const_cast<Complex*>(this->x().data()); // pointer to the array of non-zero elements
        AStore.colind = const_cast<int*>(this->i().data());     // row indices
        AStore.rowptr = const_cast<int*>(this->p().data());     // column pointers
        
        SuperMatrix A;
        A.Stype = SLU_NR;       // storage type: compressed sparse, row-major
        A.Dtype = SLU_Z;        // data type: double complex
        A.Mtype = SLU_GE;       // mathematical type: general
        A.nrow  = this->rows(); // number of rows
        A.ncol  = this->cols(); // number of columns
        A.Store = &AStore;      // data structure pointer
    
    //
    // Create the (empty) right hand side and solution matrix.
    //
    
        DNformat BStore, XStore;
        SuperMatrix L, U, B, X;
        B.ncol = X.ncol = 0;
        B.Store = &BStore;
        X.Store = &XStore;
    
    //
    // Prepare SuperLU environment.
    //
    
        // calculation options
        superlu_options_t options;
        set_default_options(&options);
        options.ColPerm = MMD_AT_PLUS_A;
        options.SymPattern = YES;
        options.ILU_DropRule = DROP_BASIC;
        options.ILU_DropTol = droptol;
        
        // calculation diagnostic information
        SuperLUStat_t stat;
        StatInit(&stat);
        
        // memory usage
        mem_usage_t mem_usage;
    
    //
    // Compute the factorization.
    //
        
        // permutation arrays, elimination tree
        iArray perm_r(A.nrow), perm_c(A.ncol), etree(A.nrow);
        
        // row and column scale factors, reciprocal condition number, reciprocal pivot growth factor
        rArray R(A.nrow), C(A.ncol); double rcond, rpg;
        
        // forward and backward errors (one element per one right hand side)
        double ferr, berr;
        
        // equilibration done
        char equed;
        
        // status indicator
        int info;
        
        // LU factorization
        zgssvx
        (
            &options,       // calculation options
            &A,             // matrix data structure
            &perm_c[0],     // column permutation
            &perm_r[0],     // row permutation
            &etree[0],      // elimination tree
            &equed,         // equilibration done
            &R[0],          // row scale factors
            &C[0],          // column scale factors
            &L,             // L-factor
            &U,             // U-factor
            nullptr,        // workspace (not used)
            0,              // size of the workspace (0 = automatic allocation)
            &B,             // right-hand sides (empty)
            &X,             // solution matrix (empty)
            &rpg,           // reciprocal pivot growth factor
            &rcond,         // reciprocal condition number
            &ferr,          // forward error
            &berr,          // backward error
            &mem_usage,     // memory usage
            &stat,          // diagnostic infomation
            &info           // result status
        );
        
        if (info < 0)
            HexException("SuperLU/zgssvx: Parameter %d has illegal value.", -info);
        if (info > A.ncol + 1)
            HexException("SuperLU/zgssvx: Memory allocation failure after %d bytes.", info);
        if (info == A.ncol + 1)
            HexException("SuperLU/zgssvx: Badly conditioned system.");
        if (info > 0)
            HexException("SuperLU/zgssvx: Singular factor.");
        
    // create a new LU factorization container
    LUft<int,Complex> * lu_ptr = new LUft_SUPERLU<int,Complex>
    (
        this, perm_c, perm_r, etree, equed,
        R, C, L, U, mem_usage.for_lu, droptol
    );
    
    // wrap the pointer into smart pointer
    return std::shared_ptr<LUft<int,Complex>>(lu_ptr);
}

#endif // WITH_SUPERLU


// -------------------------------------------------------------------------------------
// SUPERLU-DIST-dependent functions.
//

#ifdef WITH_SUPERLU_DIST
#include <superlu_zdefs.h>

template<>
std::shared_ptr<LUft<LU_int_t,Complex>> CsrMatrix<LU_int_t,Complex>::factorize_superlu_dist (double droptol, void * data) const
{
    //
    // Create matrix of the system.
    //

        cArray xdata (this->x());
        NumberArray<LU_int_t> idata (this->i());
        NumberArray<LU_int_t> pdata (this->p());
        
        NCformat AStore;
        AStore.nnz    = xdata.size();                           // number of non-zero elements
        AStore.nzval  = &xdata[0];                              // pointer to the array of non-zero elements
        AStore.rowind = reinterpret_cast<int_t*>(&idata[0]);    // row indices
        AStore.colptr = reinterpret_cast<int_t*>(&pdata[0]);    // column pointers
        
        SuperMatrix A;
        A.Stype = SLU_NC;       // storage type: compressed sparse, column-major (SuperLU-dist suports no other)
        A.Dtype = SLU_Z;        // data type: double complex
        A.Mtype = SLU_GE;       // mathematical type: general
        A.nrow  = this->rows(); // number of rows
        A.ncol  = this->cols(); // number of columns
        A.Store = &AStore;      // data structure pointer

    //
    // Prepare SuperLU environment.
    //
    
        // get process grid
        gridinfo_t * grid = (gridinfo_t *)data;
        
        // calculation options
        superlu_options_t options;
        set_default_options_dist(&options);
        options.ParSymbFact = YES;
        options.ColPerm = METIS_AT_PLUS_A;
        options.PrintStat = NO;
        options.SymPattern = YES;
//         options.ILU_DropRule = DROP_BASIC;
//         options.ILU_DropTol = droptol;
        
        // distributed scale and permutation data
        ScalePermstruct_t ScalePermstruct;
        ScalePermstructInit(A.nrow, A.ncol, &ScalePermstruct);
        
        // distributed factorization data
        LUstruct_t LUstruct;
        LUstructInit(A.nrow, &LUstruct);
        
        // calculation diagnostic information
        SuperLUStat_t stat;
        PStatInit(&stat);
    
    //
    // Compute the factorization.
    //
    
        // backward error (one element per one right hand side)
        double berr[1];
        
        // status indicator
        int info;
        
        // LU factorization
        pzgssvx_ABglobal
        (
            &options,           // calculation options
            &A,                 // matrix to factorize
            &ScalePermstruct,   // scaling and permutation data
            nullptr,            // right hand sides (not used)
            A.nrow,             // right hand sides leading dimension
            0,                  // number of right hand sides (none, only factorizing)
            grid,               // process grid
            &LUstruct,          // factorization data
            berr,               // backward error
            &stat,              // diagnostic information
            &info               // result status
        );
        
        mem_usage_t mem_usage;
        zQuerySpace_dist(A.nrow, &LUstruct, grid, &stat, &mem_usage);
        
        // TODO : Reduce memory usage over the group.
        
        if (info > A.ncol)
            HexException("SuperLU/zgssvx: Memory allocation failure after %d bytes.", info);
        if (info > 0)
            HexException("SuperLU/zgssvx: Singular factor.");
    
    // create a new LU factorization container
    LUft<LU_int_t,Complex> * lu_ptr = new LUft_SUPERLU_DIST<LU_int_t,Complex>
    (
        this, ScalePermstruct, LUstruct, grid, mem_usage.for_lu
    );
    
    // wrap the pointer into smart pointer
    return std::shared_ptr<LUft<LU_int_t,Complex>>(lu_ptr);
}

#endif // WITH_SUPERLU_DIST
