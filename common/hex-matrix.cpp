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

#include "hex-arrays.h"
#include "hex-hdffile.h"
#include "hex-densematrix.h"
#include "hex-luft.h"
#include "hex-misc.h"
#include "hex-symbandmatrix.h"

#ifdef SINGLE
extern "C" void cgetrf_ (int*, int*, Complex*, int*, int*, int*);
extern "C" void cgetri_ (int*, Complex*, int*, int*, Complex*, int*, int*);
extern "C" void cgeev_ (char*, char*, int*, Complex*, int*, Complex*, Complex*, int*, Complex*, int*, Complex*, int*, Real*, int*);
#else
extern "C" void zgetrf_ (int*, int*, Complex*, int*, int*, int*);
extern "C" void zgetri_ (int*, Complex*, int*, int*, Complex*, int*, int*);
extern "C" void zgeev_ (char*, char*, int*, Complex*, int*, Complex*, Complex*, int*, Complex*, int*, Complex*, int*, Real*, int*);
#endif

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
#ifdef SINGLE
    cgetrf_
#else
    zgetrf_
#endif
    (
        &N,             // number of rows in the matrix
        &N,             // number of columns of the matrix
        inv.begin(),    // matrix data (on return the LU factors)
        &N,             // leading dimension of A
        &IPIV[0],       // pivot indices
        &INFO           // diagnostic information
    );
    
    // query for optimal work size for inversion
#ifdef SINGLE
    cgetri_
#else
    zgetri_
#endif
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
#ifdef SINGLE
    cgetri_
#else
    zgetri_
#endif
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
#ifdef SINGLE
    cgeev_
#else
    zgeev_
#endif
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
#ifdef SINGLE
    cgeev_
#else
    zgeev_
#endif
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

// -------------------------------------------------------------------------------------
// Sparse matrix routines.
//

template<>
void SymBandMatrix<Complex>::sym_band_dot (int n, int d, const cArrayView M, Complex alpha, const cArrayView A, Complex beta, cArrayView B, MatrixSelection::Selection tri)
{
    // check dimensions
    if (A.size() % n != 0)
        HexException("Incompatible dimensions: %d (mat) × %ld (vec). You are probably mixing radial data for different grids.", n, A.size());
    
    // get number of source vectors
    int Nvec = A.size() / n;
    
    // skip if there are no data
    if (n == 0 or d == 0)
        return;
    
    // matrix data pointer
    Complex const * const restrict pM = M.data();
    
    // source and destination vector pointers
    Complex const * restrict pA = A.data();
    Complex       * restrict pB = B.data();
    
    // for all source vectors (columns)
    for (int j = 0; j < Nvec; j++)
    {
        // scale the destination vector
        for (int k = 0; k < n; k++)
            pB[k] *= beta;
        
        // for all rows/cols of the matrix (and of the target vector)
        for (int k = 0; k < n; k++)
        {
            // Diagonal element
            
                if (tri & MatrixSelection::Diagonal)
                {
                    pB[k] += pM[k * d] * alpha * pA[k];
                }
            
            // Direct pass (~ strict upper triangle)
        
                if (tri & MatrixSelection::StrictUpper)
                {
                    for (int l = 1; l < d and k + l < n; l++)
                        pB[k] += pM[k * d + l] * alpha * pA[k + l];
                }
        
            // Mirror pass (~ strict lower triangle)
                
                if (tri & MatrixSelection::StrictLower)
                {
                    for (int l = 1; l < d and k + l < n; l++)
                        pB[k + l] += pM[k * d + l] * alpha * pA[k];
                }
        }
        
        // move on to the next source vector
        pA += n;
        pB += n;
    }
    
    return;
}

template<>
CsrMatrix<LU_int_t,Complex> CooMatrix<LU_int_t,Complex>::tocsr () const
{
    // get number of structurally non-zero elements
    LU_int_t nz = x_.size();
    
    // get row lengths
    std::vector<LU_int_t> len (m_, 0);
    for (LU_int_t n = 0; n < nz; n++) if (x_[n] != 0.0_r)
        len[i_[n]]++;
    
    // create element pointer array for each matrix row
    std::vector<std::vector<LU_int_t>> elem_ptrs (m_);
    
    // reserve memory for the row data
    for (LU_int_t n = 0; n < m_; n++)
        elem_ptrs[n].reserve(len[n]);
    
    // store index of each element into the appropriate row
    for (LU_int_t n = 0; n < nz; n++) if (x_[n] != 0.0_r)
        elem_ptrs[i_[n]].push_back(n);
    
    // sort element pointers by their column index
    for (LU_int_t n = 0; n < m_; n++)
        std::sort(elem_ptrs[n].begin(), elem_ptrs[n].end(), [&](LU_int_t a, LU_int_t b) { return j_[a] < j_[b]; });
    
    // allocate output arrays
    NumberArray<LU_int_t> Ap(m_ + 1), Ai(nz);
    cArray Ax(nz);
    
    // copy & sum entries
    LU_int_t pos = 0;
    for (LU_int_t m = 0; m < m_; m++)
    {
        // for all non-zero elements
        for (LU_int_t n = 0; n < (LU_int_t)elem_ptrs[m].size(); n++)
        {
            // get global index
            LU_int_t gid = elem_ptrs[m][n];
            
            // update elements
            Ai[pos]  = j_[gid];
            Ax[pos] += x_[gid];
            
            // if the next row entry is in the same column, do not advance position counter,
            // so that the element will be added to current position
            if (n < (LU_int_t)elem_ptrs[m].size() - 1 and j_[gid] == j_[elem_ptrs[m][n + 1]])
                continue;
            
            // otherwise move on to the next position
            pos++;
        }
        
        // insert next row pointer
        Ap[m + 1] = pos;
    }
    
    return CsrMatrix<LU_int_t,Complex> (m_, n_, Ap, Ai, Ax);
}


// -------------------------------------------------------------------------------------
// UMFPACK-dependent functions.
//

#ifdef WITH_UMFPACK
#include <umfpack.h>

template<>
std::shared_ptr<LUft<LU_int_t,Complex>> CsrMatrix<LU_int_t,Complex>::factorize_umfpack (Real droptol, void * data) const
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
    
    // diagnostic information
    rArray Info (UMFPACK_INFO);
    
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
        std::cerr << "\nSymbolic factorization error " << status << std::endl;
        UMFPACK_REPORT_STATUS_F(0, status);
        std::exit(EXIT_FAILURE);
    }
    
    // do some factorizations
    status = UMFPACK_NUMERIC_F
    (
        p_.data(), i_.data(),    // column and row indices
        reinterpret_cast<const double*>(x_.data()), 0,    // matrix data
        Symbolic, &Numeric, Control, &Info[0]    // UMFPACK internals
    );
    if (status != 0)
    {
        std::cerr << "\nNumeric factorization error " << status << std::endl;
        UMFPACK_REPORT_STATUS_F(0, status);
        std::exit(EXIT_FAILURE);
    }
    
    // release unused data
    UMFPACK_FREE_SYMBOLIC_F(&Symbolic);
    
    // create a new LU factorization container
    LUft_UMFPACK<LU_int_t,Complex> * lu_ptr = new LUft_UMFPACK<LU_int_t,Complex>(*this, Numeric);
    lu_ptr->info_ = Info;
    
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
            if (x_[i] != 0.0_r)
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

#endif // WITH_UMFPACK


// -------------------------------------------------------------------------------------
// SUPERLU-dependent functions.
//

#ifdef WITH_SUPERLU
#ifdef SINGLE
    #include <slu_cdefs.h>
#else
    #include <slu_zdefs.h>
#endif

template<>
std::shared_ptr<LUft<int,Complex>> CsrMatrix<int,Complex>::factorize_superlu (Real droptol, void * data) const
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
#ifdef SINGLE
        A.Dtype = SLU_C;        // data type: single complex
#else
        A.Dtype = SLU_Z;        // data type: double complex
#endif
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
        
        // reusable information
        GlobalLU_t Glu;
        
        // memory usage
        mem_usage_t mem_usage;
    
    //
    // Compute the factorization.
    //
        
        // permutation arrays, elimination tree
        iArray perm_r(A.nrow), perm_c(A.ncol), etree(A.nrow);
        
        // row and column scale factors, reciprocal condition number, reciprocal pivot growth factor
        rArray R(A.nrow), C(A.ncol); Real rcond, rpg;
        
        // forward and backward errors (one element per one right hand side)
        Real ferr, berr;
        
        // equilibration done
        char equed;
        
        // status indicator
        int info;
        
        // LU factorization
#ifdef SINGLE
        cgssvx
#else
        zgssvx
#endif
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
//             &Glu,           // reusable information
            &mem_usage,     // memory usage
            &stat,          // diagnostic infomation
            &info           // result status
        );
        
        if (info < 0)
            HexException("SuperLU/?gssvx: Parameter %d has illegal value.", -info);
        if (info > A.ncol + 1)
            HexException("SuperLU/?gssvx: Memory allocation failure after %d bytes.", info);
        if (info == A.ncol + 1)
            HexException("SuperLU/?gssvx: Badly conditioned system.");
        if (info > 0)
            HexException("SuperLU/?gssvx: Singular factor.");
        
    // create a new LU factorization container
    LUft<int,Complex> * lu_ptr = new LUft_SUPERLU<int,Complex>
    (
        *this, perm_c, perm_r, etree, equed,
        R, C, L, U, Glu, mem_usage.for_lu, droptol
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
std::shared_ptr<LUft<LU_int_t,Complex>> CsrMatrix<LU_int_t,Complex>::factorize_superlu_dist (Real droptol, void * data) const
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
#ifdef SINGLE
        A.Dtype = SLU_C;        // data type: single complex
#else
        A.Dtype = SLU_Z;        // data type: double complex
#endif
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
        *this, ScalePermstruct, LUstruct, grid, mem_usage.for_lu
    );
    
    // wrap the pointer into smart pointer
    return std::shared_ptr<LUft<LU_int_t,Complex>>(lu_ptr);
}

#endif // WITH_SUPERLU_DIST

// -------------------------------------------------------------------------------------
// MUMPS-dependent functions.
//

#ifdef WITH_MUMPS

template<>
std::shared_ptr<LUft<LU_int_t,Complex>> CsrMatrix<LU_int_t,Complex>::factorize_mumps (Real droptol, void * data) const
{
    //
    // Extract parameters.
    //
    
        bool out_of_core = *((MUMPS_INT*)data);
        MUMPS_INT verbosity_level = *((MUMPS_INT*)data + 1);
        MUMPS_INT comm = *((MUMPS_INT*)data + 2);
    
    //
    // Create matrix of the system (i.e. the IJV triplet).
    //

        // estimate non-zero element count of the upper triangle
        LU_int_t nz = (i_.size() + n_) / 2;
        
        // data arrays
        NumberArray<MUMPS_INT> I, J;
        cArray A;
        
        // allocate memory
        I.reserve(nz);
        J.reserve(nz);
        A.reserve(nz);
        
        // for all rows
        for (LU_int_t row = 0, nz = 0; row < m_; row++)
        {
            // for all columns with structurally non-zero entries
            for (LU_int_t idx = p_[row]; idx < p_[row + 1]; idx++)
            {
                // get column index
                LU_int_t col = i_[idx];
                
                // only consider upper triangle (and the main diagonal)
                if (row <= col and x_[idx] != 0.0_z)
                {
                    // insert the element
                    I.push_back(row + 1);
                    J.push_back(col + 1);
                    A.push_back(x_[idx]);
                    
                    // update true element count
                    nz++;
                }
            }
        }
    
    //
    // Prepare MUMPS environment.
    //
    
        MUMPS_STRUC_C * settings = new MUMPS_STRUC_C;
        
        // initialize
        settings->job = MUMPS_INITIALIZE;
        settings->sym = 2;
        settings->par = 1;
#ifdef WITH_MPI
        settings->comm_fortran = comm;
#endif
        MUMPS_C(settings);
        
        // analyze
        settings->job = MUMPS_ANALYZE;
        settings->ICNTL(1) = (verbosity_level == 0 ? 0 : 6); // errors to STDOUT (default: 6)
        settings->ICNTL(2) = 0; // diagnostics to /dev/null
        settings->ICNTL(3) = (verbosity_level == 0 ? 0 : 6); // global info to STDOUT (default: 6)
        settings->ICNTL(4) = verbosity_level; // verbosity level (default: 2)
        settings->ICNTL(5) = 0; // COO format
        settings->ICNTL(22) = out_of_core; // OOC factorization
        std::strcpy(settings->ooc_tmpdir, ".");
        std::strcpy(settings->ooc_prefix, "ooc_");
        settings->n = this->n_;
        settings->nz = nz;
        settings->irn = I.data();
        settings->jcn = J.data();
        settings->a = reinterpret_cast<MUMPS_COMPLEX*>(A.data());
        MUMPS_C(settings);
    
    //
    // Compute the factorization.
    //
    
        settings->job = MUMPS_FACTORIZE;
        MUMPS_C(settings);
    
    //
    // Create a new LU factorization container
    //
    
        LUft<LU_int_t,Complex> * lu_ptr = new LUft_MUMPS<LU_int_t,Complex>
        (
            settings,
            std::move(I),
            std::move(J),
            std::move(A)
        );
    
    // wrap the pointer into smart pointer
    return std::shared_ptr<LUft<LU_int_t,Complex>>(lu_ptr);
}

#endif // WITH_MUMPS
