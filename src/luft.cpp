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

#include "luft.h"

// ------------------------------------------------------------------------------------
// UMFPACK-dependent functions.
//

#ifdef WITH_UMFPACK

template<>
std::size_t LUft_UMFPACK<LU_int_t,Complex>::size () const
{
    if (numeric_ == nullptr)
        return 0;
    
    LU_int_t lnz, unz, m, n, nz_udiag;
    LU_int_t status = UMFPACK_GET_LUNZ_F
    (
        &lnz,       // number of non-zero elements in L-factor
        &unz,       // number of non-zero elements in U-factor
        &m,         // row count
        &n,         // column count
        &nz_udiag,  // ?
        numeric_    // factorization object
    );
    
    return status == 0 ? (lnz + unz) * 16LL : 0LL; // Byte count
}

template<>
void LUft_UMFPACK<LU_int_t,Complex>::solve (const cArrayView b, cArrayView x, int eqs) const
{
    // check sizes
    assert(eqs * matrix_->n_ == (int)x.size());
    assert(eqs * matrix_->n_ == (int)b.size());
    
    // solve for all RHSs
    for (int eq = 0; eq < eqs; eq++)
    {
        // solve for current RHS
        int status = UMFPACK_SOLVE_F
        (
            UMFPACK_Aat,            // matrix orientation
            matrix_->p().data(),    // row pointers
            matrix_->i().data(),    // column indices
            
            // matrix elements (interleaved)
            reinterpret_cast<const double*>(matrix_->x().data()),
            nullptr,
            
            // solutions (interleaved)
            reinterpret_cast<double*>(&x[0] + eq * matrix_->rows()),
            nullptr,
            
            // right-hand side vectors (interleaved)
            reinterpret_cast<const double*>(&b[0] + eq * matrix_->rows()),
            nullptr,
            
            numeric_,   // factorization object
            nullptr,    // ?
            &info_[0]   // diagnostic information
        );
        
        // check output
        if (status != UMFPACK_OK)
        {
            std::cerr << "\n[CsrMatrix::LUft::solve] Exit status " << status << std::endl;
            UMFPACK_REPORT_STATUS_F(0, status);
        }
    }
}

template<>
void LUft_UMFPACK<LU_int_t,Complex>::save (std::string name) const
{
    int err = UMFPACK_SAVE_NUMERIC_F(numeric_, const_cast<char*>(name.c_str()));
    
    if (err == UMFPACK_ERROR_invalid_Numeric_object)
        HexException("[LUft::save] Invalid numeric object.");
    
    if (err == UMFPACK_ERROR_file_IO)
        HexException("[LUft::save] Failed to save LU object \"%s\" (size = %ld).", name.c_str(), size());
}

template<>
void LUft_UMFPACK<LU_int_t,Complex>::load (std::string name, bool throw_on_io_failure)
{
    int err = UMFPACK_LOAD_NUMERIC_F(&numeric_, const_cast<char*>(name.c_str()));
    
    if (err == UMFPACK_ERROR_out_of_memory)
        HexException("[LUft::load] Out of memory.");
    
    if (err == UMFPACK_ERROR_file_IO and throw_on_io_failure)
        HexException("[LUft::save] Failed to load LU object \"%s\".", name.c_str());
}

template<>
void LUft_UMFPACK<LU_int_t,Complex>::drop ()
{
    if (numeric_ != nullptr)
    {
        UMFPACK_FREE_NUMERIC_F(&numeric_);
        numeric_ = nullptr;
    }
}

#endif // WITH_UMFPACK

// ------------------------------------------------------------------------------------
// SUPERLU-dependent functions.
//

#ifdef WITH_SUPERLU
#include <slu_zdefs.h>

template<>
void LUft_SUPERLU<int,Complex>::solve (const cArrayView b, cArrayView x, int eqs) const
{
    //
    // Create matrix of the system.
    //
    
        NRformat AStore;
        AStore.nnz    = matrix_->x().size();                        // number of non-zero elements
        AStore.nzval  = const_cast<Complex*>(matrix_->x().data());  // pointer to the array of non-zero elements
        AStore.colind = const_cast<int*>(matrix_->i().data());      // row indices
        AStore.rowptr = const_cast<int*>(matrix_->p().data());      // column pointers
        
        SuperMatrix A;
        A.Stype = SLU_NR;           // storage type: compressed sparse, row-major
        A.Dtype = SLU_Z;            // data type: double complex
        A.Mtype = SLU_GE;           // mathematical type: general
        A.nrow  = matrix_->rows();  // number of rows
        A.ncol  = matrix_->cols();  // number of columns
        A.Store = &AStore;          // data structure pointer
    
    //
    // Create the right hand side.
    //
    
        DNformat BStore;
        BStore.lda = matrix_->cols();                   // leading dimension
        BStore.nzval = const_cast<Complex*>(b.data());  // data pointer
        
        SuperMatrix B;
        B.Stype = SLU_DN;               // storage type: Fortran dense matrix
        B.Dtype = SLU_Z;                // data type: double
        B.Mtype = SLU_GE;               // mathematical type: general
        B.nrow = matrix_->cols();       // number of rows (= number of equations)
        B.ncol = eqs;                   // number of columns (= number of right hand sides)
        B.Store = &BStore;              // data structure pointer
        
    //
    // Create the solution matrix.
    //
        
        DNformat XStore;
        XStore.lda = matrix_->cols();   // leading dimension
        XStore.nzval = x.data();        // data pointer
        
        SuperMatrix X;
        X.Stype = SLU_DN;               // storage type: Fortran dense matrix
        X.Dtype = SLU_Z;                // data type: double
        X.Mtype = SLU_GE;               // mathematical type: general
        X.nrow = matrix_->cols();       // number of rows (= number of equations)
        X.ncol = eqs;                   // number of columns (= number of right hand sides)
        X.Store = &XStore;              // data structure pointer
        
    //
    // Prepare SuperLU environment.
    //
    
        // calculation options
        superlu_options_t options;
        set_default_options(&options);
        options.ColPerm = MMD_AT_PLUS_A;
        options.SymPattern = YES;
        options.ILU_DropRule = DROP_BASIC;
        options.ILU_DropTol = droptol_;
        options.Fact = FACTORED;
        
        // calculation diagnostic information
        SuperLUStat_t stat;
        StatInit(&stat);
        
        // memory usage
        mem_usage_t mem_usage;
    
    //
    // Solve the right hand sides.
    //
        
        // forward nad backward errors (one element per one right hand side)
        rArray ferr(eqs), berr(eqs);
        
        // reciprocal condition number, reciprocal pivot growth factor
        double rcond, rpg;
        
        // back-substitution
        int info;
        zgssvx
        (
            &options,                       // calculation options
            &A,                             // matrix data structure
            const_cast<int*>(&perm_c_[0]),  // column permutation
            const_cast<int*>(&perm_r_[0]),  // row permutation
            const_cast<int*>(&etree_[0]),   // elimination tree
            const_cast<char*>(&equed_),     // equilibration done
            const_cast<double*>(&R_[0]),    // row scale factors
            const_cast<double*>(&C_[0]),    // column scale factors
            const_cast<SuperMatrix*>(&L_),  // L-factor
            const_cast<SuperMatrix*>(&U_),  // U-factor
            nullptr,                        // workspace (not used)
            0,              // size of the workspace (0 = dynamic allocation)
            &B,             // right-hand sides
            &X,             // solution matrix
            &rpg,           // reciprocal pivot growth factor
            &rcond,         // reciprocal condition number
            &ferr[0],       // forward error
            &berr[0],       // backward error
            &mem_usage,     // memory usage
            &stat,          // diagnostic infomation
            &info           // result status
        );
    
    // check status
    if (info != 0)
        HexException("SuperLU/zgssvx failed with status %d.", info);
}

#endif // WITH_SUPERLU

// ------------------------------------------------------------------------------------
// SUPERLU-DIST-dependent functions.
//

#ifdef WITH_SUPERLU_DIST
#include <superlu_zdefs.h>

template<>
void LUft_SUPERLU_DIST<LU_int_t,Complex>::solve (const cArrayView b, cArrayView x, int eqs) const
{
    //
    // Create matrix of the system.
    //
    
        cArray xdata (matrix_->x());
        NumberArray<LU_int_t> idata (matrix_->i());
        NumberArray<LU_int_t> pdata (matrix_->p());
        
        NCformat AStore;
        AStore.nnz    = xdata.size();                           // number of non-zero elements
        AStore.nzval  = &xdata[0];                              // pointer to the array of non-zero elements
        AStore.rowind = reinterpret_cast<int_t*>(&idata[0]);    // row indices
        AStore.colptr = reinterpret_cast<int_t*>(&pdata[0]);    // column pointers
        
        SuperMatrix A;
        A.Stype = SLU_NC;           // storage type: compressed sparse, column-major (SuperLU-dist suports no other)
        A.Dtype = SLU_Z;            // data type: double complex
        A.Mtype = SLU_GE;           // mathematical type: general
        A.nrow  = matrix_->rows();  // number of rows
        A.ncol  = matrix_->cols();  // number of columns
        A.Store = &AStore;          // data structure pointer
    
    //
    // Prepare SuperLU environment.
    //
    
        // calculation options
        superlu_options_t options;
        set_default_options_dist(&options);
        options.ParSymbFact = YES;
        options.ColPerm = METIS_AT_PLUS_A;
        options.PrintStat = NO;
        options.SymPattern = YES;
        options.Fact = FACTORED;
        options.IterRefine = NOREFINE;
//         options.ILU_DropRule = DROP_BASIC;
//         options.ILU_DropTol = droptol;
        
        // calculation diagnostic information
        SuperLUStat_t stat;
        PStatInit(&stat);
    
    //
    // Solve the system.
    //
        
        // backward error (one element per one right hand side)
        double berr[eqs];
        
        // status indicator
        int info;
        
        // LU factorization
        x = b;
        pzgssvx_ABglobal
        (
            &options,                                           // calculation options
            &A,                                                 // matrix to factorize
            const_cast<ScalePermstruct_t*>(&ScalePermstruct_),  // scaling and permutation data
            reinterpret_cast<doublecomplex*>(&x[0]),            // right hand sides (assume packed)
            A.nrow,                                             // right hand sides leading dimension
            eqs,                                                // number of right hand sides
            grid_,                                              // process grid
            const_cast<LUstruct_t*>(&LUstruct_),                // factorization data
            berr,                                               // backward error
            &stat,                                              // diagnostic information
            &info                                               // result status
        );
        
        if (info > A.ncol)
            HexException("SuperLU/zgssvx: Memory allocation failure after %d bytes.", info);
        if (info > 0)
            HexException("SuperLU/zgssvx: Singular factor.");
    
    //
    // Clean up.
    //
    
        PStatFree(&stat);
}

#endif // WITH_SUPERLU_DIST
