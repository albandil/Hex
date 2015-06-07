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
// UMFPACK-dependent functions (LP64)
//

#ifdef WITH_UMFPACK

template<>
std::size_t LUft_UMFPACK<int,Complex>::size () const
{
    if (numeric_ == nullptr)
        return 0;
    
    int lnz, unz, m, n, nz_udiag;
    int status = umfpack_zi_get_lunz
    (
        &lnz,       // number of non-zero elements in L-factor
        &unz,       // number of non-zero elements in U-factor
        &m,         // row count
        &n,         // column count
        &nz_udiag,  // ?
        numeric_    // factorization object
    );
    
    return status == 0 ? (lnz + unz) * 16 : 0; // Byte count
}

template<>
void LUft_UMFPACK<int,Complex>::solve (const cArrayView b, cArrayView x, int eqs) const
{
    // check sizes
    assert(eqs * matrix_->n_ == (int)x.size());
    assert(eqs * matrix_->n_ == (int)b.size());
    
    // solve for all RHSs
    for (int eq = 0; eq < eqs; eq++)
    {
        // solve for current RHS
        int status = umfpack_zi_solve
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
            umfpack_zi_report_status(0, status);
        }
    }
}

template<>
void LUft_UMFPACK<int,Complex>::save (std::string name) const
{
    int err = umfpack_zi_save_numeric(numeric_, const_cast<char*>(name.c_str()));
    
    if (err == UMFPACK_ERROR_invalid_Numeric_object)
        HexException("[LUft::save] Invalid numeric object.");
    
    if (err == UMFPACK_ERROR_file_IO)
        HexException("[LUft::save] Failed to save LU object \"%s\" (size = %ld).", name.c_str(), size());
}

template<>
void LUft_UMFPACK<int,Complex>::load (std::string name, bool throw_on_io_failure)
{
    int err = umfpack_zi_load_numeric(&numeric_, const_cast<char*>(name.c_str()));
    
    if (err == UMFPACK_ERROR_out_of_memory)
        HexException("[LUft::load] Out of memory.");
    
    if (err == UMFPACK_ERROR_file_IO and throw_on_io_failure)
        HexException("[LUft::save] Failed to load LU object \"%s\".", name.c_str());
}

template<>
void LUft_UMFPACK<int,Complex>::drop ()
{
    if (numeric_ != nullptr)
    {
        umfpack_zi_free_numeric(&numeric_);
        numeric_ = nullptr;
    }
}

#endif // WITH_UMFPACK


// ------------------------------------------------------------------------------------
// UMFPACK-dependent functions (ILP64)
//

#ifdef WITH_UMFPACK

template<>
std::size_t LUft_UMFPACK<std::int64_t,Complex>::size () const
{
    if (numeric_ == nullptr)
        return 0;
    
    std::int64_t lnz, unz, m, n, nz_udiag;
    std::int64_t status = umfpack_zl_get_lunz
    (
        &lnz,       // number of non-zero elements in L-factor
        &unz,       // number of non-zero elements in U-factor
        &m,         // row count
        &n,         // column count
        &nz_udiag,  // ?
        numeric_    // factorization object
    );
    
    return status == 0 ? (lnz + unz) * 16 : 0; // Byte count
}

template<>
void LUft_UMFPACK<std::int64_t,Complex>::solve (const cArrayView b, cArrayView x, int eqs) const
{
    // check sizes
    assert(eqs * matrix_->n_ == (int)x.size());
    assert(eqs * matrix_->n_ == (int)b.size());
    
    // solve for all RHSs
    for (int eq = 0; eq < eqs; eq++)
    {
        // solve for current RHS
        std::int64_t status = umfpack_zl_solve
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
            umfpack_zl_report_status(0, status);
        }
    }
}

template<>
void LUft_UMFPACK<std::int64_t,Complex>::save (std::string name) const
{
    std::int64_t err = umfpack_zl_save_numeric(numeric_, const_cast<char*>(name.c_str()));
    
    if (err == UMFPACK_ERROR_invalid_Numeric_object)
        HexException("[LUft::save] Invalid numeric object.");
    
    if (err == UMFPACK_ERROR_file_IO)
        HexException("[LUft::save] Failed to save LU object \"%s\" (size = %ld).", name.c_str(), size());
}

template<>
void LUft_UMFPACK<std::int64_t,Complex>::load (std::string name, bool throw_on_io_failure)
{
    std::int64_t err = umfpack_zl_load_numeric(&numeric_, const_cast<char*>(name.c_str()));
    
    if (err == UMFPACK_ERROR_out_of_memory)
        HexException("[LUft::load] Out of memory.");
    
    if (err == UMFPACK_ERROR_file_IO and throw_on_io_failure)
        HexException("[LUft::save] Failed to load LU object \"%s\".", name.c_str());
}

template<>
void LUft_UMFPACK<std::int64_t,Complex>::drop ()
{
    if (numeric_ != nullptr)
    {
        umfpack_zl_free_numeric(&numeric_);
        numeric_ = nullptr;
    }
}

#endif // WITH_UMFPACK


// ------------------------------------------------------------------------------------
// SUPERLU-dependent functions (LP64).
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
