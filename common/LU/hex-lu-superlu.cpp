//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2016, Jakub Benda, Charles University in Prague                    //
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

#ifdef WITH_SUPERLU

// --------------------------------------------------------------------------------- //

#include "hex-lu-superlu.h"

// --------------------------------------------------------------------------------- //

template<>
LUft_SUPERLU<LU_int_t,Complex>::LUft_SUPERLU ()
    : LUft<LU_int_t,Complex>(), size_(0)
{
    // nothing
}

template<>
void LUft_SUPERLU<LU_int_t,Complex>::factorize (CsrMatrix<LU_int_t,Complex> const & matrix, LUftData data)
{
    //
    // Create matrix of the system.
    //

        P_.resize(matrix.p().size());
        for (std::size_t i = 0; i < P_.size(); i++)
            P_[i] = matrix.p()[i];
        
        I_.resize(matrix.i().size());
        for (std::size_t i = 0; i < I_.size(); i++)
            I_[i] = matrix.i()[i];
        
        X_.resize(matrix.x().size());
        for (std::size_t i = 0; i < X_.size(); i++)
            X_[i] = Complex(matrix.x()[i].real(), matrix.x()[i].imag());
        
        NRformat AStore;
        AStore.nnz    = X_.size();  // number of non-zero elements
        AStore.nzval  = X_.data();  // pointer to the array of non-zero elements
        AStore.colind = I_.data();  // row indices
        AStore.rowptr = P_.data();  // column pointers
        
        SuperMatrix A;
        A.Stype = SLU_NR;           // storage type: compressed sparse, row-major
#ifdef SINGLE
        A.Dtype = SLU_C;            // data type: single complex
#else
        A.Dtype = SLU_Z;            // data type: double complex
#endif
        A.Mtype = SLU_GE;           // mathematical type: general
        A.nrow  = P_.size() - 1;    // number of rows
        A.ncol  = P_.size() - 1;    // number of columns
        A.Store = &AStore;          // data structure pointer
    
    //
    // Create the (empty) right hand side and solution matrix.
    //
    
        DNformat BStore, XStore;
        SuperMatrix B, X;
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
        options.ILU_DropTol = data.drop_tolerance;
        
        // calculation diagnostic information
        SuperLUStat_t stat;
        StatInit(&stat);
        
        // memory usage
        mem_usage_t mem_usage;
    
    //
    // Compute the factorization.
    //
        
        // permutation arrays, elimination tree
        perm_r_.resize(A.nrow);
        perm_c_.resize(A.ncol);
        etree_.resize(A.nrow);
        
        // row and column scale factors, reciprocal condition number, reciprocal pivot growth factor
        R_.resize(A.nrow);
        C_.resize(A.ncol);
        Real rcond, rpg;
        
        // forward and backward errors (one element per one right hand side)
        Real ferr, berr;
        
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
            &perm_c_[0],    // column permutation
            &perm_r_[0],    // row permutation
            &etree_[0],     // elimination tree
            &equed_,        // equilibration done
            &R_[0],         // row scale factors
            &C_[0],         // column scale factors
            &L_,            // L-factor
            &U_,            // U-factor
            nullptr,        // workspace (not used)
            0,              // size of the workspace (0 = automatic allocation)
            &B,             // right-hand sides (empty)
            &X,             // solution matrix (empty)
            &rpg,           // reciprocal pivot growth factor
            &rcond,         // reciprocal condition number
            &ferr,          // forward error
            &berr,          // backward error
            &Glu_,          // reusable information
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
        
        size_ = mem_usage.for_lu;
}

template<>
void LUft_SUPERLU<LU_int_t,Complex>::solve (const cArrayView b, cArrayView x, int eqs) const
{
    //
    // Create matrix of the system.
    //
    
        NRformat AStore;
        AStore.nnz    = X_.size();                         // number of non-zero elements
        AStore.nzval  = const_cast<Complex*>(X_.data());   // pointer to the array of non-zero elements
        AStore.colind = const_cast<int*>(I_.data());       // row indices
        AStore.rowptr = const_cast<int*>(P_.data());       // column pointers
        
        SuperMatrix A;
        A.Stype = SLU_NR;           // storage type: compressed sparse, row-major
#ifdef SINGLE
        A.Dtype = SLU_C;            // data type: single complex
#else
        A.Dtype = SLU_Z;            // data type: double complex
#endif
        A.Mtype = SLU_GE;           // mathematical type: general
        A.nrow  = P_.size() - 1;    // number of rows
        A.ncol  = P_.size() - 1;    // number of columns
        A.Store = &AStore;          // data structure pointer
    
    //
    // Create the right hand side.
    //
    
        DNformat BStore;
        BStore.lda = P_.size() - 1;                     // leading dimension
        BStore.nzval = const_cast<Complex*>(b.data());  // data pointer
        
        SuperMatrix B;
        B.Stype = SLU_DN;               // storage type: Fortran dense matrix
#ifdef SINGLE
        B.Dtype = SLU_C;                // data type: single complex
#else
        B.Dtype = SLU_Z;                // data type: double complex
#endif
        B.Mtype = SLU_GE;               // mathematical type: general
        B.nrow  = P_.size() - 1;        // number of rows (= number of equations)
        B.ncol  = eqs;                  // number of columns (= number of right hand sides)
        B.Store = &BStore;              // data structure pointer
        
    //
    // Create the solution matrix.
    //
        
        DNformat XStore;
        XStore.lda = P_.size() - 1;     // leading dimension
        XStore.nzval = x.data();        // data pointer
        
        SuperMatrix X;
        X.Stype = SLU_DN;               // storage type: Fortran dense matrix
#ifdef SINGLE
        X.Dtype = SLU_C;                // data type: single complex
#else
        X.Dtype = SLU_Z;                // data type: double complex
#endif
        X.Mtype = SLU_GE;               // mathematical type: general
        X.nrow  = P_.size() - 1;        // number of rows (= number of equations)
        X.ncol  = eqs;                  // number of columns (= number of right hand sides)
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
        Real rcond, rpg;
        
        // back-substitution
        int info;
#ifdef SINGLE
        cgssvx
#else
        zgssvx
#endif
        (
            &options,                       // calculation options
            &A,                             // matrix data structure
            const_cast<int*>(&perm_c_[0]),  // column permutation
            const_cast<int*>(&perm_r_[0]),  // row permutation
            const_cast<int*>(&etree_[0]),   // elimination tree
            const_cast<char*>(&equed_),     // equilibration done
            const_cast<Real*>(&R_[0]),      // row scale factors
            const_cast<Real*>(&C_[0]),      // column scale factors
            const_cast<SuperMatrix*>(&L_),  // L-factor
            const_cast<SuperMatrix*>(&U_),  // U-factor
            nullptr,                        // workspace (not used)
            0,                              // size of the workspace (0 = dynamic allocation)
            &B,                             // right-hand sides
            &X,                             // solution matrix
            &rpg,                           // reciprocal pivot growth factor
            &rcond,                         // reciprocal condition number
            &ferr[0],                       // forward error
            &berr[0],                       // backward error
            const_cast<GlobalLU_t*>(&Glu_), // reusable information
            &mem_usage,                     // memory usage
            &stat,                          // diagnostic infomation
            &info                           // result status
        );
    
    // check status
    if (info != 0)
        HexException("SuperLU/?gssvx failed with status %d.", info);
}

template<>
bool LUft_SUPERLU<LU_int_t,Complex>::valid () const
{
    return size_ != 0;
}

template<>
void LUft_SUPERLU<LU_int_t,Complex>::save (std::string name) const
{
    /*HexException("SuperLU factorizer does not yet support --out-of-core option.");*/
}

template<>
void LUft_SUPERLU<LU_int_t,Complex>::load (std::string name, bool throw_on_io_failure)
{
    if (throw_on_io_failure)
        HexException("SuperLU factorizer does not yet support --out-of-core option.");
}

template<>
void LUft_SUPERLU<LU_int_t,Complex>::drop ()
{
    if (size_ != 0)
    {
        Destroy_SuperNode_Matrix(&L_);
        Destroy_CompCol_Matrix(&U_);
        size_ = 0;
    }
}

// --------------------------------------------------------------------------------- //

addFactorizerToRuntimeSelectionTable(SUPERLU, LU_int_t, Complex)

// --------------------------------------------------------------------------------- //

#endif // WITH_SUPERLU
