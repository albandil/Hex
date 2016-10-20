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

#ifdef WITH_SUPERLU_DIST

// --------------------------------------------------------------------------------- //

#include "hex-lu-superlu-dist.h"

// --------------------------------------------------------------------------------- //

#include <superlu_zdefs.h>

// --------------------------------------------------------------------------------- //

template<>
LUft_SUPERLU_DIST<LU_int_t,Complex>::LUft_SUPERLU_DIST ()
    : LUft<LU_int_t,Complex>(), size_(0)
{
    // nothing
}

template<>
void LUft_SUPERLU_DIST<LU_int_t,Complex>::drop ()
{
    if (size_ != 0)
    {
        Destroy_LU(P_.size() - 1, grid_, &LUstruct_);
        ScalePermstructFree(&ScalePermstruct_);
        LUstructFree(&LUstruct_);
        P_.drop();
        I_.drop();
        X_.drop();
        size_ = 0;
    }
}

template<>
LUft_SUPERLU_DIST<LU_int_t,Complex>::~LUft_SUPERLU_DIST ()
{
    drop ();
}

template<>
void LUft_SUPERLU_DIST<LU_int_t,Complex>::factorize (CsrMatrix<LU_int_t,Complex> const & matrix, LUftData data)
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
        
        NCformat AStore;
        AStore.nnz    = X_.size();       // number of non-zero elements
        AStore.nzval  = &X_[0];          // pointer to the array of non-zero elements
        AStore.rowind = &I_[0];          // row indices
        AStore.colptr = &P_[0];          // column pointers
        
        SuperMatrix A;
        A.Stype = SLU_NC;        // storage type: compressed sparse, column-major (SuperLU-dist suports no other)
#ifdef SINGLE
        A.Dtype = SLU_C;         // data type: single complex
#else
        A.Dtype = SLU_Z;         // data type: double complex
#endif
        A.Mtype = SLU_GE;        // mathematical type: general
        A.nrow  = P_.size() - 1; // number of rows
        A.ncol  = P_.size() - 1; // number of columns
        A.Store = &AStore;       // data structure pointer

    //
    // Prepare SuperLU environment.
    //
    
        // get process grid
        gridinfo_t * grid = (gridinfo_t *)data.superlu_dist_grid;
        
        // calculation options
        superlu_dist_options_t options;
        set_default_options_dist(&options);
        options.ParSymbFact = YES;
        options.ColPerm = METIS_AT_PLUS_A;
        options.PrintStat = NO;
        options.SymPattern = YES;
//         options.ILU_DropRule = DROP_BASIC;
//         options.ILU_DropTol = droptol;
        
        // distributed scale and permutation data
        ScalePermstructInit(A.nrow, A.ncol, &ScalePermstruct_);
        
        // distributed factorization data
        LUstructInit(A.nrow, &LUstruct_);
        
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
            &ScalePermstruct_,  // scaling and permutation data
            nullptr,            // right hand sides (not used)
            A.nrow,             // right hand sides leading dimension
            0,                  // number of right hand sides (none, only factorizing)
            grid,               // process grid
            &LUstruct_,         // factorization data
            berr,               // backward error
            &stat,              // diagnostic information
            &info               // result status
        );
        
        superlu_dist_mem_usage_t mem_usage;
        zQuerySpace_dist(A.nrow, &LUstruct_, grid, &stat, &mem_usage);
        
        if (info > A.ncol)
            HexException("SuperLU/zgssvx: Memory allocation failure after %d bytes.", info);
        if (info > 0)
            HexException("SuperLU/zgssvx: Singular factor.");
        
        size_ = mem_usage.for_lu;
}

template<>
void LUft_SUPERLU_DIST<LU_int_t,Complex>::solve (const cArrayView b, cArrayView x, int eqs) const
{
    //
    // Create matrix of the system.
    //
    
        NCformat AStore;
        AStore.nnz    = X_.size();                      // number of non-zero elements
        AStore.nzval  = const_cast<Complex*>(&X_[0]);   // pointer to the array of non-zero elements
        AStore.rowind = const_cast<int_t*>(&I_[0]);     // row indices
        AStore.colptr = const_cast<int_t*>(&P_[0]);     // column pointers
        
        SuperMatrix A;
        A.Stype = SLU_NC;           // storage type: compressed sparse, column-major (SuperLU-dist suports no other)
        A.Dtype = SLU_Z;            // data type: double complex
        A.Mtype = SLU_GE;           // mathematical type: general
        A.nrow  = P_.size() - 1;    // number of rows
        A.ncol  = P_.size() - 1;    // number of columns
        A.Store = &AStore;          // data structure pointer
    
    //
    // Prepare SuperLU environment.
    //
    
        // calculation options
        superlu_dist_options_t options;
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
        std::vector<double> berr(eqs);
        
        // status indicator
        int info;
        
        // LU factorization
        x = b;
        pzgssvx_ABglobal
        (
            &options,                                           // calculation options
            &A,                                                 // matrix to factorize
            const_cast<ScalePermstruct_t*>(&ScalePermstruct_),  // scaling and permutation data
            const_cast<doublecomplex*>(reinterpret_cast<const doublecomplex*>(&x[0])), // right hand sides (assume packed)
            A.nrow,                                             // right hand sides leading dimension
            eqs,                                                // number of right hand sides
            grid_,                                              // process grid
            const_cast<LUstruct_t*>(&LUstruct_),                // factorization data
            &berr[0],                                           // backward error
            &stat,                                              // diagnostic information
            &info                                               // result status
        );
        
        if (info > A.ncol)
            HexException("SuperLU/?gssvx: Memory allocation failure after %d bytes.", info);
        if (info > 0)
            HexException("SuperLU/?gssvx: Singular factor.");
    
    //
    // Clean up.
    //
    
        PStatFree(&stat);
}

template<>
void LUft_SUPERLU_DIST<LU_int_t,Complex>::save (std::string name) const
{ 
    HexException("SuperLU_dist factorizer does not yet support --out-of-core option.");
}

template<>
void LUft_SUPERLU_DIST<LU_int_t,Complex>::load (std::string name, bool throw_on_io_failure)
{
    if (throw_on_io_failure)
        HexException("SuperLU_dist factorizer does not yet support --out-of-core option.");
}

// --------------------------------------------------------------------------------- //

addFactorizerToRuntimeSelectionTable(SUPERLU_DIST, LU_int_t, Complex)

// --------------------------------------------------------------------------------- //

#endif // WITH_SUPERLU_DIST
