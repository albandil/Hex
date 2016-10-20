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

#ifdef WITH_UMFPACK

// --------------------------------------------------------------------------------- //

#include "hex-lu-umfpack.h"
#include "hex-csrmatrix.h"

// --------------------------------------------------------------------------------- //

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
bool LUft_UMFPACK<LU_int_t,Complex>::valid () const
{
    return numeric_ != nullptr and size() > 0;
}

template<>
double LUft_UMFPACK<LU_int_t,Complex>::cond () const
{
    return info_[UMFPACK_RCOND];
}

template<>
void LUft_UMFPACK<LU_int_t,Complex>::factorize (CsrMatrix<LU_int_t,Complex> const & matrix, LUftData data)
{
    // Use standard UMFPACK sequence
    void *Symbolic, *Numeric;
    LU_int_t status;
    
    // get default setting
    double Control[UMFPACK_CONTROL];
    UMFPACK_DEFAULTS_F(Control);
    
    // modify the drop tolerance
    Control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_SYMMETRIC;
    Control[UMFPACK_DROPTOL] = data.drop_tolerance;
    
    // diagnostic information
    rArray Info (UMFPACK_INFO);
    
    // matrix data
    LU_int_t m = matrix.rows();
    LU_int_t n = matrix.cols();
    p_ = matrix.p();
    i_ = matrix.i();
    x_ = matrix.x();
    
    // analyze the sparse structure
    status = UMFPACK_SYMBOLIC_F
    (
        m, n,                       // matrix dimensions
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
    
    // release symbolic data
    UMFPACK_FREE_SYMBOLIC_F(&Symbolic);
    
    // store numeric data
    numeric_ = Numeric;
}

template<>
void LUft_UMFPACK<LU_int_t,Complex>::solve (const cArrayView b, cArrayView x, int eqs) const
{
    // number of unknowns
    std::size_t N = p_.size() - 1;
    
    // check sizes
    assert(eqs * N == x.size());
    assert(eqs * N == b.size());
    
    // solve for all RHSs
    for (int eq = 0; eq < eqs; eq++)
    {
        // solve for current RHS
        int status = UMFPACK_SOLVE_F
        (
            UMFPACK_Aat,    // matrix orientation
            p_.data(),      // row pointers
            i_.data(),      // column indices
            
            // matrix elements (interleaved)
            reinterpret_cast<const double*>(x_.data()),
            nullptr,
            
            // solutions (interleaved)
            reinterpret_cast<double*>(&x[0] + eq * N),
            nullptr,
            
            // right-hand side vectors (interleaved)
            reinterpret_cast<const double*>(&b[0] + eq * N),
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
    p_.hdfsave("csr-p-" + name);
    i_.hdfsave("csr-i-" + name);
    x_.hdfsave("csr-x-" + name);
    
    int err = UMFPACK_SAVE_NUMERIC_F(numeric_, const_cast<char*>(name.c_str()));
    
    if (err == UMFPACK_ERROR_invalid_Numeric_object)
        HexException("[LUft::save] Invalid numeric object.");
    
    if (err == UMFPACK_ERROR_file_IO)
        HexException("[LUft::save] Failed to save LU object \"%s\" (size = %ld).", name.c_str(), size());
}

template<>
void LUft_UMFPACK<LU_int_t,Complex>::load (std::string name, bool throw_on_io_failure)
{
    if (not p_.hdfload("csr-p-" + name) or
        not i_.hdfload("csr-i-" + name) or
        not x_.hdfload("csr-x-" + name))
    {
        if (throw_on_io_failure)
            HexException("[LUft::load] Failed to load the matrix data from files \"csr-*-%s\".", name.c_str());
        
        return;
    }
    
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
        p_.drop();
        i_.drop();
        x_.drop();
    }
}

// --------------------------------------------------------------------------------- //

addFactorizerToRuntimeSelectionTable(UMFPACK, LU_int_t, Complex)

// --------------------------------------------------------------------------------- //

#endif // WITH_UMFPACK
