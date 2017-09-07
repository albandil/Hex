//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2017, Jakub Benda, Charles University in Prague                    //
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

#include "lu-umfpack.h"

// --------------------------------------------------------------------------------- //

LUft_UMFPACK::LUft_UMFPACK ()
    : LUft(), numeric_(nullptr), info_(UMFPACK_INFO)
{
    rdata_["drop_tolerance"] = 1e-8;
}

std::size_t LUft_UMFPACK::size () const
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

bool LUft_UMFPACK::valid () const
{
    return numeric_ != nullptr and size() > 0;
}

Real LUft_UMFPACK::cond () const
{
    return info_[UMFPACK_RCOND];
}

void LUft_UMFPACK::factorize (CsrMatrix<LU_int_t,Complex> const & matrix)
{
    // Use standard UMFPACK sequence
    void *Symbolic, *Numeric;
    LU_int_t status;
    
    // get default setting
    double Control[UMFPACK_CONTROL];
    UMFPACK_DEFAULTS_F(Control);
    
    // modify the drop tolerance
    Control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_SYMMETRIC;
    Control[UMFPACK_DROPTOL] = rdata_["drop_tolerance"];
    
    // diagnostic information
    double Info[UMFPACK_INFO];
    
    // matrix data
    LU_int_t m = matrix.rows();
    LU_int_t n = matrix.cols();
#ifndef SINGLE
    x_ = matrix.x();
#else
    x_.resize(0);
    x_.reserve(matrix.x().size());
    for (Complex x : matrix.x())
	x_.push_back(std::complex<double>(x.real(), x.imag()));
#endif
    p_ = matrix.p();
    i_ = matrix.i();
#ifndef SINGLE
    x_ = matrix.x();
#else
    x_.resize(matrix.x().size());
    for (std::size_t i = 0; i < x_.size(); i++) x_[i] = Complex(matrix.x()[i].real(), matrix.x()[i].imag());
#endif
    
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
        Symbolic, &Numeric, Control, Info    // UMFPACK internals
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

void LUft_UMFPACK::solve (const cArrayView b, cArrayView x, int eqs) const
{
    // number of unknowns
    std::size_t N = p_.size() - 1;
    
    // check sizes
    assert(eqs * N == x.size());
    assert(eqs * N == b.size());
    
#ifdef SINGLE
    NumberArray<std::complex<double>> B(b.size()), X(x.size());
    for (std::size_t i = 0; i < b.size(); i++)
    {
        B[i].real(b[i].real());
        B[i].imag(b[i].imag());
        X[i].real(x[i].real());
        X[i].imag(x[i].imag());
    }
#else
    cArrayView B(b), X(x);
#endif
    
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
            reinterpret_cast<double*>(&X[0] + eq * N),
            nullptr,
            
            // right-hand side vectors (interleaved)
            reinterpret_cast<const double*>(&B[0] + eq * N),
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
    
#ifdef SINGLE
    for (std::size_t i = 0; i < x.size(); i++)
    {
        x[i].real(X[i].real());
        x[i].imag(X[i].imag());
    }
#endif
}

void LUft_UMFPACK::save (std::string name) const
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

void LUft_UMFPACK::load (std::string name, bool throw_on_io_failure)
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

void LUft_UMFPACK::drop ()
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

addClassToParentRunTimeSelectionTable(LUft, LUft_UMFPACK)

// --------------------------------------------------------------------------------- //

#endif // WITH_UMFPACK
