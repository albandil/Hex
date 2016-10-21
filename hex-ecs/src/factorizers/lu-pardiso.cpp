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

#ifdef WITH_PARDISO

#ifdef _OPENMP
#include <omp.h>
#endif

// --------------------------------------------------------------------------------- //

#include "lu-pardiso.h"

// --------------------------------------------------------------------------------- //

template<>
LUft_Pardiso<LU_int_t,Complex>::LUft_Pardiso ()
    : LUft<LU_int_t,Complex>()
{
    std::memset(pt_, 0, sizeof(pt_));
    std::memset(iparm_, 0, sizeof(iparm_));
    std::memset(dparm_, 0, sizeof(dparm_));
}

template<>
LUft_Pardiso<LU_int_t,Complex>::~LUft_Pardiso ()
{
    // nothing
}

template<>
void LUft_Pardiso<LU_int_t,Complex>::factorize (CsrMatrix<LU_int_t,Complex> const & matrix, LUftData data)
{
    int mtype = 6;  // complex symmetric
    int solver = 0; // sparse direct solver
    int error;
    
    //
    // Cast matrix data to the data types expected by Pardiso (32-bit Int and 8-byte Real).
    //
    
        // copy row pointers
        P_.resize(matrix.p().size());
        for (std::size_t i = 0; i < P_.size(); i++)
            P_[i] = matrix.p()[i];
        
        // copy column indices
        I_.resize(matrix.i().size());
        for (std::size_t i = 0; i < I_.size(); i++)
            I_[i] = matrix.i()[i];
        
        // copy elements
        X_.resize(matrix.x().size());
        for (std::size_t i = 0; i < X_.size(); i++)
            X_[i] = std::complex<double>(matrix.x()[i].real(), matrix.x()[i].imag());
    
    //
    // Initialize Pardiso.
    //
    
        std::memset(pt_, 0, sizeof(pt_));
        std::memset(iparm_, 0, sizeof(iparm_));
        std::memset(dparm_, 0, sizeof(dparm_));
        
        pardisoinit (pt_, &mtype, &solver, iparm_, dparm_, &error);
    
    //
    // Calculate the LU factorization.
    //
    
        int maxfct = 1;                 // maximal number of factorizations
        int mnum = 1;                   // index of the matrix to factorize
        int phase = 12;                 // analysis & numerical factorization
        int n = P_.size() - 1;          // rank of the matrix
        std::vector<double> perm (n);   // permutation
        int msglvl = 1;                 // verbosity
        
        DPARM(5) = data.drop_tolerance;     // ILU drop tolerance
        
#ifdef _OPENMP
        IPARM(3) = omp_get_num_threads();   // OpenMP threads
#else
        IPARM(3) = 1;
#endif
#ifdef WITH_MPI
        MPI_Comm_size(MPI_Comm_f2c(data.fortran_comm), &IPARM(52)); // number of MPI processes used for this factorization
#else
        IPARM(52) = 1;
#endif
        
        pardiso
        (
            pt_,
            &maxfct,
            &mnum,
            &mtype,
            &phase,
            &n,
            reinterpret_cast<double*>(X_.data()),
            const_cast<int*>(P_.data()),
            const_cast<int*>(I_.data()),
            nullptr,
            nullptr,
            iparm_,
            &msglvl,
            nullptr,
            nullptr,
            &error,
            dparm_
        );
        
        if (error != 0)
            HexException("Pardiso failed to factorize the matrix (code = %d).", error);
}

template <>
void LUft_Pardiso<LU_int_t,Complex>::solve (const cArrayView b, cArrayView x, int eqs) const
{
    int maxfct = 1;
    int mnum = 1;
    int mtype = 6;
    int phase = -1;
    int n = b.size();
    int msglvl = 1;
    int error;
    
    pardiso
    (
        const_cast<void**>(pt_),
        &maxfct,
        &mnum,
        &mtype,
        &phase,
        &n,
        const_cast<double*>(reinterpret_cast<const double*>(X_.data())),
        const_cast<int*>(P_.data()),
        const_cast<int*>(I_.data()),
        nullptr,
        &eqs,
        const_cast<int*>(iparm_),
        &msglvl,
        const_cast<double*>(reinterpret_cast<const double*>(b.data())),
        const_cast<double*>(reinterpret_cast<const double*>(x.data())),
        &error,
        const_cast<double*>(dparm_)
    );
    
    if (error != 0)
        HexException("Pardiso failed to solve.");
}

template <>
void LUft_Pardiso<LU_int_t,Complex>::drop ()
{
    int maxfct = 1;
    int mnum = 1;
    int mtype = 6;
    int phase = -1;
    int n = P_.size() - 1;
    int msglvl = 1;
    int error = 1;
    
    pardiso
    (
        pt_,
        &maxfct,
        &mnum,
        &mtype,
        &phase,
        &n,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        iparm_,
        &msglvl,
        nullptr,
        nullptr,
        &error,
        dparm_
    );
    
    if (error != 0)
        HexException("Pardiso failed to release memory.");
}


template <>
std::size_t LUft_Pardiso<LU_int_t,Complex>::size () const
{
    return IPARM(18) * std::size_t(16);
}

template<>
void LUft_Pardiso<LU_int_t,Complex>::save (std::string name) const
{ 
    HexException("Pardiso factorizer does not yet support --out-of-core option.");
}

template<>
void LUft_Pardiso<LU_int_t,Complex>::load (std::string name, bool throw_on_io_failure)
{
    if (throw_on_io_failure)
        HexException("Pardiso factorizer does not yet support --out-of-core option.");
}

// --------------------------------------------------------------------------------- //

addFactorizerToRuntimeSelectionTable(Pardiso, LU_int_t, Complex)

// --------------------------------------------------------------------------------- //

#endif // WITH_PARDISO
