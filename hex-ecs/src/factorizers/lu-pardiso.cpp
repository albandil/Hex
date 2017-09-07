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

#if (defined(WITH_PARDISO) || defined(WITH_MKL))

// --------------------------------------------------------------------------------- //

#ifdef _OPENMP
    #include <omp.h>
#endif

// --------------------------------------------------------------------------------- //

#include "lu-pardiso.h"

// --------------------------------------------------------------------------------- //

void LUft_Pardiso::pardisoerror (int error) const
{
    switch (error)
    {
        case 0:
            return;
        case -1:
            HexException("Pardiso: Input inconsistent.");
        case -2:
            HexException("Pardiso: Not enough memory.");
        case -3:
            HexException("Pardiso: Reordering problem.");
        case -4:
            HexException("Pardiso: Zero pivot, numerical fact. or iterative refinement problem.");
        case -5:
            HexException("Pardiso: Unclassified (internal) error.");
        case -6:
            HexException("Pardiso: Preordering failed.");
        case -7:
            HexException("Pardiso: Diagonal matrix problem.");
        case -8:
            HexException("Pardiso: 32-bit integer overflow problem.");
        case -10:
            HexException("Pardiso: No license file \"pardiso.lic\" found.");
        case -11:
            HexException("Pardiso: License is expired.");
        case -12:
            HexException("Pardiso: Wrong username or hostname.");
        case -100:
            HexException("Pardiso: Reached maximum number of Krylov-subspace iteration in iterative solver.");
        case -101:
            HexException("Pardiso: No sufficient convergence in Krylov-subspace iteration within 25 iterations.");
        case -102:
            HexException("Pardiso: Error in Krylov-subspace iteration.");
        case -103:
            HexException("Pardiso: Break-down in Krylov-subspace iteration.");
        default:
            HexException("Pardiso: Unknown error.");
    }
}

LUft_Pardiso::LUft_Pardiso () : LUft()
{
    std::memset(pt_,    0, sizeof(pt_));
    std::memset(iparm_, 0, sizeof(iparm_));
    std::memset(dparm_, 0, sizeof(dparm_));
    
    int mtype  = 13;     // complex symmetric
    int solver = 0;     // sparse direct solver
    int error  = 0;     // success indicator
    
    idata_["groupsize"] = 1;
    rdata_["drop_tolerance"] = 1e-8;
    
#ifdef WITH_MKL
    pardisoinit(pt_, &mtype, iparm_);
#else
    pardisoinit(pt_, &mtype, &solver, iparm_, dparm_, &error);
    pardisoerror(error);
#endif
}

void LUft_Pardiso::drop ()
{
    int maxfct = 1;     // maximal number of numerical factorizations
    int mtype  = 13;    // matrix type: complex symmetric
    int phase  = -1;    // release all internal memory for all matrices
    int msglvl = 0;     // verbosity
    int error  = 0;     // success indicator
    
    int O      = 0;     // dummy integer
    double D   = 0;     // dummy double
    
#ifdef WITH_MKL
    pardiso(pt_, &maxfct, &O, &mtype, &phase, &O, &D, &O, &O, &O, &O, iparm_, &msglvl, &D, &D, &error);
#else
    pardiso(pt_, &maxfct, &O, &mtype, &phase, &O, &D, &O, &O, &O, &O, iparm_, &msglvl, &D, &D, &error, dparm_);
#endif
    
    pardisoerror(error);
}

LUft_Pardiso::~LUft_Pardiso ()
{
    drop();
}

void LUft_Pardiso::factorize (CsrMatrix<LU_int_t,Complex> const & matrix)
{
    //
    // Cast matrix data to the data types expected by Pardiso (32-bit Int and 8-byte Real).
    //
    
        // copy row pointers
        P_.resize(matrix.p().size());
        for (std::size_t i = 0; i < P_.size(); i++)
            P_[i] = matrix.p()[i] + 1;
        
        // copy column indices
        I_.resize(matrix.i().size());
        for (std::size_t i = 0; i < I_.size(); i++)
            I_[i] = matrix.i()[i] + 1;
        
        // copy elements
        X_.resize(matrix.x().size());
        for (std::size_t i = 0; i < X_.size(); i++)
            X_[i] = std::complex<double>(matrix.x()[i].real(), matrix.x()[i].imag());
    
    //
    // Calculate the LU factorization.
    //
    
        int maxfct = 1;                 // maximal number of factorizations
        int mnum   = 1;                 // index of the matrix to factorize
        int mtype  = 13;                // complex symmetric
        int phase  = 12;                // analysis & numerical factorization
        int n      = P_.size() - 1;     // rank of the matrix
        int msglvl = 0;                 // verbosity
        int error  = 0;                 // success indicator
        int nrhs   = 1;                 // number of right-hand sides
        
#ifdef WITH_PARDISO
    #ifdef _OPENMP
        IPARM(3) = omp_get_max_threads();
    #else
        IPARM(3) = 1;
    #endif
        IPARM(52) = idata_["groupsize"];
#endif
        
        perm_.resize(n);
        
        DPARM(5) = rdata_["drop_tolerance"];     // ILU drop tolerance
        
        pardiso
        (
            pt_,
            &maxfct,
            &mnum,
            &mtype,
            &phase,
            &n,
            reinterpret_cast<double*>(X_.data()),
            P_.data(),
            I_.data(),
            perm_.data(),
            &nrhs,
            iparm_,
            &msglvl,
            nullptr,
            nullptr,
            &error
#ifdef WITH_PARDISO
            , dparm_
#endif
        );
        
        pardisoerror(error);
}

void LUft_Pardiso::solve (const cArrayView b, cArrayView x, int eqs) const
{
    int maxfct = 1;
    int mnum = 1;
    int mtype = 13;
    int phase = 33;
    int n = b.size() / eqs;
    int msglvl = 0;
    int error = 0;
    
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
        const_cast<int*>(perm_.data()),
        &eqs,
        const_cast<int*>(iparm_),
        &msglvl,
        const_cast<double*>(reinterpret_cast<const double*>(b.data())),
        const_cast<double*>(reinterpret_cast<const double*>(x.data())),
        &error
#ifdef WITH_PARDISO
        , const_cast<double*>(dparm_)
#endif
    );
    
    pardisoerror(error);
}

std::size_t LUft_Pardiso::size () const
{
    // IPARM(18) is negative for invalid / uninitialized setups,
    // but Hex-ecs is expecting zero
    
    return IPARM(18) > 0 ? IPARM(18) * std::size_t(16) : 0;
}

void LUft_Pardiso::save (std::string name) const
{ 
    HexException("Pardiso factorizer does not yet support --out-of-core option.");
}

void LUft_Pardiso::load (std::string name, bool throw_on_io_failure)
{
    if (throw_on_io_failure)
        HexException("Pardiso factorizer does not yet support --out-of-core option.");
}


// --------------------------------------------------------------------------------- //

addClassToParentRunTimeSelectionTable(LUft, LUft_Pardiso)

// --------------------------------------------------------------------------------- //

#endif // WITH_PARDISO
