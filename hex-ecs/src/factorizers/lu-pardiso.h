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

#if (defined(WITH_PARDISO) || defined(WITH_MKL))

// --------------------------------------------------------------------------------- //

#include "hex-csrmatrix.h"

#include "luft.h"

// --------------------------------------------------------------------------------- //

#define DPARM(x) dparm_[x-1]
#define IPARM(x) iparm_[x-1]

// --------------------------------------------------------------------------------- //

#ifdef WITH_MKL

// old interfaces still in Intel MKL

extern "C" void pardisoinit
(
    void* pt,               // internal data structure
    const int* mtype,       // matrix type
    int* iparm              // integer parameters
);

extern "C" void pardiso
(
    void* pt,               // internal data structure
    const int* maxfct,      // maximal number of factorizations
    const int* mnum,        // factorization index
    const int* mtype,       // matrix type
    const int* phase,       // solver phase
    const int* n,           // matrix rank
    const void* a,          // matrix non-zero elements
    const int* ia,          // matrix column pointers
    const int* ja,          // matrix row indices
    int* perm,              // permutation array
    const int* nrhs,        // number of right-hand sides
    int* iparm,             // integer parameters
    const int* msglvl,      // verbosity level
    void* b,                // right-hand sides
    void* x,                // solutions
    int* error              // status indicator
);

#else

// new interfaces (as of Pardiso 5.0.0)

extern "C" void pardisoinit
(
    void* pt,               // internal data structure
    const int* mtype,       // matrix type
    const int* solver,      // solver
    int* iparm,             // integer parameters
    double* dparm,          // real parameters
    int* error              // status indicator
);

extern "C" void pardiso
(
    void* pt,               // internal data structure
    const int* maxfct,      // maximal number of factorizations
    const int* mnum,        // factorization index
    const int* mtype,       // matrix type
    const int* phase,       // solver phase
    const int* n,           // matrix rank
    const void* a,          // matrix non-zero elements
    const int* ia,          // matrix column pointers
    const int* ja,          // matrix row indices
    int* perm,              // permutation array
    const int* nrhs,        // number of right-hand sides
    int* iparm,             // integer parameters
    const int* msglvl,      // verbosity level
    void* b,                // right-hand sides
    void* x,                // solutions
    int* error,             // status indicator
    double* dparm           // double parameters
);

#endif

// --------------------------------------------------------------------------------- //

/**
 * @brief LU factorization using Pardiso.
 * 
 * Uses direct sparse solver Pardiso.
 */
template <class IdxT, class DataT>
class LUft_Pardiso : public LUft<IdxT,DataT>
{
    public:
    
        /// Default constructor.
        LUft_Pardiso ();
        
        /// Destructor.
        virtual ~LUft_Pardiso();
        
        // Disable bitwise copy
        LUft_Pardiso const & operator= (LUft_Pardiso const &) = delete;
        
        /// New instance of the factorizer.
        virtual LUft<IdxT,DataT> * New () const { return new LUft_Pardiso<IdxT,DataT>(); }
        
        /// Get name of the factorizer.
        virtual std::string name () const { return "pardiso"; }
        
        /// Factorize.
        virtual void factorize (CsrMatrix<IdxT,DataT> const & matrix, LUftData data);
        
        /// Validity indicator.
        virtual bool valid () const { return size() != 0; }
        
        /// Return LU byte size.
        virtual std::size_t size () const;
        
        /// Solve equations.
        virtual void solve (const ArrayView<DataT> b, ArrayView<DataT> x, int eqs) const;
        
        /// Save to disk.
        virtual void save (std::string name) const;
        
        /// Load from disk.
        virtual void load (std::string name, bool throw_on_io_failure = true);
        
        /// Release memory.
        virtual void drop ();
    
    private:
        
        /// Inspect the returned success indicator.
        void pardisoerror (int error) const;
        
        /// Matrix that has been factorized.
        NumberArray<int> P_;
        NumberArray<int> I_;
        NumberArray<std::complex<double>> X_;
        
        /// Permutation.
        NumberArray<int> perm_;
        
        /// Internal data of Pardiso.
        void* pt_[64];
        int iparm_[64];
        double dparm_[64];
};

// --------------------------------------------------------------------------------- //

#endif // WITH_MKL or WITH_PARDISO
