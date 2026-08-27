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

#if (!defined(HEX_LU_PARDISO_H) && (defined(WITH_PARDISO) || defined(WITH_MKL)))
#define HEX_LU_PARDISO_H

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
class LUft_Pardiso : public LUft
{
    public:

        // run-time selection mechanism
        factorizerRunTimeSelectionDefinitions(LUft_Pardiso, "pardiso")

        /// Default constructor.
        LUft_Pardiso ();

        /// Destructor.
        ~LUft_Pardiso() override;

        // Disable bitwise copy
        LUft_Pardiso const & operator= (LUft_Pardiso const &) = delete;

        /// Factorize.
        void factorize (CsrMatrix<LU_int_t,Complex> const & matrix) override;

        /// Validity indicator.
        bool valid () const override { return size() != 0; }

        /// Return LU byte size.
        std::size_t size () const override;

        /// Solve equations.
        void solve (const cArrayView b, cArrayView x, int eqs) const override;

        /// Save to disk.
        void save (std::string name) const override;

        /// Load from disk.
        void load (std::string name, bool throw_on_io_failure = true) override;

        /// Release memory.
        void drop () override;

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
