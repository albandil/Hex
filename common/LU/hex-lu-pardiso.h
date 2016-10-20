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

// --------------------------------------------------------------------------------- //

#include "hex-luft.h"
#include "hex-csrmatrix.h"

// --------------------------------------------------------------------------------- //

#define DPARM(x) dparm_[x-1]
#define IPARM(x) iparm_[x-1]

// --------------------------------------------------------------------------------- //

extern "C" void pardisoinit
(
    void *pt[64],
    int *mtype,
    int *solver,
    int iparm[64],
    double dparm[64],
    int *error
);

extern "C" void pardiso
(
    void *pt[64],
    int *maxfct,
    int *mnum,
    int *mtype,
    int *phase,
    int *n,
    double a[],
    int ia[],
    int ja[],
    int perm[],
    int *nrhs,
    int iparm[64],
    int *msglvl,
    double b[],
    double x[],
    int *error,
    double dparm[64]
);

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
        
        /// Matrix that has been factorized.
        NumberArray<int> P_;
        NumberArray<int> I_;
        NumberArray<std::complex<double>> X_;
        
        /// Internal data of Pardiso.
        void* pt_[64];
        int iparm_[64];
        double dparm_[64];
};

// --------------------------------------------------------------------------------- //

#endif // WITH_PARDISO
