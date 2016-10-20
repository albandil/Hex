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

#include "hex-luft.h"
#include "hex-csrmatrix.h"

// --------------------------------------------------------------------------------- //

#include <superlu_zdefs.h>

// --------------------------------------------------------------------------------- //

/**
 * @brief LU factorization object - SuperLU-dist specialization.
 * 
 * This class holds information on LU factorization as computed by the free
 * library SuperLU (distributed version). It is derived from LUft and
 * shares interface with that class.
 */
template <class IdxT, class DataT>
class LUft_SUPERLU_DIST : public LUft<IdxT,DataT>
{
    public:
        
        /// Default constructor.
        LUft_SUPERLU_DIST ();
        
        /// Destructor.
        virtual ~LUft_SUPERLU_DIST ();
        
        // Disable bitwise copy
        LUft_SUPERLU_DIST const & operator= (LUft_SUPERLU_DIST const &) = delete;
        
        /// New instance of the factorizer.
        virtual LUft<IdxT,DataT> * New () const { return new LUft_SUPERLU_DIST<IdxT,DataT>(); }
        
        /// Get name of the factorizer.
        virtual std::string name () const { return "superlu_dist"; }
        
        /// Factorize.
        virtual void factorize (CsrMatrix<IdxT,DataT> const & matrix, LUftData data);
        
        /// Validity indicator.
        virtual bool valid () const { return size_ != 0; }
        
        /// Return LU byte size.
        virtual std::size_t size () const { return size_; }
        
        /// Solve equations.
        virtual void solve (const ArrayView<DataT> b, ArrayView<DataT> x, int eqs) const;
        
        /// Save factorization data to disk.
        virtual void save (std::string name) const;
        
        /// Load factorization data from disk.
        virtual void load (std::string name, bool throw_on_io_failure = true);
        
        /// Release memory.
        virtual void drop ();
        
    private:
        
        /// Matrix that has been factorized.
        NumberArray<int_t> P_;
        NumberArray<int_t> I_;
        NumberArray<std::complex<double>> X_;
        
        // scaling and permutation data
        ScalePermstruct_t ScalePermstruct_;
        
        // factorization data
        LUstruct_t LUstruct_;
        
        // process grid
        gridinfo_t * grid_;
        
        /// Memory size.
        std::size_t size_;
};

// --------------------------------------------------------------------------------- //

#endif // WITH_SUPERLU_DIST
