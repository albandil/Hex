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

#include <umfpack.h>

// --------------------------------------------------------------------------------- //

#include "hex-csrmatrix.h"
#include "luft.h"

// --------------------------------------------------------------------------------- //

#ifdef _LONGINT
    #define UMFPACK_DEFAULTS_F          umfpack_zl_defaults
    #define UMFPACK_SYMBOLIC_F          umfpack_zl_symbolic
    #define UMFPACK_FREE_SYMBOLIC_F     umfpack_zl_free_symbolic
    #define UMFPACK_NUMERIC_F           umfpack_zl_numeric
    #define UMFPACK_SAVE_NUMERIC_F      umfpack_zl_save_numeric
    #define UMFPACK_LOAD_NUMERIC_F      umfpack_zl_load_numeric
    #define UMFPACK_FREE_NUMERIC_F      umfpack_zl_free_numeric
    #define UMFPACK_SOLVE_F             umfpack_zl_solve
    #define UMFPACK_GET_LUNZ_F          umfpack_zl_get_lunz
    #define UMFPACK_REPORT_STATUS_F     umfpack_zl_report_status
    #define UMFPACK_COL_TO_TRIPLET_F    umfpack_zl_col_to_triplet
    #define UMFPACK_TRIPLET_TO_COL_F    umfpack_zl_triplet_to_col
#else
    #define UMFPACK_DEFAULTS_F          umfpack_zi_defaults
    #define UMFPACK_SYMBOLIC_F          umfpack_zi_symbolic
    #define UMFPACK_FREE_SYMBOLIC_F     umfpack_zi_free_symbolic
    #define UMFPACK_NUMERIC_F           umfpack_zi_numeric
    #define UMFPACK_SAVE_NUMERIC_F      umfpack_zi_save_numeric
    #define UMFPACK_LOAD_NUMERIC_F      umfpack_zi_load_numeric
    #define UMFPACK_FREE_NUMERIC_F      umfpack_zi_free_numeric
    #define UMFPACK_SOLVE_F             umfpack_zi_solve
    #define UMFPACK_GET_LUNZ_F          umfpack_zi_get_lunz
    #define UMFPACK_REPORT_STATUS_F     umfpack_zi_report_status
    #define UMFPACK_COL_TO_TRIPLET_F    umfpack_zi_col_to_triplet
    #define UMFPACK_TRIPLET_TO_COL_F    umfpack_zi_triplet_to_col
#endif

// --------------------------------------------------------------------------------- //

/**
 * @brief LU factorization object - UMFPACK specialization.
 * 
 * This class holds information on LU factorization as computed by the free
 * library UMFPACK (part of SuiteSparse toolkit). It is derived from LUft and
 * shares interface with that class.
 */
template <class IdxT, class DataT>
class LUft_UMFPACK : public LUft<IdxT,DataT>
{
    public:
    
        /// Default constructor.
        LUft_UMFPACK ()
            : LUft<IdxT,DataT>(), numeric_(nullptr), info_(UMFPACK_INFO) {}
        
        /// Destructor.
        virtual ~LUft_UMFPACK () { drop(); }
        
        // Disable bitwise copy
        LUft_UMFPACK const & operator= (LUft_UMFPACK const &) = delete;
        
        /// New instance of the factorizer.
        virtual LUft<IdxT,DataT> * New () const { return new LUft_UMFPACK<IdxT,DataT>(); }
        
        /// Get name of the factorizer.
        virtual std::string name () const { return "umfpack"; }
        
        /// Factorize.
        virtual void factorize (CsrMatrix<IdxT,DataT> const & matrix, LUftData data);
        
        /// Return factorization information.
        rArray const & info () const { return info_; }
        
        /// Validity indicator.
        virtual bool valid () const;
        
        /// Return LU byte size.
        virtual std::size_t size () const;
        
        /// Return condition number.
        virtual Real cond () const;
        
        /// Solve equations.
        virtual void solve (const ArrayView<DataT> b, ArrayView<DataT> x, int eqs) const;
        
        /// Save factorization data to disk.
        virtual void save (std::string name) const;
        
        /// Load factorization data from disk.
        virtual void load (std::string name, bool throw_on_io_failure = true);
        
        /// Release memory.
        virtual void drop ();
        
    private:
        
        /// Numeric decomposition as produced by UMFPACK.
        void * numeric_;
        
        /// Matrix data, needed for solution.
        NumberArray<LU_int_t> p_;
        NumberArray<LU_int_t> i_;
        NumberArray<std::complex<double>> x_;
        
    public:
        
        /// Set of status flags produced by UMFPACK.
        mutable rArray info_;
};

#endif // WITH_UMFPACK
