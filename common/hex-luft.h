//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2015, Jakub Benda, Charles University in Prague                    //
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

#ifndef HEX_LUFT_H
#define HEX_LUFT_H

/// LU factorization methods.
enum
{
    LUFT_ANY,
    LUFT_UMFPACK,
    LUFT_SUPERLU,
    LUFT_SUPERLU_DIST,
    LUFT_MUMPS
};

// Forward declaration of classes.
template <class IdxT, class DataT> class LUft;
template <class IdxT, class DataT> class LUft_UMFPACK;
template <class IdxT, class DataT> class LUft_SUPERLU;
template <class IdxT, class DataT> class LUft_SUPERLU_DIST;
template <class IdxT, class DataT> class LUft_MUMPS;

// Load available matrices.
#include "hex-arrays.h"
#include "hex-matrix.h"
#include "hex-csrmatrix.h"
#include "hex-densematrix.h"

/**
 * @brief LU factorization object
 * 
 * This class is returned by the function CsrMatrix::factorize() and
 * it provides some functions that can be used when solving equations with
 * that LU factorization. Also, it is possible to store the decomposition
 * to disk (link, save), destroy the data (drop) and load later when needed
 * (load). The most important function is "solve".
 */
template <class IdxT, class DataT>
class LUft
{
    public:

        /// Default constructor.
        LUft () : filename_() {}

        /// Get new empty factorization object.
        static LUft * New (int factorizer);

        /// Destructor.
        virtual ~LUft () {}

        /**
         * @brief Validity indicator.
         * 
         * Returns true when the object contains a valid LU factorization.
         */
        virtual bool valid () const { return false; }

        /**
         * @brief Free memory.
         * 
         * Release memory occupied by the LU-factorization numeric object.
         */
        virtual void drop () {}

        /**
         * @brief Size of the numerical data.
         * 
         * Return the number of bytes occupied by the stored elements
         * of the LU-factorization. This doesn't contain any other structural data.
         */
        virtual std::size_t size () const { return 0; }

        /**
         * @brief Estimation of the condition number.
         * 
         * Note: Currently implemented only for UMFPACK backend.
         */
        virtual Real cond () const { return 0; }

        /**
         * @brief Solve equations.
         * 
         * The parameter "b" is assumed to contain several right hand
         * side vectors (their count is supplied as the optional parameter
         * "eqs"). The results are stored in "x", which has the same size
         * as "b".
         */
        virtual void solve (const ArrayView<DataT> b, ArrayView<DataT> x, int eqs) const
        {
            std::cout << "Warning: I'm placeholder LUft::solve, not solving anything!" << std::endl;
        }

        /**
         * @brief Solve equations.
         * 
         * The parameter "b" is assumed to contain several right hand
         * side vectors (their count is supplied as the optional parameter
         * "eqs").
         */
        virtual NumberArray<DataT> solve (const ArrayView<DataT> b, int eqs = 1) const
        {
            // reserve space for the solution
            NumberArray<DataT> x (b.size());

            // solve
            this->solve(b, x, eqs);

            // return the result
            return x;
        }

        /**
         * @brief Link to a disk file.
         * 
         * This function will set a filename that will be used if
         * any of the functions @ref save or @ref load is used without
         * a specific filename.
         */
        virtual void link (std::string name) { filename_ = name; }

        /**
         * @brief Unlink disk file.
         */
        virtual void unlink () { filename_.clear(); }

        /**
         * @brief Name of the linked disk file.
         */
        virtual std::string name () const { return filename_; }

        /**
         * @brief Save factorization object to a disk file.
         * 
         * Stores the LU-factorization data to a disk file in the native format
         * of the library used.
         */
        virtual void save (std::string name) const
        {
            std::cout << "Warning: I'm placeholder LUft::save, not saving anything!" << std::endl;
        }
        virtual void save () const { this->save(filename_); }

        /**
         * @brief Load factorization object from a disk file.
         */
        virtual void load (std::string name, bool throw_on_io_failure = true)
        {
            if (throw_on_io_failure)
                HexException("I'm placeholder LUft::load, not loading anything!");
        }
        virtual void load () { this->load(filename_, true); }
        virtual void silent_load () { this->load(filename_, false); }

    private:

        /// Name of the disk file.
        std::string filename_;
};

#ifdef WITH_UMFPACK
#include <umfpack.h>

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

        /// Initialize the structure using the matrix and its numeric decomposition.
        LUft_UMFPACK (CsrMatrix<IdxT,DataT> const & matrix, void * numeric)
            : LUft<IdxT,DataT>(), numeric_(numeric), matrix_(matrix), info_(UMFPACK_INFO) {}

        /// Destructor.
        virtual ~LUft_UMFPACK () { drop(); }

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
        mutable CsrMatrix<IdxT,DataT> matrix_;

    public:

        /// Set of status flags produced by UMFPACK.
        mutable rArray info_;

    private:

        // Disable bitwise copy
        LUft_UMFPACK const & operator= (LUft_UMFPACK const &);
};
#endif // WITH_UMFPACK

#ifdef WITH_SUPERLU
#ifdef SINGLE
    #include <slu_cdefs.h>
#else
    #include <slu_zdefs.h>
#endif

/**
 * @brief LU factorization object - SuperLU specialization.
 * 
 * This class holds information on LU factorization as computed by the free
 * library SuperLU (sequential version). It is derived from LUft and
 * shares interface with that class.
 */
template <class IdxT, class DataT>
class LUft_SUPERLU : public LUft<IdxT,DataT>
{
    public:

        /// Default constructor.
        LUft_SUPERLU ()
            : LUft<IdxT,DataT>(), matrix_(), size_(0) {}

        /// Initialize the structure using the matrix and its numeric decomposition.
        LUft_SUPERLU
        (
            CsrMatrix<IdxT,DataT> const & matrix,
            iArray const & perm_c,
            iArray const & perm_r,
            iArray const & etree,
            char equed,
            rArray const & R,
            rArray const & C,
            SuperMatrix L,
            SuperMatrix U,
            GlobalLU_t Glu,
            int bytes,
            Real droptol
        ) : LUft<IdxT,DataT>(), matrix_(matrix), perm_c_(perm_c), perm_r_(perm_r), etree_(etree), equed_(equed),
            R_(R), C_(R), L_(L), U_(U), Glu_(Glu), size_(bytes), droptol_(droptol)
        {
            // nothing to do
        }

        /// Destructor.
        virtual ~LUft_SUPERLU () { drop(); }

        /// Validity indicator.
        virtual bool valid () const { return matrix_.size() != 0; }

        /// Return LU byte size.
        virtual std::size_t size () const { return size_; }

        /// Solve equations.
        virtual void solve (const ArrayView<DataT> b, ArrayView<DataT> x, int eqs) const;

        /// Save factorization data to disk.
        virtual void save (std::string name) const { /*HexException("SuperLU factorizer does not yet support --out-of-core option.");*/ }

        /// Load factorization data from disk.
        virtual void load (std::string name, bool throw_on_io_failure = true) { if (throw_on_io_failure) HexException("SuperLU factorizer does not yet support --out-of-core option."); }

        /// Release memory.
        virtual void drop ()
        {
            if (size_ != 0)
            {
                Destroy_SuperNode_Matrix(&L_);
                Destroy_CompCol_Matrix(&U_);
                size_ = 0;
            }
        }

    private:

        /// Matrix that has been factorized.
        CsrMatrix<IdxT,DataT> matrix_;

        /// Row permutations.
        iArray perm_c_;

        /// Column permutations.
        iArray perm_r_;

        /// Elimitation tree.
        iArray etree_;

        /// Equilibration done.
        char equed_;

        /// Row scale factors.
        rArray R_;

        /// Column scale factors.
        rArray C_;

        /// L-factor.
        SuperMatrix L_;

        /// U-factor.
        SuperMatrix U_;

        /// Reusable information.
        GlobalLU_t Glu_;

        /// Memory size.
        std::size_t size_;

        /// Drop tolerance.
        Real droptol_;

        // Disable bitwise copy
        LUft_SUPERLU const & operator= (LUft_SUPERLU const &);
};
#endif // WITH_SUPERLU

#ifdef WITH_SUPERLU_DIST
#include <superlu_zdefs.h>

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
        LUft_SUPERLU_DIST ()
            : LUft<IdxT,DataT>() {}

        /// Initialize the structure using the matrix and its numeric decomposition.
        LUft_SUPERLU_DIST
        (
            CsrMatrix<IdxT,DataT> const & matrix,
            ScalePermstruct_t ScalePermstruct,
            LUstruct_t LUstruct,
            gridinfo_t * grid,
            std::size_t bytes
        )
            : LUft<IdxT,DataT>(), matrix_(matrix), ScalePermstruct_(ScalePermstruct), LUstruct_(LUstruct),
              grid_(grid), size_(bytes)
        {
            // nothing to do
        }

        /// Destructor.
        virtual ~LUft_SUPERLU_DIST () { drop (); }

        /// Validity indicator.
        virtual bool valid () const;

        /// Return LU byte size.
        virtual std::size_t size () const { return size_; }

        /// Solve equations.
        virtual void solve (const ArrayView<DataT> b, ArrayView<DataT> x, int eqs) const;

        /// Save factorization data to disk.
        virtual void save (std::string name) const { HexException("SuperLU_dist factorizer does not yet support --out-of-core option."); }

        /// Load factorization data from disk.
        virtual void load (std::string name, bool throw_on_io_failure = true) { HexException("SuperLU_dist factorizer does not yet support --out-of-core option."); }

        /// Release memory.
        virtual void drop ()
        {
            if (size_ != 0)
            {
                Destroy_LU(matrix_->cols(), grid_, &LUstruct_);
                ScalePermstructFree(&ScalePermstruct_);
                LUstructFree(&LUstruct_);
                matrix_ = nullptr;
                size_ = 0;
            }
        }

    private:

        /// Matrix that has been factorized.
        CsrMatrix<IdxT,DataT> matrix_;

        // scaling and permutation data
        ScalePermstruct_t ScalePermstruct_;

        // factorization data
        LUstruct_t LUstruct_;

        // process grid
        gridinfo_t * grid_;

        /// Memory size.
        std::size_t size_;

        // Disable bitwise copy
        LUft_SUPERLU_DIST const & operator= (LUft_SUPERLU_DIST const &);
};
#endif // WITH_SUPERLU_DIST

#ifdef WITH_MUMPS

#define ICNTL(x) icntl[(x)-1]
#define MUMPS_INITIALIZE    (-1)
#define MUMPS_FINISH        (-2)
#define MUMPS_ANALYZE       1
#define MUMPS_FACTORIZE     2
#define MUMPS_SOLVE         3

#ifdef SINGLE
    #include <cmumps_c.h>
    #define MUMPS_STRUC_C CMUMPS_STRUC_C
    #define MUMPS_C cmumps_c
    #define MUMPS_COMPLEX mumps_float_complex
#else
    #include <zmumps_c.h>
    #define MUMPS_STRUC_C ZMUMPS_STRUC_C
    #define MUMPS_C zmumps_c
    #define MUMPS_COMPLEX mumps_double_complex
#endif

/**
 * @brief LU factorization object - MUMPS specialization.
 * 
 * This class holds information on LU factorization as computed by the free
 * library MUMPS. It is derived from LUft and shares interface with that class.
 * 
 * \warning This class expect only symmetric matrices. The COO matrix passed
 * to this class should contain only upper or lower part of the matrix (together
 * with the main diagonal).
 */
template <class IdxT, class DataT>
class LUft_MUMPS : public LUft<IdxT,DataT>
{
    public:

        /// Default constructor.
        LUft_MUMPS ()
            : LUft<IdxT,DataT> (), settings (nullptr) {}

        /// Construct from data.
        LUft_MUMPS (MUMPS_STRUC_C * s, NumberArray<MUMPS_INT> && i, NumberArray<MUMPS_INT> && j, NumberArray<DataT> && a)
            : LUft<IdxT,DataT>(), settings(s), I(std::move(i)), J(std::move(j)), A(std::move(a)) {}

        /// Destructor.
        virtual ~LUft_MUMPS ();

        /// Validity indicator.
        virtual bool valid () const
        {
            return settings != nullptr and mmin(I.size(), J.size(), A.size()) > 0;
        }

        /// Return LU byte size.
        virtual std::size_t size () const;

        /// Condition number.
        virtual Real cond () const
        {
            #define RINFO(x) rinfo[x-1]
            return settings ? settings->RINFO(11) : 0.0_r;
        }

        /// Solve equations.
        virtual void solve (const ArrayView<DataT> b, ArrayView<DataT> x, int eqs) const;

        /// Save large data to disk.
        virtual void save (std::string name) const
        {
            I.hdfsave("I-" + name);
            J.hdfsave("J-" + name);
            A.hdfsave("A-" + name);
        }

        /// Load large data from disk.
        virtual void load (std::string name, bool throw_on_io_failure = true)
        {
            if (not I.hdfload("I-" + name) or
                not J.hdfload("J-" + name) or
                not A.hdfload("A-" + name))
            {
                if (throw_on_io_failure)
                    HexException("Failed to load MUMPS IJV matrices.");
            }
        }

        /// Release memory.
        virtual void drop ()
        {
            I.drop();
            J.drop();
            A.drop();
        }

    private:

        // Internal data of the library.
        mutable MUMPS_STRUC_C * settings;

        // data arrays
        NumberArray<MUMPS_INT> I, J;
        NumberArray<DataT> A;

        // Disable bitwise copy
        LUft_MUMPS const & operator= (LUft_MUMPS const &);
};
#endif // WITH_MUMPS

#endif // HEX_LUFT_H
