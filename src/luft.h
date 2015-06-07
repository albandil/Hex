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

#ifdef WITH_UMFPACK
#include <umfpack.h>
#endif

/// LU factorization methods.
enum
{
    LUFT_ANY,
    LUFT_UMFPACK,
    LUFT_SUPERLU,
    LUFT_SUPERLU_DIST
};

// Forward declaration of classes (needed by csrmatrix.h).
template <class IdxT, class DataT> class LUft;
template <class IdxT, class DataT> class LUft_UMFPACK;
template <class IdxT, class DataT> class LUft_SUPERLU;
template <class IdxT, class DataT> class LUft_SUPERLU_DIST;

// Load available matrices.
#include "matrix.h"

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
        
        /// Destructor.
        virtual ~LUft () {}
        
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
            std::cout << "b = " << b << std::endl;
            
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
            std::cout << "Warning: I'm placeholder LUft::load, not loading anything!" << std::endl;
        }
        virtual void load () { this->load(filename_, true); }
        virtual void silent_load () { this->load(filename_, false); }
        
    private:
        
        /// Name of the disk file.
        std::string filename_;
};

#ifdef WITH_UMFPACK
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
            : LUft<IdxT,DataT>(), numeric_(nullptr), matrix_(nullptr), info_(UMFPACK_INFO) {}
        
        /// Initialize the structure using the matrix and its numeric decomposition.
        LUft_UMFPACK (CsrMatrix<IdxT,DataT> const * matrix, void * numeric)
            : LUft<IdxT,DataT>(), numeric_(numeric), matrix_(matrix), info_(UMFPACK_INFO) {}
        
        /// Destructor.
        virtual ~LUft_UMFPACK () { drop (); }
        
        /// Return factorization information.
        rArray const & info () const { return info_; }
        
        /// Return LU byte size.
        virtual std::size_t size () const;
        
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
        
        /// Pointer to the matrix that has been factorized. Necessary for validity of @ref numeric_.
        CsrMatrix<IdxT,DataT> const * matrix_;
        
        /// Set of status flags produced by UMFPACK.
        mutable rArray info_;
        
        // Disable bitwise copy
        LUft_UMFPACK const & operator= (LUft_UMFPACK const &);
};
#endif // WITH_UMFPACK

#ifdef WITH_SUPERLU
#include <slu_zdefs.h>
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
            : LUft<IdxT,DataT>(), matrix_(nullptr), size_(0) {}
        
        /// Initialize the structure using the matrix and its numeric decomposition.
        LUft_SUPERLU
        (
            CsrMatrix<IdxT,DataT> const * matrix,
            iArray const & perm_c,
            iArray const & perm_r,
            iArray const & etree,
            char equed,
            rArray const & R,
            rArray const & C,
            SuperMatrix L,
            SuperMatrix U,
            int bytes,
            double droptol
        ) : LUft<IdxT,DataT>(), matrix_(matrix), perm_c_(perm_c), perm_r_(perm_r), etree_(etree), equed_(equed),
            R_(R), C_(R), L_(L), U_(U), size_(bytes), droptol_(droptol)
        {
            // nothing to do
        }
        
        /// Destructor.
        virtual ~LUft_SUPERLU () { drop (); }
        
        /// Return LU byte size.
        virtual std::size_t size () const { return size_; }
        
        /// Solve equations.
        virtual void solve (const ArrayView<DataT> b, ArrayView<DataT> x, int eqs) const;
        
        /// Save factorization data to disk.
        virtual void save (std::string name) const { /* TODO */ }
        
        /// Load factorization data from disk.
        virtual void load (std::string name, bool throw_on_io_failure = true) { /* TODO */ }
        
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
        
        /// Pointer to the matrix that has been factorized. Necessary for validity of @ref numeric_.
        CsrMatrix<IdxT,DataT> const * matrix_;
        
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
        
        /// Memory size.
        std::size_t size_;
        
        /// Drop tolerance.
        double droptol_;
        
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
            CsrMatrix<IdxT,DataT> const * matrix,
            ScalePermstruct_t ScalePermstruct,
            LUstruct_t LUstruct,
            gridinfo_t * grid,
            std::size_t bytes
        )
            : LUft<IdxT,DataT>(), ScalePermstruct_(ScalePermstruct), LUstruct_(LUstruct),
              grid_(grid), size_(bytes)
        {
            // nothing to do
        }
        
        /// Destructor.
        virtual ~LUft_SUPERLU_DIST () { drop (); }
        
        /// Return LU byte size.
        virtual std::size_t size () const { return size_; }
        
        /// Solve equations.
        virtual void solve (const ArrayView<DataT> b, ArrayView<DataT> x, int eqs) const;
        
        /// Save factorization data to disk.
        virtual void save (std::string name) const { /* TODO */ }
        
        /// Load factorization data from disk.
        virtual void load (std::string name, bool throw_on_io_failure = true) { /* TODO */ }
        
        /// Release memory.
        virtual void drop ()
        {
            // TODO
        }
        
    private:
        
        /// Pointer to the matrix that has been factorized. Necessary for validity of @ref numeric_.
        CsrMatrix<IdxT,DataT> const * matrix_;
        
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

#endif // HEX_LUFT_H
