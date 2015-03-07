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

#ifndef HEX_CSRMATRIX_H
#define HEX_CSRMATRIX_H

#ifndef NO_PNG
#include <png++/png.hpp>
#endif

#ifndef NO_UMFPACK
#include <umfpack.h>
#endif

#include "arrays.h"
#include "coomatrix.h"
#include "densematrix.h"

/**
 * @brief Complex CSR matrix.
 * 
 * Complex compressed-storage-row-major-ordered sparse matrix.
 * The data are stored in four arrays,
 * @code
 * p, i, x, z
 * @endcode
 * Here, \c Ap is an array of N+1 elements, where N is the number of rows.
 * It contains indices of beginnings of single rows in all remaining arrays.
 * The last element contains one-past-end pointer, i.e. it should be aqual to
 * the element count of the matrix (zero and nonzero). The array \c Ai
 * holds column indices of elements, which are stored in \c Ax and \c Az, real
 * ang imaginary part being separated.
 */
class CsrMatrix
{
    
public:
    
    // Constructors
    
    CsrMatrix ()
        : m_(0), n_(0), name_() {}
    CsrMatrix (std::size_t m, std::size_t n)
        : m_(m), n_(n), name_() {}
    CsrMatrix (CsrMatrix const & A)
        : m_(A.m_), n_(A.n_), p_(A.p_), i_(A.i_), x_(A.x_), name_() {}
    CsrMatrix (std::size_t m, std::size_t n, lArrayView const & p, lArrayView const & i, cArrayView const & x)
        : m_(m), n_(n), p_(p), i_(i), x_(x), name_() {}
    
    // Destructor
    
    ~CsrMatrix () {}
    
    /**
     * @brief Clear data.
     */
    void drop ()
    {
        m_ = n_ = 0;
        p_.drop();
        i_.drop();
        x_.drop();
    }
    
    /**
     * Ordinary matrix-vector product, \f$ A\cdot b \f$.
     * @param b Vector to multiply with.
     */
    cArray dot (cArrayView const & b) const;
    
    // Getters
    
    std::size_t size () const { return i_.size(); }
    std::size_t rows () const { return m_; }
    std::size_t cols () const { return n_; }
    lArray const & p () const { return p_; }
    lArray const & i () const { return i_; }
    cArray const & x () const { return x_; }
    
    // return absolute value of the (in absolute value) largest element
    double norm () const;
    
    /**
     * Write the matrix data to a file.
     * 
     * @param filename Filename.
     * 
     * Format of the fields i,j,x,z is:
     * @code
     * %d\t%d\t%g\t%g
     * @endcode
     */
    void write (const char* filename) const;
    
#ifndef NO_PNG
    /**
     * PNG row data generator for use in \ref plot function.
     * 
     * @note This class requires PNG++. It can be disabled if the macro
     * NO_PNG is defined.
     */
    class PngGenerator : public png::generator<png::gray_pixel_1,PngGenerator>
    {
        typedef png::generator<png::gray_pixel_1,PngGenerator> base_t;
        typedef png::packed_pixel_row<png::gray_pixel_1> row;
        typedef png::row_traits<row> row_traits;
        
        public:
            
            PngGenerator (CsrMatrix const * mat, double threshold);
            ~PngGenerator ();
            
            png::byte* get_next_row (std::size_t pos);
            
        private:
            
            CsrMatrix const * M_;
            row buff_;
            double threshold_;
    };

    /**
     * Save matrix structure as a black-and-white image.
     * @param filename File name.
     * @param threshold Largest absolute value represented by white colour.
     */
    void plot (const char* filename, double threshold = 0.) const;
#endif

#ifndef NO_UMFPACK
    /**
     * @brief LU factorization object
     * 
     * This class is returned by the function CsrMatrix::factorize() and
     * it provides some functions that can be used when solving equations with
     * that LU factorization. Also, it is possible to store the decomposition
     * to disk (link,save), destro the data (drop) and load later when needed
     * (load). The most important function is "solve".
     */
    class LUft
    {
    public:
        
        /// Default constructor.
        LUft ()
            : numeric_(nullptr), matrix_(nullptr), filename_(), info_(UMFPACK_INFO) {}
        
        /// Copy constructor.
        LUft (LUft const & lu)
            : numeric_(lu.numeric_), matrix_(lu.matrix_), filename_(lu.filename_), info_(lu.info_) {}
        
        /// Move constructor.
        LUft (LUft && lu)
            : numeric_(lu.numeric_), matrix_(lu.matrix_), filename_(lu.filename_), info_(lu.info_) { lu.numeric_ = nullptr; }
        
        /// Initialize the structure using the matrix and its numeric decomposition.
        LUft (const CsrMatrix* matrix, void* numeric)
            : numeric_(numeric), matrix_(matrix), filename_(), info_(UMFPACK_INFO) {}
        
        /// Destructor.
        ~LUft ()
        {
            drop ();
        }
        
        /**
         * @brief Transfer data from another factorization object.
         * 
         * Move contents of another LUft object to this one. If there are already some
         * data in this object, delete them.
         */
        void transfer (LUft && lu)
        {
            // clean this object
            drop();
            
            // transfer data
            numeric_  = lu.numeric_;
            matrix_   = lu.matrix_;
            filename_ = lu.filename_;
            info_     = lu.info_;
            
            // clean source
            lu.numeric_  = nullptr;
        }
        
        /**
         * @brief Size of the numerical data.
         * 
         * Return the number of bytes occupied by the stored elements
         * of the LU-factorization. This doesn't contain any other structural data.
         */
        std::size_t size () const
        {
            if (numeric_ == nullptr)
                return 0;
            
            std::int64_t lnz, unz, m, n, nz_udiag;
            std::int64_t status = umfpack_zl_get_lunz
            (
                &lnz, &unz, &m, &n, &nz_udiag, numeric_
            );
            return status == 0 ? (lnz + unz) * 16 : 0; // Byte count
        }
        
        /**
         * @brief Solve equations.
         * 
         * The parameter "b" is assumed to contain several right hand
         * side vectors (their count is supplied as the optional parameter
         * "eqs"). The results are stored in "x", which has the same size
         * as "b".
         */
        //@{
        cArray solve (const cArrayView b, unsigned eqs = 1) const
        {
            // reserve space for the solution
            cArray x(b.size());
            
            // solve
            solve(b, x, eqs);
            
            // return the result
            return x;
        }
        void solve (const cArrayView b, cArrayView x, int eqs = 1) const
        {
            // check sizes
            assert (eqs * matrix_->n_ == (int)x.size());
            assert (eqs * matrix_->n_ == (int)b.size());
            
            // solve for all RHSs
            for (int eq = 0; eq < eqs; eq++)
            {
                // solve for current RHS
                std::int64_t status = umfpack_zl_solve
                (
                    UMFPACK_Aat,
                    matrix_->p_.data(), matrix_->i_.data(),
                    reinterpret_cast<const double*>(matrix_->x_.data()), nullptr,
                    reinterpret_cast<double*>(&x[0] + eq * matrix_->n_), nullptr,
                    reinterpret_cast<const double*>(&b[0] + eq * matrix_->n_), nullptr,
                    numeric_, nullptr, &info_[0]
                );
                
                // check output
                if (status != UMFPACK_OK)
                {
                    std::cerr << "\n[CsrMatrix::LUft::solve] Exit status " << status << std::endl;
                    umfpack_zl_report_status(0, status);
                }
            }
        }
        //@}
        
        /**
         * @brief Get info array.
         * 
         * Get UMFPACK "info" array.
         */
        rArray const & info () const { return info_; }
        
        /**
         * @brief Link to a disk file.
         * 
         * This function will set a filename that will be used if
         * any of the functions @ref save or @ref load is used without
         * a specific filename.
         */
        void link (std::string name) { filename_ = name; }
        void unlink () { filename_.clear(); }
        
        /**
         * @brief Name of the linked disk file.
         */
        std::string name () const { return filename_; }
        
        /**
         * @brief Save Numeric object to a disk file.
         * 
         * Stores the LU-factorization data in the native UMFPACK format
         * to a disk file.
         */
        void save (std::string name) const
        {
            std::int64_t err = umfpack_zl_save_numeric (numeric_, const_cast<char*>(filename_.c_str()));
            
            if (err == UMFPACK_ERROR_invalid_Numeric_object)
                HexException("[LUft::save] Invalid numeric object.");
            
            if (err == UMFPACK_ERROR_file_IO)
                HexException("[LUft::save] Failed to save LU object \"%s\" (size = %ld).", name.c_str(), size());
        }
        void save () const
        {
            save (filename_);
        }
        
        /**
         * @brief Load Numeric object from a disk file.
         * 
         * The expected format is the format of umfpack_zl_save_numeric.
         */
        //@{
        void load (std::string name, bool throw_on_io_failure = true)
        {
            std::int64_t err = umfpack_zl_load_numeric (&numeric_, const_cast<char*>(filename_.c_str()));
            
            if (err == UMFPACK_ERROR_out_of_memory)
                HexException("[LUft::load] Out of memory.");
            
            if (err == UMFPACK_ERROR_file_IO and throw_on_io_failure)
                HexException("[LUft::save] Failed to load LU object \"%s\".", name.c_str());
        }
        void load ()
        {
            load (filename_, true);
        }
        void silent_load ()
        {
            load (filename_, false);
        }
        //@}
        
        /**
         * @brief Free memory.
         * 
         * Release memory occupied by the LU-factorization numeric object.
         */
        void drop ()
        {
            if (numeric_ != nullptr)
            {
                umfpack_zl_free_numeric (&numeric_);
                numeric_ = nullptr;
            }
        }
        
    private:
        
        /// Numeric decomposition as produced by UMFPACK.
        void* numeric_;
        
        /// Pointer to the matrix that has been factorized. Necessary for validity of @ref numeric_.
        const CsrMatrix* matrix_;
        
        /// Linked HDF file name.
        std::string filename_;
        
        /// Set of status flags produced by UMFPACK.
        mutable rArray info_;
        
        // Disable bitwise copy
        LUft const & operator= (LUft const &);
    };

    /**
     * @brief Compute (incomplete) LU factorization.
     * 
     * This function computes the LU factorization of the matrix. It uses
     * the free UMFPACK library.
     * 
     * @param droptol Drop tolerance. If an element of the factorization matrices
     *                should be in absolute value smaller than "droptol" the
     *                library will omit it completely effectively making the result
     *                more sparse.
     * @return Special structure of the LUft type, holding opaque information about
     *         the factorization.
     */
    LUft factorize (double droptol = 0) const;
#endif
    
    /**
     * @brief Solve the Ax = b problem, where "b" can be a matrix.
     * 
     * Uses UMFPACK through the LUft::solve function.
     * 
     * @param b Complex vector containing column-major ordered data; it may be
     *          a flattened matrix.
     * @param eqs Number of columns.
     * @return Array of roots in the same shape as "b".
     */
    cArray solve (const cArrayView b, size_t eqs = 1) const;
    
    /**
     * @brief Link to a disk file.
     * 
     * Set a default I/O file that will be used if any of the functions
     * @ref hdfload or @ref hdfsave will be used without an explicit filename.
     */
    void hdflink (std::string name);
    void unlink () { name_.clear(); }
    
    /**
     * Get linked HDF name.
     */
    std::string hdfname () const { return name_; }
    
    /**
     * Save matrix to HDF file.
     */
    //@{
    bool hdfsave () const { return hdfsave(name_); }
    bool hdfsave (std::string name) const;
    //@}
    
    /**
     * Load matrix from HDF file.
     */
    //@{
    bool hdfload () { return hdfload(name_); }
    bool hdfload (std::string name);
    //@}
    
    /**
     * Element-wise access (const).
     * If a non-existing element is referenced, zero is returned.
     */
    Complex operator() (unsigned i, unsigned j) const;
    
    /**
     * Multiplication by a number.
     */
    CsrMatrix & operator *= (Complex r);
    
    /**
     * Addition of another CSR matrix.
     * The matrices MUST HAVE THE SAME SPARSE STRUCTURE, as no indices are
     * checked.
     */
    CsrMatrix & operator &= (CsrMatrix const & B);
    
    /**
     * Subtraction of another CSR matrix.
     * The matrices MUST HAVE THE SAME SPARSE STRUCTURE, as no indices are
     * checked.
     */
    CsrMatrix& operator ^= (CsrMatrix const & B);
    
    /**
     * Sets fill-in elements so that the storage structure of this matrix
     * will be identical to that of the other CSR matrix. We can the use
     * “fast” arithmetic operators & and ^.
     * @return self
     */
    CsrMatrix sparse_like (CsrMatrix const & B) const;
    
    /**
     * Return dense array with diagonal elements of the matrix.
     */
    cArray diag() const;
    
    /**
     * Convert to COO format.
     */
    CooMatrix tocoo() const;
    
    /**
     * Convert to dense matrix.
     */
    RowMatrix<Complex> torow() const;
    
    /**
     * Solves upper triangular system of equations using backsubstitution.
     * The matrix needs not be triangular, but elements under the diagonal
     * won't be used or changed.
     * @param b Right-hand side.
     */
    cArray upperSolve (cArrayView const & b) const;
    
    /**
     * Solves lower triangular system of equations using backsubstitution.
     * The matrix needs not be triangular, but elements above the diagonal
     * won't be used or changed.
     * @param b Right-hand side.
     */
    cArray lowerSolve (cArrayView const & b) const;
    
    /**
     * Applies a user transformation on <b>nonzero</b> matrix elements.
     * @param f A functor compatible with following declaration:
     * @code
     * Complex (*f) (size_t i, size_t j, Complex z);
     * @endcode
     */
    template <class Functor> CsrMatrix nzTransform (Functor f) const
    {
        // create output matrix
        CsrMatrix A = *this;
        
        // loop over rows
        for (size_t row = 0; row < rows(); row++)
        {
            // get relevant columns
            size_t idx1 = p_[row];
            size_t idx2 = p_[row + 1];
            
            // loop over columns
            for (size_t idx = idx1; idx < idx2; idx++)
            {
                // which column is this?
                size_t col = i_[idx];
                
                // transform the output matrix
                A.x_[idx] = f(row, col, x_[idx]);
            }
        }
        
        return A;
    }
    
private:
    
    // dimensions
    std::int64_t m_;
    std::int64_t n_;
    
    // representation
    lArray p_;
    lArray i_;
    cArray x_;
    
    // linked HDF file
    std::string name_;
};

/**
 * Computes a sum of two csr-matrices OF THE SAME SPARSE STRUCTURE.
 */
inline CsrMatrix operator & (CsrMatrix const & A, CsrMatrix const & B)
{
    CsrMatrix C = A;
    return C &= B;
}

/**
 * Computes a difference of two csr-matrices OF THE SAME SPARSE STRUCTURE.
 */
inline CsrMatrix operator ^ (CsrMatrix const & A, CsrMatrix const & B)
{
    CsrMatrix C = A;
    return C ^= B;
}

/**
 * Multiplication of csr-matrix by a number.
 */
inline CsrMatrix operator * (Complex z, CsrMatrix const &  B)
{
    CsrMatrix C = B;
    return C *= z;
}

/**
 * Multiplication of csr-matrix by a number.
 */
inline CsrMatrix operator * (CsrMatrix const & A, double r)
{
    CsrMatrix C = A;
    return C *= r;
}

#endif // HEX_CSRMATRIX_H
