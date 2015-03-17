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

#ifndef HEX_SYMBANDMATRIX_H
#define HEX_SYMBANDMATRIX_H

#include "arrays.h"
#include "coomatrix.h"
#include "densematrix.h"
#include "hdffile.h"

/**
 * @brief Matrix parts.
 * 
 * This enumeration contains some matrix part identificators
 * that can be used e.g. together with matrix-vector multiplication
 * routines (if implemented). The multiplication is done, then, as if
 * the other parts were zero. For example
   @code
       y = A.dot(x, strict_upper | strict_lower)
   @endcode
 * will multiply vector @c x as if the diagonal of @c A were zero.
 */
typedef enum {
    none         = 0,  // b000,
    strict_lower = 1,  // b001, strict_lower
    strict_upper = 2,  // b010,                           strict_upper
    strict_both  = 3,  // b011, strict_lower |            strict_upper
    diagonal     = 4,  // b100,                diagonal
    lower        = 5,  // b101, strict_lower | diagonal
    upper        = 6,  // b110,                diagonal | strict_upper
    both         = 7   // b111, strict_lower | diagonal | strict_upper
} MatrixTriangle;

/**
 * @brief Symmetric diagonal matrix.
 * 
 * A multi-diagonal matrix is a sparse matrix of the structure
 * @f[
 *     A = \pmatrix {
 *            \ast   & \dots  & \ast   &        &        &        \cr
 *            \vdots & \ddots &        & \ddots &        &        \cr
 *            \ast   &        & \ddots &        & \ddots &        \cr
 *                   & \ddots &        & \ddots &        & \ast   \cr
 *                   &        & \ddots &        & \ddots & \vdots \cr
 *                   &        &        & \ast   & \dots  & \ast   \cr
 *     } \ ,
 * @f]
 * i.e. it is banded, with nonzero elements only near to the diagonal.
 * Some of the diagonals may be identically zero. This class holds all
 * nonzero main and upper diagonals (lower diagonals are not necessary
 * in the symmetric case). The diagonal storage has several advantages:
 * 
 * - Matrix-vector multiplication can be vectorized, in contrast to the
 *   CSR-matrix-vector multiplication which requires a strongly irregular
 *   memory access.
 * - Matrix-matrix multiplication conserves the DIA structure, with just
 *   a simple increase of bandwidth. It holds that the bandwidth of the
 *   matrix C = AB is equal to the sum of the bandwidths of the factors
 *   decreased by one.
 */
class SymBandMatrix
{

public:
    //
    // Constructors
    //

    /// Empty constructor.
    SymBandMatrix ();
    
    /// Size constructor.
    SymBandMatrix (std::size_t n);
    
    /**
     * @brief Data constructor.
     * 
     * @param n Size of the matrix.
     * @param d Number of main and upper diagonals.
     */
    SymBandMatrix (std::size_t n, std::size_t d);
    
    /**
     * @brief Data constructor.
     * 
     * This constructor allows population of the matrix by precomputed data. The
     * data needs to be stored in the internal format, which is concatenation of
     * nonzero elements from the matrix lines (but only main and upper diagonals).
     * The several trailing rows are zero padded to make a rectangular block of
     * numbers.
     * 
     * @param n Size of the matrix.
     * @param d Number of main and upper diagonals.
     * @param v Row-padded (main- and upper-) diagonal elements.
     */
    SymBandMatrix (std::size_t n, std::size_t d, const cArrayView v);
    
    /// Copy constructor.
    SymBandMatrix (SymBandMatrix const & A);

    /// Move constructor.
    SymBandMatrix (SymBandMatrix && A);
    
    /// Constructor - HDF loader.
    SymBandMatrix (std::string filename);
    
    /**
     * @brief Plain symmetrical populator.
     *
     * Given a functor of the signature
       @code
           Complex (*) (int, int);
       @endcode
     * the function will call the functor with row and column number of every
     * element that is to be set.
     * 
     * @param d How many (main and) upper diagonals to populate.
     * @param f The functor that will compute the matrix elements.
     */
    template <class Functor> SymBandMatrix & populate (Functor f)
    {
        // throw away old data
        elems_.resize(n_ * d_);
        
        // evaluate the elements
        for (std::size_t irow = 0; irow < n_; irow++)
        for (std::size_t id = 0; id < d_; id++)
        if (irow + id < n_)
            elems_[irow * d_ + id] = f(irow, irow + id);
        
        return *this;
    }
    
    //
    // Destructor
    //

    ~SymBandMatrix () {}
    
    /// Free all fields, set dimensions to zero.
    void drop ()
    {
        n_ = 0;
        d_ = 0;
        elems_.drop();
    }
    
    //
    // Getters
    //
    
    Complex operator() (int i, int j) const
    {
        std::size_t irow = std::min(i,j);
        std::size_t icol = std::max(i,j);
        std::size_t idia = icol - irow;
        
        assert(irow < n_ and icol < n_);
        
        return idia < d_ ? elems_[irow * d_ + idia] : 0.;
    }
    Complex & operator() (int i, int j)
    {
        std::size_t irow = std::min(i,j);
        std::size_t icol = std::max(i,j);
        std::size_t idia = icol - irow;
        
        assert(irow < n_ and icol < n_ and idia < d_);
        
        return elems_[irow * d_ + idia];
    }
    
    cArray const & data () const { return elems_; }
    cArray       & data ()       { return elems_; }
    
    /**
     * @brief Matrix dimension.
     * 
     * Return row/column count. The matrix is symmetric and so both
     * counts are equal.
     */
    std::size_t   size () const { return n_; }
    std::size_t & size ()       { return n_; }
    
    /**
     * @brief Bandwidth.
     * 
     * Return the bandwidth of the matrix, i.e. number of all (upper, main an lower)
     * diagonals that would have to be stored in a full banded-matrix format.
     */
    std::size_t bandwidth () const { return 2 * d_ - 1; }
    std::size_t halfbw () const { return d_; }
    
    /**
     * @brief Check compatibility of matrices.
     * 
     * Check that the matrix B has the same dimensions as *this matrix and
     * also that they keep the same diagonals. Such matrices can be very effectively
     * summed and subtracted -- just by doing the operation on the stored element arrays.
     */
    bool is_compatible (SymBandMatrix const & B) const;
    
    //
    // Arithmetic and other operators
    //
    
    SymBandMatrix const & operator = (SymBandMatrix && A);
    SymBandMatrix const & operator = (SymBandMatrix const & A);
    
    SymBandMatrix const & operator += (SymBandMatrix const & B);
    SymBandMatrix const & operator -= (SymBandMatrix const & B);
    
    /**
     * @brief Dot product with one or more vectors.
     *
     * @param B Set of vectors to multiply as a column-major dense matrix.
     */
    cArray dot (const cArrayView B) const;
    
    /**
     * @brief Dot product with one or more vectors.
     *
     * @param n Size of the symmetrical matrix.
     * @param d Number of main and upper diagonals.
     * @param M Matrix elements (concatenated row-padded upper triangle rows).
     * @param X Concatenate source vectors.
     */
    static cArray sym_band_dot (int n, int d, const cArrayView M, const cArrayView X);
    
    /**
     * @brief Back-substitution (lower).
     * 
     * Assume the matrix is normalized lower-triangular (i.e. has unit main diagonal
     * and zero upper triangle) and do the triangular solve.
     * 
     * @param b Right hand side of the triangular system.
     */
    cArray lowerSolve (const cArrayView b) const;
    
    /**
     * @brief Back-substitution (upper).
     * 
     * Assume the matrix is normalized upper-triangular (i.e. has unit main diagonal
     * and zero lower triangle) and do the triangular solve.
     * 
     * @param b Right hand side of the triangular system.
     */
    cArray upperSolve (const cArrayView b) const;
    
    /**
     * @brief Kronecker product.
     * 
     * Compute Kronecker product with other matrix.
     */
    SymBandMatrix kron (SymBandMatrix const & B) const;

    //
    // HDF interface
    //
    
    /// Link matrix to a disk file.
    void hdflink (std::string name);
    
    /// Return the name of the linked disk file.
    std::string hdfname () const { return name_; }
    
    /// Return content of the 'name' file as a new SymBandMatrix object.
    SymBandMatrix hdfget () const { return SymBandMatrix(name_); }
    
    /**
     * @brief Load from file.
     * 
     * Load the matrix from a HDF5 file created by the routine @ref hdfsave.
     * 
     * @return True on successful read, false otherwise (mostly when doesn't exist).
     */
    //@{
    bool hdfload () { return hdfload (name_); }
    bool hdfload (std::string name);
    bool hdfload (HDFFile & hdf, std::string prefix = "");
    //@}
    
    /**
     * @brief Save data to file.
     * 
     * Save the matrix to the HDF5 file as a set of four datasets:
     * - idiag - identifiers of (positive) diagonals that contain nonzeros,
     * - n - single number containing size of the matrix (rows or columns)
     * - x - concatenated diagonal elements (interleaved complex values)
     * - zero_blocks - information on the compression of the "x" array, see @ref NumberArray::compress.
     * 
     * @return True on successful write, false otherwise.
     */
    //@{
    bool hdfsave
    (
        std::string name,
        HDFFile::FileAccess flags = HDFFile::overwrite,
        bool docompress = false,
        std::size_t consec = 10
    ) const;
    bool hdfsave
    (
        HDFFile::FileAccess flags = HDFFile::overwrite,
        bool docompress = false,
        std::size_t consec = 10
    ) const
    {
        return hdfsave(name_, flags, docompress, consec);
    }
    //@}
    
    //
    // Conversions to other formats
    //
    
    /**
     * @brief Zero-pad rows.
     * 
     * Pad rows with zeros as below:
     * @f[
     *      \left(
     *           \matrix {
     *               \ast & \ast &      &      &      \cr
     *               \ast & \ast & \ast &      &      \cr
     *               \ast & \ast & \ast & \ast &      \cr
     *               \ast & \ast & \ast & \ast & \ast \cr
     *                    & \ast & \ast & \ast & \ast \cr
     *                    &      & \ast & \ast & \ast \cr
     *                    &      &      & \ast & \ast \cr
     *           }
     *      \right)
     *      \longrightarrow
     *      \matrix {
     *           0    & 0    & 0 \cr
     *                & 0    & 0 \cr
     *                &      & 0 \cr
     *                &      & . \cr
     *                &      & . \cr
     *                &      & . \cr
     *                &      & . \cr
     *      }
     *      \left(
     *           \matrix {
     *               \ast & \ast &      &      &      \cr
     *               \ast & \ast & \ast &      &      \cr
     *               \ast & \ast & \ast & \ast &      \cr
     *               \ast & \ast & \ast & \ast & \ast \cr
     *                    & \ast & \ast & \ast & \ast \cr
     *                    &      & \ast & \ast & \ast \cr
     *                    &      &      & \ast & \ast \cr
     *           }
     *      \right)
     *      \matrix {
     *           .    &      &      \cr
     *           .    &      &      \cr
     *           .    &      &      \cr
     *           .    &      &      \cr
     *           0    &      &      \cr
     *           0    & 0    &      \cr
     *           0    & 0    & 0    \cr
     *      }
     *      \longrightarrow
     *      \matrix {
     *          0    & 0    & 0    & \ast & \ast \cr
     *          0    & 0    & \ast & \ast & \ast \cr
     *          0    & \ast & \ast & \ast & \ast \cr
     *          \ast & \ast & \ast & \ast & \ast \cr
     *          \ast & \ast & \ast & \ast & 0    \cr
     *          \ast & \ast & \ast & 0    & 0    \cr
     *          \ast & \ast & 0    & 0    & 0    \cr
     *      }
     * @f]
     * Then concatenate rows and return as a single array.
     */
    cArray toPaddedRows () const;
    
    /**
     * @brief Zero-pad columns.
     * 
     * Pad rows with zeros as in @ref toPaddedRows, then concatenate
     * columns and return as a single array.
     */
    cArray toPaddedCols () const;
    
    /// Convert matrix part to CooMatrix.
    CooMatrix tocoo (MatrixTriangle triangle = both) const;
    
    /// Convert matrix part to RowMatrix.
    RowMatrix<Complex> torow (MatrixTriangle triangle = both) const;
    
    /// Output to a text stream.
    friend std::ostream & operator << (std::ostream & out, SymBandMatrix const & A);
    
private:

    // dimension (only square matrices allowed)
    std::size_t n_;
    
    // main and upper diagonal count
    std::size_t d_;
    
    // diagonals concatenated in row-oriented way (and padded, if necessary) of size n_*d_
    cArray elems_;
    
    // name of linked HDF file
    std::string name_;
};

class BlockSymBandMatrix
{
    private:
        
        /// Name of HDF5 scratch disk file.
        std::string diskfile_;
        
        /// Whether to keep in memory.
        bool inmemory_;
        
        /// Size of a matrix block and also of the block structure.
        std::size_t size_;
        
        /// Half bandwidth.
        std::size_t halfbw_;
        
        /// Data array.
        cArray data_;
        
    public:
        
        //
        // Constructors.
        //
        
        BlockSymBandMatrix (int size = 0, int halfbw = 0)
            : diskfile_(), inmemory_(true), size_(size), halfbw_(halfbw), data_() {}
        
        /**
         * @brief Main constructor.
         * 
         * @param Nblock Number of blocks in a row (column) of the block matrix.
         * @param blocksize Number of element in a row (column) of every single block.
         * @param blockstructure Vector of block positions; only upper part of the matrix (+ main diagonal) allowed.
         * @param name Name of the optional scratch disk file.
         */
        BlockSymBandMatrix (int size, int halfbw, bool inmemory = true, std::string name = "")
            : diskfile_(name), inmemory_(inmemory), size_(size), halfbw_(halfbw), data_()
        {
            if (inmemory_)
            {
                data_.resize(size_ * size_ * halfbw_ * halfbw_);
            }
        }
        
        /// Is this object cached in memory?
        bool inmemory () const
        {
            return inmemory_;
        }
        
        /// Access individual blocks.
        cArray getBlock (int i) const;
        void setBlock (int i, const cArrayView data);
        
        //
        // Arithmetic operators.
        //
        
        /**
         * @brief Multiply (block) vector by the matrix.
         * 
         * @param v Vector to multiply. It should be of equal size to the size of the matrix.
         *          The result will be again of the same size.
         * @param parallelize Multiply by several blocks at once (OpenMP used).
         */
        cArray dot (cArrayView v, bool parallelize = false) const;
        
        //
        // Coversions to other matrix types.
        //
        
        /**
         * @brief Convert to COO format.
         * 
         * @param loadblocks Use blocks from scratch file instead of those in memory (if any).
         */
        CooMatrix tocoo () const;
        
        //
        // Scratch file I/O.
        //
        
        /// Get name of the scratch disk file.
        std::string   hdfname () const { return diskfile_; }
        std::string & hdfname ()       { return diskfile_; }
        
        /// Release data from memory, but keep HDF link.
        void drop ();
        
        /// Check that the scratch disk file exists.
        bool hdfcheck () const;
        
        /// Reset (= create empty) HDF scratch file needed when writing individual blocks.
        bool hdfinit () const;
        
        /// Write all blocks from memory to the disk file.
        bool hdfsave () const;
        
        /// Load data from disk to memory.
        bool hdfload ();
};

SymBandMatrix operator + (SymBandMatrix const & A, SymBandMatrix const & B);
SymBandMatrix operator - (SymBandMatrix const & A, SymBandMatrix const & B);
SymBandMatrix operator * (SymBandMatrix const & A, SymBandMatrix const & B);
SymBandMatrix operator * (double z, SymBandMatrix const & A);
SymBandMatrix operator * (Complex z, SymBandMatrix const & A);

BlockSymBandMatrix kron (SymBandMatrix const & A, SymBandMatrix const & B);
cArray kron_dot (SymBandMatrix const & A, SymBandMatrix const & B, const cArrayView v);

#endif // HEX_SYMBANDMATRIX_H
