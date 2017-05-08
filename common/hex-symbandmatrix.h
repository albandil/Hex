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

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"
#include "hex-coomatrix.h"
#include "hex-densematrix.h"
#include "hex-hdffile.h"
#include "hex-openmp.h"

// --------------------------------------------------------------------------------- //

/**
 * @brief Matrix parts.
 * 
 * This enumeration contains some matrix part identificators
 * that can be used e.g. together with matrix-vector multiplication
 * routines (if implemented). The multiplication is done, then, as if
 * the other parts were zero. For example
   @code
       y = A.dot(x, MatrixSelection::StrictUpper | MatrixSelection::StrictLower);
   @endcode
 * will multiply vector @c x as if the diagonal of @c A were zero.
 */
class MatrixSelection
{
    public:
        
        typedef int Selection;
        
        static const Selection None        = 0; // = 0b0000
        static const Selection StrictLower = 1; // = 0b0001
        static const Selection StrictUpper = 2; // = 0b0010
        static const Selection StrictBoth  = 3; // = 0b0011 = StrictLower | StrictUpper
        static const Selection Diagonal    = 4; // = 0b0100
        static const Selection Lower       = 5; // = 0b0101 = StrictLower | Diagonal
        static const Selection Upper       = 6; // = 0b0110 = StrictUpper | Diagonal
        static const Selection Both        = 7; // = 0b0111 = StrictLower | StrictUpper | Diagonal
        static const Selection BlockStrictLower  = StrictLower << 3; // = 0b0001000
        static const Selection BlockStrictUpper  = StrictLower << 3; // = 0b0010000
        static const Selection BlockStrictBoth   = StrictBoth  << 3; // = 0b0011000
        static const Selection BlockDiagonal     = Diagonal    << 3; // = 0b0100000
        static const Selection BlockLower        = Lower       << 3; // = 0b0101000
        static const Selection BlockUpper        = Upper       << 3; // = 0b0110000
        static const Selection BlockBoth         = Both        << 3; // = 0b0111000
};

// --------------------------------------------------------------------------------- //

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
template <class DataT> class SymBandMatrix
{

public:
    //
    // Constructors
    //

    /// Empty constructor.
    SymBandMatrix () : n_(0), d_(0), elems_(0), name_() {}
    
    /// Size constructor.
    SymBandMatrix (std::size_t n);
    
    /**
     * @brief Data constructor.
     * 
     * @param n Size of the matrix.
     * @param d Number of main and upper diagonals.
     */
    SymBandMatrix (std::size_t n, std::size_t d) : n_(n), d_(d), elems_(n * d), name_() {}
    
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
    SymBandMatrix (std::size_t n, std::size_t d, const ArrayView<DataT> v) : n_(n), d_(d), elems_(v), name_() {}
    
    /// Copy constructor.
    SymBandMatrix (SymBandMatrix const & A) : n_(A.n_), d_(A.d_), elems_(A.elems_), name_() {}

    /// Move constructor.
    SymBandMatrix (SymBandMatrix && A) : n_(std::move(A.n_)), d_(std::move(A.d_)), elems_(std::move(A.elems_)), name_(A.name_) {}
    
    /// Constructor - HDF loader.
    SymBandMatrix (std::string filename) : name_(filename)
    {
        if (not hdfload())
            HexException("Unable to load from %s.", filename.c_str());
    }
    
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
    template <class Functor> SymBandMatrix<DataT> & populate (Functor f, bool parallelize = true)
    {
        // throw away old data
        elems_.resize(n_ * d_);
        
        // evaluate the elements
        # pragma omp parallel for collapse (2) if (parallelize)
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
    
    DataT operator() (int i, int j) const
    {
        std::size_t irow = std::min(i,j);
        std::size_t icol = std::max(i,j);
        std::size_t idia = icol - irow;
        
        assert(irow < n_ and icol < n_);
        
        return idia < d_ ? elems_[irow * d_ + idia] : 0.;
    }
    DataT & operator() (int i, int j)
    {
        std::size_t irow = std::min(i,j);
        std::size_t icol = std::max(i,j);
        std::size_t idia = icol - irow;
        
        assert(irow < n_ and icol < n_ and idia < d_);
        
        return elems_[irow * d_ + idia];
    }
    
    SymBandMatrix<DataT> operator- () const
    {
        return SymBandMatrix<DataT>(n_, d_, -elems_);
    }
    
    NumberArray<DataT> const & data () const { return elems_; }
    NumberArray<DataT>       & data ()       { return elems_; }
    
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
    bool is_compatible (SymBandMatrix const & B) const
    {
        if (n_ != B.n_)
            HexException("Unequal ranks (%d != %d).", n_, B.n_);
        
        if (d_ != B.d_)
            HexException("Unequal half-bandwidths (%d != %d).", d_, B.d_);
        
        return true;
    }
    
    //
    // Arithmetic and other operators
    //
    
    SymBandMatrix const & operator = (SymBandMatrix && A)
    {
        n_ = std::move(A.n_);
        d_ = std::move(A.d_);
        elems_ = std::move(A.elems_);
        return *this;
    }

    SymBandMatrix const & operator = (SymBandMatrix const & A)
    {
        n_ = A.n_;
        d_ = A.d_;
        elems_ = A.elems_;
        return *this;
    }
    
    SymBandMatrix const & operator += (SymBandMatrix const & B)
    {
        is_compatible(B);
        elems_ += B.elems_;
        return *this;
    }
    
    SymBandMatrix const & operator -= (SymBandMatrix const & B)
    {
        is_compatible(B);
        elems_ -= B.elems_;
        return *this;
    }
    
    /**
     * @brief Dot product with one or more vectors.
     *
     * Calculates
     * \f[
     *     A M A + b B \rightarrow B
     * \f]
     * 
     * @param B Set of vectors to multiply as a column-major dense matrix.
     */
    void dot (DataT a, const ArrayView<DataT> A, DataT b, ArrayView<DataT> B, MatrixSelection::Selection tri = MatrixSelection::Both) const
    {
        sym_band_dot(n_, d_, elems_, a, A, b, B, tri);
    }
    NumberArray<DataT> dot (const ArrayView<DataT> A, MatrixSelection::Selection tri = MatrixSelection::Both) const
    {
        NumberArray<DataT> B (A.size());
        sym_band_dot(n_, d_, elems_, 1., A, 0., B, tri);
        return B;
    }
    
    /**
     * @brief Dot product with one or more vectors.
     *
     * @param n Size of the symmetrical matrix.
     * @param d Number of main and upper diagonals.
     * @param M Matrix elements (concatenated row-padded upper triangle rows).
     * @param X Concatenate source vectors.
     */
    static void sym_band_dot
    (
        int n, int d, const ArrayView<DataT> M,
        DataT a, const ArrayView<DataT> A,
        DataT b,       ArrayView<DataT> B,
        MatrixSelection::Selection tri = MatrixSelection::Both
    );
    
    /**
     * @brief Back-substitution (lower).
     * 
     * Assume the matrix is normalized lower-triangular (i.e. has unit main diagonal
     * and zero upper triangle) and do the triangular solve.
     * 
     * @param b Right hand side of the triangular system.
     */
    NumberArray<DataT> lowerSolve (const ArrayView<DataT> b) const;
    
    /**
     * @brief Back-substitution (upper).
     * 
     * Assume the matrix is normalized upper-triangular (i.e. has unit main diagonal
     * and zero lower triangle) and do the triangular solve.
     * 
     * @param b Right hand side of the triangular system.
     */
    NumberArray<DataT> upperSolve (const ArrayView<DataT> b) const;
    
    //
    // HDF interface
    //
    
    /// Link matrix to a disk file.
    void hdflink (std::string name) { name_ = name; }
    
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
    bool hdfload (std::string name)
    {
        // open the HDF file
        HDFFile hdf(name, HDFFile::readonly);
        if (not hdf.valid())
            return false;
        
        return hdfload(hdf);
    }
    bool hdfload (HDFFile & hdf, std::string prefix = "")
    {
        // check the HDF file
        if (not hdf.valid())
            return false;
        
        // set prefix
        if (not prefix.empty())
            hdf.prefix = prefix;
        
        // read dimension
        if (not hdf.read("n", &n_, 1))
            return false;
        
        // read dimension
        if (not hdf.read("d", &d_, 1))
            return false;
        
        // compressed array info
        iArray zero_blocks_re, zero_blocks_im;
        rArray elements_re, elements_im;
        
        // read data from file
        if ((zero_blocks_re.resize(hdf.size("zero_blocks_re")) and not hdf.read("zero_blocks_re", &(zero_blocks_re[0]), zero_blocks_re.size())) or
            (zero_blocks_im.resize(hdf.size("zero_blocks_im")) and not hdf.read("zero_blocks_im", &(zero_blocks_im[0]), zero_blocks_im.size())) or
            (elements_re.resize(hdf.size("re")) and not hdf.read("re", &(elements_re[0]), elements_re.size())) or
            (elements_im.resize(hdf.size("im")) and not hdf.read("im", &(elements_im[0]), elements_im.size())))
            return false;
        
        // decompress
        elems_ = std::move
        (
            interleave
            (
                elements_re.decompress(zero_blocks_re),
                elements_im.decompress(zero_blocks_im)
            )
        );
        
        return true;
    }
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
    ) const
    {
        HDFFile hdf(name, flags);
        
        if (not hdf.valid())
            return false;
        
        // write dimension and diagonal info
        if (not hdf.write("n", &n_, 1) or
            not hdf.write("d", &d_, 1))
            return false;
        
        // compress elements array
        iArray zero_blocks_re, zero_blocks_im;
        rArray elements_re, elements_im;
        
        if (docompress)
        {
            std::tie(zero_blocks_re, elements_re) = realpart(elems_).compress(consec);
            std::tie(zero_blocks_im, elements_im) = imagpart(elems_).compress(consec);
        }
        else
        {
            elements_re = std::move(realpart(elems_));
            elements_im = std::move(imagpart(elems_));
        }
        
        // write compressed elements array (if non-empty); check result
        if ((not zero_blocks_re.empty() and
            not hdf.write("zero_blocks_re", &(zero_blocks_re[0]), zero_blocks_re.size())) or
        ((not zero_blocks_im.empty() and
            not hdf.write("zero_blocks_im", &(zero_blocks_im[0]), zero_blocks_im.size()))) or
        ((not elements_re.empty() and
            not hdf.write("re", &(elements_re[0]), elements_re.size()))) or
        ((not elements_im.empty() and
            not hdf.write("im", &(elements_im[0]), elements_im.size()))))
            return false;
        
        return true;
    }
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
    NumberArray<DataT> toPaddedRows () const;
    
    /**
     * @brief Zero-pad columns.
     * 
     * Pad rows with zeros as in @ref toPaddedRows, then concatenate
     * columns and return as a single array.
     */
    NumberArray<DataT> toPaddedCols () const;
    
    /// Convert matrix part to CooMatrix.
    template <class IdxT> CooMatrix<IdxT,DataT> tocoo (MatrixSelection::Selection triangle = MatrixSelection::Both) const
    {
        NumberArray<IdxT> I, J;
        NumberArray<DataT> V;
        
        DataT const * el = elems_.data();
        
        // for all elements
        for (IdxT i = 0; i < (IdxT)n_; i++)
        for (IdxT d = 0; d < (IdxT)d_; d++)
        {
            if (i + d < (IdxT)n_)
            {
                // skip zero elements
                if (*el == 0.0_r)
                {
                    el++;
                    continue;
                }
                
                // add this element to COO (upper triangle)
                if ((d != 0) and (triangle & MatrixSelection::StrictUpper))
                {
                    I.push_back(i);
                    J.push_back(i+d);
                    V.push_back(*el);
                }
                
                // add this element to COO (lower triangle)
                if ((d != 0) and (triangle & MatrixSelection::StrictLower))
                {
                    I.push_back(i+d);
                    J.push_back(i);
                    V.push_back(*el);
                }
                
                // main diagonal
                if ((d == 0) and (triangle & MatrixSelection::Diagonal))
                {
                    I.push_back(i);
                    J.push_back(i);
                    V.push_back(*el);
                }
            }
            
            // move on to the next element
            el++;
        }
        
        return CooMatrix<IdxT,DataT> (n_, n_, I, J, V);
    }
    
    /// Convert matrix part to RowMatrix.
    RowMatrix<DataT> torow (MatrixSelection::Selection triangle = MatrixSelection::Both) const
    {
        RowMatrix<DataT> M (n_, n_);
        
        DataT const * el = elems_.data();
        
        // for all elements
        for (std::size_t i = 0; i < n_; i++)
        for (std::size_t d = 0; d < d_; d++)
        {
            if (i + d < n_)
            {
                // main diagonal
                if ((d == 0) and (triangle & MatrixSelection::Diagonal))
                    M(i,i) = *el;
                
                // upper triangle
                if ((d != 0) and (triangle & MatrixSelection::StrictUpper))
                    M(i,i+d) = *el;
                
                // lower triangle
                if ((d != 0) and (triangle & MatrixSelection::StrictLower))
                    M(i+d,i) = *el;
            }
            // move on to the next elements
            el++;
        }
        
        return M;
    }
    
    /// Output to a text stream.
//     friend std::ostream & operator << (std::ostream & out, SymBandMatrix<DataT> const & A);
    
private:

    // dimension (only square matrices allowed)
    std::size_t n_;
    
    // main and upper diagonal count
    std::size_t d_;
    
    // diagonals concatenated in row-oriented way (and padded, if necessary) of size n_*d_
    NumberArray<DataT> elems_;
    
    // name of linked HDF file
    std::string name_;
};

// --------------------------------------------------------------------------------- //

template <class DataT> NumberArray<DataT> operator | (SymBandMatrix<DataT> const & A, const ArrayView<DataT> v) { return A.dot(v); }
template <class DataT> NumberArray<DataT> operator | (const ArrayView<DataT> v, SymBandMatrix<DataT> const & A) { return A.dot(v); }

// --------------------------------------------------------------------------------- //

template <class DataT> class BlockSymBandMatrix
{
    private:
        
        /// Name of HDF5 scratch disk file.
        std::string diskfile_;
        
        /// Whether to keep in memory.
        bool inmemory_;
        
        /// Number of blocks in a block row.
        std::size_t blockcount_;
        
        /// Block structure half-bandwidth.
        std::size_t blockhalfbw_;
        
        /// Size of a matrix block.
        std::size_t size_;
        
        /// Half bandwidth.
        std::size_t halfbw_;
        
        /// Data array.
        NumberArray<DataT> data_;
        
    public:
        
        //
        // Constructors.
        //
        
        BlockSymBandMatrix ()
            : diskfile_(), inmemory_(true), blockcount_(0), blockhalfbw_(0), size_(0), halfbw_(0), data_()
        {
            // nothing
        }
        
        /**
         * @brief Main constructor.
         * 
         * @param blockcount Number of blocks in a row (column) of the block matrix.
         * @param size Number of element in a row (column) of every single block.
         * @param name Name of the optional scratch disk file.
         */
        BlockSymBandMatrix (int blockcount, int blockhalfbw, int size, int halfbw, bool inmemory = true, std::string name = "")
            : diskfile_(name), inmemory_(inmemory), blockcount_(blockcount), blockhalfbw_(blockhalfbw), size_(size), halfbw_(halfbw), data_()
        {
            if (inmemory_)
                data_.resize(blockcount_ * size_ * blockhalfbw_ * halfbw_);
        }
        
        /// Copy constructor.
        BlockSymBandMatrix (BlockSymBandMatrix const & B)
            : diskfile_(), inmemory_(true), blockcount_(B.blockcount_), blockhalfbw_(B.blockhalfbw_), size_(B.size_), halfbw_(B.halfbw_), data_(B.data_)
        {
            if (not B.inmemory())
                HexException("Copy constructor of BlockSymBandMatrix only implemented for in-memory matrices!");
        }
        
        /// Is this object cached in memory?
        bool inmemory () const
        {
            return inmemory_;
        }
        
        /// Access to the memory buffer (maay be empty if data not in memory).
        NumberArray<DataT> const & data () const { return data_; }
        
        // Structure information.
        std::size_t blockcount () const { return blockcount_; }
        std::size_t blockhalfbw () const { return blockhalfbw_; }
        std::size_t size () const { return size_; }
        std::size_t halfbw () const { return halfbw_; }
        
        /**
         * @brief Access block.
         * 
         * Return the block [i, k].
         */
        SymBandMatrix<DataT> operator() (int i, int k) const
        {
            return SymBandMatrix<Complex>
            (
                size_, halfbw_,
                cArrayView
                (
                    data_,
                    (std::min(i, k) * blockhalfbw_ + std::abs(i - k)) * size_ * halfbw_,
                    size_ * halfbw_
                )
            );
        }
        
        /**
         * @brief Access element.
         * 
         * Return the element [i * size + j, k * size + l].
         */
        DataT operator() (int i, int j, int k, int l) const
        {
            return data_[((std::min(i, k) * blockhalfbw_ + std::abs(i - k)) * size_ + std::min(j, l)) * halfbw_ + std::abs(j - l)];
        }
        DataT & operator() (int i, int j, int k, int l)
        {
            return data_[((std::min(i, k) * blockhalfbw_ + std::abs(i - k)) * size_ + std::min(j, l)) * halfbw_ + std::abs(j - l)];
        }
        
        /// Reduced arithmetic operators.
        BlockSymBandMatrix<DataT> & operator += (BlockSymBandMatrix<DataT> const & B)
        {
            if (blockcount_ != B.blockcount_ or blockhalfbw_ != B.blockhalfbw_ or size_ != B.size_ or halfbw_ != B.halfbw_ or data_.size() != B.data_.size())
                HexException("The matrices supplied to operator \"+=\" have not the same structure.");
            
            data_ += B.data_;
            return *this;
        }
        BlockSymBandMatrix<DataT> & operator -= (BlockSymBandMatrix<DataT> const & B)
        {
            if (blockcount_ != B.blockcount_ or blockhalfbw_ != B.blockhalfbw_ or size_ != B.size_ or halfbw_ != B.halfbw_ or data_.size() != B.data_.size())
                HexException("The matrices supplied to operator \"-=\" have not the same structure.");
            
            data_ -= B.data_;
            return *this;
        }
        BlockSymBandMatrix<DataT> & operator *= (DataT z)
        {
            data_ *= z;
            return *this;
        }
        BlockSymBandMatrix<DataT> & operator /= (DataT z)
        {
            data_ /= z;
            return *this;
        }
        
        /// Access individual blocks.
        NumberArray<DataT> getBlock (int i) const
        {
            if (inmemory_)
            {
                return ArrayView<DataT> (data_, i * size_ * halfbw_, size_ * halfbw_);
            }
            else
            {
                // open disk file
                HDFFile hdf (diskfile_, HDFFile::readonly);
                
                // check success
                if (not hdf.valid())
                    HexException("Cannot open file \"%s\".\nError stack:\n%s", diskfile_.c_str(), hdf.error().c_str());
                
                // create output array
                NumberArray<DataT> data (size_ * halfbw_);
                
                // read data
                if (not hdf.read("data", &data[0], size_ * halfbw_, i * size_ * halfbw_))
                    HexException("Cannot access block %d in file \"%s\".\nError stack:\n%s", i, diskfile_.c_str(), hdf.error().c_str());
                
                return data;
            }
        }
        
        /// Access individual blocks.
        void setBlock (int i, const ArrayView<DataT> data)
        {
            if (data.size() != size_ * halfbw_)
                HexException("Wrong dimensions: %ld != %ld.", data.size(), size_ * halfbw_);
            
            if (inmemory_)
            {
                ArrayView<DataT>(data_, i * size_ * halfbw_, size_ * halfbw_) = data;
            }
            else
            {
                # pragma omp critical
                {
                    // check that the file exists
                    HDFFile * phdf = new HDFFile (diskfile_, HDFFile::readwrite);
                    if (not phdf->valid())
                    {
                        // create a new file
                        phdf = new HDFFile (diskfile_, HDFFile::overwrite);
                        if (not phdf->valid())
                            HexException("Cannot open file \"%s\" for writing.\nError stack:\n%s", diskfile_.c_str(), phdf->error().c_str());
                        
                        // initialize the dataset to its full length
                        if (not phdf->write("data", (Complex*)nullptr, size_ * halfbw_ * size_ * halfbw_))
                            HexException("Failed to initialize file \"%s\" (size %.1f).\n%s", diskfile_.c_str(), double(size_ * halfbw_)/std::pow(2,30), phdf->error().c_str());
                    }
                    
                    // write data
                    if (not phdf->write("data", &data[0], size_ * halfbw_, i * size_ * halfbw_))
                        HexException("Cannot access block %d in file \"%s\".\nError stack:\n%s", i, diskfile_.c_str(), phdf->error().c_str());
                    
                    delete phdf;
                }
            }
        }
        
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
        void dot (DataT a, const ArrayView<DataT> v, DataT b, ArrayView<DataT> w, bool parallelize = false, MatrixSelection::Selection tri = MatrixSelection::Both) const
        {
            // check vector size
            if (v.size() != blockcount_ * size_)
                HexException("[BlockSymBandMatrix::dot] Different size of matrix and vector in multiplication routine: %ld (mat) != %ld (vec).", blockcount_ * size_, v.size());
            
            // check if matrix is empty
            if (halfbw_ == 0)
                return;
            
            // scale destination vector
            w *= b;
            
            // create synchronization locks
            OMP_CREATE_LOCKS(blockcount_);
            
            // parallel section start
            # pragma omp parallel if (parallelize)
            {
                // thread-private workspace
                NumberArray<DataT> product (size_);
                
                // open data file for reading
                NumberArray<DataT> diskdata;
                HDFFile * hdf = nullptr;
                if (not inmemory_)
                    hdf = new HDFFile (diskfile_, HDFFile::readonly);
                
                // for all block diagonals
                for (std::size_t d = 0; d < blockhalfbw_; d++)
                {
                    // parallel processing of blocks in this diagonal
                    # pragma omp for schedule (dynamic,1)
                    for (std::size_t i = 0; i < blockcount_; i++)
                    if (i + d < blockcount_)
                    {
                        // block volume and offset
                        std::size_t vol = size_ * halfbw_;
                        std::size_t offset = (i * blockhalfbw_ + d) * vol;
                        
                        // data view of this block
                        ArrayView<DataT> view;
                        
                        // select matrix block data
                        if (inmemory_)
                        {
                            // use data in memory
                            view.reset(vol, data_.begin() + offset);
                        }
                        else
                        {
                            // read data from the disk
                            diskdata.resize(vol);
                            if (not hdf->read("data", &diskdata[0], vol, offset))
                                HexException("Failed to read HDF file \"%s\".\nHDF error stack:\n%s", diskfile_.c_str(), hdf->error().c_str());
                            
                            // reset view to the new data
                            view.reset(vol, diskdata.data());
                        }
                        
                        // multiply by diagonal block
                        if (d == 0)
                        {
                            SymBandMatrix<DataT>::sym_band_dot
                            (
                                size_, halfbw_, view,
                                1., ArrayView<DataT>(v, i * size_, size_),
                                0., product,
                                MatrixSelection::Both
                            );
                            
                            OMP_LOCK_LOCK(i);
                            {
                                # pragma omp simd
                                for (std::size_t pos = 0; pos < size_; pos++)
                                    w[i * size_ + pos] += a * product[pos];
                            }
                            OMP_UNLOCK_LOCK(i);
                        }
                        
                        // multiply by the other diagonals (both symmetries)
                        else
                        {
                            if (tri & MatrixSelection::StrictUpper)
                            {
                                SymBandMatrix<DataT>::sym_band_dot
                                (
                                    size_, halfbw_, view,
                                    1., ArrayView<DataT>(v, (i + d) * size_, size_),
                                    0., product
                                );
                                
                                OMP_LOCK_LOCK(i);
                                {
                                    # pragma omp simd
                                    for (std::size_t pos = 0; pos < size_; pos++)
                                        w[i * size_ + pos] += a * product[pos];
                                }
                                OMP_UNLOCK_LOCK(i);
                            }
                            
                            if (tri & MatrixSelection::StrictLower)
                            {
                                SymBandMatrix<DataT>::sym_band_dot
                                (
                                    size_, halfbw_, view,
                                    1., ArrayView<DataT>(v, i * size_, size_),
                                    0., product
                                );
                                
                                OMP_LOCK_LOCK(i + d);
                                {
                                    # pragma omp simd
                                    for (std::size_t pos = 0; pos < size_; pos++)
                                        w[(i + d) * size_ + pos] += a * product[pos];
                                }
                                OMP_UNLOCK_LOCK(i + d);
                            }
                        }
                    }
                }
                
                // release disk file
                if (not inmemory_)
                    delete hdf;
            }
            
            // delete synchronization locks
            OMP_DELETE_LOCKS();
        }
        
        //
        // Coversions to other matrix types.
        //
        
        /**
         * @brief Convert to COO format.
         * 
         * @param loadblocks Use blocks from scratch file instead of those in memory (if any).
         */
        template <class IdxT> CooMatrix<IdxT,DataT> tocoo () const
        {
            // number of structurally non-zero blocks (both upper and lower)
            std::size_t nblocks = blockcount_ * (2 * blockhalfbw_ - 1);
            
            // number of structurally non-zero elements (the structure is recursive)
            std::size_t nelem = nblocks * size_ * (2 * halfbw_ - 1);
            
            // allocate the ijv arrays
            NumberArray<IdxT> I; I.reserve(nelem);
            NumberArray<IdxT> J; J.reserve(nelem); 
            NumberArray<DataT> V; V.reserve(nelem);
            
            // open data file for reading
            HDFFile * hdf = nullptr;
            if (not inmemory_)
                hdf = new HDFFile (diskfile_, HDFFile::readonly);
            
            // for all blocks
            for (std::size_t i = 0; i < blockcount_; i++)
            for (std::size_t d = 0; d < blockhalfbw_; d++)
            if (i + d < blockcount_)
            {
                // data view of this block diagonal
                cArrayView view;
                
                // it may be necessary to load the data from disk
                cArray diskdata;
                if (inmemory_)
                {
                    // use data from memory
                    view.reset(size_ * halfbw_, data_.begin() + (i * blockhalfbw_ + d) * size_ * halfbw_);
                }
                else
                {
                    // read data from the disk
                    diskdata.resize(size_ * halfbw_);
                    if (not hdf->read("data", &diskdata[0], size_ * halfbw_, (i * blockhalfbw_ + d) * size_ * halfbw_))
                        HexException("Failed to read HDF file \"%s\".\nHDF error stack:\n%s", diskfile_.c_str(), hdf->error().c_str());
                    
                    // reset view to the new data
                    view.reset(size_ * halfbw_, diskdata.data());
                }
                
                // convert the block to COO format
                CooMatrix<IdxT,DataT> coo = SymBandMatrix<DataT>(size_, halfbw_, view).template tocoo<IdxT>();
                
                // copy all elements to whole-matrix arrays
                for (std::size_t k = 0; k < coo.v().size(); k++)
                {
                    I.push_back( coo.i()[k] +  i      * size_ );
                    J.push_back( coo.j()[k] + (i + d) * size_ );
                    V.push_back( coo.v()[k]                   );
                }
                if (d != 0)
                for (std::size_t k = 0; k < coo.v().size(); k++)
                {
                    I.push_back( coo.i()[k] + (i + d) * size_ );
                    J.push_back( coo.j()[k] +  i      * size_ );
                    V.push_back( coo.v()[k]                   );
                }
            }
            
            if (not inmemory_)
                delete hdf;
            
            // compose the final matrix
            return CooMatrix<IdxT,DataT>
            (
                blockcount_ * size_,  // number of rows
                blockcount_* size_,  // number of columns
                std::move(I),   // row indices
                std::move(J),   // column indices
                std::move(V)    // structurally nonzero matrix entries
            );
        }
        
        //
        // Scratch file I/O.
        //
        
        /// Get name of the scratch disk file.
        std::string   hdfname () const { return diskfile_; }
        std::string & hdfname ()       { return diskfile_; }
        
        /// Release data from memory, but keep HDF link.
        void drop ()
        {
            // release memory
            data_.drop();
            
            // turn off the "in memory" flag
            inmemory_ = false;
        }
        
        /// Check that the scratch disk file exists.
        bool hdfcheck () const { return HDFFile(diskfile_, HDFFile::readonly).valid(); }
        
        /// Reset (= create empty) HDF scratch file needed when writing individual blocks.
        bool hdfinit () const
        {
            // create a new file
            HDFFile hdf (diskfile_, HDFFile::overwrite);
            if (not hdf.valid())
                HexException("Cannot open file \"%s\" for writing.\nHDF error stack:\n%s", diskfile_.c_str(), hdf.error().c_str());
            
            // initialize the dataset to its full length
            if (not hdf.write("data", (Complex*)nullptr, blockcount_ * blockhalfbw_ * size_ * halfbw_))
                HexException("Failed to initialize file \"%s\" (size %.1f).\n%s", diskfile_.c_str(), double(size_ * halfbw_)/std::pow(2,30), hdf.error().c_str());
            
            return true;
        }
        
        /// Write all blocks from memory to the disk file.
        bool hdfsave () const
        {
            // open disk file
            HDFFile hdf (diskfile_, HDFFile::overwrite);
            
            // check success
            if (not hdf.valid())
                return false;
            
            // save data
            if (not hdf.write("data", &data_[0], data_.size()))
                return false;
            
            return true;
        }
        
        /// Load data from disk to memory.
        bool hdfload ()
        {
            // open disk file
            HDFFile hdf (diskfile_, HDFFile::readonly);
            
            // check success
            if (not hdf.valid())
                return false;
            
            // allocate memory
            data_.resize(blockcount_ * size_ * blockhalfbw_ * halfbw_);
            
            // read data
            if (not hdf.read("data", &data_[0], data_.size()))
                return false;
            
            // turn on the "in memory" flag
            inmemory_ = true;
            
            return true;
        }
};

template <class DataT>
SymBandMatrix<DataT> operator + (SymBandMatrix<DataT> const & A, SymBandMatrix<DataT> const & B)
{
    A.is_compatible(B);
    return SymBandMatrix<DataT>(A.size(), A.halfbw(), A.data() + B.data());
}

template <class DataT>
SymBandMatrix<DataT> operator - (SymBandMatrix<DataT> const & A, SymBandMatrix<DataT> const & B)
{
    A.is_compatible(B);
    return SymBandMatrix<DataT>(A.size(), A.halfbw(), A.data() - B.data());
}

template <class DataT>
SymBandMatrix<DataT> operator * (DataT z, SymBandMatrix<DataT> const & A)
{
    return SymBandMatrix<DataT>(A.size(), A.halfbw(), z * A.data());
}

template <class DataT>
BlockSymBandMatrix<DataT> kron
(
    SymBandMatrix<DataT> const & A,
    SymBandMatrix<DataT> const & B
)
{
    BlockSymBandMatrix<DataT> C (A.size(), A.halfbw(), B.size(), B.halfbw());
    
    for (std::size_t i = 0; i < A.size(); i++)
    for (std::size_t u = 0; u < A.halfbw(); u++)
    if (i + u < A.size())
    for (std::size_t j = 0; j < B.size(); j++)
    for (std::size_t v = 0; v < B.halfbw(); v++)
    if (j + v < B.size())
    {
        C(i,j,i+u,j+v) = A(i,i+u) * B(j,j+v);
    }
    
    return C;
}

template<class DataT>
BlockSymBandMatrix<DataT> operator *
(
    DataT z,
    BlockSymBandMatrix<DataT> const & A
)
{
    BlockSymBandMatrix<DataT> C (A);
    return C *= z;
}

template<class DataT>
BlockSymBandMatrix<DataT> operator *
(
    BlockSymBandMatrix<DataT> const & A,
    DataT z
)
{
    BlockSymBandMatrix<DataT> C (A);
    return C *= z;
}

template<class DataT>
BlockSymBandMatrix<DataT> operator +
(
    BlockSymBandMatrix<DataT> const & A,
    BlockSymBandMatrix<DataT> const & B
)
{
    BlockSymBandMatrix<DataT> C (A);
    return C += B;
}

template<class DataT>
BlockSymBandMatrix<DataT> operator -
(
    BlockSymBandMatrix<DataT> const & A,
    BlockSymBandMatrix<DataT> const & B
)
{
    BlockSymBandMatrix<DataT> C (A);
    return C -= B;
}

/**
 * @brief Kronecker product.
 * 
 * Applies the following operation:
 * \f[
 *     \mathbf{w} = a \mathbf{w} + b (\mathsf{A} \otimes \mathsf{B}) \cdot \mathbf{v} \,.
 * \f]
 */
template <class DataT> void kron_dot
(
    Real a,       ArrayView<DataT> w,
    Real b, const ArrayView<DataT> v,
    SymBandMatrix<DataT> const & A,
    SymBandMatrix<DataT> const & B
)
{
    std::size_t A_size = A.size();
    std::size_t B_size = B.size();
    
    # pragma omp parallel for collapse (2)
    for (std::size_t i = 0; i < A_size; i++)
    for (std::size_t j = 0; j < B_size; j++)
    {
        // iteration bounds
        std::size_t kmin = (i >= A.halfbw() ? i - A.halfbw() : 0);
        std::size_t lmin = (j >= B.halfbw() ? j - B.halfbw() : 0);
        std::size_t kmax = std::min(A.size() - 1, i + A.halfbw() - 1);
        std::size_t lmax = std::min(B.size() - 1, j + B.halfbw() - 1);
        
        // calculate the scalar product
        DataT res = 0;
        for (std::size_t k = kmin; k <= kmax; k++)
        for (std::size_t l = lmin; l <= lmax; l++)
            res += A(i,k) * B(j,l) * v[k * B.size() + l];
        
        // save result
        w[i * B.size() + j] = a * w[i * B.size() + j] + b * res;
    }
}

// --------------------------------------------------------------------------------- //

#endif // HEX_SYMBANDMATRIX_H
