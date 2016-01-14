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

#ifndef HEX_CSCMATRIX_H
#define HEX_CSCMATRIX_H

#include "hex-arrays.h"
#include "hex-matrix.h"
#include "hex-coomatrix.h"

/**
 * @brief Complex CSC matrix.
 * 
 * Complex compressed-storage-column-major-ordered sparse matrix.
 * The data are stored in three arrays,
 * @code
 * p, i, x
 * @endcode
 * The format is that used in UMFPACK with interleaved real and comples parts;
 * for explanation see UMFPACK Uses's Guide.
 */
template <class IdxT, class DataT> class CscMatrix
{
    
public:
    
    // Constructors
    
    CscMatrix ()
        : m_(0), n_(0) {}
    CscMatrix (IdxT m, IdxT n)
        : m_(m), n_(n) {}
    CscMatrix (CscMatrix const & A)
        : m_(A.m_), n_(A.n_), p_(A.p_), i_(A.i_), x_(A.x_) {}
    CscMatrix (IdxT m, IdxT n, const ArrayView<IdxT> p, const ArrayView<IdxT> i, const ArrayView<DataT> x)
        : m_(m), n_(n), p_(p), i_(i), x_(x) {}
    
    // Destructor
    
    ~CscMatrix () {}
    
    /**
     * Convert to COO-matrix.
     */
    CooMatrix<IdxT,DataT> tocoo () const;
    
    /**
     * Matrix-vector product using a transposed matrix, \f$ A^T \cdot b \f$.
     * For ordinary matrix vector product convert this matrix first to CSR
     * format.
     * @param b Vector to multiply with.
     */
    NumberArray<DataT> dotT (const ArrayView<DataT> b) const
    {
        // create output array
        NumberArray<DataT> c (n_);
        
        // the matrix "*this" is actually transposed
        for (IdxT icol = 0; icol < n_; icol++)
        {
            IdxT idx1 = p_[icol];
            IdxT idx2 = p_[icol+1];
            
            // for all nonzero elements in this column
            for (IdxT idx = idx1; idx < idx2; idx++)
            {
                // get row number
                IdxT irow = i_[idx];
                
                // store product
                c[icol] += x_[idx] * b[irow];
            }
        }
        
        return c;
    }
    
    // Getters
    
    std::size_t size () const { return i_.size(); }
    std::size_t rows () const { return m_; }
    std::size_t cols () const { return n_; }
    lArray const & p () const { return p_; }
    lArray const & i () const { return i_; }
    cArray const & x () const { return x_; }
    
    /**
     * Multiplication by a number.
     */
    CscMatrix<IdxT,DataT> & operator *= (DataT r)
    {
        std::size_t N = i_.size();
        for (std::size_t i = 0; i < N; i++)
            x_[i] *= r;
        return *this;
    }
    
    /**
     * Addition of another CSC matrix.
     * The matrices MUST HAVE THE SAME SPARSE STRUCTURE, as no indices are
     * checked.
     */
    CscMatrix<IdxT,DataT> & operator &= (CscMatrix<IdxT,DataT> const &  B)
    {
        std::size_t N = i_.size();
        
        assert(m_ == B.m_);
        assert(n_ == B.n_);
        assert(N == B.i_.size());
        
        for (std::size_t i = 0; i < N; i++)
        {
            assert(i_[i] == B.i_[i]);
            
            x_[i] += B.x_[i];
        }
        
        return *this;
    }
    
    /**
     * Subtraction of another CSC matrix.
     * The matrices MUST HAVE THE SAME SPARSE STRUCTURE, as no indices are
     * checked.
     */
    CscMatrix<IdxT,DataT> & operator ^= (CscMatrix<IdxT,DataT> const &  B)
    {
        std::size_t N = i_.size();
        
        assert(m_ == B.m_);
        assert(n_ == B.n_);
        assert(N == B.i_.size());
        
        for (std::size_t i = 0; i < N; i++)
        {
            assert(i_[i] == B.i_[i]);
            
            x_[i] -= B.x_[i];
        }
        
        return *this;
    }
    
    /**
     * Save matrix to HDF file.
     * @param name Filename.
     */
    bool hdfsave (const char* name) const
    {
        HDFFile file(name, HDFFile::overwrite);
        
        // write dimensions
        file.write("m", &m_, 1);
        file.write("n", &n_, 1);
        
        // write indices
        if (not p_.empty())
            file.write("p", &(p_[0]), p_.size());
        if (not i_.empty())
            file.write("i", &(i_[0]), i_.size());
        
        // write complex data as a "double" array
        if (not x_.empty())
        {
            file.write
            (
                "x",
                reinterpret_cast<typename typeinfo<DataT>::cmpttype const*>( &(x_[0]) ),
                x_.size() * typeinfo<DataT>::ncmpt
            );
        }
        
        return true;
    }
    
    /**
     * Load matrix from HDF file.
     * @param name Filename.
     */
    bool hdfload (const char* name)
    {
        HDFFile hdf(name, HDFFile::readonly);
        
        // read dimensions
        hdf.read("m", &m_, 1);
        hdf.read("n", &n_, 1);
        
        // read indices
        if (p_.resize(hdf.size("p")))
            hdf.read("p", &(p_[0]), p_.size());
        if (i_.resize(hdf.size("i")))
            hdf.read("i", &(i_[0]), i_.size());
        
        // read data
        if (x_.resize(hdf.size("x") / 2))
        {
            hdf.read
            (
                "x",
                reinterpret_cast<typename typeinfo<DataT>::cmpttype*>(&(x_[0])),
                x_.size() * typeinfo<DataT>::ncmpt
            );
        }
        
        return true;
    }
    
private:
    
    // dimensions
    IdxT m_;
    IdxT n_;
    
    // representation
    NumberArray<IdxT> p_;
    NumberArray<IdxT> i_;
    NumberArray<DataT> x_;
};

/**
 * Computes a sum of two csc-matrices OF THE SAME SPARSE STRUCTURE.
 */
template <class IdxT, class DataT>
CscMatrix<IdxT,DataT> operator & (CscMatrix<IdxT,DataT> const & A, CscMatrix<IdxT,DataT> const & B)
{
    CscMatrix<IdxT,DataT> C = A;
    return C &= B;
}

/**
 * Computes a difference of two csc-matrices OF THE SAME SPARSE STRUCTURE.
 */
template <class IdxT, class DataT>
CscMatrix<IdxT,DataT> operator ^ (CscMatrix<IdxT,DataT> const & A, CscMatrix<IdxT,DataT> const & B)
{
    CscMatrix<IdxT,DataT> C = A;
    return C ^= B;
}

/**
 * Multiplication of csc-matrix by a number.
 */
template <class IdxT, class DataT>
CscMatrix<IdxT,DataT> operator * (DataT r, CscMatrix<IdxT,DataT> const & B)
{
    CscMatrix<IdxT,DataT> C = B;
    return C *= r;
}

/**
 * Multiplication of csc-matrix by a number.
 */
template <class IdxT, class DataT>
CscMatrix<IdxT,DataT> operator * (CscMatrix<IdxT,DataT> const & A, DataT r)
{
    CscMatrix<IdxT,DataT> C = A;
    return C *= r;
}

#endif // HEX_CSCMATRIX_H
