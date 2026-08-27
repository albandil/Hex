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

#ifndef HEX_LU_LAPACK_H
#define HEX_LU_LAPACK_H

// --------------------------------------------------------------------------------- //

#include "hex-csrmatrix.h"

// --------------------------------------------------------------------------------- //

#include "luft.h"

// --------------------------------------------------------------------------------- //

/**
 * @brief LU factorization using LAPACK.
 * 
 * This makes sense only for very small matrices.
 */
class LUft_LAPACK : public LUft
{
    public:

        // run-time selection mechanism
        factorizerRunTimeSelectionDefinitions(LUft_LAPACK, "lapack")

        /// Default constructor.
        LUft_LAPACK ();

        /// Destructor.
        ~LUft_LAPACK() override;

        // Disable bitwise copy
        LUft_LAPACK const & operator= (LUft_LAPACK const &) = delete;

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

        /// Cuthill-McKee ordering.
        NumberArray<LU_int_t> R_;

        /// Factorization computed by xGBTRF use in xGBTRS.
        cArray LU_;

        /// Pivot sequence.
        NumberArray<blas::Int> ipiv_;

        /// Matrix rank.
        blas::Int n_;

        /// Matrix half-bandwidth.
        blas::Int k_;
};

// --------------------------------------------------------------------------------- //

#endif // HEX_LU_LAPACK_H
