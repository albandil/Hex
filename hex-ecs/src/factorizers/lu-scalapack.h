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

#if (!defined(HEX_LU_SCALAPACK_H) && defined(WITH_SCALAPACK))
#define HEX_LU_SCALAPACK_H

// --------------------------------------------------------------------------------- //

#include "hex-csrmatrix.h"

// --------------------------------------------------------------------------------- //

#include "luft.h"

// --------------------------------------------------------------------------------- //

/**
 * @brief LU factorization using ScaLAPACK.
 * 
 * This is a really last-resort option tractable only for very small matrices.
 * They are converted to a banded format (using reverse Cuthill McKee algorithm)
 * and solved by a direct application of the banded PZGBTRF + PZGBTRS routines.
 */
class LUft_SCALAPACK : public LUft
{
    public:

        // run-time selection mechanism
        factorizerRunTimeSelectionDefinitions(LUft_SCALAPACK, "scalapack")

        /// Default constructor.
        LUft_SCALAPACK ();

        /// Destructor.
        ~LUft_SCALAPACK() override;

        // Disable bitwise copy
        LUft_SCALAPACK const & operator= (LUft_SCALAPACK const &) = delete;

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
};

// --------------------------------------------------------------------------------- //

#endif // WITH_SCALAPACK
