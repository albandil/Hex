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

#ifdef WITH_SCALAPACK

// --------------------------------------------------------------------------------- //

#include "lu-scalapack.h"

// --------------------------------------------------------------------------------- //

template<>
LUft_SCALAPACK<LU_int_t,Complex>::LUft_SCALAPACK ()
{
    // nothing
}

template<>
void LUft_SCALAPACK<LU_int_t,Complex>::drop ()
{
    
}

template<>
LUft_SCALAPACK<LU_int_t,Complex>::~LUft_SCALAPACK ()
{
    // nothing
}

template<>
void LUft_SCALAPACK<LU_int_t,Complex>::factorize (CsrMatrix<LU_int_t,Complex> const & matrix, LUftData data)
{
    
}

template<>
void LUft_SCALAPACK<LU_int_t,Complex>::solve (const ArrayView<Complex> b, ArrayView<Complex> x, int eqs) const
{
    
}

template<>
void LUft_SCALAPACK<LU_int_t,Complex>::save (std::string name) const
{
    HexException("ScaLAPACK factorizer does not yet support --out-of-core option.");
}

template<>
void LUft_SCALAPACK<LU_int_t,Complex>::load (std::string name, bool throw_on_io_failure)
{
    if (throw_on_io_failure)
        HexException("ScaLAPACK factorizer does not yet support --out-of-core option.");
}

template<>
std::size_t LUft_SCALAPACK<LU_int_t,Complex>::size () const
{
    return 0;
}

// --------------------------------------------------------------------------------- //

/* addFactorizerToRuntimeSelectionTable(SCALAPACK, LU_int_t, Complex) */

// --------------------------------------------------------------------------------- //

#endif
