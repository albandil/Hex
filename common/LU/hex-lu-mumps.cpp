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

#ifdef WITH_MUMPS

// --------------------------------------------------------------------------------- //

#include "hex-lu-mumps.h"

// --------------------------------------------------------------------------------- //

template<>
LUft_MUMPS<LU_int_t,Complex>::LUft_MUMPS ()
    : LUft<LU_int_t,Complex>()
{
    settings.job = MUMPS_NOACTION;
}

template<>
LUft_MUMPS<LU_int_t,Complex>::~LUft_MUMPS ()
{
    drop ();
    
    // destroy MUMPS data
    if (settings.job != MUMPS_NOACTION)
    {
        settings.job = MUMPS_FINISH;
        zmumps_c(&settings);
    }
}

template<>
std::size_t LUft_MUMPS<LU_int_t,Complex>::size () const
{
    #define INFO(x) info[x-1]
    std::size_t elems_size = (settings.INFO(9) > 0 ? settings.INFO(9) : std::size_t{1000000} * std::abs(settings.INFO(9)));
    std::size_t index_size = settings.INFO(10);
    
    return sizeof(MUMPS_COMPLEX) * elems_size + sizeof(MUMPS_INT) * index_size;
}

template<>
void LUft_MUMPS<LU_int_t,Complex>::factorize (CsrMatrix<LU_int_t,Complex> const & matrix, LUftData data)
{
    //
    // Create matrix of the system (i.e. the IJV triplets of the symmetric part).
    //

        n_ = matrix.rows();
        
        // estimate non-zero element count of the upper triangle
        LU_int_t nz = (matrix.i().size() + n_) / 2;
        
        // allocate memory
        I.resize(0); I.reserve(nz);
        J.resize(0); J.reserve(nz);
        A.resize(0); A.reserve(nz);
        
        // for all rows
        for (LU_int_t row = 0, nz = 0; row < n_; row++)
        {
            // for all columns with structurally non-zero entries
            for (LU_int_t idx = matrix.p()[row]; idx < matrix.p()[row + 1]; idx++)
            {
                // get column index
                LU_int_t col = matrix.i()[idx];
                
                // only consider upper triangle (and the main diagonal)
                if (row <= col and matrix.x()[idx] != 0.0_z)
                {
                    // insert the element
                    I.push_back(row + 1);
                    J.push_back(col + 1);
                    A.push_back(matrix.x()[idx]);
                    
                    // update true element count
                    nz++;
                }
            }
        }
    
    //
    // Prepare MUMPS environment.
    //
    
        // initialize
        settings.job = MUMPS_INITIALIZE;
        settings.sym = 2;
        settings.par = 1;
#ifdef WITH_MPI
        settings.comm_fortran = data.fortran_comm;
#endif
        MUMPS_C(&settings);
        
        // analyze
        settings.job = MUMPS_ANALYZE;
        settings.ICNTL(1) = (data.verbosity == 0 ? 0 : 6); // errors to STDOUT (default: 6)
        settings.ICNTL(2) = 0; // diagnostics to /dev/null
        settings.ICNTL(3) = (data.verbosity == 0 ? 0 : 6); // global info to STDOUT (default: 6)
        settings.ICNTL(4) = data.verbosity; // verbosity level (default: 2)
        settings.ICNTL(5) = 0; // COO format
        settings.ICNTL(22) = data.out_of_core; // OOC factorization
        std::strcpy(settings.ooc_tmpdir, ".");
        std::strcpy(settings.ooc_prefix, "ooc_");
        settings.n = n_;
        settings.nz = nz;
        settings.irn = I.data();
        settings.jcn = J.data();
        settings.a = reinterpret_cast<MUMPS_COMPLEX*>(A.data());
        MUMPS_C(&settings);
    
    //
    // Compute the factorization.
    //
    
        settings.job = MUMPS_FACTORIZE;
        MUMPS_C(&settings);
}

template<>
void LUft_MUMPS<LU_int_t,Complex>::solve (const cArrayView b, cArrayView x, int eqs) const
{
    // copy right-hand side to the solution vector
    if (x.data() != b.data())
        x = b;
    
    // run the back-substitution
    settings.nrhs = 1;
    settings.lrhs = x.size();
    settings.rhs = reinterpret_cast<MUMPS_COMPLEX*>(x.data());
    settings.irn = const_cast<MUMPS_INT*>(I.data());
    settings.jcn = const_cast<MUMPS_INT*>(J.data());
    settings.a   = reinterpret_cast<MUMPS_COMPLEX*>(const_cast<Complex*>(A.data()));
    settings.job = 3;
    MUMPS_C(&settings);
}

// --------------------------------------------------------------------------------- //

addFactorizerToRuntimeSelectionTable(MUMPS, LU_int_t, Complex)

// --------------------------------------------------------------------------------- //

#endif // WITH_MUMPS
