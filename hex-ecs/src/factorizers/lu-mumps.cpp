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

#ifdef WITH_MUMPS

#define ICNTL(x) icntl[(x)-1]
#define INFO(x)  info[(x)-1]
#define INFOG(x) infog[(x)-1]

// --------------------------------------------------------------------------------- //

#include "lu-mumps.h"

// --------------------------------------------------------------------------------- //

LUft_MUMPS::LUft_MUMPS () : LUft()
{
    settings.job = MUMPS_NOACTION;

    workspace_ = nullptr;
    workspace_size_ = 0;

    static const char * ooc_dir  = ".";

    idata_["verbosity"]          = 0;
    idata_["centralized_matrix"] = true;
    idata_["out_of_core"]        = false;
    idata_["fortran_comm"]       = 0;
    idata_["disk_workspace"]     = false;

    rdata_["memory_relaxation"]  = 20.0;

    pdata_["ooc_dir"]            = (void*)ooc_dir;
}

LUft_MUMPS::~LUft_MUMPS ()
{
    drop ();

    // destroy MUMPS data
    if (settings.job != MUMPS_NOACTION)
    {
        settings.job = MUMPS_FINISH;
        MUMPS_C(&settings);

        if (settings.INFOG(1) != 0)
            HexException("MUMPS finish failed, error code: %d, detail %d", settings.INFOG(1), settings.INFOG(2));
    }
}

std::size_t LUft_MUMPS::size () const
{
    std::size_t elems_size = (settings.INFO(9) > 0 ? settings.INFO(9) : std::size_t{1000000} * std::abs(settings.INFO(9)));
    std::size_t index_size = settings.INFO(10);

    return sizeof(MUMPS_COMPLEX) * elems_size + sizeof(MUMPS_INT) * index_size;
}

void LUft_MUMPS::factorize (CsrMatrix<LU_int_t,Complex> const & matrix)
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
        nz = 0;
        for (LU_int_t row = 0; row < n_; row++)
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
        settings.comm_fortran = idata_["fortran_comm"];
#endif
        MUMPS_C(&settings);

        if (settings.INFOG(1) != 0)
            HexException("MUMPS initialization failed, error code: %d, detail %d", settings.INFOG(1), settings.INFOG(2));

        // analyze
        settings.job = MUMPS_ANALYZE;
        settings.ICNTL(1) = (idata_["verbosity"] == 0 ? 0 : 6); // errors to STDOUT (default: 6)
        settings.ICNTL(2) = 0; // diagnostics to /dev/null
        settings.ICNTL(3) = (idata_["verbosity"] == 0 ? 0 : 6); // global info to STDOUT (default: 6)
        settings.ICNTL(4) = idata_["verbosity"]; // verbosity level (default: 2)
        settings.ICNTL(5) = 0; // COO format
        settings.ICNTL(14) = rdata_["memory_relaxation"];
        settings.ICNTL(18) = rdata_["centralized_matrix"] ? 0 : 3;
        settings.ICNTL(22) = idata_["out_of_core"]; // OOC factorization
        settings.ICNTL(28) = 2; // parallel ordering
        std::strcpy(settings.ooc_tmpdir, (const char*)pdata_["ooc_dir"]);
        std::strcpy(settings.ooc_prefix, "ooc_");
        settings.n = n_;
        settings.nz  = settings.nz_loc  = nz;
        settings.irn = settings.irn_loc = I.data();
        settings.jcn = settings.jcn_loc = J.data();
        settings.a   = settings.a_loc   = reinterpret_cast<MUMPS_COMPLEX*>(A.data());
        MUMPS_C(&settings);

        if (settings.INFOG(1) != 0)
            HexException("MUMPS analysis failed, error code: %d, detail %d", settings.INFOG(1), settings.INFOG(2));

        if (idata_["disk_workspace"])
        {
            int param = (settings.ICNTL(22) ? settings.ICNTL(20) : settings.ICNTL(8));

            workspace_size_ = (param >= 0 ? param : (-param) * std::size_t{1000000});
            workspace_ = vMemAllocator<MUMPS_COMPLEX>::alloc(workspace_size_, settings.ooc_tmpdir, settings.ooc_prefix);

            settings.lwk_user = param;
            settings.wk_user = workspace_;
        }

    //
    // Compute the factorization.
    //

        settings.job = MUMPS_FACTORIZE;
        MUMPS_C(&settings);

        if (settings.INFOG(1) != 0)
            HexException("MUMPS factorization failed, error code: %d, detail %d", settings.INFOG(1), settings.INFOG(2));
}

void LUft_MUMPS::solve (const cArrayView b, cArrayView x, int eqs) const
{
    // copy right-hand side to the solution vector
    if (x.data() != b.data())
        x = b;

    // run the back-substitution
    settings.job  = MUMPS_SOLVE;
    settings.nrhs = 1;
    settings.lrhs = x.size();
    settings.rhs  = reinterpret_cast<MUMPS_COMPLEX*>(x.data());
    MUMPS_C(&settings);

    if (settings.INFOG(1) != 0)
        HexException("MUMPS backsubstitution failed, error code: %d, detail %d", settings.INFOG(1), settings.INFOG(2));
}

void LUft_MUMPS::save (std::string name) const
{ 
    // TODO
}

void LUft_MUMPS::load (std::string name, bool throw_on_io_failure)
{
    // TODO

    if (throw_on_io_failure)
        HexException("Failed to load MUMPS LU decomposition from disk.");
}


// --------------------------------------------------------------------------------- //

addClassToParentRunTimeSelectionTable(LUft, LUft_MUMPS)

// --------------------------------------------------------------------------------- //

#endif // WITH_MUMPS
