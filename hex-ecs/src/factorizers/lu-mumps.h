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

#if ( !defined(HEX_LU_MUMPS_H) && defined(WITH_MUMPS) )
#define HEX_LU_MUMPS_H

// --------------------------------------------------------------------------------- //

#include "hex-csrmatrix.h"

// --------------------------------------------------------------------------------- //

#include "luft.h"

// --------------------------------------------------------------------------------- //

#define MUMPS_NOACTION      0
#define MUMPS_INITIALIZE    (-1)
#define MUMPS_FINISH        (-2)
#define MUMPS_ANALYZE       1
#define MUMPS_FACTORIZE     2
#define MUMPS_SOLVE         3

// --------------------------------------------------------------------------------- //

#ifdef SINGLE
    #include <cmumps_c.h>
    #define MUMPS_STRUC_C CMUMPS_STRUC_C
    #define MUMPS_C cmumps_c
    #define MUMPS_COMPLEX mumps_complex
#else
    #include <zmumps_c.h>
    #define MUMPS_STRUC_C ZMUMPS_STRUC_C
    #define MUMPS_C zmumps_c
    #define MUMPS_COMPLEX mumps_double_complex
#endif

// --------------------------------------------------------------------------------- //

/**
 * @brief LU factorization object - MUMPS specialization.
 * 
 * This class holds information on LU factorization as computed by the free
 * library MUMPS. It is derived from LUft and shares interface with that class.
 * 
 * @warning This class expect only symmetric matrices. The COO matrix passed
 * to this class should contain only upper or lower part of the matrix (together
 * with the main diagonal).
 */
class LUft_MUMPS : public LUft
{
    public:
        
        typedef struct
        {
            MUMPS_INT out_of_core;
            MUMPS_INT verbosity_level;
            MUMPS_INT comm;
        }
        Data;
        
        // run-time selection mechanism
        factorizerRunTimeSelectionDefinitions(LUft_MUMPS, "mumps")
        
        /// Default constructor.
        LUft_MUMPS ();
        
        /// Destructor.
        virtual ~LUft_MUMPS ();
        
        // Disable bitwise copy
        LUft_MUMPS const & operator= (LUft_MUMPS const &) = delete;
        
        /// Factorize.
        virtual void factorize (CsrMatrix<LU_int_t,Complex> const & matrix);
        
        /// Validity indicator.
        virtual bool valid () const
        {
            return mmin(I.size(), J.size(), A.size()) > 0;
        }
        
        /// Return LU byte size.
        virtual std::size_t size () const;
        
        /// Condition number.
        virtual Real cond () const
        {
            #define RINFO(x) rinfo[x-1]
            return settings.RINFO(11);
        }
        
        /// Solve equations.
        virtual void solve (const cArrayView b, cArrayView x, int eqs) const;
        
        /// Save large data to disk.
        virtual void save (std::string name) const
        {
            I.hdfsave("I-" + name);
            J.hdfsave("J-" + name);
            A.hdfsave("A-" + name);
        }
        
        /// Load large data from disk.
        virtual void load (std::string name, bool throw_on_io_failure = true)
        {
            if (not I.hdfload("I-" + name) or
                not J.hdfload("J-" + name) or
                not A.hdfload("A-" + name))
            {
                if (throw_on_io_failure)
                    HexException("Failed to load MUMPS IJV matrices.");
            }
        }
        
        /// Release memory.
        virtual void drop ()
        {
            I.drop();
            J.drop();
            A.drop();
            
            if (workspace_ != nullptr)
            {
                vMemAllocator<MUMPS_COMPLEX>::free(workspace_);
                
                workspace_ = nullptr;
                workspace_size_ = 0;
            }
        }
        
    private:
        
        // Internal data of the library.
        mutable MUMPS_STRUC_C settings;
        
        // rank
        MUMPS_INT n_;
        
        // data arrays
        NumberArray<MUMPS_INT> I, J;
        NumberArray<Complex> A;
        
        // additional workspace
        MUMPS_COMPLEX* workspace_;
        std::size_t workspace_size_;
};

#endif // HEX_LU_MUMPS_H, WITH_MUMPS
