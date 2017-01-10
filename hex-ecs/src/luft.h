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

#ifndef HEX_LUFT_H
#define HEX_LUFT_H

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"
#include "hex-rts.h"

// --------------------------------------------------------------------------------- //

template <class IdxT, class DataT> class CsrMatrix;

// --------------------------------------------------------------------------------- //

typedef struct
{
    Real drop_tolerance;
    bool out_of_core;
    int verbosity;
    int fortran_comm;
    int groupsize;
    const char * ooc_dir;
    void * superlu_dist_grid;
}
LUftData;

// --------------------------------------------------------------------------------- //

extern LUftData defaultLUftData;

// --------------------------------------------------------------------------------- //

/**
 * @brief LU factorization object
 * 
 * This class is returned by the function CsrMatrix::factorize() and
 * it provides some functions that can be used when solving equations with
 * that LU factorization. Also, it is possible to store the decomposition
 * to disk (link, save), destroy the data (drop) and load later when needed
 * (load). The most important function is "solve".
 */
class LUft
{
    public:
        
        //
        // Run-time selection mechanism (object factory).
        //
        
            baseClassRunTimeSelectionDefinitions(LUft, ())
        
        //
        // Class member functions.
        //
        
            /// Default constructor.
            LUft () : filename_() {}
            
            /// Destructor.
            virtual ~LUft () {}
            
            /// Factorize.
            virtual void factorize (CsrMatrix<LU_int_t,Complex> const & matrix, LUftData data = defaultLUftData) = 0;
            
            /**
             * @brief Validity indicator.
             * 
             * Returns true when the object contains a valid LU factorization.
             */
            virtual bool valid () const { return false; }
            
            /**
             * @brief Free memory.
             * 
             * Release memory occupied by the LU-factorization numeric object.
             */
            virtual void drop () {}
            
            /**
             * @brief Size of the numerical data.
             * 
             * Return the number of bytes occupied by the stored elements
             * of the LU-factorization. This doesn't contain any other structural data.
             */
            virtual std::size_t size () const { return 0; }
            
            /**
             * @brief Estimation of the condition number.
             * 
             * Note: Currently implemented only for UMFPACK backend.
             */
            virtual Real cond () const { return 0; }
            
            /**
             * @brief Solve equations.
             * 
             * The parameter "b" is assumed to contain several right hand
             * side vectors (their count is supplied as the optional parameter
             * "eqs"). The results are stored in "x", which has the same size
             * as "b".
             */
            virtual void solve (const cArrayView b, cArrayView x, int eqs) const = 0;
            
            /**
             * @brief Solve equations.
             * 
             * The parameter "b" is assumed to contain several right hand
             * side vectors (their count is supplied as the optional parameter
             * "eqs").
             */
            virtual cArray solve (const cArrayView b, int eqs = 1) const
            {
                // reserve space for the solution
                cArray x (b.size());
                
                // solve
                this->solve(b, x, eqs);
                
                // return the result
                return x;
            }
            
            /**
             * @brief Link to a disk file.
             * 
             * This function will set a filename that will be used if
             * any of the functions @ref save or @ref load is used without
             * a specific filename.
             */
            virtual void link (std::string name) { filename_ = name; }
            
            /**
             * @brief Unlink disk file.
             */
            virtual void unlink () { filename_.clear(); }
            
            /**
             * @brief Name of the linked disk file.
             */
            virtual std::string filename () const { return filename_; }
            
            /**
             * @brief Save factorization object to a disk file.
             * 
             * Stores the LU-factorization data to a disk file in the native format
             * of the library used.
             */
            virtual void save (std::string name) const = 0;
            virtual void save () const { this->save(filename_); }
            
            /**
             * @brief Load factorization object from a disk file.
             */
            virtual void load (std::string name, bool throw_on_io_failure = true) = 0;
            virtual void load () { this->load(filename_, true); }
            virtual void silent_load () { this->load(filename_, false); }
        
    private:
        
        /// Name of the disk file.
        std::string filename_;
};

// --------------------------------------------------------------------------------- //

#define factorizerRunTimeSelectionDefinitions(TYPE,NAME) \
    derivedClassRunTimeSelectionDefinitions(LUft, (), TYPE, (), NAME)

// --------------------------------------------------------------------------------- //

#endif // HEX_LUFT_H
