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

#ifndef HEX_ECS_PRECONDITIONERS
#define HEX_ECS_PRECONDITIONERS

// --------------------------------------------------------------------------------- //

#include <tuple>

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"
#include "hex-matrix.h"
#include "hex-rts.h"

// --------------------------------------------------------------------------------- //

#include "ang.h"
#include "radial.h"
#include "parallel.h"

// --------------------------------------------------------------------------------- //

/**
 * @brief Preconditioner template.
 * 
 * This interface class declares all necessary methods that a valid preconditioner
 * object needs to implement. The preconditioner thus needs to be able to multiply
 * a vector by the matrix of the set of equations that is to be solved and, of course,
 * to precondition the solution by solving the auxiliary preconditioner equation
 * @f$ \mathbf{M}\mathbf{z} = \mathbf{r} @f$.
 */
class PreconditionerBase
{
    public:
        
        //
        // Run-time selection mechanism (object factory).
        //
        
            baseClassRunTimeSelectionDefinitions
            (
                PreconditionerBase,
                (
                    Parallel const & par,
                    InputFile const & inp,
                    AngularBasis const & ll,
                    Bspline const & bspline_inner,
                    Bspline const & bspline_full,
                    CommandLine const & cmd
                )
            )
        
        //
        // Class member functions.
        //
        
            /**
             * @brief Dummy default constructor needed by the run-time selection.
             */
            PreconditionerBase () {}
            
            /**
             * @brief Virtual destructor.
             * 
             * To be overridden in derived classes.
             */
            virtual ~PreconditionerBase () {}
            
            /**
             * @brief Description of the preconditioner.
             * 
             * Simple documentation of the preconditioner.
             */
            virtual std::string description () const { return ""; }
            
            /**
             * @brief Initialize the preconditioner.
             * 
             * This function contains all computation intensive preparations
             * for the preconditioner, e.g. computation of radial integrals.
             * It may use only SMP environment.
             */
            virtual void setup () {}
            
            /**
             * @brief Update the preconditioner for the next energy.
             * 
             * This function updates the preconditioner for another right hand side.
             * It may use the MPI environment.
             */
            virtual void update (Real E) {}
            
            /**
             * @brief Clean up memory etc.
             * 
             * This function is called when the preconditioner will no longer be used
             * to release resources. To re-enable the same preconditioner then would require
             * a new call to @ref setup.
             */
            virtual void finish () {}
            
            /**
             * @brief Calculate the right-hand side.
             */
            virtual void rhs (BlockArray<Complex> & chi, int ie, int instate) const {}
            
            /**
             * @brief Multiply by the matrix equation.
             * 
             * This function implements matrix multiplication by the matrix of
             * the set of equations that is to be solved.
             */
            virtual void multiply (BlockArray<Complex> const & p, BlockArray<Complex> & q, MatrixSelection::Selection tri = MatrixSelection::Both) const {}
            
            /**
             * @brief Precondition the equation.
             * 
             * This function preconditions the equation, solving the preconditioner
             * equation
             * @f[
             *       \mathbf{M}\mathbf{z} = \mathbf{r} \ .
             * @f]
             * It may use the MPI environment.
             */
            virtual void precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const {}
};

// --------------------------------------------------------------------------------- //

#define preconditionerRunTimeSelectionDefinitions(TYPE,NAME) \
    derivedClassRunTimeSelectionDefinitions \
    ( \
        PreconditionerBase, \
        ( \
            Parallel const & par, \
            InputFile const & inp, \
            AngularBasis const & ll, \
            Bspline const & bspline_inner, \
            Bspline const & bspline_full, \
            CommandLine const & cmd \
        ), \
        TYPE, \
        ( \
            par, \
            inp, \
            ll, \
            bspline_inner, \
            bspline_full, \
            cmd \
        ), \
        NAME \
    )

// --------------------------------------------------------------------------------- //

#endif // HEX_ECS_PRECONDITIONERS
