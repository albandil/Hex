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
                    CommandLine  const & cmd,
                    InputFile    const & inp,
                    Parallel     const & par,
                    AngularBasis const & ang,
                    Bspline const & bspline_inner,
                    Bspline const & bspline_full,
                    Bspline const & bspline_panel_x,
                    Bspline const & bspline_panel_y
                )
            )
        
        //
        // Class member functions.
        //
        
            /**
             * @brief Dummy default constructor needed by the run-time selection.
             */
            PreconditionerBase () : verbose_(true) {}
            
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
             * It may use the MPI environment. The energy is in Ry.
             */
            virtual void update (Real E) {}
            
            /**
             * @brief Get the number of allowed bound states.
             * 
             * Returns the number of energetically allowed bound states of the first and second
             * particle with specific angular momenta. This is necessary because it determines
             * the size of the solution vector.
             * 
             * This function can be called only after a call to the @ref setup function,
             * where it is assumed that the initialization of these numbers takes place.
             * 
             * The energy is in Ry.
             */
            virtual std::pair<int,int> bstates (Real E, int l1, int l2) const
            {
                return std::make_pair(0,0);
            }
            
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
            virtual void multiply
            (
                BlockArray<Complex> const & p,
                BlockArray<Complex> & q,
                MatrixSelection::Selection tri = MatrixSelection::BlockBoth | MatrixSelection::Both
            ) const {}
            
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
            
            /**
             * @brief Set verbosity level.
             * 
             * Determine whether the preconditioner will produce a text output.
             */
            virtual void verbose (bool b) { verbose_ = b; }
            
    protected:
        
        bool verbose_;
};

// --------------------------------------------------------------------------------- //

#define preconditionerRunTimeSelectionDefinitions(TYPE,NAME) \
    derivedClassRunTimeSelectionDefinitions \
    ( \
        PreconditionerBase, \
        ( \
            CommandLine  const & cmd, \
            InputFile    const & inp, \
            Parallel     const & par, \
            AngularBasis const & ang, \
            Bspline const & bspline_inner, \
            Bspline const & bspline_full,  \
            Bspline const & bspline_panel_x, \
            Bspline const & bspline_panel_y  \
        ), \
        TYPE, \
        ( \
            cmd, \
            inp, \
            par, \
            ang, \
            bspline_inner, \
            bspline_full,  \
            bspline_panel_x, \
            bspline_panel_y  \
        ), \
        NAME \
    )

// --------------------------------------------------------------------------------- //

#endif // HEX_ECS_PRECONDITIONERS
