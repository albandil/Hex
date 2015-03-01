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

#ifndef HEX_PRECONDITIONERS
#define HEX_PRECONDITIONERS

#include <tuple>

#include "arrays.h"
#include "matrix.h"
#include "radial.h"
#include "parallel.h"

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
        
        /**
         * @brief Initialize the preconditioner.
         * 
         * This function contains all computation intensive preparations
         * for the preconditioner, e.g. computation of radial integrals.
         * It may use only SMP environment.
         */
        virtual void setup () = 0;
        
        /**
         * @brief Update the preconditioner for the next energy.
         * 
         * This function updates the preconditioner for another right hand side.
         * It may use the MPI environment.
         */
        virtual void update (double E) = 0;
        
        /**
         * @brief Calculate the right-hand side.
         */
        virtual void rhs (BlockArray<Complex> & chi, int ie, int instate, int Spin) const = 0;
        
        /**
         * @brief Multiply by the matrix equation.
         * 
         * This function implements matrix multiplication by the matrix of
         * the set of equations that is to be solved.
         */
        virtual void multiply (BlockArray<Complex> const & p, BlockArray<Complex> & q) const = 0;
        
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
        virtual void precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const = 0;
};

class NoPreconditioner;
class CGPreconditioner;
class KPACGPreconditioner;
class ILUCGPreconditioner;
class GPUCGPreconditioner;

#include "preconditioners/NoPreconditioner.h"
#include "preconditioners/CGPreconditioner.h"
#include "preconditioners/KPAPreconditioner.h"
#include "preconditioners/ILUPreconditioner.h"
#include "preconditioners/GPUPreconditioner.h"

/**
 * @brief Preconditioner traits.
 * 
 * This class is used for accessing all available preconditioners. Preconditioner types and objects
 * are never accessed in any other way than by calling the member functions of this class or through
 * the pointer retrieved from the method "choose".
 * 
 * For this reason, addition of a new preconditioner is very simple. One just needs to derive the new
 * preconditioner class from PreconditionerBase, specify a name to be used on command line as static
 * attribute (see code below) and add the class name to the tuple AvailableTypes.
 * 
 * @code
 *     class MyNewPreconditioner
 *     {
 *         public:
 *             static const std::string name;
 * 
 *             // other compulsory methods
 *             // ...
 *     };
 * @endcode
 */
class Preconditioners
{
    public:
        
        /**
         * @brief List of available preconditioners.
         * 
         * First one (ILUCGPreconditioner at the moment) is considered default. This tuple is the only place that has to be modified
         * when adding a new preconditioner to Hex-ECS (besides the obvious declaration and definition somewhere in reach, ideally
         * in "preconditioners.h" and "preconditioners.cpp").
         */
        typedef std::tuple <
            ILUCGPreconditioner         // Solve diagonal blocks by drop-tolerance incomplete LU factorization.
            , NoPreconditioner          // No preconditioner.
            , CGPreconditioner          // Solve diagonal blocks by non-preconditioned CG iterations.
#ifndef NO_LAPACK
            , KPACGPreconditioner       // Solve diagonal blocks by separate electrons preconditioned CG iterations.
#endif
#ifndef NO_OPENCL
            , GPUCGPreconditioner         // KPA implemented on GPU.
#endif
        > AvailableTypes;
        
        /**
         * @brief Number of available preconditioners.
         */
        static std::size_t size ()
        {
            return std::tuple_size<AvailableTypes>::value;
        }
        
        /**
         * @brief Return pointer to new preconditioner object.
         * 
         * The preconditioner index from CommandLine::preconditioner is used
         * to choose the correct preconditioner.
         */
        //@{
        template <int i = 0> static inline typename std::enable_if<(i < std::tuple_size<AvailableTypes>::value), PreconditionerBase*>::type choose
        (
            Parallel const & par,
            InputFile const & inp,
            std::vector<std::pair<int,int>> const & ll,
            Bspline const & bspline,
            CommandLine const & cmd
        )
        {
            // check if i-th preconditioner is the right one
            // - Yes : return new pointer of its type
            // - No : try the next preconditioner
            return   (i == cmd.preconditioner)
                   ? new typename std::tuple_element<i,AvailableTypes>::type (par, inp, ll, bspline, cmd)
                   : choose<i+1>(par, inp, ll, bspline, cmd);
        }
        template <int i = 0> static inline typename std::enable_if<(i == std::tuple_size<AvailableTypes>::value), PreconditionerBase*>::type choose
        (
            Parallel const & par,
            InputFile const & inp,
            std::vector<std::pair<int,int>> const & ll,
            Bspline const & bspline,
            CommandLine const & cmd
        )
        {
            // we visited all available preconditioners without success : return NULL pointer
            return nullptr;
        }
        //@}
        
        /**
         * @brief Find preconditioner by name.
         * 
         * Return index of the corresponding preconditioner in AvailableTypes or (-1) if
         * none such exists. "Corresponding preconditioner" is such type that has a static
         * member
         * @code
         *     std::string name;
         * @endcode
         * equal to the argument of this function.
         */
        //@{
        template <int i = 0> static inline typename std::enable_if<(i < std::tuple_size<AvailableTypes>::value),int>::type findByName (std::string name)
        {
            // check if i-th preconditioner is the right one
            // - Yes : return new pointer of its type
            // - No : try the next preconditioner
            return (name == std::tuple_element<i,AvailableTypes>::type::name) ? i : findByName<i+1>(name);
        }
        template <int i = 0> static inline typename std::enable_if<(i == std::tuple_size<AvailableTypes>::value),int>::type findByName (std::string name)
        {
            // we visited all available preconditioners without success : return invalid index (-1)
            return -1;
        }
        //@}
        
        /**
         * @brief Get name of the preconditioner.
         */
        //@{
        template <int i = 0> static inline typename std::enable_if<(i < std::tuple_size<AvailableTypes>::value),std::string>::type name (unsigned idx)
        {
            return (i == idx) ? std::tuple_element<i,AvailableTypes>::type::name : name<i+1>(idx);
        }
        template <int i = 0> static inline typename std::enable_if<(i == std::tuple_size<AvailableTypes>::value),std::string>::type name (unsigned idx)
        {
            return "";
        }
        //@}
        
        /**
         * @brief Get description of the preconditioner.
         */
        //@{
        template <int i = 0> static inline typename std::enable_if<(i < std::tuple_size<AvailableTypes>::value),std::string>::type description (unsigned idx)
        {
            return (i == idx) ? std::tuple_element<i,AvailableTypes>::type::description : description<i+1>(idx);
        }
        template <int i = 0> static inline typename std::enable_if<(i == std::tuple_size<AvailableTypes>::value),std::string>::type description (unsigned idx)
        {
            return "";
        }
        //@}
};
#endif
