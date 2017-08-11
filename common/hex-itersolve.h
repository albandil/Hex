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

#ifndef HEX_ITERSOLVE
#define HEX_ITERSOLVE

#include <chrono>
#include <iostream>

#include "hex-arrays.h"
#include "hex-matrix.h"
#include "hex-misc.h"
#include "hex-special.h"

/**
 * @brief Process the solution.
 * 
 * This function is called every time the solution is updated.
 */
template <class TArrayView> void default_process_solution (unsigned iteration, const TArrayView r)
{
    // no processing by default
}

/**
 * @brief Constrain the residual vector.
 */
template <class TArrayView> void default_constraint (const TArrayView x, TArrayView r)
{
    // no constraint by default
}

/**
 * @brief Allocate and return a new working array of length 'n'.
 */
template <class TArray>  [[noreturn]] TArray default_new_array (std::size_t n, std::string name)
{
    HexException("CG: You need to specify \"new_array\" routine.");
}

/**
 * @brief Calculate the scalar product of two arrays.
 */
template <class T, class TArrayView> [[noreturn]] T default_scalar_product (const TArrayView x, const TArrayView y)
{
    HexException("CG: You need to specify \"scalar_product\" routine.");
}

/**
 * @brief Compute the standard 2-norm of an array.
 */
template <class TArrayView> [[noreturn]] double default_compute_norm (const TArrayView x)
{
    HexException("CG: You need to specify \"compute_norm\" routine.");
}

/**
 * @brief Do the @f$ \alpha x + \beta y @f$ operation.
 */
template <class T, class TArrayView> [[noreturn]] void default_axby (T a, TArrayView x, T b, const TArrayView y)
{
    HexException("CG: You need to specify \"axby\" routine.");
}

/**
 * @brief Conjugate gradients solver.
 * 
 * Callback-based conjugate gradients. There is a variety of template arguments and
 * function parameters that aim at generality so that the function can be used with
 * various implementations of array arithmetic. The numerically intensive sections
 * are being computed by user-supplied routines. Also the array types have a template
 * character; they are
 * - TArray: A stand-alone array type that can be allocated. Has to support
 *           the methods "norm" (computation of 2-norm), "size" (length of the array)
 *           and "data" (pointer to the raw data).
 * - TArrayView: Either a reference to TArray, or a different compatible type that
 *               supports the three listed methods. It is not necessary to be allocable.
 * 
 * 
 * - NewArray: routine that allocates an array and possibly does something more, e.g.
 *             registers the array at the computing device (GPU). NewArray has to be
 *             compatible with the signature
 *   @code
 *       TArray (*NewArray)
 *   @endcode
 * - Preconditioner: routine that 
 * 
 * @param b Right hand side.
 * @param x Output vector. An initial guess may be present at beginning.
 * @param eps Relative tolerance.
 * @param min_iterations Minimal iterations count.
 * @param max_iterations Maximal iterations count.
 * @param apply_preconditioner Functor compatible with the signature
 *        @code
 *            void (*) (const TArrayView, TArrayView)
 *        @endcode
 *        Applies a custom preconditioner to the vector given as first argument
 *        and stores the result in the second argument.
 * @param matrix_multiply Functor compatible with the signature
 *        @code
 *            void (*) (const TArrayView, TArrayView)
 *        @endcode
 *        Multiplies the vector given in the first argument by the matrix of the
 *        equation set. The result will be stored in the second argument.
 * @param verbose Whether to display progress information.
 * @param new_array Functor compatible with the signature
 *        @code
 *            TArray (*) (size_t)
 *        @endcode
 *        Allocates (and preprocesses, if needed) a new array of type TArray and
 *        returns a copy of the array. This can be useful when the new arrays need
 *        to be uploaded to a different computing device (mostly GPU).
 * @param axby Functor compatible with the signature
 *        @code
 *            void (*) (Complex a, const TArrayView x, Complex b, const TArrayView y)
 *        @endcode
 *        Stores the linear combination @f$ a\mathbf{x} + b\mathbf{y} @f$ into the
 *        array @f$ \mathbf{x} @f$.
 * @param scalar_product Functor compatible with the signature
 *        @code
 *            Complex (*) (const TArrayView x, const TArrayView y)
 *        @endcode
 *        Computes the scalar product (without complex conjugation), i.e. the sum
 *        @f[
 *                 \sum_{i = 1}^N x_i y_i \ .
 *        @f]
 * @param compute_norm Functor compatible with signature
 *        @code
 *            double (*) (const TArrayView x)
 *        @endcode
 *        Computes the norm of the array.
 * 
 * @return Iteration count.
 */
template
<
    class T = double,
    class TArray = rArray,
    class TArrayView = rArray&
>
class ConjugateGradients
{
    public:
        
        ConjugateGradients ()
            : verbose(false), 
              compute_norm(default_compute_norm<TArrayView>),
              scalar_product(default_scalar_product<T,TArrayView>),
              axby(default_axby<T,TArrayView>),
              new_array(default_new_array<TArray>),
              constrain(default_constraint<TArrayView>),
              process_solution(default_process_solution<TArrayView>),
              k(1), recovered(false)
        {
            // nothing
        }
        
        // verbosity control
        bool verbose;
        
        // callback functions
        std::function<void (const TArrayView, TArrayView)> apply_preconditioner;
        std::function<void (const TArrayView, TArrayView)> matrix_multiply;
        std::function<double (TArrayView)> compute_norm;
        std::function<T (TArrayView, TArrayView)> scalar_product;
        std::function<void (T, const TArrayView, T, TArrayView)> axby;
        std::function<TArray (std::size_t,std::string)> new_array;
        std::function<void (const TArrayView, TArrayView)> constrain;
        std::function<void (unsigned, const TArrayView)> process_solution;
        
        unsigned solve
        (
            const TArrayView b,
            TArrayView x,
            double eps,
            unsigned min_iterations,
            unsigned max_iterations
        )
        {
            Timer timer (time_offset);
            ok = false;
            
            // compute norm of the right hand side
            double bnorm = compute_norm(b), rnorm;
            
            // trivial solution for zero right hand side
            if (bnorm == 0)
            {
                axby(0., x, 0., b); // x = 0 b
                ok = true;
                return 0;
            }
            
            // get size of the problem
            std::size_t N = b.size();
            
            // residual; initialized to starting residual using the initial guess
            TArray r (new_array(N, "cg-r"));
            if (not recovered)
            {
                matrix_multiply(x, r); // r = A x
                axby(-1., r, 1., b); // r = b - r
            }
//             constrain(r);
            rnorm = compute_norm(r);
            recovered = false;
            
            // terminate if the initial guess is already fine enough
            if (rnorm / bnorm < eps)
            {
                ok = true;
                return 0;
            }
            
            // if the (non-zero) initial guess seems horribly wrong,
            //    use rather the right hand side as the initial guess
            if (rnorm / bnorm > 1000)
            {
                axby(0., x, 0., b); // x = 0 b
                axby(0., r, 1., b); // r = b
//                 constrain(r);
            }
            
            // some auxiliary arrays (search directions etc.)
            TArray p (new_array(N, "cg-p"));
            TArray z (new_array(N, "cg-z"));
            
            // Iterate
            
            for (; k <= max_iterations; k++)
            {
                if (verbose)
                {
                    std::cout << '\t';
                    std::cout << std::setw(4) << std::right << k;
                    std::cout << " | ";
                    std::cout << std::setw(11) << std::left << timer.nice_time();
                    std::cout << " | ";
                    std::cout << std::setw(15) << std::left << rnorm / bnorm;
                }
                
                time_offset = timer.seconds();
                
                // apply desired preconditioner
                apply_preconditioner(r, z); // z = M⁻¹r
                
                // compute projection ρ = r·z
                rho_new = scalar_product(r, z);
                
                // setup search direction p
                if (k == 1)
                {
                    axby(0., p, 1., z); // p = z
                }
                else
                {
                    beta = rho_new / rho_old;
                    if (not std::isfinite(std::abs(beta)))
                    {
                        std::cout << "\t    Warning: Iterative method breakdown: (r|z) -> 0." << std::endl;
                        ok = false;
                        return k;
                    }
                    axby(beta, p, 1, z); // p = beta p + z
                }
                
                // move to next Krylov subspace by multiplying A·p
                matrix_multiply(p, z);
                
                // compute projection ratio α
                alpha = rho_new / scalar_product(p, z);
                if (not std::isfinite(std::abs(alpha)))
                {
                    std::cout << "\t    Warning: Iterative method breakdown: (p|z) -> 0." << std::endl;
                    ok = false;
                    return k;
                }
                
                // update the solution and the residual
                axby(1., x, alpha, p); // x = x + α p
                axby(1., r, -alpha, z); // r = r - α z
                process_solution(k, x);
                constrain(x, r);
                
                // compute and check norm
                rnorm = compute_norm(r);
                residual = rnorm / bnorm;
                if (verbose)
                    std::cout << std::endl;
                if (not std::isfinite(residual))
                {
                    std::cout << "\t     The norm of the solution is not finite. Something went wrong!" << std::endl;
                    ok = false;
                    break;
                }
                
                // check convergence, but always do at least "min_iterations" iterations
                if (k >= min_iterations and residual < eps)
                {
                    ok = true;
                    if (verbose)
                        std::cout << "\tConvergence reached, final residual " << residual << std::endl;
                    break;
                }
                
                // move to the next iteration: store previous projection
                rho_old = rho_new;
            }
            
            return k;
        }
        
        TArray r, p, z;
        
        void dump () const
        {
            std::ofstream out ("cg.dat");
            out << k << std::endl;
            out << rho_old.real() << std::endl;
            out << rho_old.imag() << std::endl;
            out << time_offset << std::endl;
        }
        
        void recover ()
        {
            std::ifstream in ("cg.dat");
            if (in.is_open())
            {
                in >> k;
                double x;
                in >> x; rho_old.real(x);
                in >> x; rho_old.imag(x);
                in >> time_offset;
                recovered = true;
            }
            else
            {
                reset();
            }
        }
        
        void reset ()
        {
            recovered = false;
            k = 1;
            rho_new = rho_old = alpha = beta = 0;
            time_offset = 0;
        }
        
        Complex rho_new;
        Complex rho_old;
        Complex alpha, beta;
        
        unsigned k;
        unsigned time_offset;
        bool recovered;
        double residual;
        bool ok;
};

#endif
