/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2014                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HEX_ITERSOLVE
#define HEX_ITERSOLVE

#include <chrono>
#include <iostream>

#include "arrays.h"
#include "complex.h"
#include "matrix.h"
#include "misc.h"

/**
 * @brief Return new complex array.
 * 
 * Create (and return a copy of) a NumberArray<Complex> object. This is used as
 * a default method of creating an array in cg_callbacks. It can be substituted
 * by a different method, if necessary; for example when the created array needs
 * to be registered at GPU first.
 */
inline cArray default_new_complex_array (size_t n)
{
    return cArray(n);
}

/**
 * @brief Compute norm of an array.
 * 
 * This routine is the default way of how to compute a norm of an object. It can be
 * overloaded by some more sophisticated implementation (e.g. BLAS or OpenCL version).
 */
inline double default_compute_norm (const cArrayView x)
{
    return x.norm();
}

/**
 * @brief Do the @f$ \alpha x + \beta y @f$ operation.
 * 
 * Computes a linear combination of two vectors and stores the output in the first
 * vector. This is a default implementation of the method. It can be substituted by
 * another, if necessary; for example one could call specialized routine from BLAS
 * or use a GPU kernel.
 */
inline void default_complex_axby (Complex a, cArrayView x, Complex b, const cArrayView y)
{
    size_t N = x.size();
    assert (N == y.size());
    
    // accelerators
    Complex       * const restrict px = x.data();
    Complex const * const restrict py = y.data();
    
    // do the axby per element
    for (size_t i = 0; i < N; i++)
        px[i] = a * px[i] + b * py[i];
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
    class TArray,
    class TArrayView,
    class Preconditioner,
    class MatrixMultiplication,
    class NewArray      = decltype(default_new_complex_array),
    class AxbyOperation = decltype(default_complex_axby),
    class ScalarProduct = decltype(operator|<Complex>),
    class ComputeNorm   = decltype(default_compute_norm)
>
unsigned cg_callbacks
(
        const TArrayView b,
              TArrayView x,
                  double eps,
                unsigned min_iterations,
                unsigned max_iterations,
          Preconditioner apply_preconditioner,
    MatrixMultiplication matrix_multiply,
                    bool verbose        = true,
                NewArray new_array      = default_new_complex_array,
           AxbyOperation axby           = default_complex_axby,
           ScalarProduct scalar_product = operator|<Complex>,
             ComputeNorm compute_norm   = default_compute_norm
)
{
    Timer timer;
    
    // compute norm of the right hand side
    double bnorm = compute_norm(b);
    
    // trivial solution for zero right hand side
    if (bnorm == 0)
    {
        x.fill(0.);
        return 0;
    }
    
    // get size of the problem
    size_t N = b.size();
    
    // residual; initialized to starting residual using the initial guess
    TArray r (std::move(new_array(N)));
    matrix_multiply(x, r); // r = A x
    axby (-1., r, 1., b); // r = b - r
    double rnorm = compute_norm(r);
    
    // terminate if the initial guess is already fine enough
    if (rnorm / bnorm < eps)
        return 0;
    
    // if the (non-zero) initial guess seems horribly wrong,
    //    use rather the right hand side as the initial guess
    if (rnorm / bnorm > 1000)
    {
        x.fill(0.);
        axby (0., r, 1., b); // r = b
    }
    
    // some auxiliary arrays (search directions etc.)
    TArray p (std::move(new_array(N)));
    TArray q (std::move(new_array(N)));
    TArray z (std::move(new_array(N)));
    
    // some other scalar variables
    Complex rho_new;        // contains inner product r_i^T · r_i
    Complex rho_old;        // contains inner product r_{i-1}^T · r_{i-1}
    Complex alpha, beta;    // contains projection ratios
    
    // Iterate
    
    unsigned k;
    for (k = 0; k < max_iterations; k++)
    {
        int sec = timer.seconds();
        
        if (verbose)
        {
            std::cout << "\t[cg] Residual relative magnitude after "
                    << k << " iterations: " << rnorm / bnorm
                    << " (" << sec / 60 << " min)\n";
        }
        
        // apply desired preconditioner
        apply_preconditioner(r, z); // z = M⁻¹r
        
        // compute projection ρ = r·z
        rho_new = scalar_product(r, z);
        
        // setup search direction p
        if (k == 0)
        {
            axby (0., p, 1., z); // p = z
        }
        else
        {
            beta = rho_new / rho_old;
            axby (beta, p, 1, z); // p = beta p + z
        }
        
        // move to next Krylov subspace by multiplying A·p
        matrix_multiply(p, q);
        
        // compute projection ratio α
        alpha = rho_new / scalar_product(p, q);
        
        // update the solution and the residual
        axby (1., x, alpha, p); // x = x + α p
        axby (1., r, -alpha, q); // r = r - α q
        
        // compute and check norm
        rnorm = compute_norm(r);
        if (not std::isfinite(rnorm))
        {
            std::cout << "\t[cg] Damn... the norm of the solution is not finite. Something went wrong!\n";
            break;
        }
        
        // check convergence, but do at least "min_iterations" iterations
        if (k >= min_iterations and rnorm / bnorm < eps)
            break;
        
        // move to the next iteration: store previous projection
        rho_old = rho_new;
    }
    
    return k;
}

/**
 * Callback-based BiCGSTAB.
 * 
 * Bi-Conjugate gradients stabilized method.
 * 
 * @param b Right hand side.
 * @param x Output vector. An initial guess may be present at beginning.
 * @param eps Relative tolerance.
 * @param min_iterations Minimal iterations count.
 * @param max_iterations Maximal iterations count.
 * @param apply_preconditioner Functor compatible with void(*)(const Array&, Array&) prototype.
 *                             Should apply a custom preconditioner to a given vector.
 * @param matrix_multiply Functor compatible with void(*)(const Array&, Array&) prototype.
 *                        Should multiply the given vector by the matrix of the equation set.
 * @param verbose Whether to comment the progress to stdout.
 * @return Iteration count.
 */
template <typename TFunctor1, typename TFunctor2>
int bicgstab_callbacks (
    const cArrayView b, cArrayView x,
    double eps,
    int min_iterations, int max_iterations,
    TFunctor1 apply_preconditioner,
    TFunctor2 matrix_multiply,
    bool verbose = false
) {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    std::chrono::duration<int> sec;
    
    int N = b.size();
    double bnorm = b.norm();
    
    cArray x_im1(N), r_im1(N), rt(N), p_i(N), p_im1(N), v_im1(N), phat(N), v_i(N), s(N), shat(N), t(N), x_i(N), r_i(N);
    Complex rho_im1, rho_im2, beta, alpha_i, alpha_im1, omega_i, omega_im1;
    
    x_im1 = x;
    matrix_multiply(x_im1,r_im1);
    rt = r_im1 = b - r_im1;
    
    int i;
    for (i = 1; i < max_iterations; i++)
    {
        sec = std::chrono::duration_cast<std::chrono::duration<int>>(std::chrono::steady_clock::now()-start);
        
        if (verbose)
        {
            std::cout << "\t[Bi-CGSTAB] Residual relative magnitude after "
                    << i << " iterations: " << r_im1.norm() / bnorm
                    << " (" << sec.count()/60 << " min)\n";
        }
        
        rho_im1 = (rt | r_im1);
        if (std::abs(rho_im1) == 0.)
            throw exception ("[Bi-CGSTAB] Failed, rho = 0.");
        
        if (i == 1)
        {
            p_i = r_im1;
        }
        else
        {
            beta = (rho_im1 / rho_im2) * (alpha_im1 / omega_im1);
            p_i = r_im1 + beta * (p_im1 - omega_im1 * v_im1);
        }
        
        apply_preconditioner(p_i, phat);
        matrix_multiply(phat, v_i);
        alpha_i = rho_im1 / (rt | v_i);
        s = r_im1 - alpha_i * v_i;
        
        if (s.norm() < eps * bnorm)
        {
            x = x_im1 + alpha_i * phat;
            break;
        }
        
        apply_preconditioner(s, shat);
        matrix_multiply(shat, t);
        omega_i = (t|s) / (t|t);
        
        x_i = x_im1 + alpha_i * phat + omega_i * s;
        r_i = s - omega_i * t;
        
        if (r_i.norm() < eps * bnorm)
        {
            x = x_i;
            break;
        }
        
        if (omega_i == 0.)
            throw exception ("[Bi-CGSTAB] Solver failed, ω = 0.");
        
        // shift vectors
        x_im1 = std::move(x_i);
        r_im1 = std::move(r_i);
        p_im1 = std::move(p_i);
        v_im1 = std::move(v_i);
        
        // shift 
        rho_im2 = rho_im1;
        alpha_im1 = alpha_i;
        omega_im1 = omega_i;
    }
    
    return i;
}

/**
 *  @brief CGS solver.
 * 
 * Conjugate gradients squared,
 * 
 * @param b Right hand side.
 * @param x Output vector. An initial guess may be present at beginning.
 * @param eps Relative tolerance.
 * @param min_iterations Minimal iterations count.
 * @param max_iterations Maximal iterations count.
 * @param apply_preconditioner Functor compatible with void(*)(const Array&, Array&) prototype.
 *                             Should apply a custom preconditioner to a given vector.
 * @param matrix_multiply Functor compatible with void(*)(const Array&, Array&) prototype.
 *                        Should multiply the given vector by the matrix of the equation set.
 * @param verbose Whether to comment the progress to stdout.
 * @return Iteration count.
 */
template <typename TFunctor1, typename TFunctor2>
int cgs_callbacks (
    const cArrayView b, cArrayView x,
    double eps,
    int min_iterations, int max_iterations,
    TFunctor1 apply_preconditioner,
    TFunctor2 matrix_multiply,
    bool verbose = false
) {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    std::chrono::duration<int> sec;
    
    int N = b.size();
    double bnorm = b.norm();
    
    Complex resid, alpha, beta, rho_1, rho_2;
    cArray r(N), rt(N), p(N), phat(N), q(N), qhat(N), vhat(N), u(N), uhat(N);
    
    matrix_multiply(x,r);
    rt = r = b - r;
    
    if (bnorm == 0.)
        bnorm = 1;
    
    if (r.norm() < eps * bnorm)
        return 0;
    
    int i;
    for (i = 1; i < max_iterations; i++)
    {
        sec = std::chrono::duration_cast<std::chrono::duration<int>>(std::chrono::steady_clock::now()-start);
        
        if (verbose)
        {
            std::cout << "\t[cgs] Residual relative magnitude after "
                    << i << " iterations: " << r.norm() / bnorm
                    << " (" << sec.count()/60 << " min)\n";
        }
        
        rho_1 = (rt | r);
        
        if (rho_1 == 0.)
        {
            throw exception ("[cgs] Solver failes, ρ = 0.");
        }
        if (i == 1)
        {
            u = r;
            p = u;
        }
        else
        {
            beta = rho_1 / rho_2;
            u = r + beta * q;
            p = u + beta * (q + beta * p);
        }
        
        apply_preconditioner(p, phat);
        matrix_multiply(phat, vhat);
        
        alpha = rho_1 / (rt | vhat);
        q = u - alpha * vhat;
        
        apply_preconditioner(u + q, uhat);
        matrix_multiply(uhat, qhat);
        
        x += alpha * uhat;
        r -= alpha * qhat;
        rho_2 = rho_1;
        
        if (r.norm() < eps * bnorm)
            break;
    }
    
    return i;
}

#endif
