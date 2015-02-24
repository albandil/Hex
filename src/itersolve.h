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

#include "arrays.h"
#include "matrix.h"
#include "misc.h"
#include "special.h"

/**
 * @brief Return new complex array.
 * 
 * Create (and return a copy of) a NumberArray<Complex> object. This is used as
 * a default method of creating an array in cg_callbacks. It can be substituted
 * by a different method, if necessary; for example when the created array needs
 * to be registered at GPU first.
 */
inline cArray default_new_complex_array (std::size_t n)
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
    std::size_t N = x.size();
    assert (N == y.size());
    
    // accelerators
    Complex       * const restrict px = x.data();
    Complex const * const restrict py = y.data();
    
    // do the axby per element
    for (std::size_t i = 0; i < N; i++)
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
    class ComputeNorm   = decltype(default_compute_norm),
    class ScalarProduct = decltype(operator|<Complex>),
    class AxbyOperation = decltype(default_complex_axby),
    class NewArray      = decltype(default_new_complex_array)
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
             ComputeNorm compute_norm   = default_compute_norm,
           ScalarProduct scalar_product = operator|<Complex>,
           AxbyOperation axby           = default_complex_axby,
                NewArray new_array      = default_new_complex_array
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
    std::size_t N = b.size();
    
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
    TArray z (std::move(new_array(N)));
    
    // some other scalar variables
    Complex rho_new;        // contains inner product r_i^T · r_i
    Complex rho_old;        // contains inner product r_{i-1}^T · r_{i-1}
    Complex alpha, beta;    // contains projection ratios
    
    // Iterate
    
    unsigned k;
    for (k = 1; k <= max_iterations; k++)
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
        
        // apply desired preconditioner
        apply_preconditioner(r, z); // z = M⁻¹r
        
        // compute projection ρ = r·z
        rho_new = scalar_product(r, z);
        
        // setup search direction p
        if (k == 1)
        {
            axby (0., p, 1., z); // p = z
        }
        else
        {
            beta = rho_new / rho_old;
            axby (beta, p, 1, z); // p = beta p + z
        }
        
        // move to next Krylov subspace by multiplying A·p
        matrix_multiply(p, z);
        
        // compute projection ratio α
        alpha = rho_new / scalar_product(p, z);
        
        // update the solution and the residual
        axby (1., x, alpha, p); // x = x + α p
        axby (1., r, -alpha, z); // r = r - α z
        
        // compute and check norm
        rnorm = compute_norm(r);
        if (verbose)
            std::cout << std::endl;
        if (not std::isfinite(rnorm))
        {
            std::cout << "\t     The norm of the solution is not finite. Something went wrong!" << std::endl;
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
template
<
    class TArray,
    class TArrayView,
    class Preconditioner,
    class MatrixMultiplication
>
unsigned bicgstab_callbacks
(
        const TArrayView b,
              TArrayView x,
                  double eps,
                unsigned min_iterations,
                unsigned max_iterations,
          Preconditioner apply_preconditioner,
    MatrixMultiplication matrix_multiply,
                    bool verbose
)
{
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    std::chrono::duration<int> sec;
    
    int N = b.size();
    double bnorm = b.norm();
    
    cArray p(N), r(N), rt(N), s(N), t(N), v(N), y(N), z(N);
    Complex alpha = 1, beta = 1, rho = 1, rho_prev = 1, omega = 1;
    
    matrix_multiply(x,r);
    rt = r = b - r;
    
    double rnorm = r.norm();
    
    Timer timer;
    
    int i;
    for (i = 1; i < max_iterations; i++)
    {
        sec = std::chrono::duration_cast<std::chrono::duration<int>>(std::chrono::steady_clock::now()-start);
        
        if (verbose)
        {
            std::cout << '\t';
            std::cout << std::setw(4) << std::right << i;
            std::cout << " | ";
            std::cout << std::setw(11) << std::left << timer.nice_time();
            std::cout << " | ";
            std::cout << std::setw(15) << std::left << rnorm / bnorm;
        }
        
        // 1.
        rho = (rt | r);
        
        // 2.
        beta = (rho / rho_prev) * (alpha / omega);
        
        // 3.
        p = r + beta * (p - omega * v);
        
        // 4.
        apply_preconditioner(p, y);
        
        // 5.
        matrix_multiply(y, v);
        
        // 6.
        alpha = rho / (rt | v);
        
        // 7.
        s = r - alpha * v;
        
        // 8.
        apply_preconditioner(s, z);
        
        // 9.
        matrix_multiply(z, t);
        
        // 10. (assumming K₁ = Id)
        omega = (t | s) / (t | t);
        
        // 11.
        x += alpha * y + omega * z;
        
        // 12. (check accuracy of x)
        matrix_multiply(x, y);
        rnorm = (y - b).norm();
        if (verbose)
            std::cout << std::endl;
        if (not std::isfinite(rnorm))
        {
            std::cout << "\t     The norm of the solution is not finite. Something went wrong!" << std::endl;
            break;
        }
        if (rnorm < eps * bnorm)
            return i;
        
        // 13.
        r = s - omega * t;
        
        // Update indices.
        rho_prev = rho;
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
int cgs_callbacks
(
    const cArrayView b, cArrayView x,
    double eps,
    int min_iterations, int max_iterations,
    TFunctor1 apply_preconditioner,
    TFunctor2 matrix_multiply,
    bool verbose = false
)
{
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
            HexException("[cgs] Solver failes, ρ = 0.");
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

template
<
    class TArray,
    class TArrayView,
    class Preconditioner,
    class MatrixMultiplication,
    class ComputeNorm   = decltype(default_compute_norm),
    class AxbyOperation = decltype(default_complex_axby),
    class ScalarProduct = decltype(operator|<Complex>),
    class NewArray      = decltype(default_new_complex_array)
>
int minres_callbacks
(
        const TArrayView b,
              TArrayView x,
                  double tol,
                unsigned min_iterations,
                unsigned max_iterations,
          Preconditioner apply_preconditioner,
    MatrixMultiplication matrix_multiply,
                    bool verbose        = true,
             ComputeNorm compute_norm   = default_compute_norm,
           AxbyOperation axby           = default_complex_axby,
           ScalarProduct scalar_product = operator|<Complex>,
                NewArray new_array      = default_new_complex_array
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
    std::size_t N = b.size();
    
    std::vector<std::string> msg(11);
    msg[0]  = " beta1 = 0.  The exact solution is  x = 0 ";
    msg[1]  = " A solution to Ax = b was found, given tol ";
    msg[2]  = " A least-squares solution was found, given tol ";
    msg[3]  = " Reasonable accuracy achieved, given eps ";
    msg[4]  = " x has converged to an eigenvector ";
    msg[5]  = " acond has exceeded 0.1/eps ";
    msg[6]  = " The iteration limit was reached ";
    msg[7]  = " A  does not define a symmetric matrix ";
    msg[8]  = " M  does not define a symmetric matrix ";
    msg[9]  = " M  does not define a pos-def preconditioner ";
    msg[10] = " beta2 = 0.  If M = I, b and x are eigenvectors ";

    int istop(0), itn(0);
    double Anorm(0.0), Acond(0.0), Arnorm(0.0);
    double rnorm(0.0), ynorm(0.0);
    bool done(false);

    // Step 1
    /*
     * Set up y and v for the first Lanczos vector v1.
     * y = beta1 P' v1, where P = C^(-1).
     * v is really P'v1
     */
    
    // residual; initialized to starting residual using the initial guess
    TArray r1 (std::move(new_array(N)));
    matrix_multiply(x, r1); // r = A x
    axby (-1., r1, 1., b); // r = b - r
    rnorm = compute_norm(r1);
    
    if (verbose)
    {
        std::cout << '\t';
        std::cout << std::setw(4) << std::right << 0;
        std::cout << " | ";
        std::cout << std::setw(11) << std::left << timer.nice_time();
        std::cout << " | ";
        std::cout << std::setw(15) << std::left << rnorm / bnorm;
    }
    
    // apply preconditioner
    TArray y (std::move(new_array(N)));
    apply_preconditioner(r1, y);
    
    if (verbose)
            std::cout << std::endl;

    double beta1 = std::abs(scalar_product(r1, y));
    
    // Test for an indefined preconditioner
    // If b = 0 exactly stop with x = x0.

    if(beta1 < 0.0)
    {
        istop = 9;
        done = true;
    }
    else
    {
        if(beta1 == 0.0)
        {
            done = true;
        }
        else
            beta1 = std::sqrt(beta1); // Normalize y to get v1 later
    }
    
    // STEP 2
    /* Initialize other quantities */
    double oldb = 0., beta = beta1, dbar = 0., epsln = 0., oldeps = 0.;
    double qrnorm = beta1, phi = 0., phibar = beta1, rhs1 = beta1;
    double rhs2 = 0., tnorm2 = 0., ynorm2 = 0.;
    double cs = -1., sn = 0.;
    double gmax = 0., gmin = std::numeric_limits<double>::max();
    double alpha = 0., gamma = 0.;
    double delta = 0., gbar = 0.;
    double z = 0.;
    
    TArray v(std::move(new_array(N)));
    TArray w(std::move(new_array(N)));
    TArray w1(std::move(new_array(N)));
    TArray w2(std::move(new_array(N)));
    TArray r2(std::move(new_array(N)));
    
    // r2 = r1
    axby(0., r2, 1., r1);
    
    /* Main Iteration */
    if (!done)
    {
        for (itn = 0; itn < max_iterations; ++itn)
        {
            // STEP 3
            /*
            -----------------------------------------------------------------
            Obtain quantities for the next Lanczos vector vk+1, k = 1, 2,...
            The general iteration is similar to the case k = 1 with v0 = 0:
            
            p1      = Operator * v1  -  beta1 * v0,
            alpha1  = v1'p1,
            q2      = p2  -  alpha1 * v1,
            beta2^2 = q2'q2,
            v2      = (1/beta2) q2.
            
            Again, y = betak P vk,  where  P = C**(-1).
            .... more description needed.
            -----------------------------------------------------------------
             */
            double s = 1./beta; //Normalize previous vector (in y)
            axby(0., v, s, y);
            
            matrix_multiply(v,y);
            
            if (itn > 0)
                axby(1., y, -beta/oldb, r1);
            
            alpha = std::abs(scalar_product(v,y));   // alphak
            axby(1., y, -alpha/beta, r2);
            axby(0., r1, 1., r2);
            axby(0., r2, 1., y);
            
            if (verbose)
            {
                std::cout << '\t';
                std::cout << std::setw(4) << std::right << itn + 2;
                std::cout << " | ";
                std::cout << std::setw(11) << std::left << timer.nice_time();
                std::cout << " | ";
                std::cout << std::setw(15) << std::left << rnorm / bnorm;
            }
            
            apply_preconditioner(r2,y);
            
            if (verbose)
                std::cout << std::endl;
            
            oldb = beta; //oldb = betak
            beta = std::abs(scalar_product(r2,y));
            
            if(beta < 0)
            {
                istop = 9;
                break;
            }
            
            beta = std::sqrt(beta);
            tnorm2 += alpha*alpha + oldb*oldb + beta*beta;
            
            if (itn == 0)    //Initialize a few things
            {
                if  (beta/beta1 < 10.0 * special::constant::eps)
                    istop = 10;
            }
            
            // Apply previous rotation Q_{k-1} to get
            // [delta_k epsln_{k+1}] = [cs sn]  [dbar_k 0]
            // [gbar_k   dbar_{k+1}]   [sn -cs] [alpha_k beta_{k+1}].
            oldeps = epsln;
            delta  = cs*dbar + sn*alpha;
            gbar   = sn*dbar - cs*alpha;
            epsln  =           sn*beta;
            dbar   =         - cs*beta;
            double root = std::sqrt(gbar * gbar + dbar * dbar);
            Arnorm = phibar * root; // ||Ar_{k-1}||
            
            // Compute next plane rotation Q_k
            gamma = std::sqrt(gbar * gbar + beta * beta); // gamma_k
            gamma = std::max(gamma, special::constant::eps);
            cs = gbar / gamma;                     // c_k
            sn = beta / gamma;                     // s_k
            phi = cs * phibar;                     // phi_k
            phibar = sn * phibar;                  // phibar_{k+1}
            
            // Update x
            double denom = 1. / gamma;
            axby(0., w1, 1., w2);
            axby(0., w2, 1, w);
            axby(0., w, -oldeps, w1);
            axby(1., w, -delta, w2);
            axby(1., w, denom, v);
            axby(1., x, phi, w);

            // go round again
            gmax    = std::max(gmax, gamma);
            gmin    = std::min(gmin, gamma);
            z       = rhs1/gamma;
            rhs1    = rhs2 - delta*z;
            rhs2    =      - epsln*z;
            
            // Estimate various norms
            
            Anorm = std::sqrt(tnorm2);
            ynorm2 = std::abs(scalar_product(x,x));
            ynorm = std::sqrt(ynorm2);
            double epsa = Anorm * special::constant::eps;
            double epsx = epsa * ynorm;
            double epsr = Anorm * ynorm * tol;
            double diag = gbar;
            if (0 == diag)
                diag = epsa;
            
            qrnorm = phibar;
            rnorm  = qrnorm;
            double test1 = 0., test2 = 0.;
            test1  = rnorm / (Anorm * ynorm); // ||r||/(||A|| ||x||)
            test2  = root / Anorm;         // ||A r_{k-1}|| / (||A|| ||r_{k-1}||)
            
            // Estimate cond(A)
            /*
             In this version we look at the diagonals of  R  in the
             factorization of the lower Hessenberg matrix,  Q * H = R,
             where H is the tridiagonal matrix from Lanczos with one
             extra row, beta(k+1) e_k^T.
             */
            Acond = gmax / gmin;
            
            //See if any of the stopping criteria is satisfied
            if (0 == istop)
            {
                double t1 = 1. + test1, t2 = 1. + test2; //This test work if tol < eps
                if (t2 <= 1.) istop = 2;
                if (t1 <= 1.) istop = 1;
                
                if (itn >= max_iterations - 1) istop = 6;
                if (Acond >= 0.1 / special::constant::eps) istop = 4;
                if (epsx >= beta1)   istop = 3;
                if (test2 <= tol)    istop = 2;
                if (test1 <= tol)    istop = 1;
            }
            
            if (0 != istop)
                break;
            
        }
    }
    
    // Display final status
    {
        std::cout << std::setfill('-') << std::setw(80) << "-" << "\n";
        std::cout << msg[istop] << "\n";
        std::cout << " Number of iterations: " << itn << "\n";
        std::cout << " Anorm = " << Anorm << "\t Acond = " << Acond << "\n";
        std::cout << " rnorm = " << rnorm << "\t ynorm = " << ynorm << "\n";
        std::cout << " Arnorm = " << Arnorm << "\n";
        std::cout << std::setfill('-') << std::setw(80) << "-" << std::endl;
        std::cout << std::setfill(' ');
    }

    return itn;
}

#endif
