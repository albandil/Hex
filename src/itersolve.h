/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2013                              *
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

inline void complex_axby (Complex a, const cArrayView x, Complex b, const cArrayView y, cArrayView z)
{
    size_t N = z.size();
    assert (N == x.size());
    assert (N == y.size());
    
    // accelerators
    Complex const * const restrict px = x.data();
    Complex const * const restrict py = y.data();
    Complex       * const restrict pz = z.data();
    
    // do the axby per element
    for (size_t i = 0; i < N; i++)
        pz[i] = a * px[i] + b * py[i];
}

/**
 * Callback-based Conjugate gradients.
 * @param b Right hand side.
 * @param x Output vector. An initial guess may be present at beginning.
 * @param eps Relative tolerance.
 * @param min_iterations Minimal iterations count.
 * @param max_iterations Maximal iterations count.
 * @param apply_preconditioner Functor compatible with void(*)(const Array&, Array&) prototype.
 *                             Should apply a custom preconditioner to a given vector.
 * @param matrix_multiply Functor compatible with void(*)(const Array&, Array&) prototype.
 *                        Should multiply the given vector by the matrix of the equation set.
 * @param verbose Whether to display progress information.
 * @return Iteration count.
 */
template <
    class TArray,
    class TArrayView,
    class Preconditioner,
    class MatrixMultiplication,
    class AxbyOperation = decltype(complex_axby),
    class ScalarProduct = decltype(operator|<Complex>)
> unsigned cg_callbacks (
    const TArrayView b,
    TArrayView x,
    double eps,
    unsigned min_iterations, unsigned max_iterations,
    Preconditioner apply_preconditioner,
    MatrixMultiplication matrix_multiply,
    bool verbose = true,
    AxbyOperation axby = complex_axby,
    ScalarProduct scalar_product = operator|<Complex>
) {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    std::chrono::duration<int> sec;
    
    // compute norm of the right hand side
    double bnorm = b.norm();
    
    // some arrays (search directions etc.)
    size_t N = b.size();
    TArray p(N), q(N), z(N);
    
    // residual; initialized to starting residual using the initial guess
    cArray r(N);
    matrix_multiply(x, r);
    axby (1., b, -1., r, r);
    
    // if the (non-zero) initial guess seems horribly wrong,
    //    use rather the right hand side as the initial guess
    if (r.norm() / bnorm > 1000)
    {
        x.fill(0.);
        r = b;
    }
    
    // some other scalar variables
    Complex rho_new;        // contains inner product r_i^T · r_i
    Complex rho_old;        // contains inner product r_{i-1}^T · r_{i-1}
    Complex alpha, beta;    // contains projection ratios
    
    // Iterate
    
    unsigned k;
    for (k = 0; k < max_iterations; k++)
    {
        sec = std::chrono::duration_cast<std::chrono::duration<int>>(std::chrono::steady_clock::now()-start);
        
        if (verbose)
        {
            std::cout << "\t[cg] Residual relative magnitude after "
                    << k << " iterations: " << r.norm() / bnorm
                    << " (" << sec.count()/60 << " min)\n";
        }
        
        // apply desired preconditioner
        apply_preconditioner(r, z);
        
        // compute projection ρ = r·z
        rho_new = scalar_product(r, z);
        
        // setup search direction p
        if (k == 0)
        {
            axby (1., z, 0., z, p);
        }
        else
        {
            beta = rho_new / rho_old;
            axby (1., z, beta, p, p);
        }
        
        // move to next Krylov subspace by multiplying A·p
        matrix_multiply(p, q);
        
        // compute projection ratio α
        alpha = rho_new / scalar_product(p, q);
        
        // update the solution and the residual
        axby (1., x, alpha, p, x);
        axby (1., r, -alpha, q, r);
        
        // compute and check norm
        double rnorm = r.norm();
        if (not finite(rnorm))
        {
            std::cout << "\t[cg] Oh my god... the norm of the solution is not finite. Something went wrong!\n";
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
