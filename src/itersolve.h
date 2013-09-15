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

/**
 * Callback-based Conjugate gradients.
 * @param bview Right hand side.
 * @param xview Output vector. An initial guess may be present at beginning.
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
template <typename Preconditioner, typename Multiplier>
unsigned cg_callbacks (
    cArrayView const bview, cArrayView xview,
    double eps,
    unsigned min_iterations, unsigned max_iterations,
    Preconditioner apply_preconditioner,
    Multiplier matrix_multiply,
    bool verbose = true
) {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    std::chrono::duration<int> sec;
    
    cArray b(bview), x(xview);
    
    // compute norm of the right hand side
    double bnorm = b.norm();
    
    // some arrays (search directions etc.)
    size_t N = b.size();
    cArray p(N), q(N), z(N);
    
    // residual; initialized to starting residual using the initial guess
    cArray r(N);
    matrix_multiply(x, r);
    r = b - r;
    
    // if the (non-zero) initial guess seems horribly wrong,
    //    use rather the right hand side as the initial guess
    if (r.norm() / bnorm > 1000)
    {
        x.clear();
        r = b;
    }
    
    // some other scalar variables
    Complex rho_new;        // contains inner product r_i^T · r_i
    Complex rho_old;        // contains inner product r_{i-1}^T · r_{i-1}
    Complex alpha, beta;    // contains projection ratios
    
    // Iterate
    
    unsigned k;
    for/*ever*/ (k = 0; ; k++)
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
        rho_new = (r|z);
        
        // setup search direction p
        if (k == 0)
        {
            p = z;
        }
        else
        {
            beta = rho_new / rho_old;
            p = z + beta * p;
        }
        
        // move to next Krylov subspace by multiplying A·p
        matrix_multiply(p, q);
        
        // compute projection ratio α
        alpha = rho_new / (p|q);
        
        // update the solution and the residual
        x += alpha * p;
        r -= alpha * q;
        
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
        
        // check iteration limit (stop at "max_iterations" iterations)
        if (k >= max_iterations)
            break;
        
        // move to the next iteration: store previous projection
        rho_old = rho_new;
    }
    
    xview = x;
    
    return k;
}

/**
 * Callback-based BiCG-STAB.
 * @param b Right hand side.
 * @param x Output vector. An initial guess may be present at beginning.
 * @param eps Relative tolerance.
 * @param min_iterations Minimal iterations count.
 * @param max_iterations Maximal iterations count.
 * @param apply_preconditioner Functor compatible with void(*)(const Array&, Array&) prototype.
 *                             Should apply a custom preconditioner to a given vector.
 * @param matrix_multiply Functor compatible with void(*)(const Array&, Array&) prototype.
 *                        Should multiply the given vector by the matrix of the equation set.
 * @return Iteration count.
 */
template <typename TFunctor1, typename TFunctor2>
unsigned bicgstab_callbacks(
    cArray const & b, cArray & x,
    double eps,
    unsigned min_iterations, unsigned max_iterations,
    TFunctor1 apply_preconditioner,
    TFunctor2 matrix_multiply
) {
    int N = b.size();
    
    cArray r_im1(N), rt(N), p_i(N), p_im1(N), v_i(N), v_im1(N), phat(N), s(N), shat(N), t(N), x_im1(N), x_i(N), r_i(N);
    Complex rho_im1, rho_im2, alpha_i, alpha_im1, omega_i, omega_im1;
    
    x_im1 = x;
    matrix_multiply(x_im1,r_im1);
    rt = r_im1 = b - r_im1;
    
    for/*ever*/ (int i = 1; ; i++)
    {
        rho_im1 = (rt | r_im1);
        if (rho_im1 == 0.)
        {
            std::cerr << "BiCG failed, ρ = 0\n";
            return i-1;
        }
        
        if (i == 1)
        {
            p_i = r_im1;
        }
        else
        {
            Complex beta_im1 = (rho_im1 / rho_im2) * (alpha_im1 / omega_im1);
            p_i = r_im1 + beta_im1 * (p_im1 - omega_im1 * v_im1);
        }
        
        apply_preconditioner(p_i, phat);
        matrix_multiply(phat, v_i);
        alpha_i = rho_im1 / (rt|v_i);
        s = r_im1 - alpha_i * v_i;
        
        std::cout << "\t[bcg] s-residual relative magnitude after " << i << " iterations: " << s.norm() / b.norm() << "\n";
        
        if (s.norm() < eps * b.norm())
        {
            x = x_im1 + alpha_i * phat;
            return i-1;
        }
        
        apply_preconditioner(s, shat);
        matrix_multiply(shat, t);
        omega_i = (t|s) / (t|t);
        x_i = x_im1 + alpha_i * phat + omega_i * shat;
        r_i = s - omega_i * t;
        
        std::cout << "\t[bcg] r-residual relative magnitude after " << i << " iterations: " << r_i.norm() / b.norm() << "\n";
        
        if (r_i.norm() < eps * b.norm())
        {
            x = x_i;
            return i;
        }
        
        if (omega_i == 0.)
        {
            std::cerr << "BiCG failed, ω = 0\n";
            return i;
        }
        
        // update
        x_im1 = x_i;
        r_im1 = r_i;
        rho_im2 = rho_im1;
        p_im1 = p_i;
        v_im1 = v_i;
        alpha_im1 = alpha_i;
        omega_im1 = omega_i;
    }
}

#endif
