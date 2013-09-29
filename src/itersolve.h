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
#include "spmatrix.h"

typedef enum {
    no_prec,
    jacobi_prec,
    ssor_prec,
    dic_prec, // diagonal incomplete Choleski
    ilu_prec , // droptol-incomplete LU
    silu_prec,
    bilu_prec
} Preconditioner;

/**
 * @brief Sparse incomplete Cholesky decomposition.
 * 
 * This routine computes the LDL-decomposition of a symmetric matrix,
 * @f[
 *     A = L D L^T \ ,
 * @f]
 * where @f$ L @f$ is a lower triangular matrix normalized so that it has
 * units on the diagonal and @f$ D @f$ is a diagonal matrix.
 * 
 * @param A Matrix elements in the form of a consecutive array @f$ \left\{a_i\right\}_{i=1}^N @f$
 *          as in
 *          @f[
 *                   \pmatrix
 *                   {
 *                       a_1 &     &     &        &       \cr
 *                       a_2 & a_3 &     &        &       \cr
 *                       a_4 & a_5 & a_6 &        &       \cr
 *                       a_7 & a_8 & a_9 & a_{10} &       \cr
 *                       \vdots &   &     &        & \ddots \cr
 *                   }
 *          @f]
 *          Whenever @f$ a_k @f$ is equal to zero, it is (or can be) omitted from the input array.
 * @param I %Array of column indices (one for every element of A).
 * @param P %Array of row pointers, i.e. starting positions of rows of A. For dense matrix
 *          it would be 0, 1, 3, 6, 10, ... The last element must be equal to the length
 *          of both A and I.
 * 
 * @return The elements of @f$ L @f$ (below diagonal) and @f$ D @f$ (at diagonal) with the
 *         exact sparse pattern as the input array A, i.e. specified by I and P arrays.
 */
cArray iChol(cArrayView const & A, lArrayView const & I, lArrayView const & P);

/**
 * @brief DIC preconditioner.
 * 
 * Setup the diagonal incomplete Cholesky preconditioner. It is a essentially the original
 * matrix with a preconditioned diagonal and the strict upper and lower triangles normalized
 * by the preconditioned diagonal, i.e. a matrix
 * @f[
 *     \matfbf{P} = \mathbf{\tilde{L}}_\mathbf{A} + \mathbf{D}^{-1} + \mathbf{\tilde{L}}_\mathbf{A}
 * @f]
 * for the preconditioner
 * @f[
 *     \mathbf{M} = (\mathbf{D} + \mathbf{L}_\mathbf{A})
 *                  \mathbf{D}^{-1}
 *                  (\mathbf{D} + \mathbf{U}_\mathbf{A})
 *                =  (1 + \mathbf{\tilde{L}}_\mathbf{A})
 *                  \mathbf{D}
 *                  (1 + \mathbf{\tilde{U}}_\mathbf{A})
 * @f]
 * The formula for the elements of @f$ \mathbf{D} @f$ is
 * @f[
 *     d_i = a_{ii} - \sum_{k < i} a_{ik} d_{kk}^{-1} a_{ki} \ ,
 * @f]
 * and is to be evaluated along the diagonal, re-using the just computed values @f$ d_i @f$.
 * Hence, the access pattern in dense matrix would be
 * @f[
 *     \pmatrix {
 *        \ast &      &      & \ast &     &     \cr
 *             & \ast &      & \ast &     &     \cr
 *             &      & \ast & \ast &     &     \cr
 *        \ast & \ast & \ast &    ? &     &     \cr
 *             &      &      &      &     &     \cr
 *             &      &      &      &     &     \cr
 *     }
 * @f]
 * In the case of the sparse SymDiaMatrix, the asterisks will occur only on the nonzero
 * diagonals.
 * 
 * @note The same preconditioner can be used for unsymmetric matrix. Then it is called
 *       DILU (diagonal incomplete LU factorization). This function, though, is implemented
 *       symmetrically (as the input is symmetrical by definition of the type).
 * 
 * @param A Matrix in SymDiaMatrix format that is to be preconditioned.
 * @return The DIC preconditioner of the symmetric matrix.
 */
SymDiaMatrix DIC_preconditioner(SymDiaMatrix const & A);

/**
 * @brief SSOR preconditioner.
 * 
 * Symmetric successive over-relaxation preconditioner for @f$ \omega = 1 @f$.
 * (Essentially symmetrized Gauss-Seidel). The resulting matrix contains
 * normalized lower (and upper) triangle and in the place of the unit diagonal
 * is the inverse diagonal of @f$ \mathbf{A} @f$. So, having the preconditioner
 */
SymDiaMatrix SSOR_preconditioner(SymDiaMatrix const & A);

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
    cCArrayView bview, cArrayView xview,
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
        
        // move to the next iteration: store previous projection
        rho_old = rho_new;
    }
    
    xview = x;
    
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
    cArrayView bview, cArrayView xview,
    double eps,
    int min_iterations, int max_iterations,
    TFunctor1 apply_preconditioner,
    TFunctor2 matrix_multiply,
    bool verbose = false
) {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    std::chrono::duration<int> sec;
    
    cArray b(bview), x(xview);
    
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
    
    xview = x;
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
    cArrayView bview, cArrayView xview,
    double eps,
    int min_iterations, int max_iterations,
    TFunctor1 apply_preconditioner,
    TFunctor2 matrix_multiply,
    bool verbose = false
) {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    std::chrono::duration<int> sec;
    
    cArray x(xview), b(bview);
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
    
    xview = x;
    return i;
}

#endif
