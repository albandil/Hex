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

#ifndef HEX_PRECONDITIONERS
#define HEX_PRECONDITIONERS

#include "arrays.h"
#include "input.h"
#include "matrix.h"
#include "radial.h"
#include "parallel.h"

typedef enum {
// default:
    ilu_prec , // droptol-incomplete LU
// other
    no_prec,
    jacobi_prec,
    ssor_prec,
    dic_prec, // diagonal incomplete Choleski
    silu_prec,
    bilu_prec,
    res_prec, // multi-resolution preconditioner
    spai_prec // sparse approximate inverse
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
SymDiaMatrix DIC(SymDiaMatrix const & A);

/**
 * @brief SSOR preconditioner.
 * 
 * Symmetric successive over-relaxation preconditioner for @f$ \omega = 1 @f$.
 * (Essentially symmetrized Gauss-Seidel). The resulting matrix contains
 * normalized lower (and upper) triangle and in the place of the unit diagonal
 * is the inverse diagonal of @f$ \mathbf{A} @f$. So, having the preconditioner
 */
SymDiaMatrix SSOR(SymDiaMatrix const & A);

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
         * @brief Get radial integrals.
         */
        virtual RadialIntegrals const & rad() const = 0;
        
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
         * @brief Return the right-hand side.
         */
        virtual void rhs (const cArrayView chi, int ienergy, int instate) const = 0;
        
        /**
         * @brief Multiply by the matrix equation.
         * 
         * This function implements matrix multiplication by the matrix of
         * the set of equations that is to be solved.
         */
        virtual void multiply (const cArrayView p, cArrayView q) const = 0;
        
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
        virtual void precondition (const cArrayView r, cArrayView z) const = 0;
};

class NoPreconditioner : public PreconditionerBase
{
    public:
        
        NoPreconditioner (
            Parallel const & par,
            InputFile const & inp,
            std::vector<std::pair<int,int>> const & ll,
            Bspline const & bspline
        ) : PreconditionerBase(), par_(par), inp_(inp), l1_l2_(ll),
            s_bspline_(bspline), s_rad_(bspline) {}
        
        virtual RadialIntegrals const & rad () const { return s_rad_; }
        
        virtual void setup ();
        virtual void update (double E);
        virtual void rhs (const cArrayView chi, int ienergy, int instate) const;
        virtual void multiply (const cArrayView p, cArrayView q) const;
        virtual void precondition (const cArrayView r, cArrayView z) const;
        
    protected:
        
        // parallel environment
        Parallel const & par_;
        
        // input parameters
        InputFile const & inp_;
        
        // coupled states
        std::vector<std::pair<int,int>> const & l1_l2_;
        
        // diagonal blocks in DIA format (these will be used in matrix multiplication)
        std::vector<SymDiaMatrix> dia_blocks_;
        
        // B-spline environment for the solution
        Bspline const & s_bspline_;
            
        // radial integrals for the solution
        RadialIntegrals s_rad_;
            
        // Kronecker products
        SymDiaMatrix S_kron_S_, S_kron_Mm1_tr_, S_kron_Mm2_, Mm1_tr_kron_S_,
            Mm2_kron_S_, half_D_minus_Mm1_tr_, half_D_minus_Mm1_tr_kron_S_,
            S_kron_half_D_minus_Mm1_tr_;
};

class JacobiPreconditioner : public NoPreconditioner
{
    public:
        
        JacobiPreconditioner (
            Parallel const & par,
            InputFile const & inp,
            std::vector<std::pair<int,int>> const & ll,
            Bspline const & bspline
        ) : NoPreconditioner(par, inp, ll, bspline) {}
        
        virtual RadialIntegrals const & rad () const { return NoPreconditioner::rad(); }
        
        virtual void setup ();
        virtual void update (double E);
        virtual void rhs (const cArrayView chi, int ienergy, int instate) const { NoPreconditioner::rhs(chi, ienergy, instate); }
        virtual void multiply (const cArrayView p, cArrayView q) const { NoPreconditioner::multiply(p, q); }
        virtual void precondition (const cArrayView r, cArrayView z) const;
        
    protected:
        
        // inverse diagonals for every block
        cArrays invd_;
};

class SSORPreconditioner : public NoPreconditioner
{
    public:
        
        SSORPreconditioner (
            Parallel const & par,
            InputFile const & inp,
            std::vector<std::pair<int,int>> const & ll,
            Bspline const & bspline
        ) : NoPreconditioner(par, inp, ll, bspline) {}
        
        virtual RadialIntegrals const & rad () const { return NoPreconditioner::rad(); }
        
        virtual void setup ();
        virtual void update (double E);
        virtual void rhs (const cArrayView chi, int ienergy, int instate) const { NoPreconditioner::rhs(chi, ienergy, instate); }
        virtual void multiply (const cArrayView p, cArrayView q) const { NoPreconditioner::multiply(p, q); }
        virtual void precondition (const cArrayView r, cArrayView z) const;
        
    protected:
        
        // inverse diagonals for every block
        std::vector<SymDiaMatrix> SSOR_;
};

class ILUPreconditioner : public NoPreconditioner
{
    public:
        
        ILUPreconditioner (
            Parallel const & par,
            InputFile const & inp,
            std::vector<std::pair<int,int>> const & ll,
            Bspline const & bspline
        ) : NoPreconditioner(par, inp, ll, bspline) {}
        
        virtual RadialIntegrals const & rad () const { return NoPreconditioner::rad(); }
        
        virtual void setup ();
        virtual void update (double E);
        virtual void rhs (const cArrayView chi, int ienergy, int instate) const;
        virtual void multiply (const cArrayView p, cArrayView q) const;
        virtual void precondition (const cArrayView r, cArrayView z) const;
        
    protected:
        
        // diagonal CSR block for every coupled state
        std::vector<CsrMatrix> csr_blocks_;
        
        // LU decompositions of the CSR blocks
        std::vector<CsrMatrix::LUft> lu_;
};

/**
 * @brief Multi-resolution preconditioner.
 * 
 * This preconditioner class is inspired by the method of multi-resolution in the finite
 * difference computations. A high-order differential operator that gives rise to the
 * multi-diagonal matrix is initially solved in low order. The result is then used as
 * a preconditioner for the full, high-order method.
 * 
 * In the case of B-spline expansion, that is being used in Hex, the bandwidth of the
 * resulting matrix is given by the order of the B-splines. Analogously to the just
 * metioned method, we compute the solution in a sharper base (lower order), which is
 * computationally less demanding. The algorithm works as follows:
 * 
 * - In every step we are solving the equation @f$ \mathbf{M}^{(s)}\mathbf{z}^{(s)} = \mathbf{r}^{(s)} @f$.
 * - The matrix @f$ \mathbf{M}^{(s)} @f$ is the diagonal block
 *   @f$ E\mathbf{S}^{(s)}\otimes\mathbf{S}^{(s)} - \mathbf{H}^{(s)} @f$. The symbol "s"
 *   stands for solution B-spline basis.
 * - But we want just to compute the simplified system @f$ \mathbf{M}^{(p)}\mathbf{z}^{(p)} = \mathbf{r}^{(p)} @f$.
 *   The symbol "p" stands for preconditioner B-spline basis.
 * - The matrices contained in @f$ \mathbf{M}^{(p)} @f$ are computed in a low-order B-spline
 *   basis; these are the matrices @f$ \mathbf{S}^{(p)} @f$, @f$ \mathbf{D}^{(p)} @f$
 *   etc., whereas fhe original full-order matrices are @f$ \mathbf{S}^{(s)} @f$,
 *   @f$ \mathbf{D}^{(s)} @f$ etc.
 * - Before and after the solution of the simplified system we need to map the solutions
 *   between the bases. This is done by the equations
 *   @f[
 *          \mathbf{S}^{(p)}\otimes\mathbf{S}^{(p)} \mathbf{r}^{(p)}
 *          = \mathbb{\Sigma}^{(p,q)} \mathbf{r}^{(q)}
 *   @f]
 * - Together we have the sequence
 *   @f[
 *          \mathbf{S}^{(p)}\otimes\mathbf{S}^{(p)} \mathbf{r}^{(p)}
 *          = \mathbb{\Sigma}^{(p,s)} \mathbf{r}^{(s)}
 *   @f]
 *   @f[
 *          \mathbf{M}^{(p)}\mathbf{z}^{(p)} = \mathbf{r}^{(p)}
 *   @f]
 *   @f[
 *          \mathbf{S}^{(s)}\otimes\mathbf{S}^{(s)} \mathbf{z}^{(s)}
 *          = \mathbb{\Sigma}^{(s,p)} \mathbf{z}^{(p)}
 *   @f]
 *   instead of the single equation
 *   @f[
 *          \mathbf{M}^{(s)}\mathbf{z}^{(s)} = \mathbf{r}^{(s)} \ .
 *   @f]
 * 
 * The transition overlap matrix @f$ \mathbb{\Sigma}^{(p,q)} @f$ is defined in components
 * as
 * @f[
 *        \Sigma_{ij}^{(p,q)} = \int B_i^{(p)}(r) B_j^{(q)}(r) \mathrm{d}r \ .
 * @f]
 * The solution of transition equations can be recast into a matrix form of Kronecker product,
 * @f[
 *        \mathbf{S}^{(p)} \mathbf{R}^{(p)} \mathbf{S}^{(p)T} = \mathrm{Matrix}(\mathbb{\Sigma}^{(p,q)} \mathbf{r}^{(q)})
 *        \qquad\Rightarrow\qquad
 *        \mathbf{R}^{(p)} = \mathbf{S}^{(p)^{-1}} \mathrm{Matrix}(\mathbb{\Sigma}^{(p,q)} \mathbf{r}^{(q)}) \mathbf{S}^{(p)^{-T}}
 * @f]
 * where the word "Matrix" indicates that we want to split long vector into a square matrix. The resulting square matrix
 * @f$ \mathbf{R} @f$ will then have the same structure.
 */
class MultiLevelPreconditioner : public SSORPreconditioner
{
    public:
        
        // constructor
        MultiLevelPreconditioner (
            Parallel const & par,
            InputFile const & inp,
            std::vector<std::pair<int,int>> const & ll,
            Bspline const & s_bspline
        ) : SSORPreconditioner(par,inp,ll,s_bspline),
            p_bspline_(1, sorted_unique(s_bspline.rknots()), s_bspline.ECStheta(), sorted_unique(s_bspline.cknots())),
            p_rad_(p_bspline_), g_(p_bspline_)
        {}
        
        virtual RadialIntegrals const & rad () const { return NoPreconditioner::rad(); }
        
        void setup();
        void update (double E);
        void rhs (const cArrayView chi, int ienergy, int instate) const;
        void multiply (const cArrayView p, cArrayView q) const;
        void precondition (const cArrayView r, cArrayView z) const;
        
    protected:
        
        // compute transfer matrix
        void computeSigma_();
        
        // compute contribution to an element of the transfer matrix
        Complex computeSigma_iknot_ (int qord, int is, int iknots, int ip, int iknotp) const;
        
        // transfer overlap matrix (s|p) multiplied from left by S⁻¹ (in s-basis)
        RowMatrix<Complex> spSigma_;
        
        // transfer overlap matrix (p|s) multiplied from left by S⁻¹ (in p-basis)
        RowMatrix<Complex> psSigma_;
        
        /// DEBUG
        ColMatrix<Complex> spSigma, SspSigma;
        ColMatrix<Complex> psSigma, SpsSigma;
        
        // preconditioner objects
            
            // - B-spline environment for the preconditioning
            Bspline p_bspline_;
            
            // - radial integrals for the preconditioning
            RadialIntegrals p_rad_;
            
            // diagonal blocks in the preconditioner basis
            std::vector<CsrMatrix> p_csr_;
            
            // factorization of p_csr
            std::vector<CsrMatrix::LUft> p_lu_;
            
            // - Kronecker products
            SymDiaMatrix p_half_D_minus_Mm1_tr_kron_S_, p_S_kron_half_D_minus_Mm1_tr_,
                         p_Mm2_kron_S_, p_S_kron_Mm2_, p_S_kron_S_;
        
        // integrator
        GaussLegendre g_;
};



#endif
