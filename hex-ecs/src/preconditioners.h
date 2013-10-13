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
        virtual void rhs (cArrayView chi, int ienergy, int instate) const = 0;
        
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
        virtual void rhs (cArrayView chi, int ienergy, int instate) const;
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
        
        virtual void setup () { NoPreconditioner::setup(); }
        virtual void update (double E);
        virtual void rhs (cArrayView chi, int ienergy, int instate) const { NoPreconditioner::rhs(chi, ienergy, instate); }
        virtual void multiply (const cArrayView p, cArrayView q) const { NoPreconditioner::multiply(p, q); }
        virtual void precondition (const cArrayView r, cArrayView z) const;
        
    protected:
        
        // inverse diagonal
        cArray invd_;
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
        virtual void rhs (cArrayView chi, int ienergy, int instate) const;
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
class MultiLevelPreconditioner : public NoPreconditioner
{
    public:
        
        // constructor
        MultiLevelPreconditioner (
            Parallel const & par,
            InputFile const & inp,
            std::vector<std::pair<int,int>> const & ll,
            Bspline const & s_bspline
        ) : NoPreconditioner(par,inp,ll,s_bspline),
            p_bspline_(1, sorted_unique(s_bspline.rknots()), s_bspline.ECStheta(), sorted_unique(s_bspline.cknots())),
            p_rad_(p_bspline_)
        {}
        
        virtual RadialIntegrals const & rad () const { return NoPreconditioner::rad(); }
        
        void setup();
        void update (double E);
        void rhs (cArrayView chi, int ienergy, int instate) const;
        void multiply (const cArrayView p, cArrayView q) const;
        void precondition (const cArrayView r, cArrayView z) const;
        
    protected:
        
        // compute transfer matrix
        void computeSigma_();
        
        // transfer overlap matrix (s|p) multiplied from left by S⁻¹ (in s-basis)
        RowMatrix spSigma_;
        
        // transfer overlap matrix (p|s) multiplied from left by S⁻¹ (in p-basis)
        RowMatrix psSigma_;
        
        // solution objects
        
            // - solution overlap matrix in CSR form (needed for LU)
            CsrMatrix s_csrS_;
            
            // - LU decomposition of the solution overlap matrix
            CsrMatrix::LUft s_luS_;
        
        // preconditioner objects
            
            // - B-spline environment for the preconditioning
            Bspline p_bspline_;
            
            // - radial integrals for the preconditioning
            RadialIntegrals p_rad_;
            
            // - solution overlap matrix in CSR form (needed for LU)
            CsrMatrix p_csrS_;
            
            // - LU decomposition of the solution overlap matrix
            CsrMatrix::LUft p_luS_;
            
            // diagonal blocks in the preconditioner basis
            std::vector<CsrMatrix> p_csr_;
            
            // factorization of p_csr
            std::vector<CsrMatrix::LUft> p_lu_;
            
            // - Kronecker products
            SymDiaMatrix p_half_D_minus_Mm1_tr_kron_S_, p_S_kron_half_D_minus_Mm1_tr_,
                         p_Mm2_kron_S_, p_S_kron_Mm2_;
};

#endif
