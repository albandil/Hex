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

#include "angular.h"
#include "arrays.h"
#include "io.h"
#include "matrix.h"
#include "radial.h"
#include "parallel.h"

#ifndef NO_OPENCL
    #include "opencl.h"
#endif

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
cArray iChol (cArrayView const & A, lArrayView const & I, lArrayView const & P);

/**
 * @brief DIC preconditioner.
 * 
 * Setup the diagonal incomplete Cholesky preconditioner. It is a essentially the original
 * matrix with a preconditioned diagonal and the strict upper and lower triangles normalized
 * by the preconditioned diagonal, i.e. a matrix
 * @f[
 *     \mathbf{P} = \mathbf{\tilde{L}}_\mathbf{A} + \mathbf{D}^{-1} + \mathbf{\tilde{L}}_\mathbf{A}^T
 * @f]
 * for the preconditioner
 * @f[
 *     \mathbf{M} = (\mathbf{D} + \mathbf{L}_\mathbf{A})
 *                  \mathbf{D}^{-1}
 *                  (\mathbf{D} + \mathbf{U}_\mathbf{A})
 *                =  (1 + \mathbf{\tilde{U}}_\mathbf{A}^T)
 *                  \mathbf{D}
 *                  (1 + \mathbf{\tilde{U}}_\mathbf{A})
 * @f]
 * The formula for the elements of @f$ \mathbf{D} @f$ is
 * @f[
 *     d_i = a_{ii} - \sum_{k < i} a_{ik} d_{k}^{-1} a_{ki} \ ,
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
SymDiaMatrix DIC (SymDiaMatrix const & A);

/**
 * @brief SSOR preconditioner.
 * 
 * Symmetric successive over-relaxation preconditioner for @f$ \omega = 1 @f$.
 * (Essentially symmetrized Gauss-Seidel). The resulting matrix contains
 * normalized lower (and upper) triangle and in the place of the unit diagonal
 * is the inverse diagonal of @f$ \mathbf{A} @f$. So, having the preconditioner
 */
SymDiaMatrix SSOR (SymDiaMatrix const & A);

/**
 * @brief SPAI preconditioner.
 * 
 * Compute sparse aproximate inverse of a given symmetrix diagonal matrix A. The sparse
 * structure of the SPAI is set by the second parameter that contains list of non-lower
 * diagonal indices (greater than or equal to zero).
 * 
 * This function uses Lapack routine ZGELSD.
 */
CooMatrix SPAI (SymDiaMatrix const & A, const iArrayView diagonals);

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

/**
 * @brief Solution driver without actual preconditioner.
 * 
 * This class "preconditions" by identity matrix, but implements all other important
 * routines, that can be used by derived classes, namely:
 * - setup : Loads / computed radial integrals for construction of the matrix of the set
 *           and for the construction of the right-hand side.
 * - update : Creates the diagonal blocks.
 * - rhs : Composes the right-hand side.
 * - multiply : Multiplies a vector by the matrix of the set of equations.
 */
class NoPreconditioner : public PreconditionerBase
{
    public:
        
        static const std::string name;
        static const std::string description;
        
        NoPreconditioner
        (
            Parallel const & par,
            InputFile const & inp,
            AngularBasis const & ang,
            Bspline const & bspline,
            CommandLine const & cmd
        ) : PreconditionerBase(), cmd_(cmd), par_(par), inp_(inp), ang_(ang),
            s_bspline_(bspline), s_rad_(s_bspline_)
        {
            // nothing to do
        }
        
        virtual RadialIntegrals const & rad () const { return s_rad_; }
        
        virtual void setup ();
        virtual void update (double E);
        virtual void rhs (cArrayView chi, int ienergy, int instate) const;
        virtual void multiply (const cArrayView p, cArrayView q) const;
        virtual void precondition (const cArrayView r, cArrayView z) const;
        
    protected:
        
        // energy
        double E_;
        
        // command line switches
        CommandLine const & cmd_;
        
        // parallel environment
        Parallel const & par_;
        
        // input parameters
        InputFile const & inp_;
        
        // coupled states
        AngularBasis const & ang_;
        
        // diagonal blocks in DIA format (these will be used in matrix multiplication)
        mutable std::vector<SymDiaMatrix> dia_blocks_;
        
        // B-spline environment for the solution
        Bspline s_bspline_;
            
        // radial integrals for the solution
        RadialIntegrals s_rad_;
            
        // Kronecker products
        SymDiaMatrix S_kron_S_, S_kron_Mm1_tr_, S_kron_Mm2_, Mm1_tr_kron_S_,
            Mm2_kron_S_, half_D_minus_Mm1_tr_, half_D_minus_Mm1_tr_kron_S_,
            S_kron_half_D_minus_Mm1_tr_, S_kron_Mm3_tr_, Mm3_tr_kron_S_;
};

/**
 * @brief CG iteration-based preconditioner.
 * 
 * This class adds some preconditioning capabilities to its base class
 * NoPreconditioner. The preconditioning is done by diagonal block solution
 * using the conjugate gradients solver (which itself is non-preconditioned).
 */
class CGPreconditioner : public NoPreconditioner
{
    public:
        
        static const std::string name;
        static const std::string description;
        
        CGPreconditioner
        (
            Parallel const & par,
            InputFile const & inp,
            AngularBasis const & ll,
            Bspline const & bspline,
            CommandLine const & cmd
        ) : NoPreconditioner(par, inp, ll, bspline, cmd) {}
        
        // reuse parent definitions
        virtual RadialIntegrals const & rad () const { return NoPreconditioner::rad(); }
        virtual void setup () { return NoPreconditioner::setup(); }
        virtual void update (double E) { return NoPreconditioner::update(E); }
        virtual void rhs (cArrayView chi, int ienergy, int instate) const { NoPreconditioner::rhs(chi, ienergy, instate); }
        virtual void multiply (const cArrayView p, cArrayView q) const { NoPreconditioner::multiply(p, q); }
        
        // declare own definitions
        virtual void precondition (const cArrayView r, cArrayView z) const;
        
        // inner CG callbacks
        virtual void CG_mmul (int iblock, const cArrayView p, cArrayView q) const;
        virtual void CG_prec (int iblock, const cArrayView r, cArrayView z) const;
};

/**
 * @brief Jacobi-preconditioned CG-based preconditioner.
 * 
 * Enhances CGPreconditioner conjugate gradients solver by Jacobi (diagonal) preconditioning.
 * This is done by redefining virtual function CG_prec.
 */
class JacobiCGPreconditioner : public CGPreconditioner
{
    public:
        
        static const std::string name;
        static const std::string description;
        
        JacobiCGPreconditioner
        (
            Parallel const & par,
            InputFile const & inp,
            AngularBasis const & ll,
            Bspline const & bspline,
            CommandLine const & cmd
        ) : CGPreconditioner(par, inp, ll, bspline, cmd) {}
        
        // reuse parent definitions
        virtual RadialIntegrals const & rad () const { return CGPreconditioner::rad(); }
        virtual void rhs (cArrayView chi, int ienergy, int instate) const { CGPreconditioner::rhs(chi, ienergy, instate); }
        virtual void multiply (const cArrayView p, cArrayView q) const { CGPreconditioner::multiply(p, q); }
        virtual void precondition (const cArrayView r, cArrayView z) const { CGPreconditioner::precondition(r, z); }
        
        // declare own definitions
        virtual void setup ();
        virtual void update (double E);
        
        // inner CG callback (needed by parent)
        virtual void CG_prec (int iblock, const cArrayView r, cArrayView z) const { z = invd_[iblock] * r; }
        
    protected:
        
        // inverse diagonals for every block
        cArrays invd_;
};

/**
 * @brief SSOR-preconditioned CG-based preconditioner.
 * 
 * Enhances CGPreconditioner conjugate gradients solver by SSOR preconditioning.
 * This is done by redefining virtual function CG_prec.
 */
class SSORCGPreconditioner : public CGPreconditioner
{
    public:
        
        static const std::string name;
        static const std::string description;
        
        SSORCGPreconditioner
        (
            Parallel const & par,
            InputFile const & inp,
            AngularBasis const & ll,
            Bspline const & bspline,
            CommandLine const & cmd
        ) : CGPreconditioner(par, inp, ll, bspline, cmd) {}
        
        // reuse parent definitions
        virtual RadialIntegrals const & rad () const { return CGPreconditioner::rad(); }
        virtual void rhs (cArrayView chi, int ienergy, int instate) const { CGPreconditioner::rhs(chi, ienergy, instate); }
        virtual void multiply (const cArrayView p, cArrayView q) const { CGPreconditioner::multiply(p, q); }
        virtual void precondition (const cArrayView r, cArrayView z) const { CGPreconditioner::precondition(r, z); }
        
        // declare own definitions
        virtual void setup ();
        virtual void update (double E);
        
        // inner CG callback (needed by parent)
        virtual void CG_prec (int iblock, const cArrayView r, cArrayView z) const;
        
    protected:
        
        // inverse diagonals for every block
        mutable std::vector<SymDiaMatrix> SSOR_;
};

/**
 * 
 */
class SepCGPreconditioner : public CGPreconditioner
{
    public:
        
        static const std::string name;
        static const std::string description;
        
        SepCGPreconditioner
        (
            Parallel const & par,
            InputFile const & inp,
            AngularBasis const & ll,
            Bspline const & bspline,
            CommandLine const & cmd
        ) : CGPreconditioner(par, inp, ll, bspline, cmd),
            invCl_invsqrtS_(inp.maxell+1), invsqrtS_Cl_(inp.maxell+1), Dl_(inp.maxell+1)
        {
            // nothing more to do
        }
        
        // reuse parent definitions
        virtual RadialIntegrals const & rad () const { return CGPreconditioner::rad(); }
        virtual void multiply (const cArrayView p, cArrayView q) const { CGPreconditioner::multiply(p,q); }
        virtual void rhs (cArrayView chi, int ienergy, int instate) const { CGPreconditioner::rhs(chi,ienergy,instate); }
        virtual void precondition (const cArrayView r, cArrayView z) const { CGPreconditioner::precondition(r,z); }
        virtual void update (double E) { CGPreconditioner::update(E); }
        
        // declare own definitions
        virtual void setup ();
        
        // inner CG callback (needed by parent)
        virtual void CG_prec (int iblock, const cArrayView r, cArrayView z) const;
        
    protected:
        
        // one-electron hamiltonian eigen-diagonalization
        std::vector<RowMatrix<Complex>> invCl_invsqrtS_, invsqrtS_Cl_;
        
        // diagonal parts of one-electron hamiltonians
        std::vector<cArray> Dl_;
};

/**
 * @brief ILU-preconditioned CG-based preconditioner.
 * 
 * Enhances CGPreconditioner conjugate gradients solver by incomplete LU factorization
 * preconditioning. This is done by redefining virtual function CG_prec. The factorization
 * is drop tolerance based and is computed by UMFPACK.
 */
class ILUCGPreconditioner : public CGPreconditioner
{
    public:
        
        static const std::string name;
        static const std::string description;
        
        ILUCGPreconditioner
        (
            Parallel const & par,
            InputFile const & inp,
            AngularBasis const & ll,
            Bspline const & bspline,
            CommandLine const & cmd
        ) : CGPreconditioner(par, inp, ll, bspline, cmd), droptol_(cmd.droptol), csr_blocks_(ll.size()), lu_(ll.size()) {}
        
        // reuse parent definitions
        virtual RadialIntegrals const & rad () const { return CGPreconditioner::rad(); }
        virtual void multiply (const cArrayView p, cArrayView q) const { CGPreconditioner::multiply(p,q); }
        virtual void rhs (cArrayView chi, int ienergy, int instate) const { CGPreconditioner::rhs(chi,ienergy,instate); }
        virtual void precondition (const cArrayView r, cArrayView z) const { CGPreconditioner::precondition(r,z); }
        virtual void setup () { CGPreconditioner::setup(); }
        
        // declare own definitions
        virtual void update (double E);
        
        // inner CG callback (needed by parent)
        virtual void CG_prec (int iblock, const cArrayView r, cArrayView z) const;
        
    protected:
        
        // drop tolarance for the factorizations
        double droptol_;
        
        // diagonal CSR block for every coupled state
        mutable std::vector<CsrMatrix> csr_blocks_;
        
        // LU decompositions of the CSR blocks
        mutable std::vector<CsrMatrix::LUft> lu_;
};

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
            ILUCGPreconditioner,        // Solve diagonal blocks by drop-tolerance incomplete LU factorization.
            NoPreconditioner,           // No preconditioner.
            CGPreconditioner,           // Solve diagonal blocks by non-preconditioned CG iterations.
            JacobiCGPreconditioner,     // Solve diagonal blocks by Jacobi-preconditioned CG iterations.
            SSORCGPreconditioner,       // Solve diagonal blocks by SSOR-preconditioned CG iterations.
            SepCGPreconditioner         // Solve diagonal blocks by separate electrons preconditioned CG iterations.
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
            AngularBasis const & ll,
            Bspline const & bspline,
            CommandLine const & cmd
        )
        {
            // check if i-th preconditioner is the right one
            if (i == cmd.preconditioner)
                // - Yes : return new pointer of its type
                return new typename std::tuple_element<i,AvailableTypes>::type(par, inp, ll, bspline, cmd);
            else
                // - No : try the next preconditioner
                return choose<i+1>(par, inp, ll, bspline, cmd);
        }
        template <int i = 0> static inline typename std::enable_if<(i == std::tuple_size<AvailableTypes>::value), PreconditionerBase*>::type choose
        (
            Parallel const & par,
            InputFile const & inp,
            AngularBasis const & ll,
            Bspline const & bspline,
            CommandLine const & cmd
        ){
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
