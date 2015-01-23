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
        virtual cArray rhs (int ienergy, int instate, int Spin) const = 0;
        
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
            std::vector<std::pair<int,int>> const & ll,
            Bspline const & bspline,
            CommandLine const & cmd
        ) : PreconditionerBase(), cmd_(cmd), par_(par), inp_(inp), l1_l2_(ll),
            dia_blocks_(l1_l2_.size()), s_bspline_(bspline), s_rad_(s_bspline_)
        {
            // nothing to do
        }
        
        virtual RadialIntegrals const & rad () const { return s_rad_; }
        
        virtual void setup ();
        virtual void update (double E);
        virtual cArray rhs (int ienergy, int instate, int Spin) const;
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
        std::vector<std::pair<int,int>> const & l1_l2_;
        
        // diagonal blocks in DIA format (these will be used in matrix multiplication)
        mutable std::vector<BlockSymDiaMatrix> dia_blocks_;
        
        // B-spline environment for the solution
        Bspline s_bspline_;
            
        // radial integrals for the solution
        RadialIntegrals s_rad_;
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
            std::vector<std::pair<int,int>> const & ll,
            Bspline const & bspline,
            CommandLine const & cmd
        ) : NoPreconditioner(par, inp, ll, bspline, cmd) {}
        
        // reuse parent definitions
        virtual RadialIntegrals const & rad () const { return NoPreconditioner::rad(); }
        virtual void setup () { return NoPreconditioner::setup(); }
        virtual void update (double E) { return NoPreconditioner::update(E); }
        virtual cArray rhs (int ienergy, int instate, int Spin) const { return NoPreconditioner::rhs(ienergy, instate, Spin); }
        virtual void multiply (const cArrayView p, cArrayView q) const { NoPreconditioner::multiply(p, q); }
        
        // declare own definitions
        virtual void precondition (const cArrayView r, cArrayView z) const;
        
        // inner CG callbacks
        virtual void CG_mmul (int iblock, const cArrayView p, cArrayView q) const;
        virtual void CG_prec (int iblock, const cArrayView r, cArrayView z) const;
};

/**
 * @brief KPA-preconditioned CG-preconditioner.
 * 
 * This nested preconditioner simplifies the Hamiltonian matrix by omitting electron-electron
 * interaction. The block of the matrix can then be written as a sum of simple Kronecker products
 * @f[
 *     \mathsf{A} = E\mathsf{S}\otimes\mathsf{S}
 *      - \mathsf{H}_1\otimes\mathsf{S}
 *      - \mathsf{S}\otimes\mathsf{H}_2 \,,
 * @f]
 * which can be easily diagonalized (only diagonalization of the small matrices are needed)
 * and so also inverted, which we need for solution of the equations.
 */
class KPACGPreconditioner : public CGPreconditioner
{
    public:
        
        static const std::string name;
        static const std::string description;
        
        KPACGPreconditioner
        (
            Parallel const & par,
            InputFile const & inp,
            std::vector<std::pair<int,int>> const & ll,
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
        virtual cArray rhs (int ienergy, int instate, int Spin) const { return CGPreconditioner::rhs(ienergy,instate,Spin); }
        virtual void precondition (const cArrayView r, cArrayView z) const { CGPreconditioner::precondition(r,z); }
        virtual void update (double E) { CGPreconditioner::update(E); }
        
        // declare own definitions
        virtual void setup ();
        
        // inner CG callback (needed by parent)
        virtual void CG_prec (int iblock, const cArrayView r, cArrayView z) const;
        virtual void CG_mmul (int iblock, const cArrayView r, cArrayView z) const;
        
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
            std::vector<std::pair<int,int>> const & ll,
            Bspline const & bspline,
            CommandLine const & cmd
        ) : CGPreconditioner(par, inp, ll, bspline, cmd), droptol_(cmd.droptol),
            csr_blocks_(ll.size()), lu_(ll.size())
        {
            // nothing more to do
        }
        
        // reuse parent definitions
        virtual RadialIntegrals const & rad () const { return CGPreconditioner::rad(); }
        virtual void multiply (const cArrayView p, cArrayView q) const { CGPreconditioner::multiply(p,q); }
        virtual cArray rhs (int ienergy, int instate, int Spin) const { return CGPreconditioner::rhs(ienergy,instate,Spin); }
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
            ILUCGPreconditioner         // Solve diagonal blocks by drop-tolerance incomplete LU factorization.
            , NoPreconditioner          // No preconditioner.
            , CGPreconditioner          // Solve diagonal blocks by non-preconditioned CG iterations.
#ifndef NO_LAPACK
            , KPACGPreconditioner       // Solve diagonal blocks by separate electrons preconditioned CG iterations.
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
