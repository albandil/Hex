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

#ifndef HEX_CCC_BASIS
#define HEX_CCC_BASIS

#include "arrays.h"
#include "matrix.h"
#include "symbolic.h"

/**
 * @brief Lapack DSYEV.
 * 
 * Compute eigenvalues and optionally also eigenvectors of a symmetric real
 * matrix. The eigenvectors will overwrite the contents of the matrix. It is
 * possible to call the functions with "lwork = -1" to query for the optimal
 * workspace size. That will be stored in the first element of "work".
 * 
 * @param jobz Compute eigenvalues ("N") or also eigenvectors ("V").
 * @param uplo Use upper triangle ("U") or lower triangle ("L").
 * @param n Order of the matrix.
 * @param a Matrix elements.
 * @param lda Leading dimension of the matrix "a".
 * @param w Eigenvalues.
 * @param work Workspace.
 * @param lwork Workspace size.
 * @param info "0" on success.
 */
extern "C" void dsyev_
(
    char* jobz,
    char* uplo,
    int* n,
    double* a,
    int* lda,
    double* w,
    double* work,
    int* lwork,
    int* info
);

/**
 * @brief Hydrogen pseudostate basis.
 * 
 * This objects manages the basis of the hydrogen atom. The primary basis consists
 * of a set of "Laguerre states" @f$ \xi_i^{(\ell)}(r) @f$ that are defined by
 * @f[
 *     \xi_i^{(\ell)}(r) = \sqrt{\frac{\lambda_\ell (k-1)!}{(2\ell + 1 + k)!}}
 *     (\lambda_\ell r)^{\ell+1} \mathrm{e}^{-\lambda_\ell r/2}
 *     L_{k-1}^{2\ell+2}(\lambda_\ell r) \ ,
 * @f]
 * where the screening constants @f$ \lambda_\ell @f$ are defined per angular
 * momentum @f$ \ell @f$ and the indices @f$ i @f$ run from 1 to the size of
 * the basis for a particular angular momentum @f$ N_\ell @f$. These states are
 * orthonormal and of finite range proportional to @f$ \lambda_\ell @f$.
 * 
 * By diagonalization of the hydrogen hamiltonian the class creates a second
 * set of states (eigenstate basis), where every state has a definite energy
 * (eigenvalue of the hamiltonian). This eigenstate basis is used in the computation
 * for expansion of the Green function and also for specification of the initial
 * and final atomic states.
 */
class LaguerreBasis
{
    public:
        
        /**
         * @brief Constructor.
         * 
         * This will create the basis for every angular momentum in the given range.
         * @param maxell Maximal single-electron angular momentum (angular momentum
         *               of the basis).
         * @param Nl For every allowed angular momentum ("maxell+1" numbers total)
         *           the size of the basis @f$ N_\ell @f$.
         * @param lambda For every angular momentum ("maxell+1" numbers total)
         *               the screening constant @f$ \lambda_\ell @f$.
         */
        LaguerreBasis (int maxell, const iArrayView Nl, const ArrayView<symbolic::rational> lambda);
        
        /**
         * @brief Screening constants.
         * Get screening constants.
         */
        //@{
        ArrayView<symbolic::rational> rat_lambdas () const { return rat_lambda_; }
        rArrayView lambdas () const { return lambda_; }
        //@}
        
        /**
         * @brief Screening constant.
         * Get screening constant for a specific angular momentum.
         */
        //@{
        symbolic::rational rat_lambda (int ell) const { return rat_lambda_[ell]; }
        double lambda (int ell) const { return lambda_[ell]; }
        //@}
        
        /**
         * @brief Evaluate Laguerre state.
         * 
         * Returns the value of the Laguerre state,
         * @f[
         *     \xi_i^{(\ell)}(r) = \sqrt{\frac{\lambda_\ell (k-1)!}{(2\ell + 1 + k)!}}
         *     (\lambda_\ell r)^{\ell+1} \mathrm{e}^{-\lambda_\ell r/2}
         *     L_{k-1}^{2\ell+2}(\lambda_\ell r) \ ,
         * @f]
         * where the screening constant @f$ \lambda_\ell @f$ is taken from the stored
         * array of screening constants that have been previously used to initialize
         * the object.
         * @param ell Angular momentum of the state (0, 1, 2, ...).
         * @param i Index if the state (1, 2, 3, ...).
         * @param r Radius (equal to or greater than zero).
         */
        double basestate (int ell, int i, double r) const;
        
        /**
         * @brief Retrieve eigenenergy of an orbital.
         * 
         * The function returns the eigenvalue of Hamiltonian for a chosen
         * eigenstate obtained by diagonalization in the Laguerre basis.
         * @param ell Angular momentum.
         * @param n Index of the eigenstate.
         */
        double energy (int ell, int n) const;
        
        /**
         * @brief Eigenstate matrix.
         * 
         * Return the whole matrix of eigenstates for a particular angular momentum.
         */
        RowMatrix<double> const & matrix (int ell) const;
        
        /**
         * @brief Retrieve expansion of an orbital.
         * 
         * This method returns the real expansion of the orbital specified by
         * the angular momentum of the basis (and the state) and by the index
         * of the state in the basis (which is something like the principal
         * quantum number lessened by the angular momentum).
         * @param ell Angular momentum.
         * @param n Index of the eigenstate.
         */
        const rArrayView orbital (int ell, int n) const;
        
        /**
         * @brief Size of the basis.
         * 
         * Report back the size of the basis for the specified angular momentum.
         * If called with negative angular momentum, it will report the number of
         * available angular momenta ("maxell + 1").
         */
        int size (int ell = -1) const;
        
    private:
        
        /// Maximal angular number of the bases.
        int maxell_;
        
        /// Screening constant of the bases.
        //@{
        Array<symbolic::rational> rat_lambda_;
        rArray lambda_;
        //@}
        
        /**
         * @brief Eigenenergies.
         * 
         * Contains a vector of arrays. Every array belongs to the basis with
         * corresponding angular momentum. Every element of such array is the
         * eigenenergy of a certain state that is stored in @ref expansions_.
         */
        std::vector<rArray> energies_;
        
        /**
         * @brief Eigenstates.
         * 
         * These are the expansions of the eigenstates of the Hamiltonian in the
         * Laguerre basis. It is a vector with an element for every available angular
         * momentum ("maxell + 1" in total). Each of these elements contains an
         * object @ref RowMatrix that stores a matrix of eigenstates. Each row of the
         * matrix corresponds to one of the eigenstates.
         */
        std::vector<RowMatrix<double>> expansions_;
};

#endif /* HEX_CCC_BASIS */
