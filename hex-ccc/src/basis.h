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

#ifndef HEX_CCC_BASIS
#define HEX_CCC_BASIS

#include "arrays.h"
#include "matrix.h"

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
         *               the screening constant @f$ \lamnda_\ell @f$.
         */
        LaguerreBasis (int maxell, iArray Nl, rArray lambda);
        
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
        rArray lambda_;
        
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
