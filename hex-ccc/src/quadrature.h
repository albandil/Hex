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

#ifndef HEX_CCC_QUADRATURE
#define HEX_CCC_QUADRATURE

#include <map>
#include <tuple>

#include <gsl/gsl_integration.h>

#include "basis.h"

class QuadratureRule
{
    public:
        
        /**
         * @brief Initialize.
         * 
         * This constructor will initialize the QuadratureRule object.
         * It will precompute the positions of the singularity and
         * setup the integration nodes for every @f$ L @f$ (atomic angular
         * momentum = angular momentum of the specific basis), @f$ N @f$
         * (specific eigenstate of the basis = atomic orbital) and
         * @f$ l @f$ (the angular momentum of the partial wave that is
         * to be integrated.
         * @param basis the LaguerreBasis object containing a valid set of bases.
         * @param E the total energy of the system.
         */
        QuadratureRule (LaguerreBasis const & basis, double E);
        
        /**
         * @brief Get shallow copy of the quadrature nodes.
         * 
         * If this function is called without arguments, it will return a view
         * of the whole array of the quadrature nodes. Otherwise it will return
         * a section of the nodes array that corresponds to the specification.
         * @param ell Angular momentum of the atomic electron (basis).
         * @param l Angular momentum of the projectile (or propagator).
         * @param n Index of the atomic state (eigenstate).
         */
        const rArrayView nodes (int ell = -1, int l = -1, int n = -1) const;
        
        /**
         * @brief Get shallow copy of the quadrature weights.
         * 
         * If this function is called without arguments, it will return a view
         * of the whole array of the quadrature weights. Otherwise it will return
         * a section of the weights array that corresponds to the specification.
         * @param ell Angular momentum of the atomic electron (basis).
         * @param l Angular momentum of the projectile (or propagator).
         * @param n Index of the atomic state (eigenstate).
         */
        const rArrayView weights (int ell = -1, int l = -1, int n = -1) const;
        
    private:
        
        /**
         * @brief Auxiliary function.
         * 
         * This function will fill the arrays "nodes" and "weights" with quadrature
         * nodes and weights, respectively, for the integration over linear momentum
         * of the Green's function. The linear momentum @f$ k @f$ is taken from
         * the linear interval (a,b). The simple Gauss-Legendre quadrature weights
         * will be further modified by the Green function
         * @f[
         *     G = \frac{1}{E - \epsilon_{n\ell} - \epsilon_{k}}\ ,
         * @f]
         * where @f$ \epsilon_{n\ell} @f$ is the energy of the atom (= basis eigenstate)
         * and @f$ \epsilon_{k} @f$ is the energy of the propagator.
         */
        void linearWeightsAndNodes (
            double E, int ell, int n,
            gsl_integration_glfixed_table const * t,
            rArrayView nodes, rArrayView weights,
            double a, double b
        );
        
        /// Reference to the Laguerre basis.
        LaguerreBasis const & basis_;
        
        /// Maximal projectile (or propagator) angular momentum. FIXME Move somewhere else.
        int maxpell_;
        
        /**
         * @brief Index array.
         * 
         * This is an auxiliary dataset that contains indices of first quadrature
         * node (or weight) in the storage array "nodes_" (or "weights_") for
         * a particular set of quantum numbers @f$ (L, l, N) @f$. For example, if
         * one needed to retrieve quadrature nodes for the combination of the
         * quantum numbers (0, 0, 1), one would write
           @code
               nodes_.slice (
                   indices_[std::make_tuple(0, 0, 1)],
                   indices_[std::make_tuple(0, 0, 2)],
               )
           @endcode
         * This particular slice is the start section of "nodes_", because no
         * smaller quantum numbers can be specified.
         */
        std::map<std::tuple<int,int,int>,int> indices_;
        
        /**
         * @brief Quadrature point counts.
         * 
         * This array hold the number of quadrature points for every combination
         * of the quantum numbers @f$ \ell @f$ (atomic angular momentum, or the
         * angular momentum of the basis), @f$ n @f$ (atomic principal quantum
         * number, or the eigenstate index within a basis with fixed @f$ \ell @f$)
         * and @f$ l @f$ (the angular momentum of the projectile or propagator).
         */
        iArray npoints_;
        
        /**
         * @brief Quadrature nodes.
         * 
         * This array contains quadrature nodes in the same amount and order as there
         * are unknowns in the final matrix equation. That means that there is
         * a number for every angular momentum of the atomic electron (= angular 
         * momentum of the basis), index of the atomic orbital (= eigenstate),
         * angular momentum of the projectile (or propagator) and finally for every
         * quadrature node. The number itself is the linear momentum of the projectile
         * or propagator.
         */
        rArray nodes_;
        
        /**
         * @brief Quadrature weights.
         * 
         * This array contains quadrature weights in the same amount and order as there
         * are unknowns in the final matrix equation. That means that there is
         * a number for every angular momentum of the atomic electron (= angular 
         * momentum of the basis), index of the atomic orbital (= eigenstate),
         * angular momentum of the projectile (or propagator) and finally for every
         * quadrature node.
         */
        rArray weights_;
};

#endif /* HEX_CCC_QUADRATURE */
