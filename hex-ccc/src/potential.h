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

#ifndef HEX_CCC_POTENTIAL
#define HEX_CCC_POTENTIAL

#include "arrays.h"
#include "basis.h"
#include "gauss.h"
#include "quadrature.h"

class PotentialMatrix
{
    public:
        
        /**
         * @brief Construct the matrix.
         * 
         * This function will precompute entries of the potential matrix based
         * on the chosen basis and quadrature rule. This is one of the most
         * time expensive operations because for every element of the matrix
         * a double integral has to be evaluated.
         * @param basis LaguerreBasis object containing the bases.
         * @param quadrature QuadratureRule object containing the quadrature settings.
         * @param J Total angular momentum.
         * @param S Total spin.
         * @param Pi Total parity.
         */
        PotentialMatrix
        (
            LaguerreBasis const & basis,
            QuadratureRule const & quadrature,
            int J, int S, int Pi
        );
        
        /**
         * @brief Return reference to the matrix of the potential.
         */
        RowMatrix<double> const & matrix() const;
        
    private:
        
        /// Basis.
        LaguerreBasis const & basis_;
        
        /// Quadrature rule.
        QuadratureRule const & quadrature_;
        
        /// The (symmetrical) matrix of the potential.
        RowMatrix<double> matrix_;
        
        
        GaussLegendre<double> g_;
        
        /**
         * @brief Double integral @f$ I_{\mathrm{dir}} @f$.
         */
        double ComputeIdir
        (
            int lambda,
            int L, int i, int l, double k,
            int Lp, int ip, int lp, double kp
        ) const;
        
        double ComputeJdir
        (
            int lambda,
            symbolic::poly const & xi, int l, double k,
            symbolic::poly const & xip, int lp, double kp
        ) const;
};

class MatrixEquation
{
    public:
        
        // constructor
        MatrixEquation
        (
            QuadratureRule const & quadrature,
            PotentialMatrix const & potential
        );
        
        // solve the equation
        rArray solve () const;
};

#endif /* HEX_CCC_POTENTIAL */

