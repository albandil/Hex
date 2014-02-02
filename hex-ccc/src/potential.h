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
        
        /// The (symmetrical) matrix of the potential.
        RowMatrix<double> matrix_;
};

class MatrixEquation
{
    public:
        
        // constructor
        MatrixEquation (
            QuadratureRule const & quadrature,
            PotentialMatrix const & potential
        );
        
        // solve the equation
        rArray solve () const;
};

/**
 * @brief Compute the direct potential term.
 * 
 * This function will compute the direct potential by evaluating the double
 * integral
 * @f[
 *     I_{\mathrm{dir},NLkl,N'L'k'l'}^\lambda
 *     =
 *     \int\limits_0^\infty \int\limits_0^\infty
 *     P_{NL}(r_1) P_{N'L'}(r_1)
 *     \left(
 *         \frac{r_<^\lambda}{r_>^\lambda}
 *         -
 *         \frac{\delta_\lambda^0}{r_2}
 *     \right)
 *     \hat{j}_l(kr_2) \hat{j}_{l'}(k'r_2)
 *     \mathrm{d}r_1 \mathrm{d}r_2 \ .
 * @f]
 * Here @f$ \hat{j}_l(kr) @f$ is the Riccati-Bessel function and @f$ P_{NL}(r) @f$
 * is the atomic orbital obtained by diagonalization of the hydrogenic Hamiltonian.
 * The orbital is supplied as an array of expansion coefficients in the Laguerre basis.
 * 
 * @param basis The objects managing Laguerre basis for all angular momenta.
 * @param k Linear momentum of the projectile or propagator (initial).
 * @param kp Linear momentum of the projectile or propagator (final).
 * @param l Angular momentum of the projectile or propagator (initial).
 * @param lp Angular momentum of the projectile or propagator (final).
 * @param L Angular momentum of Laguerre basis to use (initial).
 * @param Lp Angular momentum of Laguerre basis to use (final).
 * @param PNL Expansion coefficients of the (initial) atomic orbital in the
 *            Laguerre basis with angular momentum L.
 * @param PNLp Expansion coefficients of the (initial) atomic orbital in the
 *            Laguerre basis with angular momentum Lp.
 */
double compute_Vdir (
    LaguerreBasis const & basis, int lambda,
    double k, int l, rArray const & PNL, int L,
    double kp, int lp, rArray const & PNLp, int Lp
);

/**
 * @brief Compute the exchange potential term.
 * 
 * This function will compute the exchange potential by evaluating the double
 * integral
 * @f[
 *     I_{\mathrm{dir},NLkl,N'L'k'l'}^\lambda
 *     =
 *     \int\limits_0^\infty \int\limits_0^\infty
 *     \hat{j}_l(kr_1) P_{N'L'}(r_1)
 *     \left(
 *         \frac{r_<^\lambda}{r_>^\lambda}
 *         -
 *         \frac{\delta_\lambda^0}{r_1}
 *     \right)
 *     P_{NL}(r_2) \hat{j}_{l'}(k'r_2)
 *     \mathrm{d}r_1 \mathrm{d}r_2 \ .
 * @f]
 * Here @f$ \hat{j}_l(kr) @f$ is the Riccati-Bessel function and @f$ P_{NL}(r) @f$
 * is the atomic orbital obtained by diagonalization of the hydrogenic Hamiltonian.
 * The orbital is supplied as an array of expansion coefficients in the Laguerre basis.
 * 
 * @param basis The objects managing Laguerre basis for all angular momenta.
 * @param k Linear momentum of the projectile or propagator (initial).
 * @param kp Linear momentum of the projectile or propagator (final).
 * @param l Angular momentum of the projectile or propagator (initial).
 * @param lp Angular momentum of the projectile or propagator (final).
 * @param L Angular momentum of Laguerre basis to use (initial).
 * @param Lp Angular momentum of Laguerre basis to use (final).
 * @param PNL Expansion coefficients of the (initial) atomic orbital in the
 *            Laguerre basis with angular momentum L.
 * @param PNLp Expansion coefficients of the (initial) atomic orbital in the
 *            Laguerre basis with angular momentum Lp.
 */
double compute_Vexc (
    LaguerreBasis const & basis, int lambda,
    double k, int l, rArray const & PNL, int L,
    double kp, int lp, rArray const & PNLp, int Lp
);

#endif /* HEX_CCC_POTENTIAL */

