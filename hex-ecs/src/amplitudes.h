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

#ifndef HEX_AMPLITUDES
#define HEX_AMPLITUDES

#include <vector>

#include "arrays.h"

/**
 * Compute radial integrals for evaluation of the discrete T-matrices,
 * \f[
 *    \Lambda_l^{(1)LMS} = 
 *      \int_0^{R_0}
 *        P_{n_f l_f}(r_1)
 *        \mathcal{W}\left[
 *           \psi_{\ell_1 \ell_2}^{LMS}(r_1,\bullet),
 *           \hat{j}_l(k_f\bullet)
 *        \right]_{R_0}
 *        \mathrm{d}r_1 \ .
 * \f]
 * 
 * \param kf Outgoing projectile momenta.
 * \param ki incoming projectile momenta.
 * \param maxell Maximal angular momentum of the projectile.
 * \param L Total angular momentum (partial wave).
 * \param Spin Conserved total spin.
 * \param ni Initial atomic state - principal quantum number.
 * \param li Initial atomic state - orbital quantum number.
 * \param mi Initial atomic state - magnetic quantum number.
 * \param Ei Initial projectile energies.
 * \param lf Final atomic orbital momentum.
 * \param Pf_overlaps Hydrogenic function B-spline overlap integrals.
 * \return Vector of radial integrals.
 */
cArray computeLambda (
	rArray const & kf,
	rArray const & ki,
	int maxell,
	int L, int Spin,
	int ni, int li, int mi,
	rArray const & Ei,
	int lf,
	cArray const & Pf_overlaps,
	std::vector<std::pair<int,int>> const & coupled_states
);

/**
 * Compute hyperangular integrals for evaluation of the ionization T-matrices,
 * \f[
 *     \Xi_{\ell_1 \ell_2}^{LS}(k_1,k_2) =
 *       \int_0^{\pi/2}
 *         \left(
 *           F_{\ell_1}(k_1,r_1) F_{\ell_2}(k_2,r_2)
 *           \frac{\partial}{\partial\rho}
 *           \psi^{LS}_{\ell_1\ell_2}(r_1,r_2)
 *           -
 *           \psi^{LS}_{\ell_1\ell_2}(r_1,r_2)
 *           \frac{\partial}{\partial\rho}
 *           F_{\ell_1}(k_1,r_1) F_{\ell_2}(k_2,r_2)
 *         \right)
 *         \rho\mathrm{d}\alpha \ ,
 * \f]
 * where \f$ r_1 = \rho \cos\alpha \f$ and \f$ r_2 = \rho \sin\alpha \f$.
 * 
 * \param maxell Maximal angular momentum of the free electrons.
 * \param L Total angular momentum (partial wave).
 * \param Spin Conserved total spin.
 * \param ni Initial atomic state - principal quantum number.
 * \param li Initial atomic state - orbital quantum number.
 * \param mi Initial atomic state - magnetic quantum number.
 * \param Ei Initial projectile energies.
 * \param ics Ionization cross section (on return).
 * \return Vector of radial integrals.
 */
cArrays computeXi (
	int maxell,
	int L, int Spin,
	int ni, int li, int mi,
	rArray const & Ei,
	rArray & ics,
	std::vector<std::pair<int,int>> const & coupled_states
);

#endif
