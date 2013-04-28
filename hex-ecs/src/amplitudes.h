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
 * Compute the cross sections for all energies.
 * \param kf Outgoing projectile momenta.
 * \param ki incoming projectile momenta.
 * \param maxL2 Maximal angular momentum of the projectile.
 * \param solution Scattered waves.
 * \param L Total angular momentum (partial wave).
 * \param Spin Conserved total spin.
 * \param ni
 * \param li
 * \param mi
 * \param Ei
 * \param lf Final atomic orbital momentum.
 * \param Pf_overlaps Hydrogenic function B-spline overlap integrals.
 * \param jf_overlaps Riccati-Bessel function B-spline overlap integrals.
 * \return Vector of radial integrals.
 */
cArray computeLambda(
	const rArray& kf, const rArray& ki, int maxell,
	int L, int Spin, int ni, int li, int mi, const rArray& Ei,
	int lf, cArray Pf_overlaps, cArray jf_overlaps
);

#endif
