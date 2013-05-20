/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2012                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HEX_ANGS
#define HEX_ANGS

/**
 * \return Value of \f$ f(\lambda,l_1,l_2,l_1',l_2',L) \f$.
 */
double computef(int lambda, int l1, int l2, int l1p, int l2p, int L);

/**
 * Clebsch-Gordan coefficient. In present implementation valid only for
 * integer (not half-integer) angular momenta. [One needs to correct the signs!]
 */
double ClebschGordan(int l1, int m1, int l2, int m2, int L, int M);

/**
 * Gaunt's integral.
 * \f[
 * \int_{4\pi} Y_{l_1m_1} Y_{l_2m_2} Y^{\ast}_{lm} \mathrm{d}\Omega \ .
 * \f]
 */
double Gaunt(int l1, int m1, int l2, int m2, int l, int m);

/**
 * Compute number of angular momenta pairs \f$ \ell_1 \f$ and \f$ \ell_2 \f$ that
 * are less than or equal to "maxell" and compose the total angular momentum \f$ L \f$.
 */
int triangle_count(int L, int maxell);

#endif
