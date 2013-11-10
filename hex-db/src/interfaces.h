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

#ifndef HEX_DB_INTERFACES_H
#define HEX_DB_INTERFACES_H

/**
 * Create new database
 */
void create_new_database();

/**
 * Initialize
 */
void initialize(const char* dbname);

/**
 * Import SQL batch file.
 * Equivalent to
 * \code
 *     sqlite3 "hex.db" < "batchfile.sql"
 * \endcode
 */
void import(const char* sqlname);

/**
 * Update integral cross section.
 * Used after import of T-matrices.
 */
void update();

/**
 * \brief Optimize the SQLite database file.
 * 
 * Reduces the occupied space.
 */
void optimize();

/**
 * \brief Print data summary.
 * 
 * Prints highest partial waves for consecutive energy ranges per initial state.
 * Example:
 * \code
   1   0   0    0.05     0.6    3
   1   0   0    0.65     0.85   4
   1   0   0    0.8501   0.89   3
   1   0   0    0.8905   0.97   0
   1   0   0    1.2      5.0    6
   1   0   0    5.2     40.0    9
   1   0   0   40.0     40.0    9
 * \endcode
 */
void avail();

/**
 * \brief Dump contents of the T-matrix table.
 * 
 * The output can be used to construct an equivalent table.
 * The corresponding code is
 * \code
 * sqlite> .mode insert
 * sqlite> select * from tmat
 * \endcode
 */
void dump(const char* dumpname);

#if 0

/**
 * Return the scattering amplitude.
 * \f[
 *      f^S = -\frac{1}{2\pi} \sum_{\ell L} T^{LMS}_\ell Y_{l_f M-m_f}
 * \f]
 * \param ni Initial atomic principal quantum number.
 * \param li Initial atomic orbital quantum number.
 * \param mi Initial atomic magnetic quantum number.
 * \param nf Final atomic principal quantum number.
 * \param lf Final atomic principal quantum number.
 * \param mf Final atomic magnetic quantum number.
 * \param S Total spin (only needed for proper anti/symmetrization
 *          of the initial wave function).
 * \param Ei Projectile energy in Rydbergs (\f$ = k_i^2 \f$)
 * \param N Array length.
 * \param theta Polar angle of scattered electron.
 * \param amplitudes Results (N complex numbers as 2N real numbers).
 */
extern "C" void scattering_amplitude (
    int const * ni, int const * li, int const * mi,
    int const * nf, int const * lf, int const * mf,
    int const * S, double const * Ei,
    int const * N, double const * theta, double * amplitudes
);

/**
 * Compute the differential cross section.
 * \f[
 *       \frac{\mathrm{d}\sigma^S}{\mathrm{d}\Omega} = \frac{k_f}{k_i} \frac{2S+1}{4} |f^S|^2
 * \f]
 * \param ni Initial atomic principal quantum number.
 * \param li Initial atomic orbital quantum number.
 * \param mi Initial atomic magnetic quantum number.
 * \param nf Final atomic principal quantum number.
 * \param lf Final atomic principal quantum number.
 * \param mf Final atomic magnetic quantum number.
 * \param S Total spin (only needed for proper anti/symmetrization
 *          of the initial wave function).
 * \param N Array length.
 * \param Ei Projectile energy in Rydbergs (\f$ = k_i^2 \f$)
 * \param theta Polar angle of scattered electron.
 * \param dsigma Results.
 */
extern "C" void differential_cross_section (
    int const * ni, int const * li, int const * mi,
    int const * nf, int const * lf, int const * mf,
    int const * S, double const * Ei,
    int const * N, double const * theta, double * dsigma
);

/**
 * Compute the momentum transfer.
 * \f[
 *      \eta = \int \frac{\mathrm{d}\sigma}{\mathrm{d}\Omega} (1 - \cos\theta) \mathrm{d}\Omega(\hat{\vec{k}}_f)
 *      = ...
 * \f]
 * \param ni Initial atomic principal quantum number.
 * \param li Initial atomic orbital quantum number.
 * \param mi Initial atomic magnetic quantum number.
 * \param nf Final atomic principal quantum number.
 * \param lf Final atomic principal quantum number.
 * \param mf Final atomic magnetic quantum number.
 * \param L Total angular momentum of the electrons.
 * \param S Total spin (only needed for proper anti/symmetrization
 *          of the initial wave function).
 * \param N Array length.
 * \param Ei Projectile energy in Rydbergs (\f$ = k_i^2 \f$)
 * \param eta Results.
 */
extern "C" void momentum_transfer (
    int const * ni, int const * li, int const * mi,
    int const * nf, int const * lf, int const * mf,
    int const * L, int const * S,
    int const * N, double const * Ei, double * eta
);

/**
 * Compute integral cross section \f$ \sigma^{LS} \f$,
 * \f[
 *      \sigma^{LS}(E) = \int_{4\pi} \frac{\mathrm{d}\sigma^S}{\mathrm{d}\Omega}
 *      \mathrm{d}\Omega(\hat{\vec{k}}_f)
 *      \Rightarrow
 *      \sigma^{LS}(E) = \frac{k_f}{k_i} \frac{1}{4\pi^2} \frac{2S+1}{4} \sum_\ell \left|T_\ell^{LMS}\right|^2
 * \f]
 * \param ni Initial atomic principal quantum number.
 * \param li Initial atomic orbital quantum number.
 * \param mi Initial atomic magnetic quantum number.
 * \param nf Final atomic principal quantum number.
 * \param lf Final atomic principal quantum number.
 * \param mf Final atomic magnetic quantum number.
 * \param L Total angular momentum of the electrons.
 * \param S Total spin (only needed for proper anti/symmetrization
 *          of the initial wave function).
 * \param N Array length.
 * \param Ei Projectile energy in Rydbergs (\f$ = k_i^2 \f$)
 * \param sigma Results.
 */
extern "C" void integral_cross_section (
    int const * ni, int const * li, int const * mi,
    int const * nf, int const * lf, int const * mf,
    int const * L, int const * S,
    int const * N, double const * Ei, double * sigma
);

/**
 * Sum integral cross sections.
 * \f[
 *      \sigma_{i \rightarrow f}(E) = \sum_{LS} \sigma_{i \rightarrow f}^{LS}(E)
 * \f]
 * \param ni Initial atomic principal quantum number.
 * \param li Initial atomic orbital quantum number.
 * \param mi Initial atomic magnetic quantum number.
 * \param nf Final atomic principal quantum number.
 * \param lf Final atomic principal quantum number.
 * \param mf Final atomic magnetic quantum number.
 * \param N Array length.
 * \param Ei Projectile energy in Rydbergs (\f$ = k_i^2 \f$).
 * \param sigma Results.
 */
extern "C" void complete_cross_section (
    int const * ni, int const * li, int const * mi,
    int const * nf, int const * lf, int const * mf,
    int const * N, double const * Ei, double * sigma
);

/**
 * Sum integral cross sections and extrapolate \f$ L \rightarrow \infty \f$.
 * \f[
 *      \sigma_{i \rightarrow f}(E) = \sum_{LS} \sigma_{i \rightarrow f}^{LS}(E)
 * \f]
 * The extrapolation is done using Aitken Δ²-process. 
 * \param ni Initial atomic principal quantum number.
 * \param li Initial atomic orbital quantum number.
 * \param mi Initial atomic magnetic quantum number.
 * \param nf Final atomic principal quantum number.
 * \param lf Final atomic principal quantum number.
 * \param mf Final atomic magnetic quantum number.
 * \param N Array length.
 * \param Ei Projectile energy in Rydbergs (\f$ = k_i^2 \f$)
 * \param sigma Results.
 */
extern "C" void extrapolate_cross_section (
    int const * ni, int const * li, int const * mi,
    int const * nf, int const * lf, int const * mf,
    int const * N, double const * Ei, double * sigma
);

/**
 * Compute the collision strength.
 * \f[
 *     \Omega_{i \rightarrow f}(E) = k_i^2 (2L+1) (2S+1) \sigma_{i \rightarrow f}(E)
 * \f]
 * \param ni Initial atomic principal quantum number.
 * \param li Initial atomic orbital quantum number.
 * \param mi Initial atomic magnetic quantum number.
 * \param nf Final atomic principal quantum number.
 * \param lf Final atomic principal quantum number.
 * \param mf Final atomic magnetic quantum number.
 * \param L Total angular momentum of the electrons.
 * \param S Total spin (only needed for proper anti/symmetrization
 *          of the initial wave function).
 * \param N Array length.
 * \param Ei Projectile energy in Rydbergs (\f$ = k_i^2 \f$)
 * \param omega Results.
 */
extern "C" void collision_strength (
    int const * ni, int const * li, int const * mi,
    int const * nf, int const * lf, int const * mf,
    int const * L, int const * S,
    int const * N, double const * Ei, double * omega
);

/**
 * Compute total cross section.
 * \f[
 *      \sigma_{i \rightarrow *}(E) = \sum_{n_f = 0}^\infty \sum_{l_f = 0}^{n_f - 1}
 *             \sigma_{i \rightarrow f}(E)
 * \f]
 * \param ni Initial atomic principal quantum number.
 * \param li Initial atomic orbital quantum number.
 * \param mi Initial atomic magnetic quantum number.
 * \param N Array length.
 * \param Ei Projectile energy in Rydbergs (\f$ = k_i^2 \f$).
 * \param sigma Results.
 */
extern "C" void total_cross_section (
    int const * ni, int const * li, int const * mi,
    int const * N, double const * Ei, double * sigma
);

#endif

#endif
