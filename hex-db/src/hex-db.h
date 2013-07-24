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

#ifndef HEX_HEX_DB
#define HEX_HEX_DB

#include <map>
#include <string>
#include <vector>

#include "arrays.h"
#include "complex.h"
#include "interfaces.h"

/**
 * \mainpage
 * 
 * \section theory Theory
 * 
 * The program <b>Hex</b> produces output in the form of radial part
 * of the T-matrix. The full T-matrix for transition \f$ (n_i,l_i,m_i) \rightarrow
 * (n_f,l_f,m_f) \f$ at energy \f$ E = k_i^2 - n_i^{-2} = k_f^2 - n_f^{-2} \f$
 * (energy is in Rydbergs) is
 * \f[
 *        T_{n_i l_i m_i\rightarrow n_f l_f m_f}^{L M S} (E) =
 *        \delta_{Mm_f}
 *        \sum_\ell T_\ell^{LMS}(E) Y_{l_f m_f}(\hat{\vec k}_f) ,
 * \f]
 * where the radial part depends on theory. It will be much different for exterior
 * complex scaling and for Born approximation. This utility code (and associated
 * standalone interface library) can then be used to produce following information:
 *
 * - the scattering amplitude ( \ref scattering_amplitude )
 *   \f[
 *       f_{i \rightarrow f}^{S}(E) = -\frac{1}{2\pi} T_{i \rightarrow f}^{S}(E) =
 *       -\frac{1}{2\pi} \sum_{L} f_{i \rightarrow f}^{LS}(E) = 
 *       -\frac{1}{2\pi} \sum_{\ell L} Y_{\ell m_i-m_f}(\hat{\vec{k}}_f) T_{fi,\ell}^{LS}(E) \ ,
 *   \f]
 * - the differential cross section ( \ref differential_cross_section )
 *   \f[
 *       \frac{\mathrm{d}\sigma_{i \rightarrow f}^{S}}{\mathrm{d}\Omega}(E) = 
 *       \frac{k_f}{k_i} \frac{2S+1}{4} |f_{i \rightarrow j}^{S}|^2 \ ,
 *   \f]
 * - the partial integral cross section ( \ref integral_cross_section )
 *   \f[
 *       \sigma_{i \rightarrow f}^{LS}(E) = \frac{k_f}{k_i} \frac{2S + 1}{4}
 *       \int_{4\pi} |f_{i \rightarrow f}^{LS}(E)|^2 \mathrm{d}\Omega(\hat{\vec{k}}_f)
 *       = \frac{k_f}{k_i} \frac{2S+1}{16\pi^2} \sum_\ell |T_{fi,\ell}^{LS}(E)|^2 \ ,
 *   \f]
 * - the “complete” integral cross section ( \ref complete_cross_section, \ref extrapolate_cross_section )
 *   \f[
 *       \sigma_{i \rightarrow f}(E) = \sum_{LS} \sigma_{i \rightarrow f}^{LS}(E) \ ,
 *   \f]
 * - the total cross section ( \ref total_cross_section )
 *   \f[
 *       \sigma_i(E) = \sum_{n_f = 1}^\infty \sum_{l_f = 0}^{n_f - 1}
 *       \sum_{m_f = -l_f}^{l_f} \sigma_{i \rightarrow f}(E) \ ,
 *   \f]
 * - and also the dimensionless, $i \leftrightarrow j$ symmetrical, collision strength ( \ref collision_strength )
 *   \f[
 *       \Omega_{i \rightarrow f}^{LS} (E) = k_i^2 \sigma_{i \rightarrow f}^{LS}(E)
 *   \f]
 * - or the momentum transfer ( \ref momentum_transfer )
 *   \f[
 *       \eta_{i \rightarrow f}^{LS} = \int \frac{\mathrm{d}\sigma_{i \rightarrow f}^{LS}}{\mathrm{d}\Omega}(E) (1 - \cos\vartheta)
 *       \mathrm{d}\Omega(\hat{\vec{k}}_f) \ .
 *   \f]
 * 
 * \section database The database
 * 
 * The file "hex.db" contains SQL database in the form of one table and following columns:
 * <center><table>
 * <tr><th>Field name</th> <th>Data type</th> <th>Comment</th></tr>
 * <tr><td>ni</td> <td>INTEGER</td> <td>Initial principal quantum number.</td></tr>
 * <tr><td>li</td> <td>INTEGER</td> <td>Initial orbital quantum number.</td></tr>
 * <tr><td>mi</td> <td>INTEGER</td> <td>Initial magnetic quantum number.</td></tr>
 * <tr><td>nf</td> <td>INTEGER</td> <td>Final principal quantum number.</td></tr>
 * <tr><td>lf</td> <td>INTEGER</td> <td>Final orbital quantum number.</td></tr>
 * <tr><td>mf</td> <td>INTEGER</td> <td>Final magnetic quantum number.</td></tr>
 * <tr><td>L</td> <td>INTEGER</td> <td>Conserved total angular momentum of a partial wave.</td></tr>
 * <tr><td>S</td> <td>INTEGER</td> <td>Conserved total spin of the system.</td></tr>
 * <tr><td>Ei</td> <td>DOUBLE PRECISION</td> <td>Projectile impact energy.</td></tr>
 * <tr><td>Re_T_ell</td> <td>DOUBLE PRECISION</td> <td>%Real part of \f$ T_\ell \f$.</td></tr>
 * <tr><td>Im_T_ell</td> <td>DOUBLE PRECISION</td> <td>Imaginary part of \f$ T_\ell \f$.</td></tr>
 * </table></center>
 * 
 * The initial database file can be created by the shell command
 * <pre>
 * sqlite3 hex.db 'create table "hex" (
 *    ni integer, li integer, mi integer,
 *    nf integer, lf integer, mf integer,
 *    L integer, S integer, Ei double precision, ell integer,
 *    Re_T_ell double precision, Im_T_ell double precision,
 *    primary key (ni,li,mi,nf,lf,mf,L,S,Ei,ell)
 * )'
 * </pre>
 * which of course requires 'sqlite' being installed on the system. The insertion of
 * scattering T-matrices, as produced by computational units, is then very easy:
 * <pre>
 * sqlite3 hex.db < T-matrices.sql
 * </pre>
 * The raw database file allows standard queriyng, which can be employed to
 * retrieve some low level data, like
 * <pre>
 * sqlite3 -column hex.db 'select Ei,Re_T_ell,Im_T_ell from "hex" where
 *    ni = 1 and li = 0 and mi = 0 and
 *    nf = 1 and lf = 0 and mf = 0 and
 *    Ei > 0.75 and Ei <= 0.88' > output.txt
 * </pre>
 * which will print real and imaginary parts of H(1s) \f$ \rightarrow \f$ H(1s)
 * T-matrices for energies between \f$ E = 0.75 \f$&nbsp;Ry and \f$ E = 0.88 \f$&nbsp;Ry into
 * the file "output.txt".
 * 
 * \section usage Usage
 * 
 * The code can be easily embedded into other programs, both C and C++ headers
 * are provided. The former can be used also with Fortran. Apart from direct
 * program access, one can use a CLI (command line interface) which has the
 * following switches:
 * <center><table>
 * <tr><th>Short option</th> <th>Long option</th> <th>Usage</th> </tr>
 * <tr><td>-h</td> <td>--help</td> <td>Print basic usage info and exit.</td></tr>
 * <tr><td>-D</td> <td>--database</td> <td>Use given database file.</td></tr>
 * <tr><td>-1</td> <td>--scatamp</td> <td>Write out scattering amplitudes.</td></tr>
 * <tr><td>-2</td> <td>--dcs</td> <td>Write out differential cross sections.</td></tr>
 * <tr><td>-3</td> <td>--momtransf</td> <td>Write out Momentum transfer.</td></tr>
 * <tr><td>-4</td> <td>--integcs</td> <td>Write out integral cross sections.</td></tr>
 * <tr><td>-5</td> <td>--sumintegcs</td> <td>Write out (summed) complete cross sections.</td></tr>
 * <tr><td>-6</td> <td>--collstr</td> <td>Write out collision strength.</td></tr>
 * <tr><td>-7</td> <td>--tcs</td> <td>Write out total cross sections.</td></tr>
 * <tr><td>-a</td> <td>--ni</td> <td>Initial principal quantum number.</td></tr>
 * <tr><td>-b</td> <td>--li</td> <td>Initial orbital quantum number.</td></tr>
 * <tr><td>-c</td> <td>--mi</td> <td>Initial magnetic quantum number.</td></tr>
 * <tr><td>-d</td> <td>--nf</td> <td>Final principal quantum number.</td></tr>
 * <tr><td>-e</td> <td>--lf</td> <td>Final orbital quantum number.</td></tr>
 * <tr><td>-f</td> <td>--mf</td> <td>Final magnetic quantum number.</td></tr>
 * <tr><td>-L</td> <td>--L</td> <td>Conserved total angular momentum of a partial wave.</td></tr>
 * <tr><td>-S</td> <td>--S</td> <td>Conserved total spin of the system.</td></tr>
 * <tr><td>-E</td> <td>--E</td> <td>Projectile impact energy in Rydbergs.</td></tr>
 * </table></center>
 * 
 * <br/>
 * 
 * <center><table>
 * <tr><th>Quantity</th><th>Compulsory quantum numbers</th> <th>Optional quantum numbers*</th> <th>STDIN contains**</th></tr>
 * <tr><td>Scattering amplitude</td> <td>ni, li, mi, nf, lf, mf, S, E</td> <td>---</td> <td>angles [rad]</td></tr>
 * <tr><td>Differential cross section</td> <td>ni, nf, E</td> <td>li, lf, mi, mf, S</td> <td>angles [rad]</td></tr>
 * <tr><td>Momentum transfer</td> <td>ni, li, mi, nf, lf, mf, L, S</td> <td>---</td> <td>energies [Ry]</td></tr>
 * <tr><td>Integral cross section</td> <td>ni, li, mi, nf, lf, mf, L, S</td> <td>---</td> <td>energies [Ry]</td></tr>
 * <tr><td>Complete cross section</td> <td>ni, nf</td> <td>li, lf, mi, mf</td> <td>energies [Ry]</td></tr>
 * <tr><td>Collision strength</td> <td>ni, li, mi, nf, lf, mf, L, S</td> <td>---</td> <td>energies [Ry]</td></tr>
 * <tr><td>Total cross section</td> <td>ni, li, mi</td> <td>---</td> <td>energies [Ry]</td></tr>
 * </table></center>
 * 
 * *) Will be summed over if not specified.
 * 
 * **) Standard input should contain values of angles or energies separated by white characters (spaces, newlines etc.).
 *     The list ends with EOF (Ctrl+D if manually entering data; inserted automatically when using standard "seq" invocation).
 * 
 * <br/>
 * 
 * For every STDIN entry the program will respond with one number computed (and interpolated if necessary) 
 * from the database. If there are no relevant data in database, the result will be zero. If the available
 * energy interval doesn't contain some of the required energies, appropriate error message will be printed.
 * Here are some examples of usage:
 * <pre>
 * # scattering amplitude
 * seq 0.01 0.01 3.14    | hex-db --database="hex.db" --scatamp    --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --S=0
 *
 * # differential cross section
 * seq 0.01 0.01 3.14    | hex-db --database="hex.db" --dcs        --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --S=0
 *
 * # momentum transfer
 * seq 0.650 0.001 0.850 | hex-db --database="hex.db" --momtransf  --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --L=0 --S=0
 *
 * # integral cross section
 * seq 0.650 0.001 0.850 | hex-db --database="hex.db" --integcs    --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --L=0 --S=0
 *
 * # sum integral cross section
 * seq 0.650 0.001 0.850 | hex-db --database="hex.db" --sumintegcs --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0
 *
 * # collision strength
 * seq 0.650 0.001 0.850 | hex-db --database="hex.db" --collstr    --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --L=0 --S=0
 *
 * # total cross section
 * seq 0.650 0.001 0.850 | hex-db --database="hex.db" --tcs        --ni=1 --li=0 --mi=0
 * </pre>
 * You may need to set the LC_NUMERIC environment variable to en_US.utf8 (or en_GB.utf8 or any other
 * localitation that uses decimal point instead of decimal comma). An alternative is the workaround with sed.
 * <pre>
 * # scattering amplitude in Czech environment (1)
 * export LC_NUMERIC=en_GB.utf8
 * seq 0.01 0.01 3.14    | hex-db --database="hex.db" --scatamp    --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --S=0
 * 
 * # scattering amplitude in Czech environment (2)
 * seq 0.01 0.01 3.14 | sed -e "s/,/\./g" | hex-db --database="hex.db" --scatamp    --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --S=0
 * </pre>
 */

/// Energy units
enum eUnit {
	eUnit_Ry,	// Rydberg (13.605692 eV, default)
	eUnit_au,	// Hartree (2 Ry)
	eUnit_eV	// electron-Volt
};

/// Output (length) units
enum lUnit {
	lUnit_au,	// atomic units (Bohr radius a₀=5.29x10⁻⁹ cm)
	lUnit_cgs	// centimeters (1 cm = (1cm/a₀) a₀)
};

// Angular units
enum aUnit {
	aUnit_deg,	// degrees
	aUnit_rad	// radians
};

/**
 * Run the computations.
 */
int run (
	eUnit Eunits, lUnit lUnits, aUnit aUnits,
	std::vector<std::string> const & vars,
	std::map<std::string,std::string> const & sdata
);

#endif
