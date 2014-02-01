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

#ifndef HEX_DB_DOC
#define HEX_DB_DOC

/**
 * @mainpage
 * @author Jakub Benda, MFF UK
 * @date 1. 2. 2014
 * @section db Hex-db
 * 
 * Hex-db is the interface program that can be used to access data precomputed
 * by the computational modules. Also, it's built-in algorithms can, on demand,
 * compute various derived scattering quantities. The full list of available
 * variables is given below.
 * 
 * @subsection theory Theory
 * 
 * The computational modules of the package <b>Hex</b> produce output in the form
 * of radial part of the T-matrix. The full T-matrix for transition
 * @f$ (n_i,l_i,m_i) \rightarrow (n_f,l_f,m_f) @f$ at energy
 * @f$ E = k_i^2 - n_i^{-2} = k_f^2 - n_f^{-2} @f$ (energy is in Rydbergs) is
 * @f[
 *        T_{n_i l_i m_i\rightarrow n_f l_f m_f}^{L M S} (E) =
 *        \delta_{Mm_f}
 *        \sum_\ell T_\ell^{LMS}(E) Y_{l_f m_f}(\hat{\vec k}_f) ,
 * @f]
 * where the radial part depends on theory. It will be much different for exterior
 * complex scaling and for Born approximation. This utility code (and associated
 * standalone interface library) can then be used to produce following information:
 *
 * - the scattering amplitude
 *   @f[
 *       f_{i \rightarrow f}^{S}(E) = -\frac{1}{2\pi} T_{i \rightarrow f}^{S}(E) =
 *       -\frac{1}{2\pi} \sum_{L} f_{i \rightarrow f}^{LS}(E) = 
 *       -\frac{1}{2\pi} \sum_{\ell L} Y_{\ell m_i-m_f}(\hat{\vec{k}}_f) T_{fi,\ell}^{LS}(E) \ ,
 *   @f]
 * - the differential cross section
 *   @f[
 *       \frac{\mathrm{d}\sigma_{i \rightarrow f}^{S}}{\mathrm{d}\Omega}(E) = 
 *       \frac{k_f}{k_i} \frac{2S+1}{4} |f_{i \rightarrow j}^{S}|^2 \ ,
 *   @f]
 * - the partial integral cross section
 *   @f[
 *       \sigma_{i \rightarrow f}^{LS}(E) = \frac{k_f}{k_i} \frac{2S + 1}{4}
 *       \int_{4\pi} |f_{i \rightarrow f}^{LS}(E)|^2 \mathrm{d}\Omega(\hat{\vec{k}}_f)
 *       = \frac{k_f}{k_i} \frac{2S+1}{16\pi^2} \sum_\ell |T_{fi,\ell}^{LS}(E)|^2 \ ,
 *   @f]
 * - the “complete” integral cross section
 *   @f[
 *       \sigma_{i \rightarrow f}(E) = \sum_{LS} \sigma_{i \rightarrow f}^{LS}(E) \ ,
 *   @f]
 * - the total cross section
 *   @f[
 *       \sigma_i(E) = \sum_{n_f = 1}^\infty \sum_{l_f = 0}^{n_f - 1}
 *       \sum_{m_f = -l_f}^{l_f} \sigma_{i \rightarrow f}(E) \ ,
 *   @f]
 * - and also the dimensionless, @f$ i \leftrightarrow j @f$ symmetrical, collision strength
 *   @f[
 *       \Omega_{i \rightarrow f}^{LS} (E) = k_i^2 \sigma_{i \rightarrow f}^{LS}(E)
 *   @f]
 * - or the momentum transfer
 *   @f[
 *       \eta_{i \rightarrow f}^{LS} = \int \frac{\mathrm{d}\sigma_{i \rightarrow f}^{LS}}{\mathrm{d}\Omega}(E)
 *        (1 - \cos\vartheta)
 *       \mathrm{d}\Omega(\hat{\vec{k}}_f) \ .
 *   @f]
 * 
 * @subsection database The database
 * 
 * The file "hex.db" contains SQL database in the form of several tables. The most low-level
 * data (the T-matrices) are stored in the table "tmat" that has the following columns:
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
 * <tr><td>Re_T_ell</td> <td>DOUBLE PRECISION</td> <td>%Real part of @f$ T_\ell @f$.</td></tr>
 * <tr><td>Im_T_ell</td> <td>DOUBLE PRECISION</td> <td>Imaginary part of @f$ T_\ell @f$.</td></tr>
 * </table></center>
 * However, there are more tables (currently also "ics" for integral cross section, "ccs"
 * for complete cross section, and "ionf" for radial part of the ionization amplitude).
 * 
 * The initial database file can be created by the shell command
 * <pre>
 * hex-db --new
 * </pre>
 * which is equivalent to (if "sqlite3" program is available in PATH)
 * <pre>
 * sqlite3 hex.db 'create table "tmat" (
 *    ni integer, li integer, mi integer,
 *    nf integer, lf integer, mf integer,
 *    L integer, S integer, Ei double precision, ell integer,
 *    Re_T_ell double precision, Im_T_ell double precision,
 *    primary key (ni,li,mi,nf,lf,mf,L,S,Ei,ell)
 * )'
 * </pre>
 * and several similar statements more, for each variable that has to be precomputed.
 * Whether a variable has a dedicated table within the database can be inquired
 * by looking at its member function SQL_CreateTable. If there are similar commands,
 * one has to update the database after data insertion to be able to use the new data:
 * <pre>
 * hex-db --update
 * </pre>
 * The insertion of scattering T-matrices, as produced by computational units, is very easy:
 * <pre>
 * hex-db --import T-matrices.sql
 * </pre>
 * which is equivalent to
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
 * which will print real and imaginary parts of H(1s) @f$ \rightarrow @f$ H(1s)
 * T-matrices for energies between @f$ E = 0.75 @f$&nbsp;Ry and @f$ E = 0.88 @f$&nbsp;Ry into
 * the file "output.txt".
 * 
 * @subsection usage Usage
 * 
 * Apart from direct program access (when the sources were linked to a custom code),
 * one can use a CLI (command line interface) which has the following switches:
 * <center><table>
 * <tr><th>Short option</th> <th>Long option</th> <th>Usage</th> </tr>
 * <tr><td>-h</td> <td>--help</td> <td>Print basic usage info and exit.</td></tr>
 * <tr><td>-D</td> <td>--database</td> <td>Use given database file. This is optional; "hex.db" is used if omitted.</td></tr>
 * <tr><td>-n</td> <td>--new</td> <td>Create database file, either "hex.db" or the given name if specified by --database.</td></tr>
 * <tr><td>-i</td> <td>--import</td> <td>Import SQL data produced by some of the modules.</td></tr>
 * <tr><td>-u</td> <td>--update</td> <td>Update cache tables (variables dependent on lower-level data) after insertion.</td></tr>
 * <tr><td>-o</td> <td>--optimize</td> <td>Use VACUUM command to shrink the databse file.</td></tr>
 * <tr><td>-d</td> <td>--dump</td> <td>Export all data in the form of SQL statements.</td></tr>
 * <tr><td>-v</td> <td>--vars</td> <td>Print available variables with short info.</td></tr>
 * <tr><td>-p</td> <td>--params</td> <td>Show what quantum numbers need to be given for a specific variable.</td></tr>
 * <tr><td></td> <td>--Eunits</td> <td>Change energy units (Ry/a.u./eV).</td></tr>
 * <tr><td></td> <td>--Tunits</td> <td>Length units -- units of output (a.u./cgs).</td></tr>
 * <tr><td></td> <td>--Aunits</td> <td>Angular units (deg/rad).</td></tr>
 * <tr><td></td> <td>--&lt;variable&gt;</td> <td>Scattering variable to compute (scatamp, dcs, ics, ccs, tcs, asy, ...).</td></tr>
 * <tr><td></td> <td>--&lt;Q-number&gt;</td> <td>Quantum number specification (ni,li,mi,nf,lf,mf,L,S,Ei,...).</td></tr>
 * </table></center>
 * 
 * <br/>
 * 
 * <center><table>
 * <tr><th>Quantity</th><th>Compulsory quantum numbers</th> <th>STDIN contains*</th></tr>
 * <tr><td>Scattering amplitude</td> <td>ni, li, mi, nf, lf, mf, S, E</td> <td>angles [deg]</td></tr>
 * <tr><td>Differential cross section</td> <td>ni, li, mi, nf, lf, mf, S, E</td> <td>angles [deg]</td></tr>
 * <tr><td>Spin asymmetry</td> <td>ni, li, mi, nf, lf, mf, E</td> <td>angles [deg]</td></tr>
 * <tr><td>Momentum transfer</td> <td>ni, li, mi, nf, lf, mf, L, S</td> <td>energies [Ry]</td></tr>
 * <tr><td>Integral cross section</td> <td>ni, li, mi, nf, lf, mf, L, S</td> <td>energies [Ry]</td></tr>
 * <tr><td>Complete cross section</td> <td>ni, li, mi, nf, lf, mf</td> <td>energies [Ry]</td></tr>
 * <tr><td>Extrapolated cross section</td> <td>ni, li, mi, nf, lf, mf</td> <td>energies [Ry]</td></tr>
 * <tr><td>Collision strength</td> <td>ni, li, mi, nf, lf, mf, L, S</td> <td>energies [Ry]</td></tr>
 * <tr><td>Total cross section</td> <td>ni, li, mi</td> <td>energies [Ry]</td></tr>
 * </table></center>
 * 
 * *) Standard input should contain values of angles or energies separated by white characters (spaces, newlines etc.).
 *     The list ends with EOF (Ctrl+D if manually entering data; inserted automatically when using standard "seq" invocation).
 * 
 * <br/>
 * 
 * For every STDIN entry the program will respond with one number computed (and interpolated if necessary) 
 * from the database. If there are no relevant data in database, the result will be zero. If the available
 * energy interval doesn't contain some of the required energies, appropriate error message will be printed.
 * Here are some examples of usage:
   @verbatim
   # scattering amplitude
   seq 0.01 0.01 3.14    | hex-db --database="hex.db" --scatamp    --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --S=0
  
   # differential cross section
   seq 0.01 0.01 3.14    | hex-db --database="hex.db" --dcs        --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --S=0
  
   # momentum transfer
   seq 0.650 0.001 0.850 | hex-db --database="hex.db" --momtransf  --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --L=0 --S=0
  
   # integral cross section
   seq 0.650 0.001 0.850 | hex-db --database="hex.db" --integcs    --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --L=0 --S=0
  
   # sum integral cross section
   seq 0.650 0.001 0.850 | hex-db --database="hex.db" --sumintegcs --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0
  
   # collision strength
   seq 0.650 0.001 0.850 | hex-db --database="hex.db" --collstr    --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --L=0 --S=0
  
   # total cross section
   seq 0.650 0.001 0.850 | hex-db --database="hex.db" --tcs        --ni=1 --li=0 --mi=0
   @endverbatim
 * You may need to set the LC_NUMERIC environment variable to en_US.utf8 (or en_GB.utf8 or any other
 * localitation that uses decimal point instead of decimal comma). An alternative is the workaround with sed.
   @verbatim
   # scattering amplitude in Czech environment (1)
   export LC_NUMERIC=en_GB.utf8
   seq 0.01 0.01 3.14    | hex-db --database="hex.db" --scatamp    --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --S=0
   
   # scattering amplitude in Czech environment (2)
   seq 0.01 0.01 3.14 | sed -e "s/,/\./g" | hex-db --database="hex.db" --scatamp    --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --S=0
   @endverbatim
 */

#endif /* HEX_DB_DOC */
