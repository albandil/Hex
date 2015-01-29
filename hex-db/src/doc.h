//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2015, Jakub Benda, Charles University in Prague                    //
//                                                                                   //
// MIT License:                                                                      //
//                                                                                   //
//  Permission is hereby granted, free of charge, to any person obtaining a          //
// copy of this software and associated documentation files (the "Software"),        //
// to deal in the Software without restriction, including without limitation         //
// the rights to use, copy, modify, merge, publish, distribute, sublicense,          //
// and/or sell copies of the Software, and to permit persons to whom the             //
// Software is furnished to do so, subject to the following conditions:              //
//                                                                                   //
//  The above copyright notice and this permission notice shall be included          //
// in all copies or substantial portions of the Software.                            //
//                                                                                   //
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS          //
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       //
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE       //
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, //
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF         //
// OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.  //
//                                                                                   //
//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //

#ifndef HEX_DB_DOC
#define HEX_DB_DOC

/**
  @mainpage
  @author Jakub Benda, MFF UK
  @date 17. 5. 2014
  @section db Hex-db
  
  Hex-db is the interface program that can be used to access data precomputed
  by the computational modules. Also, its built-in algorithms can, on demand,
  compute various derived scattering quantities. The full list of available
  variables is given in the section @ref theory.
  
  @subsection install Installation
  
  The whole package is written in C++11, which is supported by most of the
  up-to-date compilers. The following have been tested:
  
  - GCC 4.8.1
  - Intel C++ Composer XE 14.0 SP1 Update 1 (with GCC 4.8.1 headers)
  
  There are several external libraries that are being used for some partial tasks.
  Here is the list of the libraries with versions that were used:
 
  - GSL 1.15<br/>
           Obtainable from http://www.gnu.org/software/gsl/<br/>
           Used for special functions and some others.
  - FFTW 3.3.3<br/>
           Obtainable from http://www.fftw.org/download.html<br/>
           Used for %Chebyshev expansion and Gauss-Chebyshev quadrature.
  - SQLite 3.7.14<br/>
           Obtainable from http://www.sqlite.org/download.html<br/>
           Used for access to the database files.
 
  Some or all of the libraries may be present in precompiled form in the repositories
  of the Linux distribution used.
 
  @subsubsection installLinux Build and run
  
  Once the needed libraries are present in the system include and library paths,
  all that is necessary to build the program is running "make" from this directory.
 
  <pre>
  make
  </pre>
 
  This should create an executable "hex-db" in the subdirectory "bin/". One can
  then create an empty database by the command
 
  <pre>
  bin/hex-db --new
  </pre>
 
  The tables present in the new database can be listed using "sqlite3" program,
 
  <pre>
  sqlite3 hex.db .tables    # prints "ccs ics ionf tmat"
  </pre>
 
  The SQL statements used to create the tables can be retrieved similarly,
 
  <pre>
  sqlite3 hex.db .schema
  </pre>
 
  As soon as there are some data to import, one will run
 
  <pre>
  bin/hex-db --import &lt;SQL_batch_file&gt; --update
  </pre>
 
  If the external libraries are installed in non-standard locations, one needs to
  edit the Makefile accordingly.
  
  @subsubsection installWindows Building on Windows
 
  The program has been successfully tested on Windows 8 (and most probably will work
  exactly the same in older versions). The compilation has been done in the
  Code::Blocks IDE (the project file *.cbp is included with
  the source code). It is necessary to install a full-featured MinGW compiler, though;
  the compiler bundled with Code::Blocks lacks all necessary features.
 
  -# Download MinGW installer from SourceForge:<br/>
    http://sourceforge.net/projects/mingw/files/Installer/mingw-get-setup.exe/download
  -# Install C++ and Fortran compilers, OpenMP and Posix threads implementations.
  -# Download and install latest Code::Blocks (tested version was 13.12).
  -# In Code::Blocks IDE, set the full MinGW installation as the compiler to use.
  
  @subsection theory Theory
  
  The computational modules of the package <b>Hex</b> (e.g. <b>hex-dwba</b>,
  <b>hex-pwba2</b> or <b>hex-ecs</b>) produce output in the form
  of radial part of the T-matrix. The full T-matrix for transition
  @f$ (n_i,l_i,m_i) \rightarrow (n_f,l_f,m_f) @f$ at energy
  @f$ E = k_i^2 - n_i^{-2} = k_f^2 - n_f^{-2} @f$ (energy is in Rydbergs) is
  @f[
         T_{n_i l_i m_i\rightarrow n_f l_f m_f}^{L M S} =
         \delta_{M}^{m_i} \sum_\ell T_{fi,\ell}^{LMS} Y_{\ell}^{m_i-m_f} ,
  @f]
  where the radial part depends on theory. It will be much different for exterior
  complex scaling and for Born approximation. This utility code (and associated
  standalone interface library) can then be used to produce following information:
 
  - the scattering amplitude
    @f[
        f_{i \rightarrow f}^{S} = -\frac{1}{2\pi} T_{i \rightarrow f}^{S} =
        -\frac{1}{2\pi} \sum_{\ell L} T_{fi,\ell}^{Lm_iS} Y_{\ell}^{m_i-m_f} \ ,
    @f]
  - the differential cross section
    @f[
        \frac{\mathrm{d}\sigma_{i \rightarrow f}^{S}}{\mathrm{d}\Omega} = 
        \frac{k_f}{k_i} \frac{2S+1}{4} |f_{i \rightarrow j}^{S}|^2 \ ,
    @f]
  - the partial integral cross section
    @f[
        \sigma_{i \rightarrow f}^{LS} = \frac{k_f}{k_i} \frac{2S + 1}{4}
        \int_{4\pi} |f_{i \rightarrow f}^{LS}|^2 \mathrm{d}\Omega
        = \frac{k_f}{k_i} \frac{2S+1}{16\pi^2} \sum_{\ell L'} T_{fi,\ell}^{Lm_iS}T_{fi,\ell}^{L'm_iS*} \ ,
    @f]
  - the “complete” integral cross section
    @f[
        \sigma_{i \rightarrow f} = \sum_{LS} \sigma_{i \rightarrow f}^{LS} \ ,
    @f]
  - the total cross section
    @f[
        \sigma_i = \sum_{n_f = 1}^\infty \sum_{l_f = 0}^{n_f - 1}
        \sum_{m_f = -l_f}^{l_f} \sigma_{i \rightarrow f} \ ,
    @f]
  - and also the dimensionless, @f$ i \leftrightarrow j @f$ symmetrical, collision strength
    @f[
        \Omega_{i \rightarrow f}^{LS} = k_i^2 \sigma_{i \rightarrow f}^{LS}
    @f]
  - or the momentum transfer
    @f[
        \eta_{i \rightarrow f}^{LS} = \int \frac{\mathrm{d}\sigma_{i \rightarrow f}^{LS}}{\mathrm{d}\Omega}
         (1 - \cos\vartheta) \, \mathrm{d}\Omega \ .
    @f]
  
  Ultimately, all these variables are computed from the partial wave expansion of the T-matrix.
  This series can converge very slow, particularly in the case of large energies. One can use the
  Born subtraction method to overcome this difficulty. Assuming that for some large total angular
  momentum the Born approximation can provide accurate partial T-matrices, it is
  @f[
      T \approx \sum_{L = 0}^{L_0} T_{\ell} Y_\ell^{m_i-m_f} + \sum_{L = L_0 + 1}^\infty T_{\mathrm{Born},\ell}
      Y_\ell^{m_i-m_f} = T_{\mathrm{Born}} + \sum_{L = 0}^{L_0} \left(T_\ell - T_{\mathrm{Born},\ell}\right)
      Y_\ell^{m_i-m_f} \ ,
  @f]
  i.e. the resulting T-matrix is the Born T-matrix corrected by the first @f$ L_0 + 1 @f$ partial waves.
  The angle-dependent full Born T-matrix @f$ T_{\mathrm{Born}} @f$ is stored in the form of %Chebyshev
  expansion. The Born partial T-matrices @f$ T_{\mathrm{Born},\ell} @f$ are stored along the exact
  partial T-matrices @f$ T_\ell @f$.
  
  @subsection database The database
  
  The file "hex.db" contains SQL database in the form of several tables. The most low-level
  data (the T-matrices) are stored in the table "tmat" that has the following columns:
  <center><table>
  <tr><th>Field name</th> <th>Data type</th> <th>Comment</th></tr>
  <tr><td>ni</td> <td>INTEGER</td> <td>Initial principal quantum number.</td></tr>
  <tr><td>li</td> <td>INTEGER</td> <td>Initial orbital quantum number.</td></tr>
  <tr><td>mi</td> <td>INTEGER</td> <td>Initial magnetic quantum number.</td></tr>
  <tr><td>nf</td> <td>INTEGER</td> <td>Final principal quantum number.</td></tr>
  <tr><td>lf</td> <td>INTEGER</td> <td>Final orbital quantum number.</td></tr>
  <tr><td>mf</td> <td>INTEGER</td> <td>Final magnetic quantum number.</td></tr>
  <tr><td>L</td> <td>INTEGER</td> <td>Conserved total angular momentum of a partial wave.</td></tr>
  <tr><td>S</td> <td>INTEGER</td> <td>Conserved total spin of the system.</td></tr>
  <tr><td>Ei</td> <td>DOUBLE PRECISION</td> <td>Projectile impact energy.</td></tr>
  <tr><td>Re_T_ell</td> <td>DOUBLE PRECISION</td> <td>%Real part of @f$ T_\ell @f$.</td></tr>
  <tr><td>Im_T_ell</td> <td>DOUBLE PRECISION</td> <td>Imaginary part of @f$ T_\ell @f$.</td></tr>
  <tr><td>Re_TBorn_ell</td> <td>DOUBLE PRECISION</td> <td>%Real part of @f$ T_{\mathrm{Born},\ell} @f$.</td></tr>
  <tr><td>Im_TBorn_ell</td> <td>DOUBLE PRECISION</td> <td>Imaginary part of @f$ T_{\mathrm{Born},\ell} @f$.</td></tr>
  </table></center>
  However, there are more tables (currently also "bornf" for angle-dependent Born T-matrices,
  "ics" for integral cross section, "ccs" for complete cross section, and "ionf" for radial part
  of the ionization amplitude).
  
  The initial database file can be created by the shell command
  <pre>
  hex-db --new
  </pre>
  which is equivalent to (if "sqlite3" program is available in PATH)
  <pre>
  sqlite3 hex.db 'CREATE TABLE "tmat" (
     ni INTEGER, li INTEGER, mi INTEGER,
     nf INTEGER, lf INTEGER, mf INTEGER,
     L INTEGER, S INTEGER, Ei DOUBLE PRECISION, ell INTEGER,
     Re_T_ell DOUBLE PRECISION, Im_T_ell DOUBLE PRECISION,
     Re_TBorn_ell DOUBLE PRECISION, Im_TBorn_ell DOUBLE PRECISION,
     PRIMARY KEY (ni,li,mi,nf,lf,mf,L,S,Ei,ell)
  )'
  </pre>
  and several similar statements more, for each variable that has to be precomputed.
  Whether a variable has a dedicated table within the database can be inquired
  by looking at its member function SQL_CreateTable. If there are similar commands,
  one has to update the database after data insertion to be able to use the new data:
  <pre>
  hex-db --update
  </pre>
  The insertion of scattering T-matrices, as produced by computational units, is very easy:
  <pre>
  hex-db --import T-matrices.sql
  </pre>
  which is equivalent to
  <pre>
  sqlite3 hex.db < T-matrices.sql
  </pre>
  The raw database file allows standard queriyng, which can be employed to
  retrieve some low level data, like
  <pre>
  sqlite3 -column hex.db 'SELECT Ei,Re_T_ell,Im_T_ell FROM "tmat" WHERE
     ni = 1 AND li = 0 AND mi = 0 AND
     nf = 1 AND lf = 0 AND mf = 0 AND
     Ei > 0.75 AND Ei <= 0.88' > output.txt
  </pre>
  which will print real and imaginary parts of H(1s) @f$ \rightarrow @f$ H(1s)
  T-matrices for energies between @f$ E = 0.75 @f$&nbsp;Ry and @f$ E = 0.88 @f$&nbsp;Ry into
  the file "output.txt".
  
  @subsection usage Usage
  
  Apart from direct program access (when the sources were linked to a custom code),
  one can use a CLI (command line interface) which has the following switches:
  <center><table>
  <tr><th>Short option</th> <th>Long option</th> <th>Usage</th> </tr>
  <tr><td><code>-h</code></td> <td><code>\--help</code></td> <td>Print basic usage info and exit.</td></tr>
  <tr><td><code>-D</code></td> <td><code>\--database</code></td> <td>Use given database file. This is optional; "hex.db" is used if omitted.</td></tr>
  <tr><td><code>-n</code></td> <td><code>\--new</code></td> <td>Create database file, either "hex.db" or the given name if specified by --database.</td></tr>
  <tr><td><code>-i</code></td> <td><code>\--import</code></td> <td>Import SQL data produced by some of the modules.</td></tr>
  <tr><td><code>-u</code></td> <td><code>\--update</code></td> <td>Update cache tables (variables dependent on lower-level data) after insertion.</td></tr>
  <tr><td><code>-o</code></td> <td><code>\--optimize</code></td> <td>Use VACUUM command to shrink the databse file.</td></tr>
  <tr><td><code>-d</code></td> <td><code>\--dump</code></td> <td>Export all data in the form of SQL statements.</td></tr>
  <tr><td><code>-v</code></td> <td><code>\--vars</code></td> <td>Print available variables with short info.</td></tr>
  <tr><td><code>-p</code></td> <td><code>\--params</code></td> <td>Show what quantum numbers need to be given for a specific variable.</td></tr>
  <tr><td></td> <td><code>\--Eunits</code></td> <td>Change energy units (Ry/a.u./eV).</td></tr>
  <tr><td></td> <td><code>\--Tunits</code></td> <td>Length units -- units of output (a.u./cgs).</td></tr>
  <tr><td></td> <td><code>\--Aunits</code></td> <td>Angular units (deg/rad).</td></tr>
  <tr><td></td> <td><code>\--&lt;variable&gt;</code></td> <td>Scattering variable to compute (scatamp, dcs, ics, ccs, tcs, asy, ...).</td></tr>
  <tr><td></td> <td><code>\--&lt;Q-number&gt;</code></td> <td>Quantum number specification (ni,li,mi,nf,lf,mf,L,S,Ei,...).</td></tr>
  </table></center>
  
  <br/>
  
  <center><table>
  <tr><th>Quantity</th><th>Compulsory quantum numbers</th> <th>STDIN contains*</th></tr>
  <tr><td>Scattering amplitude</td> <td>@f$ n_i, l_i, m_i, n_f, l_f, m_f, S, E_i @f$</td> <td>angles [deg]</td></tr>
  <tr><td>Differential cross section</td> <td>@f$ n_i, l_i, m_i, n_f, l_f, m_f, S, E_i @f$</td> <td>angles [deg]</td></tr>
  <tr><td>Spin asymmetry</td> <td>@f$ n_i, l_i, m_i, n_f, l_f, m_f, E_i @f$</td> <td>angles [deg]</td></tr>
  <tr><td>Momentum transfer</td> <td>@f$ n_i, l_i, m_i, n_f, l_f, m_f, L, S @f$</td> <td>energies [Ry]</td></tr>
  <tr><td>Integral cross section</td> <td>@f$ n_i, l_i, m_i, n_f, l_f, m_f, L, S @f$</td> <td>energies [Ry]</td></tr>
  <tr><td>Complete cross section</td> <td>@f$ n_i, l_i, m_i, n_f, l_f, m_f @f$</td> <td>energies [Ry]</td></tr>
  <tr><td>Extrapolated cross section</td> <td>@f$ n_i, l_i, m_i, n_f, l_f, m_f @f$</td> <td>energies [Ry]</td></tr>
  <tr><td>Collision strength</td> <td>@f$ n_i, l_i, m_i, n_f, l_f, m_f, L, S @f$</td> <td>energies [Ry]</td></tr>
  <tr><td>Total cross section</td> <td>@f$ n_i, l_i, m_i @f$</td> <td>energies [Ry]</td></tr>
  </table></center>
  
  \*) Standard input should contain values of angles or energies separated by white characters (spaces, newlines etc.).
      The list ends with EOF (Ctrl+D if manually entering data; inserted automatically when using standard "seq" invocation).
  
  <br/>
  
  For every STDIN entry the program will respond with one number computed (and interpolated if necessary) 
  from the database. If there are no relevant data in database, the result will be zero. If the available
  energy interval doesn't contain some of the required energies, appropriate error message will be printed.
  Here are some examples of usage:
  
  <pre>
  \# scattering amplitude
  seq 0.01 0.01 3.14    | hex-db --database="hex.db" --scatamp    --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --S=0
  
  \# differential cross section
  seq 0.01 0.01 3.14    | hex-db --database="hex.db" --dcs        --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --S=0
  
  \# momentum transfer
  seq 0.650 0.001 0.850 | hex-db --database="hex.db" --momtransf  --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --L=0 --S=0
  
  \# integral cross section
  seq 0.650 0.001 0.850 | hex-db --database="hex.db" --integcs    --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --L=0 --S=0
  
  \# sum integral cross section
  seq 0.650 0.001 0.850 | hex-db --database="hex.db" --sumintegcs --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0
  
  \# collision strength
  seq 0.650 0.001 0.850 | hex-db --database="hex.db" --collstr    --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --L=0 --S=0
  
  \# total cross section
  seq 0.650 0.001 0.850 | hex-db --database="hex.db" --tcs        --ni=1 --li=0 --mi=0
  </pre>
  
  You may need to set the LC_NUMERIC environment variable to en_US.utf8 (or en_GB.utf8 or any other
  localitation that uses decimal point instead of decimal comma). An alternative is the workaround with sed.
  
  <pre>
  \# [Czech locale] avoid commas by changing locale
  export LC_NUMERIC=en_GB.utf8
  seq 0.01 0.01 3.14    | hex-db --database="hex.db" --scatamp    --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --S=0
   
  \# [Czech locale] avoid commas by substitution
  seq 0.01 0.01 3.14 | sed -e "s/,/\./g" | hex-db --database="hex.db" --scatamp    --ni=1 --li=0 --mi=0 --nf=3 --lf=0 --mf=0 --S=0
  </pre>
  
 */

#endif /* HEX_DB_DOC */
