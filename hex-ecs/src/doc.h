//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2016, Jakub Benda, Charles University in Prague                    //
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

#ifndef HEX_ECS_DOC
#define HEX_ECS_DOC

/**
  @mainpage
  @author Jakub Benda, MFF UK, jakub.benda&at;seznam.cz
  @date 10. 3. 2016
  @section ecs Hex-ecs

  <b>Hex-ecs</b> computes partial T-matrices for elastic, excitation and ionization
  electron-hydrogen scattering. Together with the interface program <b>hex-db</b>
  it offers the possibility of generating many scatterign quantities like
  differential or integral cross sections.

  @subsection libs Language and libraries

  Hex is written in C++11 to make use of comfort of the modern C++ extensions,
  so one may need a newer compiler. Tested compilers are:

  - GCC 5.3.0
  - Intel C++ Composer XE 16.0.0

  Both worked with the same Makefile, just by setting the variable CPP to "g++"
  or "icpc". The program also uses following external packages (tested versions
  are given in parentheses):

  - <a href="http://www.gnu.org/software/gsl/">GNU Scientific Library</a> (2.1):
    for Wigner coupling coefficients and some other special functions.
  - <a href="http://www.cise.ufl.edu/research/sparse/SuiteSparse/">SuiteSparse/UMFPACK</a> (4.5.1/5.7.4):
    for sparse matrix manipulation and for a direct sparse system solver.
  - OpenMP (2.1): for parallelization at single machine.
  - MPI (OpenMPI 1.8, MSMPI 7): for parallelization at cluster.

  The next libraries are optional:
  - <a href="http://www.openblas.net/">OpenBLAS</a> (0.2.15):
    Free BLAS implementation that can be compiled for a specific
    CPU. OpenBLAS is able to run in parallel using pthreads or OpenMP.
    OpenBLAS is optional because SuiteSparse can be configured to use
    a different BLAS implementation.

  All listed libraries are open-source and easily obtainable on the internet
  and/or in the repositories of some Linux distributions.

  Equations in this documentation use MathJax, which should work in every
  up-to-date JavaScript-enabled web browser. Tested browsers are:

  - Mozilla Firefox 44.0.2
  - Konqueror 4.14.17

  @subsection usage Usage

  The program can be launched simply by
  <pre>
  hex-ecs
  </pre>
  and then it expects to find the input file "ecs.inp" in the working directory and
  after the computation is done it will produce SQL batch file for use in <b>hex-db</b>. If you run
  <pre>
  hex-ecs --help
  </pre>
  you will get full list of available options with short descriptions. For example, for computation
  of high partial waves a lot of memory is necessary, so one often uses the out-of-core functionality
  (i.e. possibility of storing temporary data on disc), together with parallelization
  <pre>
  mpiexec -n 8 hex-ecs --mpi --input ecs-L30.inp --out-of-core --preconditioner ILU --drop-tolerance 1e-7
  </pre>
  The last option will weaken the ILU preconditioner (entries smaller than 1e-7 will be discarded),
  so that less memory (and disk space) is consumed.

  The input file is expected to be something like
  @verbatim
  # B-spline order.
    4

  # ECS rotation angle in radians.
    0.63

  # B-spline knots.
  # a) Real knots of the basis that is common to atomic and projectile electron.
    L  0.0  0.0   4
    G  0.1 10.0  0.1  1.1
    L   11   60  50
   -1
  # b) Real knots of the panel overlap, if any. (Starts from zero.)
   -1
  # c) Complex region knots (Starts from zero.)
    G    0   50   1  1.02
   -1
  # d) Knots of other panels (propagator projectile basis). (Starts from zero.)
   -1

  # Initial atomic states (ni, li, mi).
    1 -1
    *
    *

  # Final atomic states (nf, lf).
    1  -1
    *

  # L  Pi limit
    0  0  4

  # Atom + projectile total energies in Rydbergs.
    E  -0.35  -0.20  -0.05   -1
   -1

  # Weak magnetic field in atomic units.
    0
  @endverbatim
  The file is not formatted (only the order matters), hash-introduced lines are comments. Initial and final states
  are given in the form of a set of columns of two or three numbers, terminated by -1 on the first line.
  The asterisk symbol stands for all possible values of the specific quantum number. The other (knot and energy)
  sequances are either linear (L), geometric (G) or explicitely listed (E). A linear sequence is specified by first value,
  last value and number of values. The geometric sequence is specified by first value, last value, length of the interval
  between the first and second value and the quotient for the interval expansion. The explicitely listed sequence is given
  by the list of values terminated by -1. B-spline knot specification consists of several sequences. The last one is terminated by -1.

  @subsection theory Theory

  Unknown scattering state, eigenfunction of the full system hamiltonian, is
  split into two parts, asymptotic “incoming particle” state and the scattered
  part, which is a solution of the driven Schrödinger equation,

  <table width = "100%" corder = "0"><tr><td width = "100%" align = "center">
  @f[
      \left(E - H\right) \Psi_{\mathrm{sc}} =
  	    H_{\mathrm{int}} \Psi_{\mathrm{inc}} \ .
  @f]
  </td><td width="0px">(1)</td></tr></table>

  The state @f$ \Psi_{\mathrm{inc}} @f$ is a product of initial atomic state
  and a projectile plane wave,

  <table width = "100%" corder = "0"><tr><td width = "100%" align = "center">
  @f[
      \Psi_{\mathrm{inc}}(\mathbf{r}_1, \mathbf{r}_2) =
       \frac{1}{k_i r_1 r_2} P_{n_i l_i}(r_1) Y_{l_i m_i}(\mathbf{{r}})
       \cdot (2\pi)^{-3/2}
       \sum_{lm} \sqrt{4\pi} \mathrm{i}^l {j}_l(k_i r_2) Y_{lm}(\mathbf{{k}}_i)
       Y_{lm}^\ast(\mathbf{{r}}_2) \ ,
  @f]
  </td><td width="0px">(2)</td></tr></table>

  and anti/symmetrized with respect to electron exchange, according to
  total spin,

  <table width = "100%" corder = "0"><tr><td width = "100%" align = "center">
  @f[
      \Psi^S_{\mathrm{inc}}(\mathbf{r}_1, \mathbf{r}_2) = 
       \Psi_{\mathrm{inc}}(\mathbf{r}_1, \mathbf{r}_2) + (-1)^S
       \Psi_{\mathrm{inc}}(\mathbf{r}_2, \mathbf{r}_1) \ .
  @f]
  </td><td width="0px">(3)</td></tr></table>

  Hamiltonian @f$ {H} @f$ is a sum of free hamiltonian and interaction hamiltonian.
  Free hamiltonian is the hamilton operator for electron and atom separated 
  far away,

  <table width = "100%" corder = "0"><tr><td width = "100%" align = "center">
  @f[
      {H}_{\mathrm{free}} = -\frac{\nabla_1^2}{2} - \frac{\nabla_2^2}{2}
        -\frac{1}{r_1} \ ,
  @f]
  </td><td width="0px">(4)</td></tr></table>

  whereas the interaction hamiltonian contains the interaction between projectile
  and the atomic constituents, that is

  <table width = "100%" corder = "0"><tr><td width = "100%" align = "center">
  @f[
      {H}_{\mathrm{int}} = \frac{1}{r_{12}} - \frac{1}{r_2}
  @f]
  </td><td width="0px">(5a)</td></tr></table>

  for direct case and

  <table width = "100%" corder = "0"><tr><td width = "100%" align = "center">
  @f[
      {H}_{\mathrm{int}} = \frac{1}{r_{12}} - \frac{1}{r_1}
  @f]
  </td><td width="0px">(5b)</td></tr></table>

  for exchange (anti/symmetrized) case. Having the solution of equation (1),
  we can extract the scattering amplitude, which is a matrix element of
  @f$ {H}_{\mathrm{int}} @f$ between the solution
  @f$ \Psi = \Psi_{\mathrm{sc}} + \Psi_{\mathrm{inc}} @f$ and some asymptotic
  state @f$ \Psi_{\mathrm{out}} @f$ that is experimentally inquired.
  Conventionally, @f$ \Psi_{\mathrm{out}} @f$ is also a product of (now final)
  atomic state, possibly
  different from the original, with an outgoing plane wave of the projectile.
  Any exchange effect are already included in the original anti/symmetrization,
  so they need not be considered now. The projection

  <table width = "100%" corder = "0"><tr><td width = "100%" align = "center">
  @f[
      f = -4\pi^2 \left<\Psi_{\mathrm{out}}\right|{H}_{\mathrm{int}}\left|\Psi\right>
  @f]
  </td><td width="0px">(6)</td></tr></table>

  can be also written as

  <table width = "100%" corder = "0"><tr><td width = "100%" align = "center">
  @f[
      f = -4\pi^2 \left<\Psi_{\mathrm{out}}\right|E - {H}_{\mathrm{free}}\left|\Psi_{\mathrm{sc}}\right> \ ,
  @f]
  </td><td width="0px">(7)</td></tr></table>

  where the following was used: @f$ {H}_{\mathrm{int}} = {H} - {H}_{\mathrm{free}} @f$,
  @f$ \Psi = \Psi_{\mathrm{sc}} + \Psi_{\mathrm{inc}} @f$ and @f$ {H}\Psi = E\Psi @f$.


  @subsection method Method

  All computations are done in time-independent way, and 
  the exterior complex scaling is used instead of boundary condition fitting.
  Radial part of sought wave-functions is expanded in a B-spline basis
  @f$ \left\{B_i\right\}_{i = 0}^{\mathrm{Nspline}-1} @f$ of a given order,

  <table width = "100%" corder = "0"><tr><td width = "100%" align = "center">
  @f[
      \Psi_{\mathrm{sc}}^{LM}(\mathbf{r}_1, \mathbf{r}_2) = \sum_{l_1 l_2}
       \psi_{l_1 l_2}^{LM}(\mathrm{r}_1, \mathrm{r}_2)
       \mathcal{Y}_{l_1 l_2}^{LM}(\mathbf{{r}}_1, \mathbf{{r}}_2) \ ,
  @f]
  </td><td width="0px">(8)</td></tr></table>

  <table width = "100%" corder = "0"><tr><td width = "100%" align = "center">
  @f[
      \psi^{LM}_{\mathrm{sc},l_1 l_2}(\mathrm{r}_1, \mathrm{r}_2) = \frac{1}{r_1 r_2}
       \sum_{ij} \psi_{l_1 l_2,ij}^{LM} B_i(r_1) B_j(r_2) \ ,
  @f]
  </td><td width="0px">(9)</td></tr></table>

  and when projecting the equation (1) on a bipolar spherical function
  @f$ \left<\mathcal{Y}_{l_1 l_2}^{LM}\right| @f$ to get rid of angular
  dependence, and on a pair of B-splines to get rid of radial dependence and
  keep only matrix elements, one arrives at a
  matrix equation for components of @f$ \psi_{l_1 l_2}^{LM}(\mathrm{r}_1,\mathrm{r}_2) @f$,

  <table width = "100%" corder = "0"><tr><td width = "100%" align = "center">
  @f[
      \left[\left(ES_{ik}S_{jl} - \frac{1}{2} D_{ik}S_{jl} - \frac{1}{2} S_{ik} D_{jl}
       - \frac{1}{2}l_1(l_1+1) M_{ik}^{(-2)}S_{jl} - \frac{1}{2}l_2(l_2+1) S_{ik}M_{jl}^{(-2)}
       + M_{ik}^{(-1)}S_{jl} + S_{ik}M_{jl}^{(-1)}\right)
       \delta_{l_1}^{l_1'}\delta_{l_2}^{l_2'} - \sum_{\lambda} f_{l_1l_2l_1'l_2';L}^\lambda
       R_{ijkl}^\lambda
       \right] \psi_{l_1' l_2',kl}^{LMS}  = \chi_{l_1l_2,kl}^{LMS} \ ,
  @f]
  </td><td width="0px">(10)</td></tr></table>

  or symbolically

  <table width = "100%" corder = "0"><tr><td width = "100%" align = "center">
  @f[
      \left[ \mathsf{Id}_1 \otimes \mathsf{Id}_2 \otimes
       \left(E\mathsf{S}\otimes\mathsf{S}
       - \frac{1}{2} \mathsf{D}\otimes\mathsf{S}
       - \frac{1}{2} \mathsf{S}\otimes\mathsf{D}
       - \frac{1}{2} l_1 (l_1+1) \mathsf{M}^{(-2)}\otimes\mathsf{S}
       - \frac{1}{2} l_2 (l_2+1) \mathsf{S}\otimes\mathsf{M}^{(-2)}
       + \mathsf{M}^{(-1)} \otimes \mathsf{S}
       + \mathsf{S} \otimes \mathsf{M}^{(-1)}\right)
        + \sum_{\lambda} \mathsf{f}_L^\lambda \otimes \mathsf{R}^\lambda
       \right] \mathsf{\psi}^{LMS}  = \mathsf{\chi}^{LMS} \ .
  @f]
  </td><td width="0px">(10*)</td></tr></table>

  Symbol @f$ \otimes @f$ stands for Kronecker product (“flattened tensor product”)
  and matrices have following meanings:
  - Matrix @f$ \mathsf{Id}_1 @f$ is identity of rank equal to maximal allowed
    angular momentum @f$ l_1 @f$. Analogically for @f$ \mathsf{Id}_2 @f$ .
  - Matrix @f$ \mathsf{D} @f$ is matrix of derivative overlaps of B-splines.
    It can be shown for B-splines basis which is zero at boundaries that in such case
  @f[
      \left<B_i\right| \left(-\frac{\mathrm{d}^2}{\mathrm{d} x^2}\right) \left|B_k\right> =
      + \int_a^b \frac{\mathrm{d}B_i}{\mathrm{d}x} \frac{\mathrm{d}B_k}{\mathrm{d}x}
      \equiv + D_{ik} \ .
  @f]
  - Matrix @f$ \mathsf{S} @f$ is just standard overlap matrix of the B-spline basis,
  @f[
      S_{ik} = \left<B_i\right|\left.\!B_k\right> \ .
  @f]
  - Matrix @f$ \mathsf{M}^{(\alpha)} @f$ is matrix element of power of coordinate
    (also called integral moment in the source code).
  - Matrix @f$ \mathsf{R}^{\lambda} @f$ is matrix of four-B-spline multipole
    integrals for multipole @f$ \lambda @f$,
  @f[
      R_{ijkl}^\lambda = \int_a^b\int_a^b B_i(r_1) B_j(r_2) \frac{r_<^\lambda}{r_>^{\lambda + 1}}
       B_k(r_1) B_l(r_2) dr_1 dr_2 \ ,
  @f]
    flattened so that “i” and “j” form one multi-index [ij] and the other indices the
    multi-index [kl].
  - Matrix @f$ f_L^\lambda @f$ is angular part of reduced matrix element,
  @f[
      \left<l_1 l_2 \right|\!| \frac{1}{r_{12}} |\!\left|l_1' l_2'\right>_L =
       \sum_{\lambda} f_{l_1 l_2 l_1' l_2'; L}^\lambda \frac{r_<^\lambda}{r_>^{\lambda + 1}} \ .
  @f]
    flattened so that “l₁” and “l₂” form one multiindex [l₁l₂] and the other two indices
    the second.

  Finally, the symbol @f$ \chi^{LMS} @f$ stands for projection of the right hand side,
  which is

  <table width = "100%" corder = "0"><tr><td width = "100%" align = "center">
  @f[
      \chi_{l_1 l_2, ij}^{LMS} = 
      \frac{1}{k_f} \sum_\ell \mathrm{i}^\ell \sqrt{2\pi(2\ell+1)} \left\{
       \left(\sum_\lambda f_{l_1 l_2 l_i \ell; L}^\lambda R_{ijkl}^\lambda
       - \delta_{l_1}^{l_i} \delta_{l_2}^{\ell} S_{ik}M_{jl}^{(-1)}\right)
       \left[P_{n_il_i}(r_1)\right]_k \left[{j}_\ell(k_i r_2)\right]_l
       + (-1)^{S+\Pi}
       \left(\sum_\lambda f_{l_1 l_2 \ell l_i; L}^\lambda R_{ijkl}^\lambda
       - \delta_{l_1}^{\ell} \delta_{l_2}^{l_i} M_{ik}^{(-1)} S_{jl}\right)
       \left[{j}_\ell(k_i r_1)\right]_k \left[P_{n_il_i}(r_2)\right]_l
       \right\}
  @f]
  </td><td width="0px">(11)</td></tr></table>

  or symbolically

  <table width = "100%" corder = "0"><tr><td width = "100%" align = "center">
  @f[
      \chi^{LMS} = 
      \frac{1}{k_f} \sum_\ell \mathrm{i}^\ell \sqrt{2\pi(2\ell+1)} C_{l_i m_i \ell 0}^{L M} \left\{
       \left(\sum_\lambda \mathsf{f}_L^\lambda \otimes \mathsf{R}^\lambda
       - \mathsf{Id}_1 \otimes \mathsf{Id}_2 \otimes \mathsf{S}
       \otimes \mathsf{M}^{(-1)}\right) \cdot
       \Delta_{l_1}^{l_i} \otimes \Delta_{l_2}^{\ell} \otimes \mathsf{P}_{n_i l_i} \otimes \mathsf{j}_{\ell,k_i}
       + (-1)^{S+\Pi}
       \left(\sum_\lambda \mathsf{f}_L^\lambda \otimes \mathsf{R}^\lambda
       - \mathsf{Id}_1 \otimes \mathsf{Id}_2 \otimes \mathsf{M}^{(-1)}
       \otimes \mathsf{S}\right) \cdot
       \Delta_{l_1}^{\ell} \otimes \Delta_{l_2}^{l_i} \otimes \mathsf{j}_{\ell,k_i} \otimes \mathsf{P}_{n_i l_i}
       \right\} \ .
  @f]
  </td><td width="0px">(11*)</td></tr></table>

  Here, the one-dimensional (!) vectors @f$ \Delta @f$ are zero vectors with
  only one element equal to one at position @f$ l_1 = l_i @f$ (etc.).
  P- and j- vectors are components of respective function in chosen B-spline
  basis.

  Factor C in the expression (11*) is the Clebsch-Gordan coefficient and the
  zero projection of @f$ \ell @f$-momentum reflect the deliberate choice
  of scattering axis along the projectile moemntum, so that the angular
  momentum projection is zero.

  @subsection restrict ECS restrictions on potential

  Exterior complex scaling of right hand side poses a serious problem for typical
  (not exponentially decreasing) potentials. One of the factors in the right hand
  side is the Riccati-Bessel functions @f$ \hat{j} @f$, which exponentially
  diverges under ECS transformation. To avoid this, potential @f$ {H}_{\mathrm{int}} @f$
  is artificially truncated at (or before) the turning point @f$ R_0 @f$, which
  has several consequences for the numerical construction:
  - In equations (6)-(7), the radial integration is done only up to @f$ r_1,r_2 = R_0 @f$
    and not further.
  - Matrices @f$ \mathsf{S}, \mathsf{M}^{(-1)}, \mathsf{R}^\lambda @f$ in (11), (11*)
    are to be computed, again, for @f$ r_1,r_2 \le R_0 @f$. These are referenced
    as “<i>truncated</i> overlap matrices” in the source code.


  @subsection amplitude Cross section

  As was said above, the scattering amplitude is

  <table width = "100%" corder = "0"><tr><td width = "100%" align = "center">
  @f[
  f = -4\pi^2 \left<\Psi_{\mathrm{out}}\right|E - {H}_{\mathrm{free}}\left|\Psi_{\mathrm{sc}}\right>_{R_0}\ ,
  @f]
  </td><td width="0px">(12)</td></tr></table>

  where the subscript @f$ R_0 @f$ means, that the radial integration is done only
  for radii less than @f$ R_0 @f$ (because the original matrix, @f$ {H}_{\mathrm{int}} @f$,
  has to be truncated at such distance to avoid far-region divergence. The outgoing – detected –
  wavefunction @f$ \Psi_{\mathrm{out}} @f$ has form simillar to (2), with replaced
  initial to final quantum numbers. Substituting such expansion into the equation (12)
  and once again using zero boundary trait of chosen B-spline basis when doing
  per parts integration one easily arrives at the formula for cross section

  <table width = "100%" corder = "0"><tr><td width = "100%" align = "center">
  @f[
      \sigma = \frac{4}{k_i k_f} \sum_{\ell L L'} C_{l_f m_f \ell 0}^{L m_f} C_{l_f m_f \ell 0}^{L' m_f}
       \left| \psi_{l_f \ell, ij}^{LMS} W[P]_i S[j]_j + \psi_{l_f \ell, ij}^{L'MS} S[P]_i W[j]_j \right|^2 \ ,
  @f]
  </td><td width="0px">(13)</td></tr></table>

  where @f$ W[P]_i @f$ and @f$ W[j]_j @f$ stand for wronskian (evaluated at @f$ R_0 - \varepsilon @f$)
  of the i-th B-spline and the (final) hydrogenic or Riccatti-Bessel function, and @f$ S[P]_i @f$
  and @f$ S[j]_i @f$ stand for (truncated) overlap integrals of the i-th B-spline and
  the (final) hydrogenic or Riccatti-Bessel function.

  @subsection code Implementation in the code

  The code runs along the following outline:
  - The program is initialized. Parameters from the command line are stored for later
    use in the class CommandLine.
  - The optional multi-process parallel environment (MPI) is set up. All parallel
    communication is mediated by the class Parallel. The parallelization by MPI
    is used whenever the switch --mpi is present on the command line. The number of
    processes (= "communicator rank") is given as an argument to the MPI launcher. E.g.
     @verbatim
         mpiexec -n 4 hex-ecs --mpi
     @endverbatim
  - Note that the parallelization by OpenMP is used always and does not use any classes,
    nor does it depend on MPI in any way. The number of OpenMP threads is completely independent
    on the communicator size and can be specified by the environment variable OMP_NUM_THREADS.
    The number of threads on every machine used in the MPI network can be different.
  - From here on, all the steps are executed by every machine in the MPI scope.
  - The input file specified on the command line (--input/-i) is found and read. If the name
    was not given on command line, the default "ecs.inp" is used. If the file doesn't exist,
    the program aborts. The values from the input file are contained in the object InputFile.
  - A B-spline basis is created according to the input (class Bspline).
  - A list of available angular states is assembled. The number of these states is restricted
    by @f$ L @f$, @f$ \Pi @f$ and @f$ n_L @f$.
  - A preconditioner is selected in accord with the user's choice. It is stored only in the
    form of a pointer to the base class PreconditionerBase, whose methods are reimplemented
    in all derived classes (specific preconditioners). The constructor is called, which will
    in all cases trigger also the computation of the radial integrals that are needed by
    all preconditioners. The radial integrals are managed by the class RadialIntegrals.
  - For all impact energies:
      - The preconditioner for the set of equations is updated for this energy (method "update").
      - For all initial states and spins:
          - Check that at least some allowed angular states contribute to this combination
            of @f$ L @f$, @f$ \Pi @f$, @f$ n_L @f$ and @f$ l_i @f$. If not, skip this inital state.
          - Check that the solution for this initial state and this impact energy hasn't been
            already computed. If it has been, skip this initial state.
          - The B-spline expansion of the right hand side is computed (method "rhs" of the chosen
            preconditioner class).
          - Run the preconditioned conjugate gradients solver. If there already is a solution for the
            previous energy, use it as an initial guess.
          - Save the solution to disk.
  - For all initial and final states:
      - Compute T-matrices for all available partial waves.
      - Save the T-matrices to SQL files.
      - Evaluate and print the integral cross sections.
 */

#endif
