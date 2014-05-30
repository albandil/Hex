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

#ifndef HEX_DWBA_DOC
#define HEX_DWBA_DOC

/**
 * @mainpage
 * @author Jakub Benda, MFF UK
 * @date 17. 5. 2014
 * @section dwba Hex-dwba
 * 
 * Hex-dwba is a high-energy computational module of the Hex package. It
 * implements the distorted wave Born approximation of the first order and
 * is expected to give trustworthy results for ehergies above 1 keV. However,
 * it can be also run in the plane wave mode (if the option <code>\--nodistort</code>
 * is used). The full usage scheme is
 * <pre>
 * hex-dwba [--nodistort] [--nodirect] [--noexchange] &lt;ni&gt; &lt;li&gt; &lt;nf&gt; &lt;lf&gt; &lt;Ei&gt; &lt;L&gt; [&lt;rmax&gt;]
 * </pre>
 * where the options allow to skip some contributions (direct/exchange) and the number
 * parameters are: the initial and final atomic quantum numbers, impact energy,
 * total angular momentum and (optional) radial grid length used in distorted wave
 * method. The result is an SQL batch file for use in <b>hex-db</b>.
 * 
 * @subsection pwba Plane wave approximation
 * 
 * Plane wave Born approximation of the first order (PWBA-1) is a perturbation
 * approach to computation of electron-atom scattering variables. The most
 * important for the whole Hex package is the scattering T-matrix, which in
 * PWBA-1 reads
 * @f[
 *         T = \left<\Psi_f\right|V\left|\Psi_i\right> \ ,
 * @f]
 * where both the final and initial wave function of the system are approximated
 * by a product of atomic orbital and projectile plane wave,
 * @f[
 *         \Psi_i(\mathbf{r}_1,\mathbf{r}_2) = \frac{1}{r_1}P_{N_iL_i}(r_1)
 *         Y_{L_iM_i}(\hat{\mathbf{r}}_1) \cdot \frac{4\pi}{k_ir_2} \sum_{l_im_i} \mathrm{i}^{l_i}
 *         \hat{j}_{l_i}(k_ir_2) Y_{l_im_i}(\hat{\mathbf{r}}_2) Y_{l_im_i}^*
 *         (\hat{\mathbf{k}}_i) \ .
 * @f]
 * The potential is the sum of the proton and atomic electron potential,
 * @f[
 *         V(\mathbf{r}_1,\mathbf{r}_2) = \frac{1}{|\mathbf{r}_1 - \mathbf{r}_2|}
 *         - \frac{1}{\mathbf{r}_2} \ .
 * @f]
 * The exchange amplitude can be computed as well if the coordinates in \f$ \Psi_i \f$
 * are exchanged and also the coordinates in \f$ V \f$.
 * 
 * After some simplifications, the resulting direct and exchange T-matrices are
 * @f[
 *         T_{\mathrm{dir},\ell} = \frac{8}{k_ik_f} \sum_{l \lambda}
 *         \frac{\mathrm{i}^{l-\ell}}{2\lambda+1}
 *         \sqrt{\frac{2l+1}{4\pi}}
 *         G_{L_i,M_i,\lambda,M_f-M_i}^{L_f,M_f}
 *         G_{\lambda,M_f-M_i,\ell,M_i-M_f}^{l,0}
 *         I_{\mathrm{dir},\ell\lambda l}^{N_i,L_i,N_f,L_f}(k_f,k_i) \ ,
 * @f]
 * @f[
 *         T_{\mathrm{exc},\ell} = \frac{8}{k_ik_f} \sum_{l \lambda}
 *         \frac{\mathrm{i}^{l-\ell}}{2\lambda+1}
 *         \sqrt{\frac{2l+1}{4\pi}}
 *         G_{l,0,\lambda,M_f}^{L_f,M_f}
 *         G_{L_i,M_i,\ell,M_i-M_f}^{L_i,M_i}
 *         I_{\mathrm{exc},\ell\lambda l}^{N_i,L_i,N_f,L_f}(k_f,k_i) \ ,
 * @f]
 * with the direct and exchange integral
 * @f[
 *         I_{\mathrm{dir},\ell\lambda l}^{N_i,L_i,N_f,L_f}(k_f,k_i) =
 *         \int_0^\infty \int_0^\infty \hat{j}_{\ell}(k_fr_2) P_{N_fL_f}(r_1)
 * 		   \left(\frac{r_<^\lambda}{r_>^{\lambda+1}} - \frac{\delta_\lambda^0}{r_2}\right)
 *         P_{N_iL_i}(r_1) \hat{j}_{l}(k_ir_2) \mathrm{d}r_1 \mathrm{d}r_2 \ ,
 * @f]
 * @f[
 *         I_{\mathrm{exc},\ell\lambda l}^{N_i,L_i,N_f,L_f}(k_f,k_i) =
 *         \int_0^\infty \int_0^\infty \hat{j}_{\ell}(k_fr_2) P_{N_fL_f}(r_1)
 * 		   \left(\frac{r_<^\lambda}{r_>^{\lambda+1}} - \frac{\delta_\lambda^0}{r_1}\right)
 *         P_{N_iL_i}(r_2) \hat{j}_{l}(k_ir_1) \mathrm{d}r_1 \mathrm{d}r_2 \ .
 * @f]
 * Here, \f$ G \f$ is the Gaunt's integral, \f$ \hat{j}_\ell \f$ the Riccati-Bessel
 * function of the \f$ \ell \f$-th order.
 * 
 * @subsubsection pwbaMethod Method
 * 
 * The twodimensional integrals could be in principle directly numerically
 * evaluated, but a fasted alternative is implemented here, which is semi-analytical.
 * All integrands can be expressed as a product of a multiplicative constant,
 * power of the coordinate, trigonometric function of some multiple of the coordinate
 * and damped exponential,
 * @f[
 *         F(r) = \sum_i k_i x^{a_i} \sin (b_ir)\, \mathrm{e}^{-c_ir}
 * @f]
 * or
 * @f[
 *         F(r) = \sum_i k_i x^{a_i} \cos (b_ir)\, \mathrm{e}^{-c_ir}
 * @f]
 * Such a function can be easily integrated. It is
 * @f[
 *         \int_0^\infty x^a \sin (bx)\, \mathrm{e}^{-cx} \mathrm{d}x =
 *         -\mathrm{Im}\,\frac{a!}{(c-\mathrm{i}b)^{a+1}} \ ,
 * @f]
 * @f[
 *         \int_r^\infty x^a \sin (bx)\, \mathrm{e}^{-cx} \mathrm{d}x =
 *         -\mathrm{Im}\,\left\{
 *            \frac{a!}{(c-\mathrm{i}b)^{a+1}} \mathrm{e}^{-(c-\mathrm{i}b)x}
 *            \sum_{k=0}^a \frac{[(c-\mathrm{i}b)x]^k}{k!}
 *         \right\} 
 * @f]
 * and analogously for cosine integral (will contain real parts instead of imaginary
 * parts). Evaluating these analytical expressions is and order faster than
 * evaluating the inner integral.
 * 
 * Necessary input are the formulas for associated Laguerre polynomials,
 * @f[
 *         L_k^s(x) = \sum_{j=s}^k \frac{ (-1)^j (k!)^2 x^{j-s} }{ (k-j)! j! (j-s)! } \ ,
 * @f]
 * hydrogen radial functions
 * @f[
 *         P_{nl}(r) = r \sqrt{\left(\frac{2}{n}\right)^3 \frac{(n-l-1)!}{2n((n+l)!)^3}}
 *           \left(\frac{2r}{n}\right)^l L_{n+l}^{2l+1}(2r/n) \mathrm{e}^{-r/n}
 * @f]
 * and Riccati-Bessel function
 * @f[
 *         \hat{j}_l(x) = \sum_{i=0}^{\lfloor\frac{l+1}{2}\rfloor} \left\{
 *                       \left(\frac{2}{x}\right)^{l-2i} \frac{(l-i)!}{i!} {-1/2-i \choose l-2i} \sin x
 *                     - \left(\frac{2}{x}\right)^{l-2i} \frac{(l-i)!}{i!} i {-1/2-i \choose l-2i+1 } \cos x \right\}\ .
 * @f]
 * 
 * @subsection dwbaTheory Distorted wave mode
 * 
 * The distorted wave Born approximation of the first order (DWBA-1) is
 * similar to the well-known plane wave Born approximation. However, to speed up
 * the convergence of the Born series (of which we are using just a single term
 * in the first-order computation), part of the potential is excluded from the
 * perurbative description and solved precisely. The spliting is the following
 * <table width = "100%"><tr><td width = "100%">
 * @f[
 *     V = W + U \,
 * @f]
 * </td><td width = "0">(1)</td></tr></table>
 * where @f$ V @f$ is the full potential felt by the projectile, @f$ U @f$ is
 * the so called "distorting potential" that will be solved precisely as
 * mentioned above and explained below. This splitting has two consequences:
 * - First of all, the free states are no longer free, so we cannot use
 *   the plane waves. Instead, one receives a distorted waves as a solution
 *   of the modified equation for partial waves
 *   <table width = "100%"><tr><td width = "100%">
 *   @f[
 *        \left(
 *            -\frac{1}{2}\frac{\mathrm{d}^2}{\mathrm{d}r^2}
 *            +\frac{\ell(\ell+1)}{2r^2}
 *            +U(r)
 *        \right)\chi_{k\ell}(r) = E_{k}\chi_{k\ell}(r) \ .
 *   @f]
 *   </td><td width = "0">(2)</td></tr></table>
 * - Secondly, the T-matrix receives two contributions now -- one from the
 *   scattering term and one as a correction due to the redefinition of the
 *   potentials. The result is the two-potential formula
 *   <table width = "100%"><tr><td width = "100%">
 *   @f[
 *        T = (N+1) \left<\chi_f^{(-)}\psi_f\right|V-U_f\left|\mathcal{A}\Psi_i^{(+)}\right>
 *      + \left<\chi_f^{(-)}\psi_f\right|U_f\left|\psi_i\beta_i\right> \ .
 *   @f]
 *   </td><td width = "0">(3)</td></tr></table>
 * 
 * @subsubsection dwbaImplementation Implementation
 * 
 * The program iterates over partial waves of the asymptotic final state.
 * For every iterations is will loop over all initial states that can
 * contribute to the requested scattering process and still result in the
 * current asymptotic final state.
 * 
 * The radial functions @f$ \chi_{k\ell}(r) @f$ are computed using numerical solution
 * of the equation (2) by the adaptive Cash-Karp Runge-Kutta method provided
 * by GSL (see below).
 * 
 * @subsection requirements Requirements
 * 
 * You will need a modern C++ compiler to compile the program as for example
 * 
 * - GCC 4.8.1 or 
 * - Intel C++ Composer XE 14.0 SP1 Update 1 (tested with GCC 4.8.1 headers)
 * 
 * DWBA uses some free external libraries
 * 
 * - <a href="http://www.gnu.org/software/gsl/">GNU Scientific Library</a> (1.16):
 *   for Wigner coupling coefficients, numerical integration and some other special functions.
 * - <a href="http://www.hdfgroup.org/HDF5/">HDF5</a> (1.8.11):
 *   for standardized binary output of intermediate results
 *   (used to speed up further simillar calculations); optional.
 * - <a href="http://www.fftw.org/">FFTW</a> (3.3.3):
 *   for fast evaluation of %Chebyshev expansions by Fast Fourier Transform.
 * - CLN (<a href="http://www.ginac.de/CLN/">http://www.ginac.de/CLN/</a>): for arbitrary-precision
 *   rational arithmetic used in the plane wave method.
 *
 * The libraries are mostly present in the repositories of all major Linux distributions. They
 * can be also compiled from source code in Windows -- this has been tested in the Code::Blocks
 * IDE. However, the compilation on Windows is not trivial and needed some manual aid.
 */

#endif /* HEX_DWBA_DOC */
