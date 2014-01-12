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

#ifndef HEX_PWBA_DOC
#define HEX_PWBA_DOC

/**
 * @mainpage
 * 
 * @author Jakub Benda, MFF UK
 * @date 12. 1. 2014
 * 
 * @section Theory Theory
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
 * @section Method Method
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
 * @section Usage Usage
 * 
 * The program can be called from the command line as
 * \verbatim # ./pwba <ni> <li> <nf> <lf> <maxL> <Ei> <eps> \endverbatim
 * Magnetic quantum numbers are not specified, program calculates all allowed
 * magnetic transitions at once (just one radial integral is needed).
 * Partial wave count can be restricted by maxL (or set the variable to -1 if you
 * intend no restriction). Otherwise the program stops when the contribution
 * of last partial wave to the direct cross section has been smaller than eps-part
 * of current sum.
 * 
 * @section Requirements
 * 
 * You will need a modern C++ compiler to compile the program as for example
 * 
 * - GCC 4.8.1 or 
 * - Intel C++ compiler 14.0
 * 
 * PWBA uses some free external libraries
 * 
 * - GSL (<a href="http://www.gnu.org/software/gsl/">http://www.gnu.org/software/gsl/</a>) and 
 * - CLN (<a href="???">???</a>)
 */

#endif
