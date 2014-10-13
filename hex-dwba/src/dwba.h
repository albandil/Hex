//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2014, Jakub Benda, Charles University in Prague                    //
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

#ifndef HEX_DWBA_DWBA1
#define HEX_DWBA_DWBA1

#include "arrays.h"
#include "potential.h"

/**
 * @brief Distorted-wave Born approximation.
 * 
 * This function computes all contributions to the T-matrix for the specified
 * total angular momentum L in the Born approximation of the first-order.
 * The contributions depend on all such "li" and "lf" that
 * the following conditions must be satisfied:
 *    |li - Li| <= L <= li + Li
 *    |lf - Lf| <= L <= lf + Lf
 * This results in the following bounds on "li" and "lf":
 *    |Li - L| <= li <= Li + L
 *    |Lf - L| <= lf <= Lf + L
 * 
 * @param Ni Initial atomic state.
 * @param Li Initial atomic state.
 * @param Nf Final atomic state.
 * @param Lf Final atomic state.
 * @param ki Initial projectile momentum.
 * @param kf Final projectile momentum.
 * @param L Total angular momentum.
 * @param Tdir List of direct scattering T-matrices for all allowed "lf".
 * @param Texc List of exchange scattering T-matrices for all allowed "lf".
 * @param rmax Maximal grid radius.
 * @param direct Whether to compute direct scattering contributions.
 * @param exchange Whether to compute exhange scattering contributions.
 */
void dwba
(
    int Ni, int Li, double ki,
    int Nf, int Lf, double kf,
    int L,
    cArrays & Tdir, cArrays & Texc,
    double rmax,
    bool direct = true, bool exchange = true
);

/**
 * Namespace members compute contributions to the first order of the distorted
 * wave Born approximation.
 */
namespace DWBA1
{
    
/**
 * Computes the one-electron part of direct T-matrix,
 * \f[
 *         D_\ell^{(1)} = \frac{1}{k^2} \sqrt{\frac{2\ell+1}{4\pi}}
 *         \int \chi_\ell(k,r) U(r) (\hat{j}_\ell(kr) - \chi_\ell(k,r))
 *         \mathrm{d}r
 * \f]
 * \param U DistortingPotential class containing information of the distorting potential.
 * \param l Partial wave (outgoing angular momentum).
 * \param k Initial (= final) wavenumber.
 */
Complex computeDirect1e
(
    DistortingPotential const& U,
    int l, double k
);
    
/**
 * Computes the two-electron part of the direct T-matrix,
 * \f[
 *            D_{l_f}^{(2)} = \frac{1}{k_i k_f} \sum_{l_i \lambda}
 *         \mathrm{i}^{l_i - l_f} \frac{\sqrt{4\pi(2l_i+1)}}{2\lambda + 1}
 *         G_{L_i M_i \lambda M_f-M_i}^{L_f M_f}
 *         G_{l_f M_i-M_f \lambda M_f-M_i}^{l_i 0} I_{\mathrm{dir}} \ ,
 * \f]
 * \f[
 *         I_{\mathrm{dir}} = \int\int \chi_{l_f}(k_f,r_2) \psi_{L_f}(E_f,r_1)
 *         \left(\frac{r_<^\lambda}{r_>^{\lambda+1}}-\frac{\delta_\lambda^0}{r_2}\right)
 *         \psi_{L_i}(E_i,r_1) \chi_{l_i}(k_i,r_2) \mathrm{d}r_1 \mathrm{d}r_2 \ .
 * \f]
 * \param U DistortingPotential class containing information of the distorting potential.
 * \param lambda Multipole moment.
 * \param Nf Principal quantum number of the final state.
 * \param Lf Orbital quantum number of the final state.
 * \param kf Final projectile wave number.
 * \param lf Final projectile partial wave.
 * \param Ni Principal quantum number of the initial state.
 * \param Li Orbital quantum number of the initial state.
 * \param ki Initial projectile wave number.
 * \param li Initial projectile partial wave.
 */
Complex computeDirect2e
(
    const DistortingPotential& U, int lambda,
    int Nf, int Lf, double kf, int lf,
    int Ni, int Li, double ki, int li
);
    
/**
 * Computes the one-electron part of exchange T-matrix,
 * \f[
 *         E^{(1)} = -\frac{1}{k_i k_f} \mathrm{i}^{L_f-L_i}
 *         \sqrt{\frac{2L_f+1}{4\pi}} \int \chi_{L_i}(k_f,r) U_f(r)
 *         \psi_{L_i}(E_i,r) \mathrm{d}r \int \psi_{L_f}(E_f,r)
 *         \chi_{L_f}(k_i,r) \mathrm{d}r
 * \f]
 * \param U DistortingPotential class containing information of the distorting potential.
 * \param Ni Initial atomic principal quantum number.
 * \param Li Initial atomic orbital quantum number.
 * \param ki Initial projectile wavenumber.
 * \param Nf Final atomic principal quantum number.
 * \param Lf Final atomic orbital quantum number.
 * \param kf Final projectile wavenumber.
 */
Complex computeExchange1e
(
    DistortingPotential const& U,
    int Ni, int Li, double ki,
    int Nf, int Lf, double kf
);

/**
 * Computes the two-electron part of the exchange T-matrix,
 * \f[
 *            E_{l_f}^{(2)} = \frac{1}{k_i k_f} \sum_{l_i \lambda}
 *         \mathrm{i}^{l_i - l_f} \frac{\sqrt{4\pi(2l_i+1)}}{2\lambda + 1}
 *         G_{l_i 0 \lambda M_f}^{L_f M_f}
 *         G_{l_f M_i-M_f \lambda M_f}^{L_i M_i} I_{\mathrm{exc}} \ ,
 * \f]
 * \f[
 *         I_{\mathrm{dir}} = \int\int \chi_{l_f}(k_f,r_2) \psi_{L_f}(E_f,r_1)
 *         \left(\frac{r_<^\lambda}{r_>^{\lambda+1}}-\frac{\delta_\lambda^0}{r_2}\right)
 *         \psi_{L_i}(E_i,r_2) \chi_{l_i}(k_i,r_1) \mathrm{d}r_1 \mathrm{d}r_2 \ .
 * \f]
 * \param U DistortingPotential class containing information of the distorting potential.
 * \param lambda Multipole moment.
 * \param Nf Principal quantum number of the final state.
 * \param Lf Orbital quantum number of the final state.
 * \param kf Final projectile wave number.
 * \param lf Final projectile partial wave.
 * \param Ni Principal quantum number of the initial state.
 * \param Li Orbital quantum number of the initial state.
 * \param ki Initial projectile wave number.
 * \param li Initial projectile partial wave.
 */
Complex computeExchange2e
(
    const DistortingPotential& U, int lambda,
    int Nf, int Lf, double kf, int lf,
    int Ni, int Li, double ki, int li
);

} // end of namespace DWBA1

#endif
