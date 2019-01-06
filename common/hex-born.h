//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2017, Jakub Benda, Charles University in Prague                    //
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

#ifndef HEX_BORN
#define HEX_BORN

// --------------------------------------------------------------------------------- //

#include <map>

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"

// --------------------------------------------------------------------------------- //

/**
 * @brief Plane-wave Born approximation.
 * 
 * This function computes all partial wave contributions to the T-matrix for the
 * specified total angular momentum @f( L @f) in the Born approximation of the
 * first order. The contributions depend on all such @f( l_i @f) and @f( lf @f) that
 * the following conditions must be satisfied:
 * @f[
 *    |l_i - L_i| <= L <= l_i + L_i
 *    |l_f - L_f| <= L <= l_f + L_f
 * @f]
 * This results in the following bounds on @f( l_i @f) and @f( l_f @f):
 * @f[
 *    |L_i - L| <= l_i <= L_i + L
 *    |L_f - L| <= l_f <= L_f + L
 * @f]
 * The function will populate the supplied two-dimensional arrays @c Tdir and @c Texc
 * in such a way that the direct and exchange T-matrix for transition
 * @f( (N_i,L_i,M_i) \rightarrow (N_f,L_f,M_f) @f) and out-going partial wave @f( \ell @)
 * is found at position @code Tdir[(Mi + Li)*(2*Lf + 1) + Mf + Lf][lf - abs(Lf-L)] @endcode
 * and @f( Texc[(M_i + L_i)*(2*L_f + 1) + M_f + L_f][\ell - \abs(L_f-L)] @f), respectively.
 * 
 * @param Ni Initial atomic state.
 * @param Li Initial atomic state.
 * @param Nf Final atomic state.
 * @param Lf Final atomic state.
 * @param ki Initial projectile momentum.
 * @param kf Final projectile momentum.
 * @param L Total angular momentum.
 * @param Tdir List of direct scattering T-matrices for all allowed @f( lf @f). Not allocated at entry.
 * @param Texc List of exchange scattering T-matrices for all allowed @f( lf @f). Not allocated at entry.
 * @param direct Whether to compute direct scattering contributions.
 * @param exchange Whether to compute exhange scattering contributions.
 */
void pwba
(
    int Ni, int Li, double ki,
    int Nf, int Lf, double kf,
    int L,
    cArrays & Tdir, cArrays & Texc,
    bool direct = true, bool exchange = true
);

// --------------------------------------------------------------------------------- //

#ifdef WITH_GINAC

    #include <ginac/ginac.h>

    #include "hex-vec3d.h"

    /**
     * @brief Spherical harmonic.
     * 
     * Spherical harmonic @f$ Y_{lm} @f$ in a symbolic form, as a function
     * of the supplied symbols "rp", "rm" and "rz", where
     * @f[
     *     rp \equiv r_+ = r_x + \mathrm{i} r_y \,,
     * @f]
     * @f[
     *     rm \equiv r_- = r_x - \mathrm{i} r_y \,,
     * @f]
     * @f[
     *     rz \equiv r_z \,.
     * @f]
     */ 
    GiNaC::ex rl_Y
    (
        int l, int m,
        GiNaC::symbol const & rp, GiNaC::symbol const & rm,
        GiNaC::realsymbol const & rz
    );

    /**
     * @brief Hydrogen function.
     * 
     * Hydrogen function @f$ \psi_{nlm} @f$ in a symbolic form, as a function
     * of the supplied symbols "r", "rp", "rm" and "rz", where
     * @f[
     *     r = \sqrt{r_x^2 + r_y^2 + r_z^2} \,,
     * @f]
     * @f[
     *     rp \equiv r_+ = r_x + \mathrm{i} r_y \,,
     * @f]
     * @f[
     *     rm \equiv r_- = r_x - \mathrm{i} r_y \,,
     * @f]
     * @f[
     *     rz \equiv r_z \,.
     * @f]
     */
    GiNaC::ex psi_nlm_poly
    (
        int n, int l, int m,
        GiNaC::possymbol const & r,
        GiNaC::symbol const & rp, GiNaC::symbol const & rm,
        GiNaC::realsymbol const & rz
    );

    /**
     * @brief Symbolic integral @f$ W_b @f$
     * 
     * This function symbolically calculates the integral
     * @f[
     *     W_b = \int
     *           \psi_{n_1 l_1 m_1}^\ast(\mathbf{r})
     *           \left(
     *                \mathrm{e}^{e}^{\mathrm{i}\mathbf{k}\cdot\mathbf{r}} - 1
     *           \right)
     *           \psi_{n_2 l_2 m_2}(\mathbf{r}) \,\mathrm{d}^3\mathbf{r} \,.
     * @f]
     * The returned tuples \f( (n_i,a_i,b_i,c_i,d_i) \f) and their coefficients @f( \alpha_i @f)
     * correspond to terms
     * @f[
     *     \alpha_i (\nu^2 + k^2)^{-n_i} \nu^{a_i} k_+^{b_i} k_-^{c_i} k_z^{d_i}
     * @f]
     * of the result. Here @f( k_\pm = k_x \pm \mathrm{i} k_y @f).
     * Use the function @ref eval_Wb to evaluate the symbolic expression for a given
     * @f( \nu @f) and @f( \mathbf{k} = (k_x,k_y,k_z) @f).
     */
    std::map<std::tuple<int,int,int,int,int>,Complex> Wb_symb_in
    (
        int n1, int l1, int m1,
        int n2, int l2, int m2
    );

    /**
     * @brief Evaluator of the @f( W_b @f) integral.
     * 
     * This function serves for evaluations of the result returned
     * from the function @ref Wb_symb_in. It also adds the Fourier factor
     * @f[
     *     \frac{4\pi}{k^2}
     * @f]
     * so that the result is equal to the T-matrix in Born approximation.
     */
    Complex eval_Wb
    (
        std::map<std::tuple<int,int,int,int,int>,Complex> const & poly,
        double nu, geom::vec3d const & k
    );

    std::map<std::tuple<int,int,int,int,int,int,int,int,int,int>,Complex> W_symb_in
    (
        int n, int l, int m
    );

    Complex eval_W
    (
        std::map<std::tuple<int,int,int,int,int,int,int,int,int,int>,Complex> const & poly,
        double nu,
        geom::vec3d const & vk,
        geom::vec3d const & vq
    );

#endif // WITH_GINAC

// --------------------------------------------------------------------------------- //

#endif // HEX_BORN
