//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2026, Jakub Benda, Charles University in Prague                    //
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

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>

#include "hex-arrays.h"
#include "hex-gausskronrod.h"
#include "hex-hydrogen.h"
#include "hex-misc.h"
#include "hex-special.h"

// --------------------------------------------------------------------------------- //

//
// Total plane-wave Born cross section (no exchange) for the transition
//
//     e(k_i) + H(1s)  ->  e(k_f) + H(nl) ,
//
// summed over the final magnetic sublevels "m". Both the excitations and the elastic
// scattering (n = 1, l = 0) are covered. Rydberg units are used throughout, as
// everywhere else in Hex: energies in Ry, momenta in a.u. (E = k^2), lengths in a_0
// and the resulting cross section in a_0^2.
//
// The direct first-order amplitude for the interaction  V = 1/|R - r| - 1/R  is
//
//     f(q) = -(2/q^2) [ <nlm| exp(i q.r) |1s> - <nlm|1s> ] ,     q = k_i - k_f ,
//
// the second term being the screening of the projectile-electron repulsion by the
// attraction to the nucleus. Expanding the plane wave into partial waves,
//
//     exp(i q.r) = 4 pi sum_{L M} i^L j_L(qr) Y*_{LM}(q) Y_{LM}(r) ,
//
// and using the sphericity of the initial state, only the single multipole L = l
// survives the angular integration. Because j_0(0) = 1, the overlap is precisely the
// q = 0 value of that multipole, so both terms can be kept together,
//
//     <nlm| exp(i q.r) |1s> - <nlm|1s> = sqrt(4 pi) i^l Y*_{lm}(q) M_l(q) ,
//
// with the screened multipole Born form factor
//
//     M_l(q) = int_0^infty P_nl(r) [ j_l(qr) - delta_{l0} ] P_1s(r) dr .    (*)
//
// The overlap is nonzero only for the elastic channel, where it is the contribution
// of the nucleus, without which the cross section would diverge; every excitation has
// <ns|1s> = 0 and is left untouched.
//
// The sum over "m" is then trivial, sum_m |Y_lm|^2 = (2l+1)/(4 pi), so that
//
//     dsigma/dOmega = (k_f/k_i) (4/q^4) (2l+1) M_l(q)^2 .
//
// Finally, dOmega = 2 pi q dq / (k_i k_f) turns the angular integration into an
// integration over the momentum transfer,
//
//     sigma = 8 pi (2l+1) / k_i^2  int_{|k_i - k_f|}^{k_i + k_f} M_l(q)^2 dq / q^3 .
//
// For the elastic channel the lower bound is q = 0, where the integrand vanishes as
// q^1 because M_0(q) = O(q^2) there.
//

// --------------------------------------------------------------------------------- //

/**
 * @brief Coefficients of the closed-form moment of the spherical Bessel function.
 * 
 * The integral
 * @f[
 *     \int_0^\infty r^m \mathrm{e}^{-cr} j_l(qr) \mathrm{d}r \,, \qquad m \ge l + 1 \,,
 * @f]
 * is known in closed form. The lowest power has the elementary value
 * @f[
 *     \int_0^\infty r^{l+1} \mathrm{e}^{-cr} j_l(qr) \mathrm{d}r
 *     = \frac{2^l l! q^l}{(c^2 + q^2)^{l+1}} \,,
 * @f]
 * and every further power of @f$ r @f$ follows from
 * @f$ r^{m+1} \mathrm{e}^{-cr} = -\partial_c\, r^m \mathrm{e}^{-cr} @f$. The set of
 * functions
 * @f[
 *     \sum_i a_i \frac{c^{e_i}}{(c^2 + q^2)^{l+1+i}}
 * @f]
 * is closed under that differentiation, and the exponent @f$ e_i = 2i - D @f$ is fixed
 * by the index and by the number @f$ D = m - l - 1 @f$ of derivatives taken, so the
 * whole expansion is carried in the single array of coefficients returned here.
 * 
 * @param m Power of the radial variable.
 * @param l Order of the spherical Bessel function.
 */
std::vector<double> bessel_moment_coefficients (int m, int l)
{
    // number of derivatives with respect to "c" that separate "m" from the base case
    int D = m - l - 1;

    if (D < 0)
        HexException("The moment r^%d of j_%d is not implemented (needs m >= l + 1).", m, l);

    std::vector<double> a (D + 1, 0.);
    a[0] = 1.;

    for (int d = 0; d < D; d++)
    {
        std::vector<double> b (D + 1, 0.);

        for (int i = 0; i <= d; i++)
        {
            //  -d/dc [ a c^e (c^2+q^2)^-p ] = -a e c^(e-1) (c^2+q^2)^-p + 2 a p c^(e+1) (c^2+q^2)^(-p-1)
            int p = l + 1 + i, e = 2*i - d;

            if (e > 0)
                b[i] -= a[i] * e;

            b[i + 1] += a[i] * 2. * p;
        }

        a = std::move(b);
    }

    return a;
}

/**
 * @brief Moment of the spherical Bessel function.
 * 
 * Evaluates @f$ \int_0^\infty r^m \mathrm{e}^{-cr} j_l(qr) \mathrm{d}r @f$ by summing
 * the expansion set up in @ref bessel_moment_coefficients.
 * 
 * @param m Power of the radial variable.
 * @param l Order of the spherical Bessel function.
 * @param c Exponential decay constant.
 * @param q Momentum transfer.
 */
double bessel_moment (int m, int l, double c, double q)
{
    std::vector<double> a = bessel_moment_coefficients(m, l);

    int D = m - l - 1;
    double u = c * c + q * q, sum = 0;

    for (int i = 0; i <= D; i++)
        sum += a[i] * gsl_sf_pow_int(c, 2*i - D) / gsl_sf_pow_int(u, l + 1 + i);

    return gsl_sf_pow_int(2.,l) * gsl_sf_fact(l) * gsl_sf_pow_int(q,l) * sum;
}

/**
 * @brief Screened moment of the spherical Bessel function.
 * 
 * Evaluates the monopole moment with its @f$ q = 0 @f$ value removed,
 * @f[
 *     \int_0^\infty r^m \mathrm{e}^{-cr} \left[ j_0(qr) - 1 \right] \mathrm{d}r \,.
 * @f]
 * Forming the difference term by term would be catastrophic at small @f$ q @f$, where
 * the two moments agree to many digits; that matters, because @f$ q = 0 @f$ is the very
 * end of the integration range of the elastic channel. Each term is therefore
 * subtracted analytically using
 * @f[
 *     \frac{1}{u^p} - \frac{1}{v^p} = -q^2 \sum_{s=0}^{p-1} \frac{1}{u^{p-s} v^{s+1}} \,,
 *     \qquad u = c^2 + q^2 \,, \quad v = c^2 \,,
 * @f]
 * which is exact and carries the factor @f$ q^2 @f$ explicitly, so that no cancellation
 * takes place at all.
 * 
 * @param m Power of the radial variable.
 * @param c Exponential decay constant.
 * @param q Momentum transfer.
 */
double bessel_moment_screened (int m, double c, double q)
{
    std::vector<double> a = bessel_moment_coefficients(m, 0);

    int D = m - 1;
    double u = c * c + q * q, v = c * c, sum = 0;

    for (int i = 0; i <= D; i++)
    {
        if (a[i] == 0.)
            continue;

        int p = i + 1;

        double geom = 0;
        for (int s = 0; s < p; s++)
            geom += 1. / (gsl_sf_pow_int(u, p - s) * gsl_sf_pow_int(v, s + 1));

        sum += a[i] * gsl_sf_pow_int(c, 2*i - D) * geom;
    }

    return -q * q * sum;
}

/**
 * @brief Screened multipole Born form factor for the 1s -> nl transition.
 * 
 * Evaluates the radial integral (*) above in closed form. The bound radial functions
 * @f[
 *     P_{nl}(r) = \sqrt{\left(\frac{2}{n}\right)^3 \frac{(n-l-1)!}{2n(n+l)!}}
 *                 \left(\frac{2r}{n}\right)^l \mathrm{e}^{-r/n}
 *                 L_{n-l-1}^{2l+1}\!\left(\frac{2r}{n}\right) r \,,
 *     \qquad P_{1s}(r) = 2 r \mathrm{e}^{-r} \,,
 * @f]
 * multiply to a polynomial times @f$ \exp(-cr) @f$ with @f$ c = 1 + 1/n @f$, so the
 * integral reduces to a combination of the moments evaluated above.
 * 
 * The monopole is taken in whichever of its two equivalent forms is better
 * conditioned, see the comment in the loop below. Neither the elastic channel, where
 * the screening by the nucleus is physically required, nor the excitations, where the
 * bare monopole would cancel to zero as @f$ q \to 0 @f$ by the orthogonality of the two
 * bound states, then suffer any loss of accuracy.
 * 
 * @param n Principal quantum number of the final state.
 * @param l Orbital quantum number of the final state.
 * @param q Momentum transfer.
 */
double formfactor (int n, int l, double q)
{
    if (n <= l)
        HexException("Invalid final state (%d,%d).", n, l);

    // normalization of the final bound state
    double norm = std::sqrt(gsl_sf_pow_int(2./n,3) * gsl_sf_fact(n-l-1) / (2. * n * gsl_sf_fact(n+l)));

    // decay constant of the product of the two bound states
    double c = 1. + 1./n;

    // sum over the terms of the generalized Laguerre polynomial
    double sum = 0;
    for (int j = 0; j <= n - l - 1; j++)
    {
        // coefficient of x^j in L_{n-l-1}^{2l+1}(x), rescaled from "x" to "2r/n"
        double cj = (j % 2 == 0 ? +1. : -1.) * gsl_sf_choose(n+l, n-l-1-j) / gsl_sf_fact(j);

        // The monopole is evaluated with the overlap <nl|1s> subtracted whenever that
        // is the better conditioned of the two forms. For the elastic channel the
        // subtraction is mandatory anyway. For an excitation the overlap vanishes, so
        // both forms are the same integral, and the subtracted one is used only below
        // q ~ c, where the bare terms would cancel against each other; above that the
        // roles are reversed and the bare form is the accurate one.
        double moment = (l == 0 and (n == 1 or q < c) ? bessel_moment_screened(j + 2, c, q)
                                                      : bessel_moment(j + l + 2, l, c, q));

        sum += cj * std::pow(2./n, j) * moment;
    }

    // the remaining factors are the "2" of P_1s and the angular factor of P_nl
    return 2. * norm * std::pow(2./n, l) * sum;
}

/**
 * @brief Run the Gauss-Kronrod quadrature and check that it converged.
 */
template <class Functor> double integrate (Functor f, double a, double b)
{
    GaussKronrod<Functor> Q (f);
    Q.integrate(a, b);

    if (not Q.ok())
        HexException("Failed to integrate over the momentum transfer - %s.", Q.status().c_str());

    return Q.result();
}

/**
 * @brief Total Born cross section for the 1s -> nl transition.
 * 
 * Integrates the form factor over the momentum transfer as sketched above.
 * 
 * For the excitations the interval @f$ [\,k_i - k_f, k_i + k_f\,] @f$ spans many decades
 * at high impact energies, so the substitution @f$ q = \mathrm{e}^t @f$ is used; the
 * integrand is smooth there and, for the dipole transitions, nearly constant. The
 * elastic channel has @f$ k_i = k_f @f$ and starts right at @f$ q = 0 @f$, where there
 * is no dynamic range to compress: the plain variable is used instead, in which the
 * integrand rises linearly from the origin. (Taking the logarithm anyway would leave a
 * semi-infinite interval carrying an integrand that is flat and negligible over almost
 * all of it, which exhausts the subdivision limit of the quadrature.)
 * 
 * @param n Principal quantum number of the final state.
 * @param l Orbital quantum number of the final state.
 * @param Ei Impact energy (Ry).
 */
double cross_section (int n, int l, double Ei)
{
    // final energy from the conservation law
    double Ef = Ei - 1. + 1./(n*n);

    // nothing to do below the excitation threshold (never the case for n = 1)
    if (Ef <= 0)
        return 0;

    double ki = std::sqrt(Ei), kf = std::sqrt(Ef);
    double qmin = std::abs(ki - kf), qmax = ki + kf;

    // integrand of the momentum-transfer integral
    auto dsigma = [n,l](double q) -> double
    {
        double M = formfactor(n, l, q);

        return M * M / (q * q * q);
    };

    // integrate, in the logarithmic variable whenever the interval avoids the origin
    double integral;
    if (qmin > 0)
    {
        auto logarithmic = [dsigma](double t) -> double { return dsigma(std::exp(t)) * std::exp(t); };

        integral = integrate(logarithmic, std::log(qmin), std::log(qmax));
    }
    else
    {
        integral = integrate(dsigma, 0., qmax);
    }

    return 8. * special::constant::pi * (2*l + 1) / (ki * ki) * integral;
}

// --------------------------------------------------------------------------------- //

const std::string usage =
    "\nUsage:\n\n"
    "  hex-fullborn [--form-factor] <n> <l> [<x> ...]\n\n"
    "Writes the total plane-wave Born cross section (no exchange) of the 1s -> nl\n"
    "transition in hydrogen, summed over the final magnetic sublevels, for every\n"
    "impact energy <x> given. The elastic channel is obtained with n = 1, l = 0.\n"
    "With \"--form-factor\" the screened multipole Born form factor M_l(q) is written\n"
    "instead, for every momentum transfer <x> given. When no <x> is given on the\n"
    "command line, the values are read from the standard input.\n\n"
    "Energies are in Rydbergs, momenta in atomic units, cross sections in a_0^2.\n\n";

int main (int argc, char * argv[])
{
    // do not let GSL abort on underflows; they are legitimate here
    gsl_set_error_handler_off();

    bool ff = false;
    std::vector<const char*> args;

    for (int iarg = 1; iarg < argc; iarg++)
    {
        if (std::strcmp(argv[iarg],"--form-factor") == 0)
            ff = true;
        else
            args.push_back(argv[iarg]);
    }

    if (args.size() < 2)
    {
        std::cout << usage;
        return EXIT_FAILURE;
    }

    int n = std::atoi(args[0]);
    int l = std::atoi(args[1]);

    if (n < 1 or l < 0 or l >= n)
    {
        std::cerr << "Invalid final state (" << n << "," << l << ")." << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << std::scientific << std::setprecision(10);

    std::cout << "# Plane-wave Born (no exchange) cross section, 1s -> " << Hydrogen::stateName(n,l) << std::endl;

    if (n == 1)
        std::cout << "# elastic channel, screened by the attraction to the nucleus" << std::endl;
    else
        std::cout << "# excitation threshold: " << 1. - 1./(n*n) << " Ry" << std::endl;
    std::cout << (ff ? "# q [a.u.]\tM_l(q)" : "# Ei [Ry]\tsigma [a_0^2]") << std::endl;

    // collect the abscissae from the command line, or from the standard input when there are none
    rArray xs;
    if (args.size() > 2)
    {
        for (std::size_t i = 2; i < args.size(); i++)
            xs.push_back(std::atof(args[i]));
    }
    else
    {
        for (double x; std::cin >> x; )
            xs.push_back(x);
    }

    for (double x : xs)
        std::cout << x << "\t" << (ff ? formfactor(n,l,x) : cross_section(n,l,x)) << std::endl;

    return EXIT_SUCCESS;
}
