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

#ifndef HEX_COMPLEX
#define HEX_COMPLEX

#include <complex>
#include <cmath>

// shorthand for std::complex<double>
typedef std::complex<double> Complex;
typedef std::complex<long double> LComplex;

// some missing complex and mixed arithmetic
template <typename C1, typename C2> auto operator * (std::complex<C1> a, std::complex<C2> b) -> std::complex<decltype(C1(0.)*C2(0.))>
{
    typedef decltype(C1(0.) * C2(0.)) T;
    return std::complex<T>(a) * std::complex<T>(b);
}
template <typename C1, typename C2> auto operator / (std::complex<C1> a, std::complex<C2> b) -> std::complex<decltype(C1(0.)/C2(0.))>
{
    typedef decltype(C1(0.) * C2(0.)) T;
    return std::complex<T>(a) / std::complex<T>(b);
}

/// Squared modulus of a complex number.
inline double sqrabs (Complex z)
{
    return z.real() * z.real() + z.imag() * z.imag();
}

/// Complex ordering by real parts.
inline bool Complex_realpart_less (Complex const & a, Complex const & b)
{
    return a.real() < b.real();
}

/// Complex ordering by imaginary parts.
inline bool Complex_imagpart_less (Complex const & a, Complex const & b)
{
    return a.imag() < b.imag();
}

/// Finite check for complex number.
inline bool Complex_finite (Complex const & z)
{
    return std::isfinite(z.real()) and std::isfinite(z.imag());
}

#endif
