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

#ifndef HEX_COMPLEX
#define HEX_COMPLEX

#include <complex>

// shorthand for std::complex<double>
typedef std::complex<double> Complex;
typedef std::complex<long double> LComplex;

// define some missing complex and mixed arithmetic
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

#endif
