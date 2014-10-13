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
