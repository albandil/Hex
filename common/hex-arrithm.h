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

#ifndef HEX_ARRITHM
#define HEX_ARRITHM

#include <type_traits>

#include "hex-arrays.h"
#include "hex-misc.h"

//
// unary operators
//

#define declareDeepUnaryOperatorTemplates(IN,OP)                              \
                                                                              \
template <class T>                                                            \
IN<T> operator OP (IN<T> const & a)                                           \
{                                                                             \
    IN<T> b(a);                                                               \
    for (T & e : b)                                                           \
        e = OP e;                                                             \
    return b;                                                                 \
}

declareDeepUnaryOperatorTemplates(NumberArray, -)

//
// reduced arithmetic operators
//

#define declareDeepReducedArithmeticTemplates(OUT,OP,IN)                      \
                                                                              \
template <                                                                    \
    class T1,                                                                 \
    class T2,                                                                 \
    class = typename std::enable_if<is_scalar<T2>::value>::type               \
> OUT<T1> operator OP##= (OUT<T1> & a, T2 const & b)                          \
{                                                                             \
    T1 * restrict pa = a.data();                                              \
                                                                              \
    size_t N = a.size();                                                      \
    for (size_t i = 0; i < N; i++)                                            \
        pa[i] OP##= b;                                                        \
                                                                              \
    return a;                                                                 \
}                                                                             \
                                                                              \
template <class T1, class T2>                                                 \
OUT<T1> operator OP##= (OUT<T1> & a, IN<T2> const & b)                        \
{                                                                             \
    assert(a.size() == b.size());                                             \
                                                                              \
    T1       * restrict pa = a.data();                                        \
    T2 const * restrict pb = b.data();                                        \
                                                                              \
    size_t N = a.size();                                                      \
    for (size_t i = 0; i < N; i++)                                            \
        pa[i] OP##= pb[i];                                                    \
                                                                              \
    return a;                                                                 \
}

declareDeepReducedArithmeticTemplates(NumberArray, +, NumberArray)
declareDeepReducedArithmeticTemplates(NumberArray, -, NumberArray)
declareDeepReducedArithmeticTemplates(NumberArray, *, NumberArray)
declareDeepReducedArithmeticTemplates(NumberArray, /, NumberArray)

#define declareShallowReducedArithmeticTemplates(OUT,OP,IN)                   \
                                                                              \
template <                                                                    \
    class T1,                                                                 \
    class T2,                                                                 \
    class = typename std::enable_if<is_scalar<T2>::value>::type               \
> OUT<T1> operator OP##= (OUT<T1> a, T2 const & b)                            \
{                                                                             \
    T1 * restrict pa = a.data();                                              \
                                                                              \
    size_t N = a.size();                                                      \
    for (size_t i = 0; i < N; i++)                                            \
        pa[i] OP##= b;                                                        \
                                                                              \
    return a;                                                                 \
}                                                                             \
                                                                              \
template <class T1, class T2>                                                 \
OUT<T1> operator OP##= (OUT<T1> a, IN<T2> b)                                  \
{                                                                             \
    assert(a.size() == b.size());                                             \
                                                                              \
    T1       * restrict pa = a.data();                                        \
    T2 const * restrict pb = b.data();                                        \
                                                                              \
    size_t N = a.size();                                                      \
    for (size_t i = 0; i < N; i++)                                            \
        pa[i] OP##= pb[i];                                                    \
                                                                              \
    return a;                                                                 \
}

declareShallowReducedArithmeticTemplates(ArrayView, +, ArrayView)
declareShallowReducedArithmeticTemplates(ArrayView, -, ArrayView)
declareShallowReducedArithmeticTemplates(ArrayView, *, ArrayView)
declareShallowReducedArithmeticTemplates(ArrayView, /, ArrayView)

//
// binary arithmetic operators for NumberArray
//

#define declareDeepBinaryOperatorTemplates(OUT,OP,IN)                         \
                                                                              \
template <                                                                    \
    class T1,                                                                 \
    class T2,                                                                 \
    class T3 = decltype(T1(1) OP T2(1))                                       \
> OUT<T3> operator OP (IN<T1> const & a, IN<T2> const & b)                    \
{                                                                             \
    assert(a.size() == b.size());                                             \
                                                                              \
    size_t N = a.size();                                                      \
    OUT<T3> c(N);                                                             \
                                                                              \
    T1 const * restrict pa = a.data();                                        \
    T2 const * restrict pb = b.data();                                        \
    T3       * restrict pc = c.data();                                        \
                                                                              \
    for (size_t i = 0; i < N; i++)                                            \
        pc[i] = pa[i] OP pb[i];                                               \
                                                                              \
    return c;                                                                 \
}                                                                             \
                                                                              \
template <                                                                    \
    class T1,                                                                 \
    class T2,                                                                 \
    class T3 = decltype(T1(1) OP T2(1)),                                      \
    class = typename std::enable_if<is_scalar<T2>::value>::type               \
> OUT<T3> operator OP (IN<T1> const & a, T2 const & b)                        \
{                                                                             \
    size_t N = a.size();                                                      \
    OUT<T3> c(N);                                                             \
                                                                              \
    T1 const * restrict pa = a.data();                                        \
    T3       * restrict pc = c.data();                                        \
                                                                              \
    for (size_t i = 0; i < N; i++)                                            \
        pc[i] = pa[i] OP b;                                                   \
                                                                              \
    return c;                                                                 \
}                                                                             \
template <                                                                    \
    class T1,                                                                 \
    class T2,                                                                 \
    class T3 = decltype(T1(1) OP T2(1)),                                      \
    class = typename std::enable_if<is_scalar<T1>::value>::type               \
> OUT<T3> operator OP (T1 const & a, IN<T2> const & b)                        \
{                                                                             \
    size_t N = b.size();                                                      \
    OUT<T3> c(N);                                                             \
                                                                              \
    T2 const * restrict pb = b.data();                                        \
    T3       * restrict pc = c.data();                                        \
                                                                              \
    for (size_t i = 0; i < N; i++)                                            \
        pc[i] = a OP pb[i];                                                   \
                                                                              \
    return c;                                                                 \
}                                                                             \

declareDeepBinaryOperatorTemplates(NumberArray, +, NumberArray)
declareDeepBinaryOperatorTemplates(NumberArray, -, NumberArray)
declareDeepBinaryOperatorTemplates(NumberArray, *, NumberArray)
declareDeepBinaryOperatorTemplates(NumberArray, /, NumberArray)

declareDeepBinaryOperatorTemplates(Array, +, ArrayView)
declareDeepBinaryOperatorTemplates(Array, -, ArrayView)
declareDeepBinaryOperatorTemplates(Array, *, ArrayView)
declareDeepBinaryOperatorTemplates(Array, /, ArrayView)

#endif
