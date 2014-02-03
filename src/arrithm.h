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

#ifndef HEX_ARRITHM
#define HEX_ARRITHM

#include <type_traits>

#include "arrays.h"
#include "misc.h"

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
