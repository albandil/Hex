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

#ifndef HEX_ARRITHM
#define HEX_ARRITHM

#include "arrays.h"
#include "misc.h"

//
// unary operators
//

#define declareDeepUnaryOperatorTemplates(IN,OP)                              \
                                                                              \
template <class T> IN<T> operator OP (IN<T> const & a)                        \
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
> auto operator OP##= (OUT<T1> & a, T2 const & b)                             \
                                    -> OUT<decltype(T1(1) OP T2(1))>          \
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
auto operator OP##= (OUT<T1> & a, IN<T2> const & b)                           \
                                    -> OUT<decltype(T1(1) OP T2(1))>          \
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
> auto operator OP##= (OUT<T1> a, T2 const & b)                               \
                                    -> OUT<decltype(T1(1) OP T2(1))>          \
{                                                                             \
    T1 * restrict pa = a.data();                                              \
                                                                              \
    size_t N = a.size();                                                      \
    for (size_t i = 0; i < N; i++)                                            \
        pa[i] OP b;                                                           \
                                                                              \
    return a;                                                                 \
}                                                                             \
                                                                              \
template <class T1, class T2>                                                 \
auto operator OP##= (OUT<T1> a, IN<T2> b)                                     \
                                    -> OUT<decltype(T1(1) OP T2(1))>          \
{                                                                             \
    assert(a.size() == b.size());                                             \
                                                                              \
    T1       * restrict pa = a.data();                                        \
    T2 const * restrict pb = b.data();                                        \
                                                                              \
    size_t N = a.size();                                                      \
    for (size_t i = 0; i < N; i++)                                            \
        pa[i] OP pb[i];                                                       \
                                                                              \
    return a;                                                                 \
}

declareShallowReducedArithmeticTemplates(ArrayView, +, ConstArrayView)
declareShallowReducedArithmeticTemplates(ArrayView, -, ConstArrayView)
declareShallowReducedArithmeticTemplates(ArrayView, *, ConstArrayView)
declareShallowReducedArithmeticTemplates(ArrayView, /, ConstArrayView)

//
// binary arithmetic operators for NumberArray
//

#define declareDeepBinaryOperatorTemplates(OUT,OP,IN)                         \
                                                                              \
template <class T1, class T2>                                                 \
auto operator OP (IN<T1> const & a, IN<T2> const & b)                         \
                                             -> OUT<decltype(T1(1) OP T2(1))> \
{                                                                             \
    assert(a.size() == b.size());                                             \
                                                                              \
    size_t N = a.size();                                                      \
    OUT<decltype(T1(1) OP T2(1))> c(N);                                       \
                                                                              \
    T1 const * restrict pa = a.data();                                        \
    T1 const * restrict pb = b.data();                                        \
    T1       * restrict pc = c.data();                                        \
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
    class = typename std::enable_if<is_scalar<T2>::value>::type               \
> auto operator OP (IN<T1> const & a, T2 const & b)                           \
                                             -> OUT<decltype(T1(1) OP T2(1))> \
{                                                                             \
    size_t N = a.size();                                                      \
    OUT<decltype(T1(1) OP T2(1))> c(N);                                       \
                                                                              \
    T1 const * restrict pa = a.data();                                        \
    T1       * restrict pc = c.data();                                        \
                                                                              \
    for (size_t i = 0; i < N; i++)                                            \
        pc[i] = pa[i] OP b;                                                   \
                                                                              \
    return c;                                                                 \
}                                                                             \
template <                                                                    \
    class T1,                                                                 \
    class T2,                                                                 \
    class = typename std::enable_if<is_scalar<T1>::value>::type               \
> auto operator OP (T1 const & a, IN<T2> const & b)                           \
                                             -> OUT<decltype(T1(1) OP T2(1))> \
{                                                                             \
    return operator OP (b, a);                                                \
}                                                                             \

declareDeepBinaryOperatorTemplates(NumberArray, +, NumberArray)
declareDeepBinaryOperatorTemplates(NumberArray, -, NumberArray)
declareDeepBinaryOperatorTemplates(NumberArray, *, NumberArray)
declareDeepBinaryOperatorTemplates(NumberArray, /, NumberArray)

declareDeepBinaryOperatorTemplates(NumberArray, +, ConstArrayView)
declareDeepBinaryOperatorTemplates(NumberArray, -, ConstArrayView)
declareDeepBinaryOperatorTemplates(NumberArray, *, ConstArrayView)
declareDeepBinaryOperatorTemplates(NumberArray, /, ConstArrayView)

#endif
