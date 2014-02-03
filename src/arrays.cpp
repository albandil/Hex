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

#include <map>
#include <vector>

#include "arrays.h"
#include "complex.h"

rArray abs (const cArrayView u)
{
    rArray v(u.size());
    
    auto iu = u.begin();
    auto iv = v.begin();
    
    while (iu != u.end())
        *(iv++) = abs(*(iu++));
    
    return v;
}

rArrays abs (cArrays const &u)
{
    rArrays v(u.size());
    
    auto iu = u.begin();
    auto iv = v.begin();
    
    while (iu != u.end())
        *(iv++) = abs(*(iu++));
    
    return v;
}

NumberArray<double> hypot (NumberArray<double> const & A, NumberArray<double> const & B)
{
    assert(A.size() == B.size());
    
    size_t N = A.size();
    NumberArray<double> C (N);

    for (size_t i = 0; i < N; i++)
        C[i] = hypot(A[i], B[i]);

    return C;
}

NumberArray<double> atan2 (NumberArray<double> const & A, NumberArray<double> const & B)
{
    assert(A.size() == B.size());
    
    size_t N = A.size();
    NumberArray<double> C (N);

    for (size_t i = 0; i < N; i++)
        C[i] = atan2(A[i], B[i]);

    return C;
}

NumberArray<double> sqrabs (NumberArray<Complex> const & A)
{
    size_t N = A.size();
    NumberArray<double> B (N);

    for (size_t i = 0; i < N; i++)
        B[i] = sqrabs(A[i]);

    return B;
}

NumberArray<double> realpart (NumberArray<Complex> const & A)
{
    size_t N = A.size();
    NumberArray<double> B (N);

    for (size_t i = 0; i < N; i++)
        B[i] = A[i].real();

    return B;
}

NumberArray<double> imagpart (NumberArray<Complex> const & A)
{
    size_t N = A.size();
    NumberArray<double> B (N);

    for (size_t i = 0; i < N; i++)
        B[i] = A[i].imag();

    return B;
}

template<> void write_array (NumberArray<double> array, const char* filename)
{
    std::ofstream fout(filename);
    for (size_t i = 0; i < array.size(); i++)
        fout << array[i] << "\n";
    fout.close();
}

template<> void write_array (NumberArray<double> grid, NumberArray<double> array, const char* filename)
{
    std::ofstream fout(filename);
    for (size_t i = 0; i < array.size(); i++)
        fout << grid[i] << "\t" << array[i] << "\n";
    fout.close();
}

template<> void write_array (NumberArray<Complex> array, const char* filename)
{
    std::ofstream fout(filename);
    for (size_t i = 0; i < array.size(); i++)
        fout << array[i].real() << "\t" << array[i].imag() << "\n";
    fout.close();
}

template<> void write_array (NumberArray<double> grid, NumberArray<Complex> array, const char* filename)
{
    std::ofstream fout(filename);
    for (size_t i = 0; i < array.size(); i++)
        fout << grid[i] << "\t" << array[i].real() << "\t" << array[i].imag() << "\n";
    fout.close();
}

rArray threshold (const rArrayView a, double eps)
{
    rArray b(a.size());
    
    for (size_t i = 0; i < a.size(); i++)
    {
        if (std::abs(a[i]) > eps)
        {
            b[i] = a[i];
        }
    }
    
    return b;
}
