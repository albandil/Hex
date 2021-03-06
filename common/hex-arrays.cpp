//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2016, Jakub Benda, Charles University in Prague                    //
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

#include <map>
#include <vector>

#include "hex-arrays.h"

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

NumberArray<Real> hypot (NumberArray<Real> const & A, NumberArray<Real> const & B)
{
    assert(A.size() == B.size());

    std::size_t N = A.size();
    NumberArray<Real> C (N);

    for (std::size_t i = 0; i < N; i++)
        C[i] = hypot(A[i], B[i]);

    return C;
}

NumberArray<Real> sqrabs (NumberArray<Complex> const & A)
{
    std::size_t N = A.size();
    NumberArray<Real> B (N);

    for (std::size_t i = 0; i < N; i++)
        B[i] = sqrabs(A[i]);

    return B;
}

NumberArray<Real> realpart (NumberArray<Complex> const & A)
{
    std::size_t N = A.size();
    NumberArray<Real> B (N);

    for (std::size_t i = 0; i < N; i++)
        B[i] = A[i].real();

    return B;
}

NumberArray<Real> imagpart (NumberArray<Complex> const & A)
{
    std::size_t N = A.size();
    NumberArray<Real> B (N);

    for (std::size_t i = 0; i < N; i++)
        B[i] = A[i].imag();

    return B;
}

template<> void write_array (const ArrayView<Real> grid, const ArrayView<Real> array, std::string filename)
{
    std::ofstream fout(filename);
    for (std::size_t i = 0; i < array.size(); i++)
        fout << grid[i] << "\t" << array[i] << "\n";
    fout.close();
}

template<> void write_array (const ArrayView<Complex> array, std::string filename)
{
    std::ofstream fout(filename);
    for (std::size_t i = 0; i < array.size(); i++)
        fout << array[i].real() << "\t" << array[i].imag() << "\n";
    fout.close();
}

template<> void write_array (const ArrayView<Real> grid, const ArrayView<Complex> array, std::string filename)
{
    std::ofstream fout(filename);
    for (std::size_t i = 0; i < array.size(); i++)
        fout << grid[i] << "\t" << array[i].real() << "\t" << array[i].imag() << "\n";
    fout.close();
}

rArray threshold (const rArrayView a, double eps)
{
    rArray b(a.size());

    for (std::size_t i = 0; i < a.size(); i++)
    {
        if (std::abs(a[i]) > eps)
        {
            b[i] = a[i];
        }
    }

    return b;
}

cArray interleave (const rArrayView re, const rArrayView im)
{
    if (re.size() != im.size())
        HexException("Cannot interleave arrays of different sizes (%ld != %ld).", re.size(), im.size());

    cArray output (re.size());

    for (std::size_t i = 0; i < re.size(); i++)
    {
        output[i].real(re[i]);
        output[i].imag(im[i]);
    }

    return output;
}
