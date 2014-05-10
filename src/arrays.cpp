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

template<> void write_array (const ArrayView<double> array, const char* filename)
{
    std::ofstream fout(filename);
    for (size_t i = 0; i < array.size(); i++)
        fout << array[i] << "\n";
    fout.close();
}

template<> void write_array (const ArrayView<double> grid, const ArrayView<double> array, const char* filename)
{
    std::ofstream fout(filename);
    for (size_t i = 0; i < array.size(); i++)
        fout << grid[i] << "\t" << array[i] << "\n";
    fout.close();
}

template<> void write_array (const ArrayView<Complex> array, const char* filename)
{
    std::ofstream fout(filename);
    for (size_t i = 0; i < array.size(); i++)
        fout << array[i].real() << "\t" << array[i].imag() << "\n";
    fout.close();
}

template<> void write_array (const ArrayView<double> grid, const ArrayView<Complex> array, const char* filename)
{
    std::ofstream fout(filename);
    for (size_t i = 0; i < array.size(); i++)
        fout << grid[i] << "\t" << array[i].real() << "\t" << array[i].imag() << "\n";
    fout.close();
}

void writeVTK
(
    std::ofstream & out,
    const cArrayView ev,
    const rArrayView xgrid,
    const rArrayView ygrid,
    const rArrayView zgrid
)
{
    // array lengths
    int nx = xgrid.size();
    int ny = ygrid.size();
    int nz = zgrid.size();
    int N = nx * ny * nz;
    
    // write VTK header
    out << "# vtk DataFile Version 3.0\n";
    out << "Wave function\n";
    out << "ASCII\n";
    out << "DATASET RECTILINEAR_GRID\n";
    out << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";
    out << "X_COORDINATES " << nx << " float\n";
    out << to_string(xgrid) << "\n";
    out << "Y_COORDINATES " << ny << " float\n";
    out << to_string(ygrid) << "\n";
    out << "Z_COORDINATES " << nz << " float\n";
    out << to_string(zgrid) << "\n";
    out << "POINT_DATA " << N << "\n";
    out << "FIELD wavefunction 2\n";
    
    // save real part
    out << "realpart 1 " << N << " float\n";
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        for (int k = 0; k < nz; k++)
            out << ev[(i * ny + j) * nz + k].real() << " ";
        out << "\n";
    }
    
    // save imaginary part
    out << "imagpart 1 " << N << " float\n";
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        for (int k = 0; k < nz; k++)
            out << ev[(i * ny + j) * nz + k].imag() << " ";
        out << "\n";
    }
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
