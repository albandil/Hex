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
    
    std::size_t N = A.size();
    NumberArray<double> C (N);

    for (std::size_t i = 0; i < N; i++)
        C[i] = hypot(A[i], B[i]);

    return C;
}

NumberArray<double> atan2 (NumberArray<double> const & A, NumberArray<double> const & B)
{
    assert(A.size() == B.size());
    
    std::size_t N = A.size();
    NumberArray<double> C (N);

    for (std::size_t i = 0; i < N; i++)
        C[i] = atan2(A[i], B[i]);

    return C;
}

NumberArray<double> sqrabs (NumberArray<Complex> const & A)
{
    std::size_t N = A.size();
    NumberArray<double> B (N);

    for (std::size_t i = 0; i < N; i++)
        B[i] = sqrabs(A[i]);

    return B;
}

NumberArray<double> realpart (NumberArray<Complex> const & A)
{
    std::size_t N = A.size();
    NumberArray<double> B (N);

    for (std::size_t i = 0; i < N; i++)
        B[i] = A[i].real();

    return B;
}

NumberArray<double> imagpart (NumberArray<Complex> const & A)
{
    std::size_t N = A.size();
    NumberArray<double> B (N);

    for (std::size_t i = 0; i < N; i++)
        B[i] = A[i].imag();

    return B;
}

template<> void write_array (const ArrayView<double> array, const char* filename)
{
    std::ofstream fout(filename);
    for (std::size_t i = 0; i < array.size(); i++)
        fout << array[i] << "\n";
    fout.close();
}

template<> void write_array (const ArrayView<double> grid, const ArrayView<double> array, const char* filename)
{
    std::ofstream fout(filename);
    for (std::size_t i = 0; i < array.size(); i++)
        fout << grid[i] << "\t" << array[i] << "\n";
    fout.close();
}

template<> void write_array (const ArrayView<Complex> array, const char* filename)
{
    std::ofstream fout(filename);
    for (std::size_t i = 0; i < array.size(); i++)
        fout << array[i].real() << "\t" << array[i].imag() << "\n";
    fout.close();
}

template<> void write_array (const ArrayView<double> grid, const ArrayView<Complex> array, const char* filename)
{
    std::ofstream fout(filename);
    for (std::size_t i = 0; i < array.size(); i++)
        fout << grid[i] << "\t" << array[i].real() << "\t" << array[i].imag() << "\n";
    fout.close();
}

void writeVTK_points
(
    std::ofstream & out,
    const cArrayView ev,
    const rArrayView xgrid,
    const rArrayView ygrid,
    const rArrayView zgrid
)
{
    // array lengths
    unsigned nx = xgrid.size();
    unsigned ny = ygrid.size();
    unsigned nz = zgrid.size();
    unsigned N = nx * ny * nz;
    
    assert (ev.size() == N);
    
    // write VTK header
    out << "# vtk DataFile Version 3.0" << std::endl;
    out << "Wave function" << std::endl;
    out << "ASCII" << std::endl;
    out << "DATASET RECTILINEAR_GRID" << std::endl;
    out << "DIMENSIONS " << nx << " " << ny << " " << nz << std::endl;
    out << "X_COORDINATES " << nx << " float" << std::endl;
    out << to_string(xgrid) << std::endl;
    out << "Y_COORDINATES " << ny << " float" << std::endl;
    out << to_string(ygrid) << std::endl;
    out << "Z_COORDINATES " << nz << " float" << std::endl;
    out << to_string(zgrid) << std::endl;
    out << "POINT_DATA " << N << std::endl;
    out << "FIELD wavefunction 2" << std::endl;
    
    // save real part
    out << "Re 1 " << N << " float" << std::endl;
    for (unsigned i = 0; i < nx; i++)
    {
        for (unsigned j = 0; j < ny; j++)
        for (unsigned k = 0; k < nz; k++)
            out << ev[(i * ny + j) * nz + k].real() << " ";
        out << std::endl;
    }
    
    // save imaginary part
    out << "Im 1 " << N << " float" << std::endl;
    for (unsigned i = 0; i < nx; i++)
    {
        for (unsigned j = 0; j < ny; j++)
        for (unsigned k = 0; k < nz; k++)
            out << ev[(i * ny + j) * nz + k].imag() << " ";
        out << std::endl;
    }
}

void writeVTK_cells
(
    std::ofstream & out,
    const cArrayView ev,
    const rArrayView xgrid,
    const rArrayView ygrid,
    const rArrayView zgrid
)
{
    // point counts, cell counts
    unsigned px = xgrid.size(), cx = px - 1;
    unsigned py = ygrid.size(), cy = py - 1;
    unsigned pz = zgrid.size(), cz = pz - 1;
    unsigned N = cx * cy * cz;
    
    assert (ev.size() == N);
    
    // write VTK header
    out << "# vtk DataFile Version 3.0" << std::endl;
    out << "Wave function" << std::endl;
    out << "ASCII" << std::endl;
    out << "DATASET RECTILINEAR_GRID" << std::endl;
    out << "DIMENSIONS " << px << " " << py << " " << pz << std::endl;
    out << "X_COORDINATES " << px << " float" << std::endl;
    out << to_string(xgrid) << std::endl;
    out << "Y_COORDINATES " << py << " float" << std::endl;
    out << to_string(ygrid) << std::endl;
    out << "Z_COORDINATES " << pz << " float" << std::endl;
    out << to_string(zgrid) << std::endl;
    out << "CELL_DATA " << N << std::endl;
    out << "FIELD wavefunction 2" << std::endl;
    
    // save real part
    out << "Re 1 " << N << " float" << std::endl;
    for (unsigned i = 0; i < cx; i++)
    {
        for (unsigned j = 0; j < cy; j++)
        for (unsigned k = 0; k < cz; k++)
            out << ev[(i * cy + j) * cz + k].real() << " ";
        out << std::endl;
    }
    
    // save imaginary part
    out << "Im 1 " << N << " float" << std::endl;
    for (unsigned i = 0; i < cx; i++)
    {
        for (unsigned j = 0; j < cy; j++)
        for (unsigned k = 0; k < cz; k++)
            out << ev[(i * cy + j) * cz + k].imag() << " ";
        out << std::endl;
    }
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
