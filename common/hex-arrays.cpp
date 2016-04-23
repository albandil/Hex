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

NumberArray<double> hypot (NumberArray<double> const & A, NumberArray<double> const & B)
{
    assert(A.size() == B.size());
    
    std::size_t N = A.size();
    NumberArray<double> C (N);

    for (std::size_t i = 0; i < N; i++)
        C[i] = hypot(A[i], B[i]);

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

NumberArray<Complex> conjugate (NumberArray<Complex> const & A)
{
    return A.conj();
}

template<> void write_array (const ArrayView<double> array, std::string filename)
{
    std::ofstream fout(filename);
    for (std::size_t i = 0; i < array.size(); i++)
        fout << array[i] << "\n";
    fout.close();
}

template<> void write_array (const ArrayView<double> grid, const ArrayView<double> array, std::string filename)
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

template<> void write_array (const ArrayView<double> grid, const ArrayView<Complex> array, std::string filename)
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
    const rArrayView zgrid,
    bool with_header
)
{
    // array lengths
    unsigned nx = xgrid.size();
    unsigned ny = ygrid.size();
    unsigned nz = zgrid.size();
    unsigned N = nx * ny * nz;
    
    // write VTK header
    if (with_header)
    {
        out << "# vtk DataFile Version 3.0" << std::endl;
        out << "Wave function" << std::endl;
        out << "ASCII" << std::endl;
    }
    
    // write grid data
    out << "DATASET RECTILINEAR_GRID" << std::endl;
    out << "DIMENSIONS " << nx << " " << ny << " " << nz << std::endl;
    out << "X_COORDINATES " << nx << " float" << std::endl;
    out << to_string(xgrid) << std::endl;
    out << "Y_COORDINATES " << ny << " float" << std::endl;
    out << to_string(ygrid) << std::endl;
    out << "Z_COORDINATES " << nz << " float" << std::endl;
    out << to_string(zgrid) << std::endl;
    out << "POINT_DATA " << N << std::endl;
    out << "FIELD wavefunction " << 2 * ev.size() / N << std::endl;
    
    for (std::size_t iblock = 0; (iblock + 1) * N <= ev.size(); iblock++)
    {
        // save real part
        out << "block_" << iblock << "_re 1 " << N << " float" << std::endl;
        for (unsigned i = 0; i < nx; i++)
        {
            for (unsigned j = 0; j < ny; j++)
            for (unsigned k = 0; k < nz; k++)
                out << ev[((iblock * nx + i) * ny + j) * nz + k].real() << " ";
            out << std::endl;
        }
        
        // save imaginary part
        out << "block_" << iblock << "_im 1 " << N << " float" << std::endl;
        for (unsigned i = 0; i < nx; i++)
        {
            for (unsigned j = 0; j < ny; j++)
            for (unsigned k = 0; k < nz; k++)
                out << ev[((iblock * nx + i) * ny + j) * nz + k].imag() << " ";
            out << std::endl;
        }
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
