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

#include <map>
#include <vector>

#ifndef NO_HDF
#include <H5Cpp.h>
#endif

#include "arrays.h"
#include "complex.h"

rArray abs (cArray const &u)
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
    Array<double> C (N);

    for (size_t i = 0; i < N; i++)
        C[i] = hypot(A[i], B[i]);

    return C;
}

NumberArray<double> atan2 (NumberArray<double> const & A, NumberArray<double> const & B)
{
    assert(A.size() == B.size());
    
    size_t N = A.size();
    Array<double> C (N);

    for (size_t i = 0; i < N; i++)
        C[i] = atan2(A[i], B[i]);

    return C;
}

NumberArray<double> sqrabs (NumberArray<Complex> const & A)
{
    size_t N = A.size();
    Array<double> B (N);

    for (size_t i = 0; i < N; i++)
        B[i] = sqrabs(A[i]);

    return B;
}

NumberArray<double> realpart (NumberArray<Complex> const & A)
{
    size_t N = A.size();
    Array<double> B (N);

    for (size_t i = 0; i < N; i++)
        B[i] = A[i].real();

    return B;
}

NumberArray<double> imagpart (NumberArray<Complex> const & A)
{
    size_t N = A.size();
    Array<double> B (N);

    for (size_t i = 0; i < N; i++)
        B[i] = A[i].imag();

    return B;
}

#ifndef NO_HDF
bool save_array(rArray const & vec, const char* name, const double * const pdelta)
{
    try
    {
        H5::H5File h5file(name, H5F_ACC_TRUNC);
        
        int rank = 1;
        hsize_t length = vec.size();
        
        // save data
        H5::DataSpace dspc(rank, &length);
        H5::IntType dtype(H5::PredType::NATIVE_DOUBLE);
        H5::DataSet dset = h5file.createDataSet("data", dtype, dspc);
        dset.write(&vec[0], H5::PredType::NATIVE_DOUBLE);
        
        // save delta (if wanted)
        if (pdelta != 0)
        {
            length = 1;
            H5::DataSpace dspc(rank, &length);
            H5::IntType dtype(H5::PredType::NATIVE_DOUBLE);
            H5::DataSet dset = h5file.createDataSet("delta", dtype, dspc);
            dset.write(pdelta, H5::PredType::NATIVE_DOUBLE);
        }
        
        return true;
    }
    catch (...)
    {
        return false;
    }
}

bool save_array(cArray const & vec, const char* name)
{
    try
    {
        H5::H5File h5file(name, H5F_ACC_TRUNC);
        
        int rank = 1;
        hsize_t length = 2 * vec.size();
        
        // save data as an interleaved array
        H5::DataSpace dspc(rank, &length);
        H5::IntType dtype(H5::PredType::NATIVE_DOUBLE);
        H5::DataSet dset = h5file.createDataSet("data", dtype, dspc);
        dset.write(reinterpret_cast<const double*>(&vec[0]), H5::PredType::NATIVE_DOUBLE);
        
        return true;
    }
    catch (...)
    {
        return false;
    }
}

bool load_array(rArray & vec, const char* name, double* pdelta)
{
    try
    {
        H5::Exception::dontPrint();
        H5::H5File h5file(name, H5F_ACC_RDONLY);
        
        // load data 
        H5::DataSet dset = h5file.openDataSet("data");
        H5::DataSpace dspc = dset.getSpace();
        size_t N = dspc.getSimpleExtentNpoints();
        vec.resize(N);
        dset.read(&vec[0], H5::PredType::NATIVE_DOUBLE, dspc, dspc);
        
        // load delta (if wanted)
        if (pdelta != 0)
        {
            dset = h5file.openDataSet("delta");
            dspc = dset.getSpace();
            dset.read(pdelta, H5::PredType::NATIVE_DOUBLE, dspc, dspc);
        }
        
        return true;
    }
    catch (...)
    {
        return false;
    }
}

bool load_array(cArray & vec, const char* name)
{
    try
    {
        H5::Exception::dontPrint();
        H5::H5File h5file(name, H5F_ACC_RDONLY);
        
        // load data 
        H5::DataSet dset = h5file.openDataSet("data");
        H5::DataSpace dspc = dset.getSpace();
        size_t N = dspc.getSimpleExtentNpoints();
        vec.resize(N/2);
        dset.read(reinterpret_cast<Complex*>(&vec[0]), H5::PredType::NATIVE_DOUBLE, dspc, dspc);
        
        return true;
    }
    catch (...)
    {
        return false;
    }
}
#endif

template<> void write_array(ArrayView<double> array, const char* filename)
{
    std::ofstream fout(filename);
    for (size_t i = 0; i < array.size(); i++)
        fout << array[i] << "\n";
    fout.close();
}

template<> void write_array(ArrayView<double> grid, ArrayView<double> array, const char* filename)
{
    std::ofstream fout(filename);
    for (size_t i = 0; i < array.size(); i++)
        fout << grid[i] << "\t" << array[i] << "\n";
    fout.close();
}

template<> void write_array(ArrayView<Complex> array, const char* filename)
{
    std::ofstream fout(filename);
    for (size_t i = 0; i < array.size(); i++)
        fout << array[i].real() << "\t" << array[i].imag() << "\n";
    fout.close();
}

template<> void write_array(ArrayView<double> grid, ArrayView<Complex> array, const char* filename)
{
    std::ofstream fout(filename);
    for (size_t i = 0; i < array.size(); i++)
        fout << grid[i] << "\t" << array[i].real() << "\t" << array[i].imag() << "\n";
    fout.close();
}
