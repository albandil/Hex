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

#include <iostream>
#include <string>
#include <H5Cpp.h>

#include "hdffile.h"
#include "misc.h"

HDFFile::HDFFile(std::string filename, FileAccess flag)
    : file_(nullptr), valid_(true)
{
    try
    {
        switch (flag)
        {
            case readonly:
                file_ = new H5::H5File(filename.c_str(), H5F_ACC_RDONLY);
                break;
            case overwrite:
                file_ = new H5::H5File(filename.c_str(), H5F_ACC_TRUNC);
                break;
            case failifexists:
                file_ = new H5::H5File(filename.c_str(), H5F_ACC_EXCL);
                break;
            case readwrite:
                file_ = new H5::H5File(filename.c_str(), H5F_ACC_RDWR);
                break;
        }
    }
    catch (H5::FileIException e)
    {
        file_ = nullptr;
        valid_ = false;
    }
}

HDFFile::~HDFFile()
{
    if (file_ != nullptr)
    {
        file_->flush(H5F_SCOPE_GLOBAL);
        file_->close();
        delete file_;
    }
}

size_t HDFFile::size(std::string dataset) const
{
    if (not valid_)
        return 0;
    
    try
    {
        H5::DataSet dset = file_->openDataSet(dataset.c_str());
        H5::DataSpace dspc = dset.getSpace();
        size_t length = dspc.getSimpleExtentNpoints();
        return length;
    }
    catch (...)
    {
        return 0;
    }
}

template<> bool HDFFile::read<int>(std::string dataset, int* buffer, size_t length) const
{
    return read_(dataset, buffer, length, H5::IntType(H5::PredType::NATIVE_INT));
}

template<> bool HDFFile::read<long>(std::string dataset, long* buffer, size_t length) const
{
    return read_(dataset, buffer, length, H5::IntType(H5::PredType::NATIVE_LONG));
}

template<> bool HDFFile::read<unsigned long>(std::string dataset, unsigned long * buffer, size_t length) const
{
    return read_(dataset, buffer, length, H5::IntType(H5::PredType::NATIVE_ULONG));
}

template<> bool HDFFile::read<double>(std::string dataset, double* buffer, size_t length) const
{
    return read_(dataset, buffer, length, H5::FloatType(H5::PredType::NATIVE_DOUBLE));
}

template<> bool HDFFile::read<Complex>(std::string dataset, Complex* buffer, size_t length) const
{
    return read_(dataset, buffer, 2*length, H5::FloatType(H5::PredType::NATIVE_DOUBLE));
}

template<> bool HDFFile::write<int>(std::string dataset, int const * buffer, size_t length)
{
    return write_(dataset, buffer, length, H5::IntType(H5::PredType::NATIVE_INT));
}

template<> bool HDFFile::write<long>(std::string dataset, long const * buffer, size_t length)
{
    return write_(dataset, buffer, length, H5::IntType(H5::PredType::NATIVE_LONG));
}

template<> bool HDFFile::write<unsigned long>(std::string dataset, unsigned long const * buffer, size_t length)
{
    return write_(dataset, buffer, length, H5::IntType(H5::PredType::NATIVE_ULONG));
}

template<> bool HDFFile::write<double>(std::string dataset, double const * buffer, size_t length)
{
    return write_(dataset, buffer, length, H5::FloatType(H5::PredType::NATIVE_DOUBLE));
}

template<> bool HDFFile::write<Complex>(std::string dataset, Complex const * buffer, size_t length)
{
    return write_(dataset, buffer, 2*length, H5::FloatType(H5::PredType::NATIVE_DOUBLE));
}

bool HDFFile::read_(std::string dataset, void * buffer, hsize_t length, H5::AtomType dtype) const
{
//     std::cout << "[HDFFile::read_] Read " << length << " elements from \"" << dataset << "\".\n";
    
    if (not valid_)
        return false;
    
    if (length == 0)
        return true;
    try
    {
        H5::DataSet dset = file_->openDataSet(dataset.c_str());
        H5::DataSpace dspc = dset.getSpace();
        
        if (length != dspc.getSimpleExtentNpoints())
            throw exception ("Dimensions do not match, %ld != %ld.", length, dspc.getSimpleExtentNpoints());
        
        dset.read(buffer, dtype, dspc, dspc);
        
        return true;
    }
    catch (H5::FileIException e)
    {
        return false;
    }
}

bool HDFFile::write_(std::string dataset, void const * buffer, hsize_t length, H5::AtomType dtype)
{
    if (not valid_)
        return false;
    
    if (length == 0)
        return true;
    
    try
    {
        H5::DataSpace dspc(1, &length);
        H5::DataSet dset = file_->createDataSet(dataset.c_str(), dtype, dspc);
        
        dset.write(buffer, dtype);
        
        return true;
    }
    catch (H5::FileIException e)
    {
        return false;
    }
}
