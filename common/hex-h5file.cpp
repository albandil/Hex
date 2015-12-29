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

#include <iostream>
#include <string>

#include <hdf5.h>

#include "hex-h5file.h"

H5File::H5File (std::string filename, FileAccess flag)
    : prefix(), file_(-1), name_(filename), valid_(false)
{
    // separate filesystem path and dataset path (by semicolon)
    std::size_t pos = name_.find(':');
    if (pos != std::string::npos)
    {
        prefix = name_.substr(pos + 1);
        name_.resize(pos);
    }
    
    // try to open the file
    switch (flag)
    {
        case overwrite:
            file_ = H5Fcreate(name_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            break;
            
        case failifexists:
            file_ = H5Fcreate(name_.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
            break;
            
        case readonly:
            file_ = H5Fopen(name_.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
            break;
            
        case readwrite:
            file_ = H5Fopen(name_.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
            break;
    }
    
    // save success indicator
    valid_ = (file_ >= 0);
}

H5File::~H5File ()
{
    if (file_ >= 0)
    {
        // write everything to disk
        H5Fflush(file_, H5F_SCOPE_GLOBAL);
        
        // close the file
        H5Fclose(file_);
        
        // invalidate IDs
        file_ = -1;
        valid_ = false;
    }
}

void H5File::save_error () const
{
    // get number of messages in the error stack
    ssize_t n = H5Eget_num(H5E_DEFAULT);
    
    // output stream
    std::stringstream out;
    
    // for all messages
    H5E_type_t type;
    for (ssize_t i = 0; i < n; i++)
    {
        // get message size
        ssize_t size = H5Eget_msg(i, &type, nullptr, 0);
        
        // skip empty messages
        if (size <= 0)
        {
            out << "  #" << i << std::endl;
            continue;
        }
        
        // retrieve message text
        char text[size];
        size = H5Eget_msg(i, &type, text, size);
        if (size <= 0)
        {
            out << "  #" << i << std::endl;
            continue;
        }
        
        out << "  #" << i << " " << text << std::endl;
    }
    
    // get resulting string
    error_ = out.str();
}

std::size_t H5File::size (std::string dataset) const
{
    // invalid file has zero size
    if (not valid_)
        return 0;
    
    // compose full dataset path
    dataset = prefix + dataset;
    
    // HDF descriptors
    hid_t dset = -1, dspc = -1;
    
    // dataset lenght
    hssize_t length = 0;
    
    // success indicator
    bool success = true;
    
    // open requested dataset
    success = (success and (dset = H5Dopen2(file_, dataset.c_str(), 0)) >= 0);
    
    // open associated dataspace
    success = (success and (dspc = H5Dget_space(dset)) >= 0);
    
    // retrieve size of the dataset
    success = (success and (length = H5Sget_simple_extent_npoints(dspc)) >= 0);
    
    // save possible error message
    save_error();
    
    // release memory
    if (dset >= 0) H5Dclose(dset);
    if (dspc >= 0) H5Sclose(dspc);
    
    // return the size
    return length;
}

bool H5File::read_ (std::string dataset, void * buffer, hsize_t length, hsize_t offset, hid_t dtype) const
{
    // skip reading invalid files
    if (not valid_)
        return false;
    
    // skip reading empty datasets
    if (length == 0)
        return true;
    
    // skip reading into null buffer
    if (buffer == nullptr)
        return false;
    
    // HDF descriptors
    hid_t dset = -1, dspc = -1, mspc = -1;
    
    // compose full dataset path
    dataset = prefix + dataset;
    
    // success indicator
    bool success = true;
    
    // open requested dataset
    success = (success and (dset = H5Dopen2(file_, dataset.c_str(), 0)) >= 0);
    
    // open associated dataspace
    success = (success and (dspc = H5Dget_space(dset)) >= 0);
    
    // select data in the file
    success = (success and H5Sselect_hyperslab(dspc, H5S_SELECT_SET, &offset, nullptr, &length, nullptr) >= 0);
    
    // create simple memory space
    success = (success and (mspc = H5Screate_simple(1, &length, nullptr)) >= 0);
    
    // read data into the buffer
    success = (success and H5Dread(dset, dtype, mspc, dspc, H5P_DEFAULT, buffer) >= 0);
    
    // save possible error mesage
    save_error();
    
    // release HDF descriptors
    if (dset >= 0) H5Dclose(dset);
    if (dspc >= 0) H5Sclose(dspc);
    if (mspc >= 0) H5Sclose(mspc);
    
    return success;
}

bool H5File::write_ (std::string dataset, void const * buffer, hsize_t length, hsize_t offset, hid_t dtype)
{
    // skip writing into invalid files
    if (not valid_)
        return false;
    
    // skip writing empty datasets
    if (length == 0)
        return true;
    
    // compose full dataset path
    dataset = prefix + dataset;
    
    // HDF descriptors
    hid_t dset = -1, dspc = -1, mspc = -1;
    
    // success indicator
    bool success = true;
    
    // try to open an existing dataset
    dset = H5Dopen2(file_, dataset.c_str(), 0);
    
    // a) dataset does not exist yet
    if (dset < 0)
    {
        // create a new data space
        success = (success and (dspc = H5Screate_simple(1, &length, nullptr)) >= 0);
        
        // create a new dataset
        success = (success and (dset = H5Dcreate2(file_, dataset.c_str(), dtype, dspc, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) >= 0);
    }
    
    // b) dataset already exists
    else
    {
        // get assigned data space
        success = (success and (dspc = H5Dget_space(dset)) >= 0);
        
        // select data in the file
        success = (success and H5Sselect_hyperslab(dspc, H5S_SELECT_SET, &offset, nullptr, &length, nullptr) >= 0);
    }
    
    // memory data space
    success = (success and (mspc = H5Screate_simple(1, &length, nullptr)) >= 0);
    
    // write data from buffer if not null
    success = (success and (buffer == nullptr or H5Dwrite(dset, dtype, mspc, dspc, H5P_DEFAULT, buffer) >= 0));
    
    // save possible error message
    save_error();
    
    // release memory
    if (dset >= 0) H5Dclose(dset);
    if (dspc >= 0) H5Sclose(mspc);
    if (mspc >= 0) H5Sclose(dspc);
    
    return success;
}

template<> bool H5File::read<int> (std::string dataset, int * buffer, std::size_t length, std::size_t offset) const
{
    return read_(dataset, buffer, length, offset, H5T_NATIVE_INT);
}

template<> bool H5File::write<int> (std::string dataset, int const * buffer, std::size_t length, std::size_t offset)
{
    return write_(dataset, buffer, length, offset, H5T_NATIVE_INT);
}

template<> bool H5File::read<unsigned> (std::string dataset, unsigned * buffer, std::size_t length, std::size_t offset) const
{
    return read_(dataset, buffer, length, offset, H5T_NATIVE_UINT);
}

template<> bool H5File::write<unsigned> (std::string dataset, unsigned const * buffer, std::size_t length, std::size_t offset)
{
    return write_(dataset, buffer, length, offset, H5T_NATIVE_UINT);
}

template<> bool H5File::read<std::int64_t> (std::string dataset, std::int64_t * buffer, std::size_t length, std::size_t offset) const
{
    return read_(dataset, buffer, length, offset, H5T_NATIVE_INT64);
}

template<> bool H5File::write<std::int64_t> (std::string dataset, std::int64_t const * buffer, std::size_t length, std::size_t offset)
{
    return write_(dataset, buffer, length, offset, H5T_NATIVE_INT64);
}

template<> bool H5File::read<std::uint64_t> (std::string dataset, std::uint64_t * buffer, std::size_t length, std::size_t offset) const
{
    return read_(dataset, buffer, length, offset, H5T_NATIVE_UINT64);
}

template<> bool H5File::write<std::uint64_t> (std::string dataset, std::uint64_t const * buffer, std::size_t length, std::size_t offset)
{
    return write_(dataset, buffer, length, offset, H5T_NATIVE_UINT64);
}

template<> bool H5File::read<double> (std::string dataset, double * buffer, std::size_t length, std::size_t offset) const
{
    return read_(dataset, buffer, length, offset, H5T_NATIVE_DOUBLE);
}

template<> bool H5File::write<double> (std::string dataset, double const * buffer, std::size_t length, std::size_t offset)
{
    return write_(dataset, buffer, length, offset, H5T_NATIVE_DOUBLE);
}

template<> bool H5File::read<Complex> (std::string dataset, Complex * buffer, std::size_t length, std::size_t offset) const
{
    return read_(dataset, buffer, length * 2, offset * 2, H5T_NATIVE_DOUBLE);
}

template<> bool H5File::write<Complex> (std::string dataset, Complex const * buffer, std::size_t length, std::size_t offset)
{
    return write_(dataset, buffer, length * 2, offset * 2, H5T_NATIVE_DOUBLE);
}
