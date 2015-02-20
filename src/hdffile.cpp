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

#ifndef NO_HDF

#include <iostream>
#include <string>

#include <hdf5.h>

#include "hdffile.h"
#include "misc.h"

HDFFile::HDFFile (std::string filename, FileAccess flag)
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
    
    valid_ = (file_ >= 0);
}

HDFFile::~HDFFile ()
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

std::size_t HDFFile::size (std::string dataset) const
{
    // invalid file has zero size
    if (not valid_)
        return 0;
    
    // compose full dataset path
    dataset = prefix + dataset;
    
    // open requested dataset
    hid_t dset = H5Dopen2(file_, dataset.c_str(), 0);
    if (dset < 0)
        return 0; // non-existing dataset has zero size
    
    // open associated dataspace
    hid_t dspc = H5Dget_space(dset);
    if (dspc < 0)
        return 0; // some unknown error, assumed invalid dataset
    
    // retrieve size of the dataset
    hssize_t length = H5Sget_simple_extent_npoints(dspc);
    
    // release memory
    H5Sclose(dspc);
    H5Dclose(dset);
    
    // return the size
    return length;
}

bool HDFFile::read_ (std::string dataset, void * buffer, hsize_t length, hsize_t offset, hid_t dtype) const
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
    
    // compose full dataset path
    dataset = prefix + dataset;
    
    // open requested dataset
    hid_t dset = H5Dopen2(file_, dataset.c_str(), 0);
    if (dset < 0)
        return false;
    
    // open associated dataspace
    hid_t dspc = H5Dget_space(dset);
    if (dspc < 0)
    {
        H5Dclose(dset);
        return false;
    }
    
    // select data in the file
    herr_t res = H5Sselect_hyperslab(dspc, H5S_SELECT_SET, &offset, nullptr, &length, nullptr);
    if (res < 0)
    {
        H5Dclose(dset);
        H5Sclose(dspc);
        return false;
    }
    
    // create simple memory space
    hid_t mspc = H5Screate_simple(1, &length, nullptr);
    if (mspc < 0)
    {
        H5Dclose(dset);
        H5Sclose(dspc);
        return false;
    }
    
    // read data
    herr_t err = H5Dread(dset, dtype, mspc, dspc, H5P_DEFAULT, buffer);
    if (err < 0)
    {
        H5Dclose(dset);
        H5Sclose(dspc);
        H5Sclose(mspc);
        return false;
    }
    
    // release memory
    H5Sclose(mspc);
    H5Sclose(dspc);
    H5Dclose(dset);
    
    return true;
}

bool HDFFile::write_ (std::string dataset, void const * buffer, hsize_t length, hsize_t offset, hid_t dtype)
{
    // skip writing into invalid files
    if (not valid_)
        return false;
    
    // skip writing empty datasets
    if (length == 0)
        return true;
    
    // compose full dataset path
    dataset = prefix + dataset;
    
    // file data space
    hid_t dspc;
    
    // try to open an existing dataset
    hid_t dset = H5Dopen2(file_, dataset.c_str(), 0);
    if (dset < 0)
    {
        // create a new data space
        dspc = H5Screate_simple(1, &length, nullptr);
        if (dspc < 0)
        {
            H5Dclose(dset);
            return false;
        }
        
        // create a new dataset
        dset = H5Dcreate2(file_, dataset.c_str(), dtype, dspc, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (dset < 0)
        {
            H5Sclose(dspc);
            H5Dclose(dset);
            return false;
        }
    }
    else
    {
        // get assigned data space
        dspc = H5Dget_space(dset);
        if (dspc < 0)
        {
            H5Dclose(dset);
            return false;
        }
        
        // select data in the file
        herr_t res = H5Sselect_hyperslab(dspc, H5S_SELECT_SET, &offset, nullptr, &length, nullptr);
        if (res < 0)
        {
            H5Dclose(dset);
            H5Sclose(dspc);
            return false;
        }
    }
    
    // memory data space
    hid_t mspc = H5Screate_simple(1, &length, nullptr);
    if (mspc < 0)
    {
        H5Dclose(dset);
        H5Sclose(dspc);
        return false;
    }
    
    // write data from buffer
    herr_t err = H5Dwrite(dset, dtype, mspc, dspc, H5P_DEFAULT, buffer);
    if (err < 0)
    {
        H5Dclose(dset);
        H5Sclose(dspc);
        H5Sclose(mspc);
        return false;
    }
    
    // release memory
    H5Dclose(dset);
    H5Sclose(dspc);
    H5Sclose(mspc);
    
    return true;
}

#endif
