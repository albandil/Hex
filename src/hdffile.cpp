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
#include <H5Cpp.h>

#include "hdffile.h"
#include "misc.h"

HDFFile::HDFFile (std::string filename, FileAccess flag)
    : prefix(), file_(nullptr), name_(filename), valid_(true)
{
    // separate filesystem path and dataset path (by semicolon)
    std::size_t pos = name_.find(':');
    if (pos != std::string::npos)
    {
        prefix = name_.substr(pos + 1);
        name_.resize(pos);
    }
    
    // try to open the file
    try
    {
        switch (flag)
        {
            case readonly:
                file_ = new H5::H5File(name_.c_str(), H5F_ACC_RDONLY);
                break;
            case overwrite:
                file_ = new H5::H5File(name_.c_str(), H5F_ACC_TRUNC);
                break;
            case failifexists:
                file_ = new H5::H5File(name_.c_str(), H5F_ACC_EXCL);
                break;
            case readwrite:
                file_ = new H5::H5File(name_.c_str(), H5F_ACC_RDWR);
                break;
        }
    }
    catch (H5::FileIException & e)
    {
        file_ = nullptr;
        valid_ = false;
    }
}

HDFFile::~HDFFile ()
{
    if (file_ != nullptr)
    {
        try
        {
            file_->flush(H5F_SCOPE_GLOBAL);
            file_->close();
            delete file_;
        }
        catch (H5::Exception & e)
        {
            e.printErrorStack();
            Exception("Failed to close the HDF file \"%s\".", name_.c_str());
        }
    }
}

std::size_t HDFFile::size (std::string dataset) const
{
    if (not valid_)
        return 0;
    
    try
    {
        dataset = prefix + dataset;
        H5::DataSet dset = file_->openDataSet(dataset.c_str());
        H5::DataSpace dspc = dset.getSpace();
        std::size_t length = dspc.getSimpleExtentNpoints();
        return length;
    }
    catch (...)
    {
        // the dataset doesn't exist
        return 0; // Exception("Failed to retrieve size of HDF dataset \"%s\".", dataset.c_str());
    }
}

bool HDFFile::read_ (std::string dataset, void * buffer, hsize_t length, hsize_t offset, H5::AtomType dtype) const
{
    if (not valid_)
        return false;
    
    if (length == 0)
        return true;
    try
    {
        
        dataset = prefix + dataset;
        H5::DataSet dset = file_->openDataSet(dataset.c_str());
        H5::DataSpace dspc = dset.getSpace();
        
        dspc.selectHyperslab(H5S_SELECT_SET, &length, &offset);
        
        H5::DataSpace mspc(1, &length);
        
        if (buffer != nullptr)
            dset.read(buffer, dtype, mspc, dspc);
        
        return true;
    }
    catch (H5::Exception & e)
    {
        e.printErrorStack();
        Exception("Failed to read HDF dataset \"%s\" from file \"%s\".", dataset.c_str(), name_.c_str());
    }
}

bool HDFFile::write_ (std::string dataset, void const * buffer, hsize_t length, hsize_t offset, H5::AtomType dtype)
{
    if (not valid_)
        return false;
    
    if (length == 0)
        return true;
    
    try
    {
        dataset = prefix + dataset;
        
        H5::DataSet dset;
        H5::DataSpace mspc(1, &length);
        H5::DataSpace dspc(1, &length);
        
        try
        {
            dset = file_->openDataSet(dataset.c_str());
            dspc = dset.getSpace();
            dspc.selectHyperslab(H5S_SELECT_SET, &length, &offset); 
        }
        catch (H5::Exception)
        {
            dset = file_->createDataSet(dataset.c_str(), dtype, dspc);
        }
        
        if (buffer != nullptr)
            dset.write(buffer, dtype, mspc, dspc);
        
        return true;
    }
    catch (H5::Exception & e)
    {
        e.printErrorStack();
        Exception("Failed to write HDF dataset \"%s\" to file \"%s\".", dataset.c_str(), name_.c_str());
    }
}

#endif
