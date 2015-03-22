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

#include <cstdio>

#include "hdffile.h"
#include "misc.h"

const std::size_t HDFFile::header_byte_size = 1024;

HDFFile::HDFFile (std::string filename, FileAccess flag)
    : prefix(), file_(nullptr), name_(filename), changed_(false), valid_(false)
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
            file_ = std::fopen(name_.c_str(), "wb+");
            break;
            
        case readonly:
            file_ = std::fopen(name_.c_str(), "rb");
            break;
            
        case readwrite:
            file_ = std::fopen(name_.c_str(), "rb+");
            break;
    }
    
    // save success indicator
    valid_ = (file_ != nullptr);
    
    // read header
    if (valid_ and (flag == readonly or flag == readwrite))
    {
        // read number of datasets
        std::size_t ndt;
        if (std::fread(&ndt, sizeof(std::size_t), 1, file_) < 1)
        {
            // no data yet...
            return;
        }
        datasets_.resize(ndt);
        
        // read dataset information
        for (std::size_t idt = 0; idt < ndt; idt++)
        {
            // read length of the dataset name
            std::size_t len;
            std::fread(&len, sizeof(std::size_t), 1, file_);
            
            // read dataset name
            datasets_[idt].name.resize(len);
            std::fread(&(datasets_[idt].name[0]), sizeof(char), len, file_);
            
            // read datatype byte size and element count
            std::fread(&(datasets_[idt].bytes), sizeof(std::size_t), 1, file_);
            std::fread(&(datasets_[idt].elements), sizeof(std::size_t), 1, file_);
        }
    }
}

HDFFile::~HDFFile ()
{
    if (file_ != nullptr)
    {
        // write the modified header
        if (changed_)
        {
            // write number of datasets
            std::size_t ndt = datasets_.size();
            std::fseek(file_, 0, SEEK_SET);
            std::fwrite(&ndt, sizeof(ndt), 1, file_);
            
            // write datasets information
            for (std::size_t idt = 0; idt < ndt; idt++)
            {
                // write length of dataset name
                std::size_t len = datasets_[idt].name.size();
                std::fwrite(&len, sizeof(ndt), 1, file_);
                
                // write dataset name
                std::fwrite(&(datasets_[idt].name[0]), sizeof(char), len, file_);
                
                // write datatype byte size and element count
                std::fwrite(&(datasets_[idt].bytes), sizeof(std::size_t), 1, file_);
                std::fwrite(&(datasets_[idt].elements), sizeof(std::size_t), 1, file_);
            }
        }
        
        std::fclose(file_);
        file_ = nullptr;
        valid_ = false;
    }
}

std::size_t HDFFile::size (std::string dataset) const
{
    // try to find the dataset with this name
    for (std::size_t i = 0; i < datasets_.size(); i++)
    {
        if (datasets_[i].name == dataset)
            return datasets_[i].elements;
    }
    
    return 0;
}

bool HDFFile::read_ (std::string dataset, void * buffer, std::size_t length, std::size_t offset, std::size_t bytes) const
{
    // skip reading invalid files
    if (not valid_)
    {
        error_ = "Invalid file.";
        return false;
    }
    
    // skip reading empty datasets
    if (length == 0)
        return true;
    
    // skip reading into null buffer
    if (buffer == nullptr)
    {
        error_ = "Null buffer.";
        return false;
    }
    
    // compose full dataset path
    dataset = prefix + dataset;
    
    // try to find the dataset with this name
    int id = -1, seek = header_byte_size;
    for (std::size_t i = 0; i < datasets_.size(); i++)
    {
        if (datasets_[i].name == dataset)
        {
            id = i;
            break;
        }
        else
        {
            seek += datasets_[i].elements * datasets_[i].bytes;
        }
    }
    
    // dataset not found
    if (id == -1)
    {
        error_ = format("Can't read dataset \"%s\" in file \"%s\".", dataset.c_str(), name_.c_str());
        return false;
    }
    
    // check that the selection fits into the dataset
    if ((offset + length) * bytes > datasets_[id].elements * datasets_[id].bytes)
    {
        error_ = format("Selection does not fit into the dataset ((%ld + %ld) * %ld >= %ld * %ld).", offset, length, bytes, datasets_[id].elements, datasets_[id].bytes);
        return false;
    }
    
    // reposition stream according to the selection
    if (std::fseek(file_, seek + offset * bytes, SEEK_SET) != 0)
    {
        error_ = "Error while reading dataset: seek failure.";
        return false;
    }
    
    // read data
    std::size_t read = std::fread(buffer, bytes, length, file_);
    if (read < length)
    {
        error_ = format("Error while reading dataset: read failure (read %ld of %ld).", read, length);
        return false;
    }
    
    return true;
}

bool HDFFile::write_ (std::string dataset, void const * buffer, std::size_t length, std::size_t offset, std::size_t bytes)
{
    // skip writing into invalid files
    if (not valid_)
    {
        error_ = "Invalid file.";
        return false;
    }
    
    // skip writing empty datasets
    if (length == 0)
        return true;
    
    // compose full dataset path
    dataset = prefix + dataset;
    
    // try to find the dataset with this name
    int id = -1, seek = header_byte_size;
    for (std::size_t i = 0; i < datasets_.size(); i++)
    {
        if (datasets_[i].name == dataset)
        {
            id = i;
            break;
        }
        else
        {
            seek += datasets_[i].elements * datasets_[i].bytes;
        }
    }
    
    // create a new dataset
    if (id < 0)
    {
        // create a new dataset
        DatasetInfo dst;
        dst.name = dataset;
        dst.bytes = bytes;
        dst.elements = length + offset;
        datasets_.push_back(dst);
        id = datasets_.size() - 1;
        changed_ = true;
    }
    
    // check that the selection fits into the dataset
    if ((offset + length) * bytes > datasets_[id].elements * datasets_[id].bytes)
    {
        error_ = format("Selection does not fit into the dataset ((%ld + %ld) * %ld >= %ld * %ld).", offset, length, bytes, datasets_[id].elements, datasets_[id].bytes);
        return false;
    }
    
    // exit if nothing to write
    if (buffer == nullptr)
        return true;
    
    // reposition stream according to the selection
    if (std::fseek(file_, seek + offset * bytes, SEEK_SET) != 0)
    {
        error_ = "Error while writing dataset.";
        return false;
    }
    
    // write data
    std::size_t written = std::fwrite(buffer, bytes, length, file_);
    if (written < length)
    {
        error_ = format("Error while writing dataset (written %ld of %ld).", written, length);
        perror("error");
        return false;
    }
    
    return true;
}
