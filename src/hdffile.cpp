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
#include <fstream>

#include "hdffile.h"
#include "misc.h"

const std::size_t HDFFile::header_byte_size = 1024;

HDFFile::HDFFile (std::string filename, FileAccess flag)
    : prefix(), file_(nullptr), name_(filename), changed_(false)
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
            file_.open(name_.c_str(), std::ios::in | std::ios::out | std::ios::binary | std::ios::trunc);
            break;
            
        case readonly:
            file_.open(name_.c_str(), std::ios::in | std::ios::binary);
            break;
            
        case readwrite:
            file_.open(name_.c_str(), std::ios::in | std::ios::out | std::ios::binary);
            break;
    }
    
    // read header
    if (file_.is_open() and (flag == readonly or flag == readwrite))
    {
        // read number of datasets
        std::size_t ndt;
        file_.read(reinterpret_cast<char*>(&ndt), sizeof(std::size_t));
        
        // exit if there are no data in the header
        if (file_.fail())
            return;
        
        // read dataset information
        datasets_.resize(ndt);
        for (std::size_t idt = 0; idt < ndt; idt++)
        {
            // read length of the dataset name
            std::size_t len;
            file_.read(reinterpret_cast<char*>(&len), sizeof(std::size_t));
            
            // read dataset name
            datasets_[idt].name.resize(len);
            file_.read(&(datasets_[idt].name[0]), len);
            
            // read datatype byte size and element count
            file_.read(reinterpret_cast<char*>(&(datasets_[idt].bytes)), sizeof(std::size_t));
            file_.read(reinterpret_cast<char*>(&(datasets_[idt].elements)), sizeof(std::size_t));
        }
    }
}

HDFFile::~HDFFile ()
{
    if (file_.is_open())
    {
        // write the modified header
        if (changed_)
        {
            // write number of datasets
            std::size_t ndt = datasets_.size();
            file_.seekp(file_.beg);
            file_.write(reinterpret_cast<char*>(&ndt), sizeof(ndt));
            
            // write datasets information
            for (std::size_t idt = 0; idt < ndt; idt++)
            {
                // write length of dataset name
                std::size_t len = datasets_[idt].name.size();
                file_.write(reinterpret_cast<char*>(&len), sizeof(ndt));
                
                // write dataset name
                file_.write(&(datasets_[idt].name[0]), len);
                
                // write datatype byte size and element count
                file_.write(reinterpret_cast<char*>(&(datasets_[idt].bytes)), sizeof(std::size_t));
                file_.write(reinterpret_cast<char*>(&(datasets_[idt].elements)), sizeof(std::size_t));
            }
        }
        file_.close();
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
    if (not file_.is_open())
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
    int id = -1;
    std::size_t seek = header_byte_size;
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
        std::cout << "Can't read dataset " << dataset << " in " << name_ << std::endl;
        return false;
    }
    
    // check that the selection fits into the dataset
    if ((offset + length) * bytes > datasets_[id].elements * datasets_[id].bytes)
    {
        error_ = format("Selection does not fit into the dataset ((%ld + %ld) * %ld >= %ld * %ld).", offset, length, bytes, datasets_[id].elements, datasets_[id].bytes);
        return false;
    }
    
    // reposition stream according to the selection
    file_.seekg(seek + offset * bytes, std::ios_base::beg);
    if (file_.bad())
    {
        std::cerr << "Error while reading dataset (seek failure)." << std::endl;
        return false;
    }
    
    // read data
    file_.read(reinterpret_cast<char*>(buffer), bytes * length);
    if (file_.fail())
    {
        error_ = format("Error while reading dataset: read failure.");
        return false;
    }
    
    return true;
}

bool HDFFile::write_ (std::string dataset, void const * buffer, std::size_t length, std::size_t offset, std::size_t bytes)
{
    // skip writing into invalid files
    if (not file_.is_open())
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
    int id = -1;
    std::size_t seek = header_byte_size;
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
    file_.seekp(seek + offset * bytes, std::ios_base::beg);
    if (file_.fail())
    {
        error_ = format("Error while writing dataset (seek failure).");
        return false;
    }
    
    // write data
    file_.write(reinterpret_cast<const char*>(buffer), bytes * length);
    if (file_.bad())
    {
        error_ = format("Error while writing dataset.");
        return false;
    }
    
    return true;
}
