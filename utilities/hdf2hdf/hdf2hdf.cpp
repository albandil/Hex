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
#include <vector>
#include <string>

#include "../../src/arrays.h"
#include "../../src/arrays.cpp"
#include "../../src/hdffile.cpp"

#include "h5file.cpp"

int main (int argc, char *argv[])
{
    // check number of arguments
    if (argc == 1 or std::string(argv[0]) == std::string("-h") or std::string(argv[0]) == std::string("--help"))
    {
        std::cout << std::endl;
        std::cout << "Converts files from Hex Data Format (HDF) to Hierarchical Data Format (H5)." << std::endl;
        std::cout << std::endl;
        std::cout << "Usage:" << std::endl;
        std::cout << "\thex-hdf2hdf <srcfile> [<dstfile>]" << std::endl;
        std::cout << std::endl;
        std::cout << "If the destination file is not given, the basename of source file will be used with \".h5\" extension." << std::endl;
        std::cout << std::endl;
        return argc == 1 ? EXIT_FAILURE : EXIT_SUCCESS;
    }
    
    // source file name
    std::string srcfile(argv[1]);
    
    // destination file name
    std::string dstfile;
    if (argc > 2)
    {
        dstfile = argv[2];
    }
    else
    {
        dstfile = argv[1];
        
        // erase original extension
        while (not dstfile.empty() and dstfile.back() != '.')
            dstfile.pop_back();
        
        // append new extension
        dstfile.push_back('h');
        dstfile.push_back('5');
    }
    
    // open source file for reading
    HDFFile hdf (srcfile.c_str(), HDFFile::readonly);
    if (not hdf.valid())
    {
        std::cerr << "Cannot open file \"" << srcfile << "\"" << std::endl;
        return EXIT_FAILURE;
    }
    
    // mute HDF library
    H5Eset_auto2(0, 0, nullptr);
    
    // open destination file for writing
    H5File h5f (dstfile.c_str(), H5File::overwrite);
    if (not h5f.valid())
    {
        std::cerr << "Cannot open file \"" << dstfile << "\"" << std::endl;
        return EXIT_FAILURE;
    }
    
    // for all datasets
    for (HDFFile::DatasetInfo const & di : hdf.datasets())
    {
        // read dataset into memory
        int * buffer = new int [di.bytes * di.elements / 4];
        if (not hdf.read(di.name, buffer, di.elements * di.bytes / 4))
        {
            std::cerr << "Failed to read dataset \"" << di.name << "\": " << hdf.error() << std::endl;
            return EXIT_FAILURE;
        }
        
        // Hex data files lack type information, we may need a hint from the user
        char choice;
        switch (di.bytes)
        {
            case 4:
                std::cout << std::endl;
                std::cout << "Dataset \"" << di.name << "\" has byte size 4 (32 bits). The data type is" << std::endl;
                std::cout << "a) int32" << std::endl;
                std::cout << "b) uint32" << std::endl;
                std::cout << "c) float" << std::endl << std::endl;
                while (1)
                {
                    std::cout << "a/b/c ? " << std::flush;
                    std::cin >> choice;
                    if (choice == 'a') { h5f.write(di.name, (std::int32_t*)buffer, di.elements);    break; }
                    if (choice == 'b') { h5f.write(di.name, (std::uint32_t*)buffer, di.elements);   break; }
                    if (choice == 'c') { h5f.write(di.name, (float*)buffer, di.elements);           break; }
                }
                break;
            
            case 8:
                std::cout << std::endl;
                std::cout << "Dataset \"" << di.name << "\" has byte size 8 (64 bits). The data type is" << std::endl;
                std::cout << "a) int64" << std::endl;
                std::cout << "b) uint64" << std::endl;
                std::cout << "c) double" << std::endl << std::endl;
                while (1)
                {
                    std::cout << "a/b/c ? " << std::flush;
                    std::cin >> choice;
                    if (choice == 'a') { h5f.write(di.name, (std::int64_t*)buffer, di.elements);    break; }
                    if (choice == 'b') { h5f.write(di.name, (std::uint64_t*)buffer, di.elements);   break; }
                    if (choice == 'c') { h5f.write(di.name, (double*)buffer, di.elements);          break; }
                }
                break;
            
            default:
                std::cerr << "Unknown data type of byte size " << di.bytes << " (" << di.bytes * 8 << " bits)." << std::endl;
                break;
        }
        
        // release buffer
        delete [] buffer;
    }
    
    std::cout << std::endl;
    
    return EXIT_SUCCESS;
}
