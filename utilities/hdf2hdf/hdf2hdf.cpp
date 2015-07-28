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

#include "../../src/arrays.cpp"
#include "../../src/hdffile.cpp"
#include "../../src/h5file.cpp"
#include "../../src/misc.cpp"

std::string get_extension (std::string filename)
{
    // skip empty filenames
    if (filename.size() == 0)
        return filename;
    
    // find last dot
    std::size_t n = filename.size() - 1;
    for (; n > 0; n--)
        if (filename[n] == '.')
            break;
    
    // return substring beyond the last dot
    return filename.substr(n + 1, std::string::npos);
}

int convert_hdf_to_h5 (std::string srcfile, std::string dstfile)
{
    // open source file for reading
    HDFFile hdf (srcfile.c_str(), HDFFile::readonly);
    if (not hdf.valid())
    {
        std::cerr << "Cannot open file \"" << srcfile << "\"" << std::endl;
        return EXIT_FAILURE;
    }
    
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
        int * buffer = new int [(di.datatype % 0x100) * di.elements / 4];
        if (not hdf.read(di.name, buffer, di.elements * (di.datatype % 0x100) / 4))
        {
            std::cerr << "Failed to read dataset \"" << di.name << "\": " << hdf.error() << std::endl;
            return EXIT_FAILURE;
        }
        
        // write array with correct data type
        switch (di.datatype)
        {
            case DataInt32:
                h5f.write(di.name, (int*)buffer, di.elements);
                break;
            case DataInt64:
                h5f.write(di.name, (std::int64_t*)buffer, di.elements);
                break;
            case DataUInt32:
                h5f.write(di.name, (unsigned*)buffer, di.elements);
                break;
            case DataUInt64:
                h5f.write(di.name, (std::uint64_t*)buffer, di.elements);
                break;
            case DataFloat32:
                h5f.write(di.name, (float*)buffer, di.elements);
                break;
            case DataDouble64:
                h5f.write(di.name, (double*)buffer, di.elements);
                break;
            case DataDouble80:
                h5f.write(di.name, (long double*)buffer, di.elements);
                break;
            default:
                std::cerr << "Unknown data type of byte size " << (di.datatype % 0x100) << " (" << (di.datatype % 0x100) * 8 << " bits)." << std::endl;
                break;
        };
        
        // release buffer
        delete [] buffer;
    }
    
    return EXIT_SUCCESS;
}

herr_t write_dataset_h5_to_hdf (hid_t id, const char * name, const H5L_info_t * info, void * data)
{
    H5O_info_t infobuf;
    if (H5Oget_info_by_name(id, name, &infobuf, H5P_DEFAULT) < 0)
    {
        std::cerr << "Failed to determine type of object \"" << name << "\"" << std::endl;
        return -1;
    }
    
    if (infobuf.type == H5O_TYPE_DATASET)
        ((std::vector<std::string>*)data)->push_back(name);
    
    return 0;
}

int convert_h5_to_hdf (std::string srcfile, std::string dstfile)
{
    // open source file for reading
    H5File h5f (srcfile.c_str(), H5File::readonly);
    if (not h5f.valid())
    {
        std::cerr << "Cannot open file \"" << srcfile << "\"" << std::endl;
        return EXIT_FAILURE;
    }
    
    // open destination file for writing
    HDFFile hdf (dstfile.c_str(), HDFFile::overwrite);
    if (not hdf.valid())
    {
        std::cerr << "Cannot open file \"" << dstfile << "\"" << std::endl;
        return EXIT_FAILURE;
    }
    
    // get all datasets
    std::vector<std::string> datasets;
    if (H5Literate(h5f.file(), H5_INDEX_NAME, H5_ITER_NATIVE, nullptr, &write_dataset_h5_to_hdf, &datasets) < 0)
    {
        std::cerr << "Iteration ofer datasets failed." << std::endl;
        return EXIT_FAILURE;
    }
    
    // copy all datasets
    for (std::string const & dataset : datasets)
    {
        hid_t dset = H5Dopen(h5f.file(), dataset.c_str(), H5P_DEFAULT);
        hid_t dtyp = H5Dget_type(dset);
        hid_t dspc = H5Dget_space(dset);
        
        hssize_t nElem = H5Sget_simple_extent_npoints(dspc);
        
        if (H5Tequal(dtyp, H5T_NATIVE_INT) > 0)
        {
            iArray arr(nElem);
            h5f.read(dataset, &arr[0], nElem);
            hdf.write(dataset, &arr[0], nElem);
        }
        else if (H5Tequal(dtyp, H5T_NATIVE_INT64) > 0)
        {
            lArray arr(nElem);
            h5f.read(dataset, &arr[0], nElem);
            hdf.write(dataset, &arr[0], nElem);
        }
        else if (H5Tequal(dtyp, H5T_NATIVE_DOUBLE) > 0)
        {
            rArray arr(nElem);
            h5f.read(dataset, &arr[0], nElem);
            hdf.write(dataset, &arr[0], nElem);
        }
        else
        {
            std::cerr << "Unknown data type of dataset \"" << dataset << "\"" << std::endl;
        }
        
        H5Sclose(dspc);
        H5Tclose(dtyp);
        H5Dclose(dset);
    }
    
    return EXIT_SUCCESS;
}

int main (int argc, char *argv[])
{
    // check number of arguments
    if (argc < 3)
    {
        std::cout << std::endl;
        std::cout << "Converts files between the Hex Data Format (*.hdf) and the Hierarchical Data Format (*.h5)." << std::endl;
        std::cout << std::endl;
        std::cout << "Use" << std::endl;
        std::cout << "\thex-hdf2hdf <srcfile>.hdf <dstfile>.h5" << std::endl;
        std::cout << "or" << std::endl;
        std::cout << "\thex-hdf2hdf <srcfile>.h5  <dstfile>.hdf" << std::endl;
        std::cout << std::endl;
        std::cout << "Please note that the Hex data format is flat, and so only flat HDF5 files will be converted successfully." << std::endl;
        std::cout << "Hex data format also lacks type information. When converting from Hex data format to HDF5, user will be prompted for data types." << std::endl;
        std::cout << std::endl;
        return argc == 1 ? EXIT_FAILURE : EXIT_SUCCESS;
    }
    
    // source file name
    std::string srcfile(argv[1]);
    std::string srcext = get_extension(srcfile);
    
    // destination file name
    std::string dstfile(argv[2]);
    std::string dstext = get_extension(dstfile);
    
    // mute HDF library
    H5Eset_auto2(0, 0, nullptr);
    
    // pick the conversion direction
    if (srcext == std::string("hdf") and dstext == std::string("h5"))
        return convert_hdf_to_h5(srcfile, dstfile);
    if (srcext == std::string("h5") and dstext == std::string("hdf"))
        return convert_h5_to_hdf(srcfile, dstfile);
    
    std::cout << "The files' format is the same." << std::endl;
    return EXIT_FAILURE;
}
