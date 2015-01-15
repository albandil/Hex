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

#ifndef HEX_HDFFILE
#define HEX_HDFFILE

#ifndef NO_HDF

#include <string>
#include <H5Cpp.h>

#include "complex.h"

/**
 * @brief HDF I/O management.
 * 
 * This class abstracts the raw functions for work with HDF files, so that
 * template functions can be used and the user doesn't need to bother with
 * type specializations. The typical usage is
 * @code
 *     rArray nums = linspace (1., 10., 10);
 *     
 *     ...
 *     
 *     HDFFile hdffile("test.hdf", readonly);
 *     
 *     if (hdffile.valid())
 *         hdffile.write("nums", nums.data(), nums.size());
 *     else
 *         throw exception ("Cannot write to test.hdf!");
 *     
 * @endcode
 */
class HDFFile
{
public:
    
    /**
     * @brief Access mode.
     * 
     * This flag is used in the constructor of HDFFile to specify
     * access mode to the file (read only, overwrite etc.).
     */
    typedef enum
    {
//         create,
        readonly,
        overwrite,
        failifexists,
        readwrite
    }
    FileAccess;
    
    /// Constructor from filename and access mode.
    HDFFile(std::string filename, FileAccess flag);
    
    /// Destructor.
    ~HDFFile();
    
    /// Get size of the dataset of valid file.
    std::size_t size (std::string dataset) const;
    
    /// load data from a valid file.
    template <typename T> bool read (std::string dataset, T * buffer, std::size_t length, std::size_t offset = 0) const;
    
    /// Write data to a valid file.
    template <typename T> bool write (std::string dataset, T const * data, std::size_t length, std::size_t offset = 0);
    
    /// Check that the file is valid.
    bool valid () const { return valid_; }
    
    /// Prefix for dataset paths.
    std::string prefix;
    
    /// Internal HDF object.
    H5::H5File * file () { return file_; }
    
private:
    
    /// Pointer to the HDF structure.
    H5::H5File * file_;
    
    /// Filename.
    std::string name_;
    
    /// Whether the file is valid.
    bool valid_;
    
    /// Auxiliary read function.
    bool read_(std::string dataset, void * buffer, hsize_t length, hsize_t offset, H5::AtomType dtype) const;
    
    /// Auxiliary write function.
    bool write_(std::string dataset, void const * buffer, hsize_t length, hsize_t offset, H5::AtomType dtype);
};

#endif

#endif
