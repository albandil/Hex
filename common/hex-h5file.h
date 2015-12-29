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

#ifndef HEX_H5FILE
#define HEX_H5FILE

#include <string>

#include <hdf5.h>

#include "hex-numbers.h"
#include "hex-misc.h"

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
 *     H5File hdffile("test.hdf", readonly);
 *     
 *     if (hdffile.valid())
 *         hdffile.write("nums", nums.data(), nums.size());
 *     else
 *         throw exception ("Cannot write to test.hdf!");
 *     
 * @endcode
 */
class H5File
{
public:
    
    /**
     * @brief Access mode.
     * 
     * This flag is used in the constructor of H5File to specify
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
    H5File (std::string filename, FileAccess flag);
    
    /// Destructor.
    ~H5File ();
    
    /// Get size of the dataset of valid file.
    std::size_t size (std::string dataset) const;
    
    /// Load data from a valid file.
    template <typename T> bool read (std::string dataset, T * buffer, std::size_t length, std::size_t offset = 0) const
    {
        HexException("Unsupported data type %s.", typeid(T).name());
    }
    
    /// Write data to a valid file.
    template <typename T> bool write (std::string dataset, T const * buffer, std::size_t length, std::size_t offset = 0)
    {
        HexException("Unsupported data type %s.", typeid(T).name());
    }
    
    /// Check that the file is valid.
    bool valid () const { return valid_; }
    
    /// Prefix for dataset paths.
    std::string prefix;
    
    /// Internal HDF object.
    hid_t file () const { return file_; }
    
    /// Error stack.
    std::string const & error () const { return error_; };
    
private:
    
    /// Pointer to the HDF structure.
    hid_t file_;
    
    /// Filename.
    std::string name_;
    
    /// HDF library error stack (set on unsuccessful return).
    mutable std::string error_;
    
    /// Whether the file is valid.
    bool valid_;
    
    /// Auxiliary read function.
    bool read_(std::string dataset, void * buffer, hsize_t length, hsize_t offset, hid_t dtype) const;
    
    /// Auxiliary write function.
    bool write_(std::string dataset, void const * buffer, hsize_t length, hsize_t offset, hid_t dtype);
    
    /// Set error and return false.
    void save_error () const;
};

template<> bool H5File::read<int> (std::string dataset, int * buffer, std::size_t length, std::size_t offset) const;
template<> bool H5File::write<int> (std::string dataset, int const * buffer, std::size_t length, std::size_t offset);

template<> bool H5File::read<unsigned> (std::string dataset, unsigned * buffer, std::size_t length, std::size_t offset) const;
template<> bool H5File::write<unsigned> (std::string dataset, unsigned const * buffer, std::size_t length, std::size_t offset);

template<> bool H5File::read<std::int64_t> (std::string dataset, std::int64_t * buffer, std::size_t length, std::size_t offset) const;
template<> bool H5File::write<std::int64_t> (std::string dataset, std::int64_t const * buffer, std::size_t length, std::size_t offset);

template<> bool H5File::read<std::uint64_t> (std::string dataset, std::uint64_t * buffer, std::size_t length, std::size_t offset) const;
template<> bool H5File::write<std::uint64_t> (std::string dataset, std::uint64_t const * buffer, std::size_t length, std::size_t offset);

template<> bool H5File::read<double> (std::string dataset, double * buffer, std::size_t length, std::size_t offset) const;
template<> bool H5File::write<double> (std::string dataset, double const * buffer, std::size_t length, std::size_t offset);

template<> bool H5File::read<Complex> (std::string dataset, Complex * buffer, std::size_t length, std::size_t offset) const;
template<> bool H5File::write<Complex> (std::string dataset, Complex const * buffer, std::size_t length, std::size_t offset);


#endif
