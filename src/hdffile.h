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

#include <string>
#include <fstream>

#include "numbers.h"
#include "misc.h"

/**
 * @brief Hex data file.
 * 
 * This class offers an interface for accessing Hex data files. Hex data files contain
 * a fixed-size header with list of datasets, which is followed by the datasets themselves.
 * Concurrent write access to the same file is not controlled and will generally result
 * in broken data. Concurrent read access, on the contrary, is allowed and should not
 * have any side effects.
 * 
 * The typical usage is
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
//         failifexists,
        readwrite
    }
    FileAccess;
    
    /// Header size.
    static const std::size_t header_byte_size;
    
    /// Constructor from filename and access mode.
    HDFFile (std::string filename, FileAccess flag);
    
    /// Destructor.
    ~HDFFile ();
    
    /// Get size of the dataset of valid file.
    std::size_t size (std::string dataset) const;
    
    /// Load data from a valid file.
    template <typename T> bool read (std::string dataset, T * buffer, std::size_t length, std::size_t offset = 0) const
    {
        return read_(dataset, buffer, length * typeinfo<T>::ncmpt, offset * typeinfo<T>::ncmpt, typeinfo<T>::hdfcmpttype());
    }
    
    /// Write data to a valid file.
    template <typename T> bool write (std::string dataset, T const * buffer, std::size_t length, std::size_t offset = 0)
    {
        return write_(dataset, buffer, length * typeinfo<T>::ncmpt, offset * typeinfo<T>::ncmpt, typeinfo<T>::hdfcmpttype());
    }
    
    /// Check that the file is valid.
    bool valid () const { return file_.is_open(); }
    
    /// Prefix for dataset paths.
    std::string prefix;
    
    /// Error stack.
    std::string const & error () const { return error_; };
    
    typedef struct
    {
        /// Name of the dataset.
        std::string name;
        
        /// Data type of the elements.
        std::size_t datatype;
        
        /// Number of elements.
        std::size_t elements;
    }
    DatasetInfo;
    
    /// File name.
    std::string const & name () const { return name_; }
    
    /// Dataset information.
    std::vector<DatasetInfo> const & datasets () const { return datasets_; }
    
private:
    
    /// File stream.
    mutable std::fstream file_;
    
    /// Header
    std::vector<DatasetInfo> datasets_;
    
    /// Filename.
    std::string name_;
    
    /// Error message.
    mutable std::string error_;
    
    /// Whether the dataset layout has been changed.
    bool changed_;
    
    /// Auxiliary read function.
    bool read_(std::string dataset, void * buffer, std::size_t length, std::size_t offset, HDFDataType dtype) const;
    
    /// Auxiliary write function.
    bool write_(std::string dataset, void const * buffer, std::size_t length, std::size_t offset, HDFDataType dtype);
    
    /// Set error and return false.
    void save_error () const;
};

#endif
