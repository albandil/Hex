/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2014                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

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
    typedef enum {
        readonly,
        overwrite,
        failifexists,
        readwrite
    } FileAccess;
    
    /// Constructor from filename and access mode.
    HDFFile(std::string filename, FileAccess flag);
    
    /// Destructor.
    ~HDFFile();
    
    /// Get size of the dataset of valid file.
    size_t size (std::string dataset) const;
    
    /// load data from a valid file.
    template <typename T> bool read (std::string dataset, T * buffer, size_t length) const;
    
    /// Write data to a valid file.
    template <typename T> bool write (std::string dataset, T const * data, size_t length);
    
    /// Check that the file is valid.
    bool valid () const { return valid_; }
    
private:
    
    /// Pointer to the HDF structure.
    H5::H5File * file_;
    
    /// Whether the file is valid.
    bool valid_;
    
    /// Auxiliary read function.
    bool read_(std::string dataset, void * buffer, hsize_t length, H5::AtomType dtype) const;
    
    /// Auxiliary write function.
    bool write_(std::string dataset, void const * buffer, hsize_t length, H5::AtomType dtype);
};

#endif

#endif
