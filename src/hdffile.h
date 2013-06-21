/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2013                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HEX_HDFFILE_H
#define HEX_HDFFILE_H

#include <string>
#include <H5Cpp.h>

#include "complex.h"

class HDFFile
{
public:
	
	typedef enum {
		readonly,
		overwrite,
		failifexists,
		readwrite
	} FileAccess;
	
	HDFFile(std::string filename, FileAccess flag);
	~HDFFile();
	
	// get size of the dataset
	size_t size (std::string dataset) const;
	
	// load data
	template <typename T> bool read (std::string dataset, T * buffer, size_t length) const;
	
	// write data
	template <typename T> bool write (std::string dataset, T const * data, size_t length);
	
	bool valid() const { return valid_; }
	
private:
	
	H5::H5File * file_;
	bool valid_;
	bool read_(std::string dataset, void * buffer, hsize_t length, H5::AtomType dtype) const;
	bool write_(std::string dataset, void const * buffer, hsize_t length, H5::AtomType dtype);
};

#endif
