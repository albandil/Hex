//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2017, Jakub Benda, Charles University in Prague                    //
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

#ifndef HEX_MEMORY
#define HEX_MEMORY

// --------------------------------------------------------------------------------- //

#include <iostream>
#include <climits>
#include <cstdint>
#include <cstddef>
#include <cstdlib>

// --------------------------------------------------------------------------------- //

#include "hex-misc.h"

// --------------------------------------------------------------------------------- //

// number array alignment (256 bits ~ AVX2)
#define SIMD_VECTOR_BITS 256u
#define SIMD_VECTOR_BYTES (SIMD_VECTOR_BITS / 8u)

// --------------------------------------------------------------------------------- //

// page size
#define PAGE_SIZE 4096

// --------------------------------------------------------------------------------- //

// number of components of vector that fits into SIMD vector type
#define simd_single_vec_size (SIMD_VECTOR_BYTES / sizeof(float))
#define simd_double_vec_size (SIMD_VECTOR_BYTES / sizeof(double))
#ifdef SINGLE
    #define simd_real_vec_size simd_single_vec_size
#else
    #define simd_real_vec_size simd_double_vec_size
#endif

// --------------------------------------------------------------------------------- //

// aligned pointers
#if (defined(__GNUC__) && !defined(__clang__))
    #define assume_aligned(x,a) __builtin_assume_aligned((x),(a))
#else
    #define assume_aligned(x,a) (x)
#endif

// --------------------------------------------------------------------------------- //

// restricted pointers
#ifdef __GNUC__
    #define restrict __restrict
#else
    #define restrict
#endif

// --------------------------------------------------------------------------------- //

// low-level file and memory functions
#ifdef __linux__
    #include <sys/mman.h>
    #include <sys/stat.h>
    
    #include <fcntl.h>
    #include <unistd.h>
#endif

// --------------------------------------------------------------------------------- //

/**
 * @brief Basic memory allocator.
 * 
 * Memory allocators are used by @ref Array -like classes for memory
 * management: allocating and freeing of memory. PlainAllocator uses
 * the usual C++ "new/delete" mechanism for manipulation of the memory.
 */
template <class T> class PlainAllocator
{
    public:
        
        /**
         * @brief Get alignment of a type.
         * 
         * Return alignment of the allocated array in bytes.
         */
        static std::size_t align ()
        {
            return std::alignment_of<T>::value;
        }
        
        /**
         * @brief Allocate array.
         * 
         * This method will allocate fields of the data type T.
         * The elements will be constructed by the implicit constructor.
         * @param n Number of items to allocate.
         */
        static T * alloc (std::size_t n)
        {
            try
            {
                return new T[n]();
            }
            catch (std::bad_alloc const & err)
            {
                std::size_t bytes = n * sizeof(T);
                std::string strsize = (
                    bytes < 1204           ? format("%d B", bytes) : (
                    bytes < 1024*1204      ? format("%.2f kiB", bytes / 1024.) : (
                    bytes < 1024*1024*1024 ? format("%.2f MiB", bytes / (1024. * 1024.)) :
                      /* else */             format("%.2f GiB", bytes / (1024. * 1024 * 1024.)))));
                
                HexException("Insufficent memory (unable to allocate next %s).", strsize.c_str());
            }
            
            return nullptr;
        }
        
        /**
         * @brief Free array.
         * 
         * This method will deallocate memory pointed to by the pointer
         * "ptr". It is expected that pointer has been previously
         * retrieved from the function PlainAllocator::alloc. The pointer
         * can have the value "nullptr"; in such a case nothing will
         * be done.
         */
        static void free (T * ptr)
        {
            if (ptr != nullptr)
                delete [] ptr;
        }
};

// --------------------------------------------------------------------------------- //

/**
 * @brief Aligned memory allocator.
 * 
 * Memory allocators are used by @ref Array -like classes for memory
 * management: allocating and freeing of memory. This particular class
 * allocates aligned memory, i.e. every allocated block will begin at
 * a memory address that is multiple of the size of the allocated type.
 * Also, all items of the allocated array will be placed in memory-aligned
 * locations. The alignment helps the system load the memory efficiently
 * and the compiler can thus do some optimizations -- particularly the
 * autovectorization, which can utilize SIMD instruction of the CPU and
 * speed up the operation on the arrays.
 */
template <class T, std::size_t alignment_ = std::alignment_of<T>::value> class AlignedAllocator
{
    public:
        
        static const std::size_t alignment = alignment_;
        
        /**
         * @brief Allocate aligned memory.
         * 
         * Allocate array of "n" items of the type "T". The returned pointer address
         * will be a multiple of the alignment given as the template argument of the
         * class. There are several restriction on value of the alignment:
         * - The smallest allowed alignment is sizeof(max_align_t). If the template
         *   argument is smalled, the value sizeof(max_align_t) will be used instead.
         * - The alignment has to be an integer multiple of sizeof(void*). If this is
         *   not the case, function will raise the requested alignment to conform to
         *   this rule.
         * 
         * Here is a graphical representation of the memory:
           <pre>
           root  origin
              v       v
              |  |xxxx|yyyy|yyyy|yyyy|yyyy|
            ->   <-   <\-\-\-\-\-\- data \-\-\-\-\-\-\->
            padding
                ->    <-
              root address
           </pre>
         */
        static T * alloc (std::size_t n)
        {
            // is there anything to allocate?
            if (n < 1)
                return nullptr;
            
            // get the alignment and raise it to nearest multiple of sizeof(void*)
            std::size_t align = std::max(alignment_, alignof(max_align_t));
            if (align % sizeof(void*) != 0)
                align += sizeof(void*) - align % sizeof(void*);
            
            // count elements (make the count even) and compute upper memory bound
            std::size_t elems = n + (n % 2);
            std::size_t bytes = elems * sizeof(T) + sizeof(std::uintptr_t) + align;
            
            // allocate memory pool
            void * root;
            try
            {
                root = new char [bytes];
            }
            catch (std::bad_alloc const & err)
            {
                std::string strsize = (
                    bytes < 1204           ? format("%d B", bytes) : (
                    bytes < 1024*1204      ? format("%.2f kiB", bytes / 1024.) : (
                    bytes < 1024*1024*1024 ? format("%.2f MiB", bytes / (1024. * 1024.)) :
                      /* else */             format("%.2f GiB", bytes / (1024. * 1024 * 1024.)))));
                
                HexException("Insufficent memory (unable to allocate next %s).", strsize.c_str());
            }
            std::uintptr_t root_address = (std::uintptr_t)root;
            
            // compute padding
            unsigned padding = 0;
            if ((root_address + sizeof(root_address)) % align != 0)
                padding = align - (root_address + sizeof(root_address)) % align;
            
            // place origin
            std::uintptr_t origin_address = root_address + sizeof(root_address) + padding;
            T * origin = reinterpret_cast<T*>(origin_address);
            
            // store the address of 'root' just before the 'origin'
            *(reinterpret_cast<std::uintptr_t*>(origin_address) - 1) = root_address;
            
            // clear the last element so that we may disregard it during multiplication
            *(origin + n + (n % 2) - 1) = 0;
            
            // return the pointer
            return origin;
        }
        
        /**
         * @brief Deallocate aligned memory.
         * 
         * Free the memory pointed to by the pointer "ptr". It is assumed
         * that the pointer was obtained from the call to the function
         * alloc of the same specialization of the class AlignedAllocator.
         * Value "nullptr" is allowed; nothing will be done in that case.
         */
        static void free (T * origin)
        {
            // only free non-null pointers
            if (origin != nullptr)
            {
                // get root address from the location just before the pointer
                std::uintptr_t root_address = *(reinterpret_cast<std::uintptr_t*>(origin) - 1);
                
                // free the whole memory pool
                delete [] reinterpret_cast<char*>(root_address);
            }
        }
        
        /**
         * @brief Test aligned pointer.
         * 
         * This routine will print some diagnostic information about a pointer
         * that is assumed to have been allocated by AlignedAllocator::alloc.
         * The information is printed to the standard output.
         */
        static void test (T const * origin)
        {
            // only test non-null pointes
            if (origin == nullptr)
                return;
            
            // get root and origin addresses
            std::uintptr_t root_address = *(reinterpret_cast<std::uintptr_t const*>(origin) - 1);
            std::uintptr_t origin_address = reinterpret_cast<std::uintptr_t>(origin);
            
            // get distance between 'root' and 'origin'; also get the root pointer
            unsigned d = origin_address - root_address;
            void * root = reinterpret_cast<void*>(root_address);
            
            // compute maximal usable alignment of this aligned pointer
            std::size_t align = 1;
            for (; ; align *= 2)
            {
                if (origin_address % (2*align) != 0)
                    break;
            }
            
            // write information to stdout
            std::cout << "Aligned pointer information:" << std::endl;
            std::cout << "   root pointer : " << root << std::endl;
            std::cout << "   origin       : " << origin << std::endl;
            std::cout << "   origin-root  : " << d << " Bytes" << std::endl;
            std::cout << "   padding      : " << d - sizeof(root_address) << " Bytes" << std::endl;
            std::cout << "   alignment    : " << align << " Bytes (default: " << alignof(max_align_t) << ")" << std::endl;
        }
};

// --------------------------------------------------------------------------------- //

/**
 * @brief Virtual memory allocator.
 * 
 * This class will create a memory mapped disk file each time an allocation is
 * requested. Such file does not occupy the physical memory and can be used to
 * offload large arrays. Still, if there is enough physical memory, it will be
 * cached in memory and the access times should not dramatically suffer.
 */
template <class T, std::size_t pagesize_ = PAGE_SIZE> class vMemAllocator
{
    public:
        
        static const std::size_t pagesize = pagesize_;
        
        /**
         * @brief Allocate virtual memory.
         * 
         * Create a temporary file with given prefix, in chosen location,
         * and map it into virtual memory space. The first three words
         * are reserved for internal usage (file descriptor, low word of
         * file size, high word of file size) and pointer to the rest is
         * returned to the used.
         */
        static T * alloc (std::size_t n, std::string tmpdir = ".", std::string prefix = "")
        {
#ifdef __linux__
            // create a unique temporary file with the given prefix
            std::string filename = tmpdir + "/" + prefix + "XXXXXX";
            int fd = mkstemp64(const_cast<char*>(filename.c_str()));
            
            if (fd == -1)
                HexException("Failed to create scratch file %s: %s", filename.c_str(), std::strerror(errno));
            
            // calculate number of pages the file will span, round up
            size_t npages = (n * sizeof(T) + 3*sizeof(int) + pagesize_ - 1) / pagesize_;
            std::size_t filesize = npages * pagesize_;
            
            std::cout << "allocate " << npages << " pages" << std::endl;
            
            // fill file with zeros
            char page[pagesize];
            std::memset(page, 0, pagesize);
            for (std::size_t ipage = 0; ipage < npages; ipage++)
            if (write(fd, page, pagesize) != pagesize)
                HexException("Cannot allocate %ld-th page in file %s: %s", ipage + 1, filename.c_str(), std::strerror(errno));
            
            // map the file to memory
            int* buffer = (int*)mmap64
            (
                nullptr,                // memory location where the mapping should start - auto
                filesize,               // size of the mapping
                PROT_READ | PROT_WRITE, // access permissions
                MAP_SHARED,             // allow immediate write to disk
                fd,                     // file desccriptor
                0                       // file offset
            );
            
            if (buffer == MAP_FAILED)
                HexException("Failed to map scratch file %s to memory: %s", filename.c_str(), std::strerror(errno));
            
            buffer[0] = fd;
            buffer[1] = filesize % INT_MAX;
            buffer[2] = filesize / INT_MAX;
            
            return (T*)(buffer + 3);
#else
            return nullptr;
#endif
        }
        
        static void free (T * origin)
        {
#ifdef __linux__
            if (origin != nullptr)
            {
                // shfit pointer to the beginning of the data array
                int* buffer = ((int*)origin) - 3;
                
                // read allocation data
                int fd = buffer[0];
                int lo = buffer[1];
                int hi = buffer[2];
                
                // calculate memory size
                size_t filesize = hi;
                filesize *= INT_MAX;
                filesize += lo;
                
                // get scratch filename
                char link[PATH_MAX], filename[PATH_MAX];
                snprintf(link, PATH_MAX, "/proc/%d/fd/%d", getpid(), fd);
                if (realpath(link, filename) == nullptr)
                    HexException("Failed to retrieve filename from %s: %s", link, std::strerror(errno));
                
                // safely unmap the disk file
                if (munmap(buffer, filesize) == -1)
                    HexException("Failed to deallocate virtual memory: %s", std::strerror(errno));
                
                // close the file
                close(fd);
                
                // remove the file
                std::remove(filename);
            }
#endif
        }
};

// --------------------------------------------------------------------------------- //

#endif /* HEX_MEMORY */
