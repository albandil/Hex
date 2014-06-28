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

#ifndef HEX_MEMORY
#define HEX_MEMORY

#include <cstdint>
#include <cstddef>

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
         * @brief Allocate array.
         * 
         * This method will allocate fields of the data type T.
         * The elements will be constructed by the implicit constructor.
         * @param n Number of items to allocate.
         */
        static T * alloc (size_t n)
        {
            return new T[n]();
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
template <class T, size_t alignment = std::alignment_of<T>::value> class AlignedAllocator
{
    public:
        
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
            std::size_t align = std::max(alignment, alignof(max_align_t));
            if (align % sizeof(void*) != 0)
                align += sizeof(void*) - align % sizeof(void*);
            
            // count elements (make the count even) and compute upper memory bound
            std::size_t elems = n + (n % 2);
            std::size_t bytes = elems * sizeof(T) + sizeof(std::uintptr_t) + align;
            
            // allocate memory pool
            void * root = std::malloc(bytes);
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
                delete [] reinterpret_cast<T*>(root_address);
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
            std::uintptr_t root_address = *(reinterpret_cast<std::uintptr_t*>(origin) - 1);
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
            std::cout << "   root pointer :   " << root << std::endl;
            std::cout << "   origin :         " << origin << std::endl;
            std::cout << "   origin to root : " << d << " bits" << std::endl;
            std::cout << "   padding :        " << d - sizeof(root_address) << " bits" << std::endl;
            std::cout << "   alignment :      " << align << " bits (default: " << alignof(max_align_t) << ")" << std::endl;
        }
};

#endif /* HEX_MEMORY */