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

#ifndef HEX_OPENCL
#define HEX_OPENCL

#ifdef WITH_OPENCL

#include <CL/cl.h>

#include "arrays.h"

/**
 * @brief OpenCL array wrapper.
 *
 * This class equips ArrayView<T> with an OpenCL context handle and
 * provides several service methods for upload and download to / from
 * the compute device (mostly GPU).
 */
template <class T> class clArrayView : public ArrayView<T>
{
    protected:
        
        /// Memory handle to the data copy on the GPU.
        cl_mem cl_handle_;
        
    public:
        
        //
        // aliases
        //
        
        typedef T DataType;
        typedef T * iterator;
        typedef T const * const_iterator;
        
        //
        // basic constructors
        //
        
        clArrayView ()
            : ArrayView<T>(), cl_handle_(nullptr) { }
        clArrayView (std::size_t n, T const * ptr)
            : ArrayView<T>(n, const_cast<T*>(ptr)), cl_handle_(nullptr) { }
        clArrayView (ArrayView<T> v, std::size_t i = 0, std::size_t n = 0)
            : ArrayView<T>(v, i, n), cl_handle_(nullptr) { }
        clArrayView (clArrayView<T> const & v)
            : ArrayView<T>(v), cl_handle_(v.cl_handle_) { }
        clArrayView (ArrayView<T> && v)
            : ArrayView<T>(v), cl_handle_(nullptr) { }
        clArrayView (clArrayView<T> && v)
            : ArrayView<T>(), cl_handle_(nullptr)
        {
            std::swap (ArrayView<T>::array_, v.ArrayView<T>::array_);
            std::swap (ArrayView<T>::N_,     v.ArrayView<T>::N_    );
            std::swap (cl_handle_,           v.cl_handle_          );
        }
        
        //
        // destructor
        //
        
        virtual ~clArrayView ()
        {
        }
        
        //
        // STL interface (except changes in size -- we don't want to bother with GPU reallocation)
        //
        
        std::size_t size () const               { return ArrayView<T>::size();        }
        T const & operator [] (std::size_t i) const  { return ArrayView<T>::operator[](i); }
        T & operator [] (std::size_t i)              { return ArrayView<T>::operator[](i); }
        T const * data () const                 { return ArrayView<T>::data();        }
        T * data ()                             { return ArrayView<T>::data();        }
        const_iterator begin () const           { return ArrayView<T>::begin();       }
        iterator begin ()                       { return ArrayView<T>::begin();       }
        const_iterator end () const             { return ArrayView<T>::end();         }
        iterator end ()                         { return ArrayView<T>::end();         }
        T const & front (std::size_t i = 0) const    { return ArrayView<T>::front(i);      }
        T & front (std::size_t i = 0)                { return ArrayView<T>::front(i);      }
        T const & back (std::size_t i = 0) const     { return ArrayView<T>::back(i);       }
        T & back (std::size_t i = 0)                 { return ArrayView<T>::back(i);       }
        
        //
        // OpenCL intrinsics
        //
        
        bool is_connected () const
        {
            return cl_handle_ != nullptr;
        }
        
        void connect (cl_context context, cl_mem_flags flags)
        {
            // release previous allocation prior to creating a new buffer !
            if (cl_handle_ != nullptr)
                HexException("[clArray::connect] Array is already connected!");
            
            // allocate memory on GPU
            cl_handle_ = clCreateBuffer(context, flags, size() * sizeof(T), data(), nullptr);
        }
        
        void disconnect ()
        {
            // release previous allocations
            if (cl_handle_ != nullptr)
                clReleaseMemObject(cl_handle_);
            
            // clear pointer
            cl_handle_ = nullptr;
        }
        
        cl_int EnqueueUpload (cl_command_queue queue)
        {
            return clEnqueueWriteBuffer(queue, cl_handle_, CL_TRUE, 0, size() * sizeof(T), data(), 0, nullptr, nullptr);
        }
        
        cl_int EnqueueDownload (cl_command_queue queue)
        {
            return clEnqueueReadBuffer(queue, cl_handle_, CL_TRUE, 0, size() * sizeof(T), data(), 0, nullptr, nullptr);
        }
        
        cl_mem const & handle () const
        {
            return cl_handle_;
        }
};

/**
 * @brief OpenCL array wrapper.
 *
 * This class combines NumberArray<T> and clArrayView<T>. It is a self-standing array of
 * memory-aligned elements. The alignment is large by default so that the class can be used
 * as the host-type memory within OpenCL. However, the aimed usage is this:
 * - Create the array.
 * - Fill the array with initial data.
 * - Connect to the computing device (GPU) and upload the data.
 * - Run computation on the compute device.
 * - Download the resulting data back to RAM.
 * - Disconnect from compute device.
 * - Destruct the class.
 */
#ifndef NO_ALIGN
#define CL_ALIGNMENT 4096
template <class T, class Alloc = AlignedAllocator<T,CL_ALIGNMENT>> class clArray : public clArrayView<T>
#else
template <class T, class Alloc = PlainAllocator<T>> class clArray : public clArrayView<T>
#endif
{
    public:
        
        //
        // aliases
        //
        
        typedef T DataType;
        typedef T * iterator;
        typedef T const * const_iterator;
        
        //
        // basic constructors
        //
        
        clArray () : clArrayView<T>()
        {
            // nothing to do
        }
        clArray (std::size_t n, T x = T(0)) : clArrayView<T>(n, Alloc::alloc(n))
        {
            for (std::size_t i = 0; i < n; i++)
                (*this)[i] = x;
        }
        clArray (std::size_t n, T const * px) : clArrayView<T>(n, Alloc::alloc(n))
        {
            for (std::size_t i = 0; i < n; i++)
                (*this)[i] = px[i];
        }
        clArray (ArrayView<T> v, std::size_t i = 0, std::size_t n = 0)
            : clArrayView<T>(((n == 0) ? v.size() : n), Alloc::alloc((n == 0) ? v.size() : n))
        {
            for (std::size_t j = 0; j < size(); j++)
                (*this)[j] = v[i+j];
        }
        clArray (std::initializer_list<T> list) : clArrayView<T>(list.size(), Alloc::alloc(list.size()))
        {
            std::size_t i = 0;
            for (auto it = list.begin(); it != list.end(); it++)
                (*this)[i++] = *it;
        }
        clArray (ArrayView<T> && rvrf) : clArrayView<T>(std::move(rvrf))
        {
            // nothing to do
        }
        clArray (clArray<T,Alloc> && rvrf) : clArrayView<T>(std::move(rvrf))
        {
            // nothing to do
        }
        
        //
        // data assignment
        //
        
        clArray & operator = (const ArrayView<T> v)
        {
            // realloc memory if needed
            resize(v.size());
            
            // copy data
            for (std::size_t j = 0; j < size(); j++)
                (*this)[j] = v[j];
            
            return *this;
        }
        
        void resize (std::size_t newsize)
        {
            if (newsize != size())
            {
                if (size() != 0)
                    Alloc::free(data());
                
                ArrayView<T>::array_ = Alloc::alloc(newsize);
                ArrayView<T>::N_ = newsize;
            }
        }
        
        //
        // destructor
        //
        
        virtual ~clArray ()
        {
            // free GPU memory
            disconnect();
            
            // free RAM memory
            Alloc::free(data());
        }
        
        //
        // STL interface (except changes in size -- we don't want to bother with GPU reallocation)
        //
        
        std::size_t size () const                   { return clArrayView<T>::size();        }
        T const & operator [] (std::size_t i) const { return clArrayView<T>::operator[](i); }
        T & operator [] (std::size_t i)             { return clArrayView<T>::operator[](i); }
        T const * data () const                     { return clArrayView<T>::data();        }
        T * data ()                                 { return clArrayView<T>::data();        }
        const_iterator begin () const               { return clArrayView<T>::begin();       }
        iterator begin ()                           { return clArrayView<T>::begin();       }
        const_iterator end () const                 { return clArrayView<T>::end();         }
        iterator end ()                             { return clArrayView<T>::end();         }
        T const & front (std::size_t i = 0) const   { return clArrayView<T>::front(i);      }
        T & front (std::size_t i = 0)               { return clArrayView<T>::front(i);      }
        T const & back (std::size_t i = 0) const    { return clArrayView<T>::back(i);       }
        T & back (std::size_t i = 0)                { return clArrayView<T>::back(i);       }
        
        //
        // OpenCL intrinsics
        //
        
        bool is_connected () const
        {
            return clArrayView<T>::is_connected();
        }
        
        void connect (cl_context context, cl_mem_flags flags = CL_MEM_COPY_HOST_PTR | CL_MEM_READ_WRITE)
        {
            clArrayView<T>::connect(context, flags);
        }
        
        void disconnect ()
        {
            clArrayView<T>::disconnect();
        }
        
        cl_int EnqueueDownload (cl_command_queue queue)
        {
            return clArrayView<T>::EnqueueDownload(queue);
        }
        
        cl_int EnqueueUpload (cl_command_queue queue)
        {
            return clArrayView<T>::EnqueueUpload(queue);
        }
        
        cl_mem const & handle () const
        {
            return clArrayView<T>::handle();
        }
};

#endif // WITH_OPENCL

#endif
