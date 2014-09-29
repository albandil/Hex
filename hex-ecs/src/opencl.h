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

#ifndef NO_OPENCL

#include <CL/cl.h>

#include "arrays.h"

/**
 * @brief OpenCL array wrapper.
 *
 * This class equips ArrayView<T> with an OpenCL context handle and
 * provides several service methods for upload and download to / from
 * the compute device (mostly GPU).
 */
template <class T> class CLArrayView : public ArrayView<T>
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
        
        CLArrayView ()
            : ArrayView<T>(), cl_handle_(nullptr) { /*std::cout << "Init Default\n";*/ }
        CLArrayView (std::size_t n, T const * ptr)
            : ArrayView<T>(n, const_cast<T*>(ptr)), cl_handle_(nullptr) { /*std::cout << "Init from ptr.\n";*/ }
        CLArrayView (ArrayView<T> v, std::size_t i = 0, std::size_t n = 0)
            : ArrayView<T>(v, i, n), cl_handle_(nullptr) { /*std::cout << "Init from ArrayView\n";*/ }
        CLArrayView (CLArrayView<T> const & v)
            : ArrayView<T>(v), cl_handle_(v.cl_handle_) { /*std::cout << "Init from CLArrayView\n";*/ }
        CLArrayView (CLArrayView<T> && v)
            : ArrayView<T>(), cl_handle_(nullptr)
        {
            std::swap (ArrayView<T>::array_, v.ArrayView<T>::array_);
            std::swap (ArrayView<T>::N_,     v.ArrayView<T>::N_    );
            std::swap (cl_handle_,           v.cl_handle_          );
        }
        
        //
        // destructor
        //
        
        virtual ~CLArrayView ()
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
                throw exception ("[clArray::connect] Release the buffer before connecting!");
            
            // allocate memory on GPU
            cl_handle_ = clCreateBuffer (context, flags, size() * sizeof(T), data(), nullptr);
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
            return clEnqueueWriteBuffer (queue, cl_handle_, CL_TRUE, 0, size() * sizeof(T), data(), 0, nullptr, nullptr);
        }
        
        cl_int EnqueueDownload (cl_command_queue queue)
        {
            return clEnqueueReadBuffer (queue, cl_handle_, CL_TRUE, 0, size() * sizeof(T), data(), 0, nullptr, nullptr);
        }
        
        cl_mem const & handle () const
        {
            return cl_handle_;
        }
};

/**
 * @brief OpenCL array wrapper.
 *
 * This class combines NumberArray<T> and CLArrayView<T>. It is a self-standing array of
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
template <class T, class Alloc = AlignedAllocator<T,CL_ALIGNMENT>> class CLArray : public CLArrayView<T>
#else
template <class T, class Alloc = PlainAllocator<T>> class CLArray : public CLArrayView<T>
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
        
        CLArray () : CLArrayView<T>()
        {
            // nothing to do
        }
        CLArray (std::size_t n, T x = T(0)) : CLArrayView<T>(n, Alloc::alloc(n))
        {
            for (std::size_t i = 0; i < n; i++)
                (*this)[i] = x;
        }
        CLArray (std::size_t n, T const * px) : CLArrayView<T>(n, Alloc::alloc(n))
        {
            for (std::size_t i = 0; i < n; i++)
                (*this)[i] = px[i];
        }
        CLArray (ArrayView<T> v, std::size_t i = 0, std::size_t n = 0)
            : CLArrayView<T>(((n == 0) ? v.size() : n), Alloc::alloc((n == 0) ? v.size() : n))
        {
            for (std::size_t j = 0; j < size(); j++)
                (*this)[j] = v[i+j];
        }
        CLArray (std::initializer_list<T> list) : CLArrayView<T>(list.size(), Alloc::alloc(list.size()))
        {
            std::size_t i = 0;
            for (auto it = list.begin(); it != list.end(); it++)
                (*this)[i++] = *it;
        }
        CLArray (CLArray<T,Alloc> && rvrf) : CLArrayView<T>(std::move(rvrf))
        {
            // nothing to do
        }
        
        //
        // data assignment
        //
        
        CLArray & operator = (const ArrayView<T> v)
        {
            // realloc memory if needed
            if (size() != v.size())
            {
                Alloc::free(data());
                ArrayView<T>::array_ = Alloc::alloc(v.size());
                ArrayView<T>::N_ = v.size();
            }
            
            // copy data
            for (std::size_t j = 0; j < size(); j++)
                (*this)[j] = v[j];
            
            return *this;
        }
        
        //
        // destructor
        //
        
        virtual ~CLArray ()
        {
            // free GPU memory
            disconnect();
            
            // free RAM memory
            Alloc::free(data());
        }
        
        //
        // STL interface (except changes in size -- we don't want to bother with GPU reallocation)
        //
        
        std::size_t size () const                   { return CLArrayView<T>::size();        }
        T const & operator [] (std::size_t i) const { return CLArrayView<T>::operator[](i); }
        T & operator [] (std::size_t i)             { return CLArrayView<T>::operator[](i); }
        T const * data () const                     { return CLArrayView<T>::data();        }
        T * data ()                                 { return CLArrayView<T>::data();        }
        const_iterator begin () const               { return CLArrayView<T>::begin();       }
        iterator begin ()                           { return CLArrayView<T>::begin();       }
        const_iterator end () const                 { return CLArrayView<T>::end();         }
        iterator end ()                             { return CLArrayView<T>::end();         }
        T const & front (std::size_t i = 0) const   { return CLArrayView<T>::front(i);      }
        T & front (std::size_t i = 0)               { return CLArrayView<T>::front(i);      }
        T const & back (std::size_t i = 0) const    { return CLArrayView<T>::back(i);       }
        T & back (std::size_t i = 0)                { return CLArrayView<T>::back(i);       }
        
        //
        // OpenCL intrinsics
        //
        
        bool is_connected () const
        {
            return CLArrayView<T>::is_connected();
        }
        
        void connect (cl_context context, cl_mem_flags flags = CL_MEM_COPY_HOST_PTR | CL_MEM_READ_WRITE)
        {
            CLArrayView<T>::connect(context, flags);
        }
        
        void disconnect ()
        {
            CLArrayView<T>::disconnect();
        }
        
        cl_int EnqueueDownload (cl_command_queue queue)
        {
            return CLArrayView<T>::EnqueueDownload(queue);
        }
        
        cl_int EnqueueUpload (cl_command_queue queue)
        {
            return CLArrayView<T>::EnqueueUpload(queue);
        }
        
        cl_mem const & handle () const
        {
            return CLArrayView<T>::handle();
        }
};

#endif

#endif
