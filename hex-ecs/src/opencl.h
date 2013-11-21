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

#ifndef HEX_OPENCL
#define HEX_OPENCL

#ifndef NO_OPENCL

#include <CL/cl.h>

#include "arrays.h"

template <class T> class CLArrayView : public ArrayView<T>
{
    private:
        
        /// Memory handle to the data copy on the GPU.
        cl_mem cl_handle_;
        
        /// Whether the data are synchronized between RAM and GPU.
        bool sync_;
        
    public:
        
        //
        // aliases
        //
        
        typedef T DataType;
        typedef T * iterator;
        typedef T const * const_iterator;
        
        //
        // constructors
        //
        
        CLArrayView ()
            : ArrayView<T>(), cl_handle_(nullptr), sync_(false) {}
        CLArrayView (const_iterator i, const_iterator j)
            : ArrayView<T>(const_cast<iterator>(i),const_cast<iterator>(j)), cl_handle_(nullptr), sync_(false) {}
        CLArrayView (ArrayView<T> v)
            : ArrayView<T>(v), cl_handle_(nullptr), sync_(false) {}
        CLArrayView (CLArrayView<T> const & v)
            : ArrayView<T>(v.begin(), v.end()), cl_handle_(v.cl_handle_), sync_(v.sync_) {}
        
        //
        // destructor
        //
        
        virtual ~CLArrayView ()
        {
            // free GPU memory
            disconnect ();
        }
        
        //
        // STL interface (except changes in size -- we don't want to bother with GPU reallocation)
        //
        
        size_t size () const                    { return ArrayView<T>::size();        }
        T const & operator [] (size_t i) const  { return ArrayView<T>::operator[](i); }
        T & operator [] (size_t i)              { return ArrayView<T>::operator[](i); }
        T const * data () const                 { return ArrayView<T>::data();        }
        T * data ()                             { return ArrayView<T>::data();        }
        const_iterator begin () const           { return ArrayView<T>::begin();       }
        iterator begin ()                       { return ArrayView<T>::begin();       }
        const_iterator end () const             { return ArrayView<T>::end();         }
        iterator end ()                         { return ArrayView<T>::end();         }
        T const & front (size_t i = 0) const    { return ArrayView<T>::front(i);      }
        T & front (size_t i = 0)                { return ArrayView<T>::front(i);      }
        T const & back (size_t i = 0) const     { return ArrayView<T>::back(i);       }
        T & back (size_t i = 0)                 { return ArrayView<T>::back(i);       }
        
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
            cl_handle_ = clCreateBuffer (context, flags, size() * sizeof(T), nullptr, nullptr);
        }
        
        void disconnect ()
        {
            // release previous allocations
            if (cl_handle_ != nullptr)
                clReleaseMemObject(cl_handle_);
            
            // clear pointer
            cl_handle_ = nullptr;
        }
        
        void enqueueUpload (cl_command_queue queue)
        {
            clEnqueueReadBuffer (queue, cl_handle_, CL_TRUE, 0, size() * sizeof(T), data(), 0, nullptr, nullptr);
            clFinish(queue);
        }
        
        void enqueueDownload (cl_command_queue queue)
        {
            clEnqueueWriteBuffer (queue, cl_handle_, CL_TRUE, 0, size() * sizeof(T), data(), 0, nullptr, nullptr);
            clFinish(queue);
        }
        
        cl_mem handle () const
        {
            return cl_handle_;
        }
};

#define CL_ALIGNMENT 4096
template <class T, class Alloc = AlignedAllocator<T,CL_ALIGNMENT>> class CLArray : public CLArrayView<T>
{
    public:
        
        //
        // aliases
        //
        
        typedef T DataType;
        typedef T * iterator;
        typedef T const * const_iterator;
        
        //
        // constructors
        //
        
        CLArray ()
            : NumberArray<T>(), cl_handle_(nullptr), sync_(false) {}
        CLArray (size_t n, T x = T(0))
            : NumberArray<T>(n,x), cl_handle_(nullptr), sync_(false) {}
        CLArray (size_t n, T const * ptr)
            : NumberArray<T>(n,ptr), cl_handle_(nullptr), sync_(false) {}
        CLArray (ArrayView<T> v)
            : NumberArray<T>(v), cl_handle_(nullptr), sync_(false) {}
        CLArray (CLArrayView<T> const & v)
            : NumberArray<T>(v.size(), v.data()), cl_handle_(v.cl_handle_), sync_(v.sync_) {}
        
        //
        // destructor
        //
        
        virtual ~CLArray ()
        {
            // free GPU memory
            disconnect ();
            
            // free RAM
            NumberArray<T>::~NumberArray();
        }
        
        //
        // STL interface (except changes in size -- we don't want to bother with GPU reallocation)
        //
        
        size_t size () const
            { return NumberArray<T>::size(); }
        T const & operator [] (size_t i) const
            { return NumberArray<T>::operator[](i); }
        T & operator [] (size_t i)
            { sync_ = false; return NumberArray<T>::operator[](i); }
        T const * data () const
            { return NumberArray<T>::data(); }
        T * data ()
            { sync_ = false; return NumberArray<T>::data(); }
        const_iterator begin () const
            { return NumberArray<T>::begin(); }
        iterator begin ()
            { sync_ = false; return ArrayView<T>::begin(); }
        const_iterator end () const
            { return NumberArray<T>::end(); }
        iterator end ()
            { sync_ = false; return ArrayView<T>::end(); }
        T const & front (size_t i = 0) const
            { return NumberArray<T>::front(i); }
        T & front (size_t i = 0)
            { sync_ = false; return ArrayView<T>::front(i); }
        T const & back (size_t i = 0) const
            { return NumberArray<T>::back(i); }
        T & back (size_t i = 0)
            { sync_ = false; return NumberArray<T>::back(i); }
        
        //
        // OpenCL intrinsics
        //
        
        bool is_connected () const
        {
            return cl_handle_ != nullptr;
        }
        
        void connect (cl_mem_flags flags)
        {
            // release previous allocation prior to creating a new buffer !
            if (cl_handle_ != nullptr)
                throw exception ("[clArray::connect] Release the buffer before connecting!");
            
            // allocate memory on GPU
            cl_handle_ = clCreateBuffer(OpenCLEnvironment::context(), flags, size() * sizeof(T), nullptr, nullptr);
            
            // mark as synchronized if the data ware copied
            if (flags & CL_MEM_COPY_HOST_PTR)
                sync_ = true;
        }
        
        void disconnect ()
        {
            // release previous allocations
            if (cl_handle_ != nullptr)
                clReleaseMemObject(cl_handle_);
            
            // clear pointer
            cl_handle_ = nullptr;
        }
        
        void download ()
        {
            clEnqueueReadBuffer (OpenCLEnvironment::queue(), cl_handle_, CL_TRUE, 0, size() * sizeof(T), data(), 0, nullptr, nullptr);
            clFinish(OpenCLEnvironment::queue());
            sync_ = true;
        }
        
        void upload ()
        {
            clEnqueueWriteBuffer (OpenCLEnvironment::queue(), cl_handle_, CL_TRUE, 0, size() * sizeof(T), data(), 0, nullptr, nullptr);
            clFinish(OpenCLEnvironment::queue());
            sync_ = true;
        }
};

#endif

#endif
