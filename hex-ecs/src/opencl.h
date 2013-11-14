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

// forward declarations
class OpenCLEnvironment;
class clKernel;

extern const char * const source;

class OpenCLEnvironment
{
    public:
        
        OpenCLEnvironment ();
        ~OpenCLEnvironment ();
        
        /**
         * @brief Return OpenCLEnvironment instance.
         * 
         * Singleton interface to OpenCLEnvironment.
         */
        static OpenCLEnvironment & get ()
        {
            static OpenCLEnvironment env_;
            return env_;
        }
        
        //
        // Kernels
        //
        
        static clKernel const & axby () { return *(get().axby_); }
        static clKernel const & amul () { return *(get().amul_); }
        static clKernel const & mmul () { return *(get().mmul_); }
        
        static cl_context context () { return get().context_; }
        static cl_program program () { return get().program_; }
        static cl_command_queue queue () { return get().queue_; }
        
    private:
        
        // OpenCL handles
        cl_platform_id platform_;
        cl_device_id device_;
        cl_context context_;
        cl_command_queue queue_;
        cl_program program_;
        
        // computational kernels
        clKernel *mmul_;
        clKernel *amul_;
        clKernel *axby_;
};

class clKernel
{
    public:
        
        clKernel ()
            : kernel_(nullptr) {}
        clKernel (const char* entry_pt)
            : kernel_(clCreateKernel(OpenCLEnvironment::program(), entry_pt, nullptr)) {}
        
        ~clKernel ()
        {
            // FIXME ? release kernel
        }
        
        template <int iarg = 0> void args ()
        {
            // do nothing
        }
        
        template <
            int iarg = 0,
            class Arg,
            class ...Args
        > void args (Arg ai, Args ...an)
        {
            setArg(ai, iarg);
            args<iarg + 1>(an...);
        }
        
        template <
            class Arg,
            class = typename std::enable_if<is_scalar<Arg>::value>::type
        > void setArg (Arg a, int iarg)
        {
            clSetKernelArg (kernel_, iarg, sizeof(Arg), &a);
        }
        
        template <
            class Arg,
            class = typename std::enable_if<!is_scalar<Arg>::value>::type
        > void setArg (Arg a, int iarg)
        {
            clSetKernelArg (kernel_, iarg, sizeof(a.handle()), &(a.handle()));
        }
        
        void exec (size_t n) const
        {
            if (kernel_ == nullptr)
                throw exception ("[clKernel::enqueueExecution] Kernel not initialized.");
                
            clEnqueueNDRangeKernel (OpenCLEnvironment::queue(), kernel_, 1, nullptr, &n, nullptr, 0, nullptr, nullptr);
            clFinish(OpenCLEnvironment::queue());
        }
        
    private:
        
        // kernel handle
        cl_kernel kernel_;
};

template <class T> class clArrayView : public ArrayView<T>
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
        
        clArrayView ()
            : ArrayView<T>(), cl_handle_(nullptr), sync_(false) {}
        clArrayView (const_iterator i, const_iterator j)
            : ArrayView<T>(const_cast<iterator>(i),const_cast<iterator>(j)), cl_handle_(nullptr), sync_(false) {}
        clArrayView (ArrayView<T> v)
            : ArrayView<T>(v), cl_handle_(nullptr), sync_(false) {}
        clArrayView (clArrayView<T> const & v)
            : ArrayView<T>(v.begin(), v.end()), cl_handle_(v.cl_handle_), sync_(v.sync_) {}
        
        //
        // destructor
        //
        
        virtual ~clArrayView ()
        {
            // free GPU memory
            disconnect ();
        }
        
        //
        // STL interface (except changes in size -- we don't want to bother with GPU reallocation)
        //
        
        size_t size () const
            { return ArrayView<T>::size(); }
        T const & operator [] (size_t i) const
            { return ArrayView<T>::operator[](i); }
        T & operator [] (size_t i)
            { sync_ = false; return ArrayView<T>::operator[](i); }
        T const * data () const
            { return ArrayView<T>::data(); }
        T * data ()
            { sync_ = false; return ArrayView<T>::data(); }
        const_iterator begin () const
            { return ArrayView<T>::begin(); }
        iterator begin ()
            { sync_ = false; return ArrayView<T>::begin(); }
        const_iterator end () const
            { return ArrayView<T>::end(); }
        iterator end ()
            { sync_ = false; return ArrayView<T>::end(); }
        T const & front (size_t i = 0) const
            { return ArrayView<T>::front(i); }
        T & front (size_t i = 0)
            { sync_ = false; return ArrayView<T>::front(i); }
        T const & back (size_t i = 0) const
            { return ArrayView<T>::back(i); }
        T & back (size_t i = 0)
            { sync_ = false; return ArrayView<T>::back(i); }
        
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
            cl_handle_ = clCreateBuffer (OpenCLEnvironment::context(), flags, size() * sizeof(T), nullptr, nullptr);
            
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
        
        void enqueueUpload ()
        {
            clEnqueueReadBuffer (OpenCLEnvironment::queue(), cl_handle_, CL_TRUE, 0, size() * sizeof(T), data(), 0, nullptr, nullptr);
            clFinish(OpenCLEnvironment::queue());
            sync_ = true;
        }
        
        void enqueueDownload ()
        {
            clEnqueueWriteBuffer (OpenCLEnvironment::queue(), cl_handle_, CL_TRUE, 0, size() * sizeof(T), data(), 0, nullptr, nullptr);
            clFinish(OpenCLEnvironment::queue());
            sync_ = true;
        }
        
        cl_mem handle () const
        {
            return cl_handle_;
        }
};

template <class T> class clNumberArray : virtual public clArrayView<T>, virtual public NumberArray<T>
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
        
        clNumberArray ()
            : NumberArray<T>(), cl_handle_(nullptr), sync_(false) {}
        clNumberArray (size_t n, T x = T(0))
            : NumberArray<T>(n,x), cl_handle_(nullptr), sync_(false) {}
        clNumberArray (size_t n, T const * ptr)
            : NumberArray<T>(n,ptr), cl_handle_(nullptr), sync_(false) {}
        clNumberArray (ArrayView<T> v)
            : NumberArray<T>(v), cl_handle_(nullptr), sync_(false) {}
        clNumberArray (clArrayView<T> const & v)
            : NumberArray<T>(v.size(), v.data()), cl_handle_(v.cl_handle_), sync_(v.sync_) {}
        
        //
        // destructor
        //
        
        virtual ~clNumberArray ()
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
