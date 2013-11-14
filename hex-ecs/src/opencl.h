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

static const char * const source =
"double2 complex_multiply (double2 a, double2 b)                                                                                           \n"
"{                                                                                                                                         \n"
"    double2 c;                                                                                                                            \n"
"    c.x = a.x * b.x - a.y * b.y;                                                                                                          \n"
"    c.y = a.x * b.y + a.y * b.x;                                                                                                          \n"
"    return c;                                                                                                                             \n"
"}                                                                                                                                         \n"
"                                                                                                                                          \n"
"__kernel void a_vec_b_vec (__constant double2 a, __constant double2 *x, __constant double2 b, __constant double2 *y, __global double2 *z) \n"
"{                                                                                                                                         \n"
"    uint i = get_global_id(0);                                                                                                            \n"
"    z[i] = complex_multiply(a,x[i]) + complex_multiply(b,y[i]);                                                                           \n"
"}                                                                                                                                         \n"
"                                                                                                                                          \n"
"__kernel void vec_mul_vec (__constant double2 *a, __constant double2 *b, __global double2 *c)                                             \n"
"{                                                                                                                                         \n"
"    uint i = get_global_id(0);                                                                                                            \n"
"    z[i] = complex_multiply(a[i],b[i]);                                                                                                   \n"
"}                                                                                                                                         \n"
"                                                                                                                                          \n"
"__kernel CSR_dot_vec (__constant uint *Ap, __constant uint *Ai, __constant double2 *Ax, __constant double2 *x, __global double2 *y)       \n"
"{                                                                                                                                         \n"
"    uint i = get_global_id(0);                                                                                                            \n"
"    double2 sprod = 0.;                                                                                                                   \n"
"                                                                                                                                          \n"
"    for (int idx = Ap[i]; idx < Ap[i + 1]; idx++)                                                                                         \n"
"        sprod += complex_multiply(Ax[idx],x[Ai[idx]]);                                                                                    \n"
"                                                                                                                                          \n"
"    y[i] = sprod;                                                                                                                         \n"
"}                                                                                                                                         \n"

class OpenCLEnvironment
{
    private:
        
        // OpenCL handles
        cl_platform_id platform_;
        cl_device_id device_;
        cl_context context_;
        cl_command_queue queue_;
        cl_program program_;
        
        // computational kernels
        clKernel mmul_;
        clKernel amul_;
        clKernel axby_;
        
    public:
        
        OpenCLEnvironment ()
        {
            // setup environment
            clGetPlatformIDs (1, &platform_, nullptr);
            clGetDeviceIDs (platform_, CL_DEVICE_TYPE_GPU, 1, &device_, nullptr);
            context_ = clCreateContext (nullptr, 1, &device_, nullptr, nullptr, nullptr);
            queue_ = clCreateCommandQueue (context_, device_, 0, nullptr);
            program_ = clCreateProgramWithSource (context_, 1, &source, nullptr, nullptr);
            clBuildProgram (program_, 1, &device_, nullptr, nullptr, nullptr);
            
            // set program entry points
            mmul_ = clKernel ("CSR_dot_vec");
            amul_ = clKernel ("vec_mul_vec");
            axby_ = clKernel ("a_vec_b_vec");
        }
        
        class clKernel
        {
            public:
                
                clKernel (const char* entry_pt)
                    : kernel_(clCreateKernel(program_, entry_pt, nullptr)) {}
                
                ~clKernel ()
                {
                    // ? release kernel
                }
                
                void enqueueExecution (size_t n)
                {
                    clEnqueueNDRangeKernel (queue_, kernel_, 1, nullptr, &n, nullptr, 0, nullptr, nullptr);
                }
                
            private:
                
                cl_kernel kernel_;
        };
        
        template <class T> class clArray : public NumberArray<T>
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
                
                clArray ()
                    : NumberArray<T>(), cl_handle_(nullptr), sync_(false) {}
                clArray (size_t n, T x = T(0))
                    : NumberArray<T>(n,x), cl_handle_(nullptr), sync_(false) {}
                clArray (std::initializer_list<T> list)
                    : NumberArray<T>(list), cl_handle_(nullptr), sync_(false) {}
                
                //
                // destructor
                //
                
                virtual ~clArray ()
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
                const_iterator cbegin () const
                    { return NumberArray<T>::cbegin(); }
                iterator begin ()
                    { sync_ = false; return NumberArray<T>::begin(); }
                const_iterator end () const
                    { return NumberArray<T>::end(); }
                const_iterator cend () const
                    { return NumberArray<T>::cend(); }
                iterator end ()
                    { sync_ = false; return NumberArray<T>::end(); }
                T const & front (size_t i = 0) const
                    { return NumberArray<T>::front(i); }
                T & front (size_t i = 0)
                    { sync_ = false; return NumberArray<T>::front(i); }
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
                    cl_handle_ = clCreateBuffer (context_, flags, size() * sizeof(T), nullptr, nullptr);
                    
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
                    clEnqueueReadBuffer (queue_, cl_handle_, CL_TRUE, 0, size() * sizeof(T), data(), 0, nullptr, nullptr);
                }
                
                void enqueueDownload ()
                {
                    clEnqueueWriteBuffer (queue_, cl_handle_, CL_TRUE, 0, size() * sizeof(T), data(), 0, nullptr, nullptr);
                }
        };
        
};

#endif
