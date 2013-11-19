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

#ifndef NO_OPENCL

#include "opencl.h"

const char * const source =
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
"}                                                                                                                                         \n";

OpenCLEnvironment::OpenCLEnvironment ()
{
    // setup environment
    clGetPlatformIDs (1, &platform_, nullptr);
    clGetDeviceIDs (platform_, CL_DEVICE_TYPE_GPU, 1, &device_, nullptr);
    context_ = clCreateContext (nullptr, 1, &device_, nullptr, nullptr, nullptr);
    queue_ = clCreateCommandQueue (context_, device_, 0, nullptr);
    program_ = clCreateProgramWithSource (context_, 1, const_cast<const char**>(&source), nullptr, nullptr);
    clBuildProgram (program_, 1, &device_, nullptr, nullptr, nullptr);
    
    // OpenCL information
    std::cout << "=== OpenCL information ===\n";
    cl_ulong size;
    clGetDeviceInfo(device_, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &size, 0);
    std::cout << "Local memory size: " << size/1024 << " kiB\n";
    clGetDeviceInfo(device_, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &size, 0);
    std::cout << "Global memory size: " << std::setprecision(3) << size/pow(1024,3) << " GiB\n";
    clGetDeviceInfo(device_, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_ulong), &size, 0);
    std::cout << "Max compute units: " << size << "\n";
    clGetDeviceInfo(device_, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(cl_ulong), &size, 0);
    std::cout << "Max work group size: " << size << "\n\n";
    
    std::cout << "=== Kernel source ===\n";
    std::cout << source << "\n\n";
    
    // set program entry points
    mmul_ = new clKernel ("CSR_dot_vec");
    amul_ = new clKernel ("vec_mul_vec");
    axby_ = new clKernel ("a_vec_b_vec");
}

OpenCLEnvironment::~OpenCLEnvironment ()
{
    delete mmul_;
    delete amul_;
    delete axby_;
}

#endif
