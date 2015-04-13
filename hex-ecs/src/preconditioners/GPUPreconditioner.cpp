//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2015, Jakub Benda, Charles University in Prague                    //
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

#if (!defined(NO_OPENCL) && !defined(NO_LAPACK))

#include <iostream>
#include <set>

#include "../arrays.h"
#include "../clarrays.h"
#include "../gauss.h"
#include "../itersolve.h"
#include "../misc.h"
#include "../preconditioners.h"
#include "../radial.h"
#include "../special.h"

#include <CL/cl.h>

const std::string GPUCGPreconditioner::prec_name = "GPU";
const std::string GPUCGPreconditioner::prec_description = "Block inversion using conjugate gradients preconditioned by Kronecker product approximation (on GPU).";

// kernels' source as byte array, generated by "xxd" from the CL source
char kernels_cl[] = {
    #include "GPUPreconditioner.inc"
    , 0x00 // terminate the string by zero
};

// pointer to the source; to be used in setup
char * source = &kernels_cl[0];

void GPUCGPreconditioner::setup ()
{
    KPACGPreconditioner::setup();
    
    // shorthands
    int order = s_bspline_.order();
    int Nspline = s_bspline_.Nspline();
    int Nreknot = s_bspline_.Nreknot();
    
    std::cout << "Setting up OpenCL environment" << std::endl;
    
    // auxiliary variables
    char platform_name[1024], platform_vendor[1024], platform_version[1024], device_name[1024], device_vendor[1024];
    cl_platform_id platforms[10]; cl_uint nplatforms;
    cl_device_id devices[10]; cl_uint ndevices;
    
    // check that the chosen platform is available
    clGetPlatformIDs(10, platforms, &nplatforms);
    if (cmd_.ocl_platform >= nplatforms)
        HexException("The requested platform index (%d) does not exist on this system. Run hex-ecs --cl-list to find all existing.", cmd_.ocl_platform);
    
    platform_ = platforms[cmd_.ocl_platform];
    clGetPlatformInfo(platform_, CL_PLATFORM_NAME, sizeof(platform_name), platform_name, nullptr);
    clGetPlatformInfo(platform_, CL_PLATFORM_VENDOR, sizeof(platform_vendor), platform_vendor, nullptr);
    clGetPlatformInfo(platform_, CL_PLATFORM_VERSION, sizeof(platform_version), platform_version, nullptr);
    
    std::cout << "\t- platform: " << platform_name << " (" << platform_vendor << ")" << std::endl;
    std::cout << "\t- available version: " << platform_version << std::endl;
    
    // check that the chosen device is available
    clGetDeviceIDs(platform_, CL_DEVICE_TYPE_ALL, 10, devices, &ndevices);
    if (cmd_.ocl_device >= nplatforms)
        HexException("The requested device index (%d) does not exist for platform %d. Run hex-ecs --cl-list to find all existing.",  cmd_.ocl_device, cmd_.ocl_platform);
    
    device_ = devices[cmd_.ocl_device];
    clGetDeviceInfo(device_, CL_DEVICE_NAME, sizeof(device_name), device_name, nullptr);
    clGetDeviceInfo(device_, CL_DEVICE_VENDOR, sizeof(device_vendor), device_vendor, nullptr);
    
    std::cout << "\t- device: " << device_name << " (" << device_vendor << ")" << std::endl;
    
    cl_ulong max_compute_units, max_work_group_size, local_memory_size, global_memory_size;
    clGetDeviceInfo(device_, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_ulong), &max_compute_units, 0);
    clGetDeviceInfo(device_, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(cl_ulong), &max_work_group_size, 0);
    clGetDeviceInfo(device_, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &local_memory_size, 0);
    clGetDeviceInfo(device_, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &global_memory_size, 0);
    
    std::cout << "\t- max compute units: " << max_compute_units << std::endl;
    std::cout << "\t- max work group size: " << max_work_group_size << std::endl;
    
    // choose default workgroup size
    Nlocal_ = std::min<std::size_t>(64, max_work_group_size);
    
    // choose block size for "mul_ABt" kernel, aim at 4 concurrent threads
    blocksize_ = std::sqrt(local_memory_size / (4/*threads*/ * 16/*double*/ * 2/*arrays*/));
    blocksize_ = std::min<std::size_t>(blocksize_, std::sqrt(max_work_group_size));
    std::cout << "\t- matrix multiplication tile size: " << blocksize_ << std::endl;
    
    // compute global memory requirements of the preconditioner
    std::size_t greq = 0;
    if (not cmd_.gpu_large_data)
    {
        // - Nspline*Nspline for: four preconditioner matrices, 5 CG arrays (x,b,r,p,z) and one temporary array (tmA)
        greq += 10 * (std::size_t)Nspline * (std::size_t)Nspline;
    }
    // - padded one-electron matrices (S, D, M1, M2)
    greq += 4 * s_rad_.S().size();
    // - preconditioner eigenvalues
    greq += 2 * Nspline;
    // - full integral moments
    greq += 2 * (s_rad_.maxlambda() + 1) * s_rad_.S().size();
    // - partial integral moments
    greq += 2 * s_rad_.Mitr_L(-1).size();
    // - all these were complex numbers
    greq *= 16;
    
    std::cout << "\t- global memory size: " << format("%.2f", global_memory_size/gsl_sf_pow_int(1024,3)) << " GiB ";
    std::cout << "(apx. " << format("%.2f", greq * 100. / global_memory_size) << " % will be used)" << std::endl;
    
    if (cmd_.gpu_large_data)
    {
        std::size_t host_mem = 16 * 10 * (std::size_t)Nspline * (std::size_t)Nspline;
        std::cout << "\t- data kept in host memory: " << format("%.2f", host_mem/gsl_sf_pow_int(1024,3)) << " GiB " << std::endl;
        std::cout << "\t- WARNING: --ocl-use-host-memory will slow down the solution due to the host-device data transfers." << std::endl;
    }
    
    // compute local memory requirements of the preconditioner
    std::size_t lreq = 0;
    // - either Complex per default Nlocal_
    lreq = Nlocal_;
    // - or the requirements of mul_ABt
    lreq = std::max<std::size_t>(lreq, 2/*array*/ * 16/*double*/ * blocksize_);
    // - all these were complex numbers
    lreq *= 16;
    
    std::cout << "\t- local memory size: " << local_memory_size/1024 << " kiB ";
    std::cout << "(apx. " << format("%.2f", lreq * 100. / local_memory_size) << " % will be used)" << std::endl << std::endl;
    
    // create context and command queue
    context_ = clCreateContext(nullptr, 1, &device_, nullptr, nullptr, nullptr);
#ifdef CL_VERSION_2_0
    queue_ = clCreateCommandQueueWithProperties(context_, device_, nullptr, nullptr);
#else
    queue_ = clCreateCommandQueue(context_, device_, 0, nullptr);
#endif
    
    // setup compile flags
    std::ostringstream flags;
    flags << " -cl-fast-relaxed-math ";
    flags << " -D ORDER="      << order      << " ";
    flags << " -D NSPLINE="    << Nspline    << " ";
    flags << " -D NREKNOT="    << Nreknot    << " ";
    flags << " -D NLOCAL="     << Nlocal_    << " ";
    flags << " -D BLOCK_SIZE=" << blocksize_ << " ";
    
    // build program
    program_ = clCreateProgramWithSource(context_, 1, const_cast<const char**>(&source), nullptr, nullptr);
    clBuildProgram(program_, 1, &device_, flags.str().c_str(), nullptr, nullptr);
    
    cl_build_status status;
    clGetProgramBuildInfo(program_, device_, CL_PROGRAM_BUILD_STATUS, sizeof(status), &status, nullptr);
    if (status != CL_SUCCESS)
    {
        std::cout << std::endl << "Source:" << std::endl << source << std::endl;
        std::cout << std::endl << "Command line:" << std::endl << flags.str() << std::endl << std::endl;
        
        char log [100000];
        clGetProgramBuildInfo(program_, device_, CL_PROGRAM_BUILD_LOG, sizeof(log), log, nullptr);
        std::cout << "clGetProgramBuildInfo: log" << std::endl << log << std::endl;
        
        HexException("Failed to initialize OpenCL.");
    }
    
    // determine where to place the data (device / host)
    largeDataFlags_ = smallDataFlags_ = CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR;
    if (cmd_.gpu_large_data)
        largeDataFlags_ = CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR;
    
    // set program entry points
    mml1_ = clCreateKernel(program_, "mmul_1el",        nullptr);
    mml2_ = clCreateKernel(program_, "mmul_2el",        nullptr);
    axby_ = clCreateKernel(program_, "a_vec_b_vec",     nullptr);
    norm_ = clCreateKernel(program_, "norm",            nullptr);
    spro_ = clCreateKernel(program_, "scalar_product",  nullptr);
    mabt_ = clCreateKernel(program_, "mul_ABt",         nullptr);
    krdv_ = clCreateKernel(program_, "kron_div",        nullptr);
}

void GPUCGPreconditioner::precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const
{
    // shorthands
    std::size_t Nspline = s_rad_.bspline().Nspline();
    std::size_t Nsegsiz = Nspline * Nspline;
    
    // performance timers
    std::size_t us_prec = 0, us_mmul = 0, us_spro = 0, us_axby = 0, us_norm = 0;
    std::size_t us_mmul_1 = 0;
    
    // round 'Nsegsiz' to nearest larger multiple of Nlocal_
    std::size_t Nglobal = Nlocal_ * ((Nsegsiz + Nlocal_ - 1) / Nlocal_);
    
    // iterations
    iArray n (l1_l2_.size());
    
    // some OpenCL auxiliary storage arrays (used by kernels for temporary data)
    clArray<Complex> tmp (Nglobal / Nlocal_);
    clArray<double>  nrm (Nglobal / Nlocal_);
    clArray<Complex> tmA (Nspline * Nspline);
    tmp.connect(context_, smallDataFlags_);
    nrm.connect(context_, smallDataFlags_);
    tmA.connect(context_, largeDataFlags_);
    
    // connect B-spline knots
    clArrayView<Complex> t (s_bspline_.t().size(), s_bspline_.t().data());
    t.connect(context_, smallDataFlags_);
    
    // create OpenCL representation of the one-electron matrices + transfer data to GPU memory
    clArrayView<Complex> S_p (s_rad_.S().data().size(), s_rad_.S().data().data());
    clArrayView<Complex> D_p (s_rad_.D().data().size(), s_rad_.D().data().data());
    clArrayView<Complex> Mm1_tr_p (s_rad_.Mm1_tr().data().size(), s_rad_.Mm1_tr().data().data());
    clArrayView<Complex> Mm2_p (s_rad_.Mm2().data().size(), s_rad_.Mm2().data().data());
    S_p.connect(context_, smallDataFlags_);
    D_p.connect(context_, smallDataFlags_);
    Mm1_tr_p.connect(context_, smallDataFlags_);
    Mm2_p.connect(context_, smallDataFlags_);
    
    // create OpenCL representation of the one-electron partial integral moments + transfer data to GPU memory
    clArrayView<Complex> Mi_L[s_rad_.maxlambda() + 1], Mi_mLm1[s_rad_.maxlambda() + 1];
    clArrayView<Complex> M_L[s_rad_.maxlambda() + 1], M_mLm1[s_rad_.maxlambda() + 1];
    for (int lambda = 0; lambda <= s_rad_.maxlambda(); lambda++)
    {
        Mi_L[lambda].reset(s_rad_.Mitr_L(lambda).size(), s_rad_.Mitr_L(lambda).data());
        Mi_mLm1[lambda].reset(s_rad_.Mitr_mLm1(lambda).size(), s_rad_.Mitr_mLm1(lambda).data());
        M_L[lambda].reset(s_rad_.Mtr_L(lambda).data().size(), const_cast<Complex*>(s_rad_.Mtr_L(lambda).data().data()));
        M_mLm1[lambda].reset(s_rad_.Mtr_mLm1(lambda).data().size(), const_cast<Complex*>(s_rad_.Mtr_mLm1(lambda).data().data()));
        
        Mi_L[lambda].connect(context_, smallDataFlags_);
        Mi_mLm1[lambda].connect(context_, smallDataFlags_);
        M_L[lambda].connect(context_, smallDataFlags_);
        M_mLm1[lambda].connect(context_, smallDataFlags_);
    }
    
    // for all diagonal blocks
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyWork(ill))
    {
        // get angular momenta
        int l1 = l1_l2_[ill].first;
        int l2 = l1_l2_[ill].second;
        
        // load blocks
        if (cmd_.outofcore)
        {
            const_cast<BlockArray<Complex>&>(r).hdfload(ill);
            z.hdfload(ill);
        }
        
        // create segment views
        clArrayView<Complex> rview (r[ill], 0, Nsegsiz);
        clArrayView<Complex> zview (z[ill], 0, Nsegsiz);
        rview.connect(context_, largeDataFlags_);
        zview.connect(context_, largeDataFlags_);
        
        // initialize block preconditioner
        this->CG_init(ill);
        
        // preconditioner matrices
        clArrayView<Complex> prec1a (Nspline * Nspline, prec_[l1].invCl_invsqrtS.data().data());
        clArrayView<Complex> prec2a (Nspline * Nspline, prec_[l2].invCl_invsqrtS.data().data());
        clArrayView<Complex> Dl1 (Nspline, prec_[l1].Dl.data());
        clArrayView<Complex> Dl2 (Nspline, prec_[l2].Dl.data());
        clArrayView<Complex> prec1b (Nspline * Nspline, prec_[l1].invsqrtS_Cl.data().data());
        clArrayView<Complex> prec2b (Nspline * Nspline, prec_[l2].invsqrtS_Cl.data().data());
        prec1a.connect(context_, largeDataFlags_);
        prec2a.connect(context_, largeDataFlags_);
        Dl1.connect(context_, smallDataFlags_);
        Dl2.connect(context_, smallDataFlags_);
        prec1b.connect(context_, largeDataFlags_);
        prec2b.connect(context_, largeDataFlags_);
        
        // angular integrals
        iArray lambdas (s_rad_.maxlambda() + 1);
        rArray fs (s_rad_.maxlambda() + 1);
        int Nlambdas = 0;
        for (int lambda = 0; lambda <= s_rad_.maxlambda(); lambda++)
        {
            double f = special::computef(lambda,l1,l2,l1,l2,inp_.L);
            
            if (not std::isfinite(f))
                HexException("Failed to compute angular integral f[%d](%d,%d,%d,%d;%d).", lambda, l1, l2, l1, l2, inp_.L);
            
            if (f == 0)
                continue;
            
            lambdas[Nlambdas] = lambda;
            fs[Nlambdas] = f;
            
            Nlambdas++;
        }
        
        // allocation (and upload) of an OpenCL array
        auto new_opencl_array = [&](std::size_t n, std::string name) -> clArray<Complex>
        {
            // create array
            clArray<Complex> a(n);
            
            // connect the array to GPU
            a.connect(context_, largeDataFlags_);
            
            // use this array
            return a;
        };
        
        // multiplies vector by the (approximate) ill-th diagonal block
        auto inner_mmul = [&](/* const */ clArrayView<Complex> a, clArrayView<Complex> b) -> void
        {
            // multiply
            //      b = A · a
            
            Timer timer;
            
            // one-electron contribution
            clSetKernelArg(mml1_, 0, sizeof(double), &E_);
            clSetKernelArg(mml1_, 1, sizeof(cl_mem), &S_p.handle());
            clSetKernelArg(mml1_, 2, sizeof(cl_mem), &D_p.handle());
            clSetKernelArg(mml1_, 3, sizeof(cl_mem), &Mm1_tr_p.handle());
            clSetKernelArg(mml1_, 4, sizeof(cl_mem), &Mm2_p.handle());
            clSetKernelArg(mml1_, 5, sizeof(int),    &l1);
            clSetKernelArg(mml1_, 6, sizeof(int),    &l2);
            clSetKernelArg(mml1_, 7, sizeof(cl_mem), &a.handle());
            clSetKernelArg(mml1_, 8, sizeof(cl_mem), &b.handle());
            clEnqueueNDRangeKernel(queue_, mml1_, 1, nullptr, &Nsegsiz, nullptr, 0, nullptr, nullptr);
            clFinish(queue_);
            
            us_mmul_1 += timer.microseconds();
            
            // two-electron contribution
            for (int ilambda = 0; ilambda < Nlambdas; ilambda++)
            {
                clSetKernelArg(mml2_, 0, sizeof(cl_mem), &t.handle());
                clSetKernelArg(mml2_, 1, sizeof(int),    &(lambdas[ilambda]));
                clSetKernelArg(mml2_, 2, sizeof(double), &(fs[ilambda]));
                clSetKernelArg(mml2_, 3, sizeof(cl_mem), &(M_L[lambdas[ilambda]].handle()));
                clSetKernelArg(mml2_, 4, sizeof(cl_mem), &(M_mLm1[lambdas[ilambda]].handle()));
                clSetKernelArg(mml2_, 5, sizeof(cl_mem), &(Mi_L[lambdas[ilambda]].handle()));
                clSetKernelArg(mml2_, 6, sizeof(cl_mem), &(Mi_mLm1[lambdas[ilambda]].handle()));
                clSetKernelArg(mml2_, 7, sizeof(cl_mem), &a.handle());
                clSetKernelArg(mml2_, 8, sizeof(cl_mem), &b.handle());
                clEnqueueNDRangeKernel(queue_, mml2_, 1, nullptr, &Nsegsiz, nullptr, 0, nullptr, nullptr);
            }
            clFinish(queue_);
            
            us_mmul += timer.microseconds();
        };
        
        // applies KPA preconditioner (two "kron-dots")
        auto inner_prec = [&](const clArrayView<Complex> x, clArrayView<Complex> y) -> void
        {
            // multiply by approximate inverse block
            Timer timer;
            
            std::size_t gsize[2] = { blocksize_ * ((Nspline + blocksize_ - 1) / blocksize_), blocksize_ * ((Nspline + blocksize_ - 1) / blocksize_) };
            std::size_t lsize[2] = { blocksize_, blocksize_ };
            
            clSetKernelArg(mabt_, 0, sizeof(cl_mem), &prec2a.handle());
            clSetKernelArg(mabt_, 1, sizeof(cl_mem), &x.handle());
            clSetKernelArg(mabt_, 2, sizeof(cl_mem), &tmA.handle());
            clEnqueueNDRangeKernel(queue_, mabt_, 2, nullptr, gsize, lsize, 0, nullptr, nullptr);
            clFinish(queue_);
            
            clSetKernelArg(mabt_, 0, sizeof(cl_mem), &prec1a.handle());
            clSetKernelArg(mabt_, 1, sizeof(cl_mem), &tmA.handle());
            clSetKernelArg(mabt_, 2, sizeof(cl_mem), &y.handle());
            clEnqueueNDRangeKernel(queue_, mabt_, 2, nullptr, gsize, lsize, 0, nullptr, nullptr);
            clFinish(queue_);
            
            clSetKernelArg(krdv_, 0, sizeof(Complex), &E_);
            clSetKernelArg(krdv_, 1, sizeof(cl_mem),  &Dl1.handle());
            clSetKernelArg(krdv_, 2, sizeof(cl_mem),  &Dl2.handle());
            clSetKernelArg(krdv_, 3, sizeof(cl_mem),  &y.handle());
            clEnqueueNDRangeKernel(queue_, krdv_, 1, nullptr, &Nspline, nullptr, 0, nullptr, nullptr);
            clFinish(queue_);
            
            clSetKernelArg(mabt_, 0, sizeof(cl_mem), &prec2b.handle());
            clSetKernelArg(mabt_, 1, sizeof(cl_mem), &y.handle());
            clSetKernelArg(mabt_, 2, sizeof(cl_mem), &tmA.handle());
            clEnqueueNDRangeKernel(queue_, mabt_, 2, nullptr, gsize, lsize, 0, nullptr, nullptr);
            clFinish(queue_);
            
            clSetKernelArg(mabt_, 0, sizeof(cl_mem), &prec1b.handle());
            clSetKernelArg(mabt_, 1, sizeof(cl_mem), &tmA.handle());
            clSetKernelArg(mabt_, 2, sizeof(cl_mem), &y.handle());
            clEnqueueNDRangeKernel(queue_, mabt_, 2, nullptr, gsize, lsize, 0, nullptr, nullptr);
            clFinish(queue_);
            
            us_prec += timer.microseconds();
        };
        
        // computes norm of the vector
        auto compute_norm = [&](/* const */ clArrayView<Complex> x) -> double
        {
            // multiply
            //     |x|² -> nrm
            
            Timer timer;
            
            clSetKernelArg(norm_, 0, sizeof(cl_mem), &x.handle());
            clSetKernelArg(norm_, 1, sizeof(cl_mem), &nrm.handle());
            clEnqueueNDRangeKernel(queue_, norm_, 1, nullptr, &Nglobal, &Nlocal_, 0, nullptr, nullptr);
            nrm.EnqueueDownload(queue_);
            clFinish(queue_);
            
            us_norm += timer.microseconds();
            
            // return square root of the sum
            return std::sqrt(sum(nrm));
        };
        
        // computes scalar product of two arrays
        auto scalar_product = [&](const clArrayView<Complex> x, /* const */ clArrayView<Complex> y) -> Complex
        {
            // multiply
            //     x * y -> tmp
            
            Timer timer;
            
            clSetKernelArg(spro_, 0, sizeof(cl_mem), &x.handle());
            clSetKernelArg(spro_, 1, sizeof(cl_mem), &y.handle());
            clSetKernelArg(spro_, 2, sizeof(cl_mem), &tmp.handle());
            clEnqueueNDRangeKernel(queue_, spro_, 1, nullptr, &Nglobal, &Nlocal_, 0, nullptr, nullptr);
            tmp.EnqueueDownload(queue_);
            clFinish(queue_);
            
            us_spro += timer.microseconds();
            
            // sum the product of the arrays
            return sum(tmp);
        };
        
        // weighted sum of two arrays
        auto axby = [&](Complex a, clArrayView<Complex> x, Complex b, const clArrayView<Complex> y) -> void
        {
            // multiply
            //     a * x + b * y -> z
            
            Timer timer;
            
            clSetKernelArg(axby_, 0, sizeof(Complex), &a);
            clSetKernelArg(axby_, 1, sizeof(cl_mem),  &x.handle());
            clSetKernelArg(axby_, 2, sizeof(Complex), &b);
            clSetKernelArg(axby_, 3, sizeof(cl_mem),  &y.handle());
            clEnqueueNDRangeKernel(queue_, axby_, 1, nullptr, &Nglobal, nullptr, 0, nullptr, nullptr);
            clFinish(queue_);
            
            us_axby += timer.microseconds();
        };
        
        // solve using the CG solver
        ConjugateGradients < clArray<Complex>, clArrayView<Complex> > CG;
        n[ill] = CG.solve
        (
            rview,                  // rhs
            zview,                  // solution
            cmd_.prec_itertol,      // preconditioner tolerance
            0,                      // min. iterations
            Nsegsiz,                // max. iteration
            inner_prec,             // preconditioner
            inner_mmul,             // matrix multiplication
            false,                  // verbose output
            compute_norm,           // norm of an array
            scalar_product,         // scalar product of two arrays
            axby,                   // weighted sum of two arrays
            new_opencl_array        // allocate and connect a new array
        );
        
        // download data arrays from the GPU
        zview.EnqueueDownload(queue_);
        clFinish(queue_);
        
        // free GPU memory
        rview.disconnect();
        zview.disconnect();
        prec1a.disconnect();
        prec2a.disconnect();
        prec1b.disconnect();
        prec2b.disconnect();
        Dl1.disconnect();
        Dl2.disconnect();
        
        // unload blocks
        if (cmd_.outofcore)
        {
            const_cast<BlockArray<Complex>&>(r)[ill].drop();
            z.hdfsave(ill);
            z[ill].drop();
        }
        
        // release bock preconditioner
        this->CG_exit(ill);
    }
    
    // free GPU memory
    tmp.disconnect();
    nrm.disconnect();
    tmA.disconnect();
    D_p.disconnect();
    S_p.disconnect();
    Mm1_tr_p.disconnect();
    Mm2_p.disconnect();
    for (int lambda = 0; lambda <= s_rad_.maxlambda(); lambda++)
    {
        Mi_L[lambda].disconnect();
        Mi_mLm1[lambda].disconnect();
    }
    
    // broadcast inner preconditioner iterations
    par_.sync(n.data(), 1, l1_l2_.size());
    
    // inner preconditioner info (max and avg number of iterations)
    std::cout << " | ";
    std::cout << std::setw(5) << (*std::min_element(n.begin(), n.end()));
    std::cout << std::setw(5) << (*std::max_element(n.begin(), n.end()));
    std::cout << std::setw(5) << format("%g", std::accumulate(n.begin(), n.end(), 0) / float(n.size()));
    
    // GPU kernel timing
    std::size_t us_total = us_axby + us_mmul + us_norm + us_prec + us_spro;
    std::cout << " [prec: " << format("%2d", int(us_prec * 100. / us_total)) << "%"
              << ", mul1: " << format("%2d", int(us_mmul_1 * 100. / us_total)) << "%"
              << ", mul2: " << format("%2d", int((us_mmul-us_mmul_1) * 100. / us_total)) << "%"
              << ", axby: " << format("%2d", int(us_axby * 100. / us_total)) << "%"
              << ", norm: " << format("%2d", int(us_norm * 100. / us_total)) << "%"
              << ", spro: " << format("%2d", int(us_spro * 100. / us_total)) << "%"
              << "]";
}

#endif
