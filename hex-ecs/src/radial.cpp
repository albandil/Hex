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

#include <algorithm>
#include <complex>
#include <cmath>
#include <cstdio>

#include <omp.h>

#include "arrays.h"
#include "bspline.h"
#include "gauss.h"
#include "parallel.h"
#include "radial.h"

inline double damp (Complex r, Complex R)
{
    // if sufficiently far, return clean zero
    if (r.real() > R.real())
        return 0.;
    
    // else damp using tanh(x) distribution
    return std::tanh(0.125 * (R.real() - r.real()));
}

rArray RadialIntegrals::computeScale (int lambda, int iknotmax) const
{
    // get last knot (end of last interval)
    if (iknotmax == 0)
        iknotmax = bspline_.Nknot() - 1;
    
    // quadrature order
    // NOTE : must match that in RadialIntegrals::computeMi !
    int Npts = std::max(2, bspline_.order() + std::abs(lambda) + 1);
    
    // output arrays
    rArray data (iknotmax);
    
    // for all knots
    for (int iknot = 0; iknot < iknotmax; iknot++)
    {
        // skip zero-length intervals
        if (bspline_.t(iknot) == bspline_.t(iknot + 1))
            continue;
        
        // get (real part of) the left-most Gauss-Legendre point
        double rho = g_.p_points(Npts, bspline_.t(iknot), bspline_.t(iknot+1))[0].real();
        
        // compute logarithms of the scale factors
        data[iknot] = log(rho);
    }
    
    return data;
}

void RadialIntegrals::M_integrand
(
    int n,
    Complex * const restrict in,
    Complex * const restrict out,
    int i, int j, int a,
    int iknot, int iknotmax,
    double & logscale
) const
{
    // extract data
    Complex R = bspline_.t(iknotmax);
    
    // evaluate B-splines
    Complex values_i[n], values_j[n];
    bspline_.B(i, iknot, n, in, values_i);
    bspline_.B(j, iknot, n, in, values_j);
    
    // all evaluations produced finite results
    bool all_finite = true;
    
    // fill output array
    if (R != 0.)
    {
        for (int k = 0; k < n; k++)
        {
            out[k] = values_i[k] * values_j[k] * pow(in[k],a) * damp(in[k],R);
            
            if (not (all_finite = Complex_finite(out[k])))
                break;
        }
    }
    else
    {
        for (int k = 0; k < n; k++)
        {
            out[k] = values_i[k] * values_j[k] * pow(in[k],a);
            
            if (not (all_finite = Complex_finite(out[k])))
                break;
        }
    }
    
    //
    // check that all elements are finite
    //
    
    if (not all_finite)
    {
        // compute logarithms of the integrand
        if (R != 0.)
        {
            for (int k = 0; k < n; k++)
            {
                out[k] = std::log(values_i[k]) + std::log(values_j[k]) + double(a) * std::log(in[k]) + std::log(damp(in[k],R));
                
                // use newly computed value as scale, if larger than the current value
                if (out[k].real() > logscale)
                    logscale = out[k].real();
            }
        }
        else
        {
            for (int k = 0; k < n; k++)
            {
                out[k] = std::log(values_i[k]) + std::log(values_j[k]) + double(a) * std::log(in[k]);
                
                // use newly computed value as scale, if larger than the current value
                if (out[k].real() > logscale)
                    logscale = out[k].real();
            }
        }
        
        // scale values by subtracting logarithms and exponentialize, so that they can be integrated by weighed summing
        for (int k = 0; k < n; k++)
            out[k] = std::exp(Complex(out[k].real() - logscale, out[k].imag()));
    }
}

cArray RadialIntegrals::computeMi (int a, int iknotmax) const
{
    int Nspline = bspline_.Nspline();
    int order = bspline_.order();
    
    int i, j, iknot;
    std::size_t size = Nspline * (2 * order + 1) * (order + 1);
    
    // (logarithms of) partial integral moments
    cArray m (size, Complex(0.,special::constant::Inf));
    
    // for all B-splines
    for (i = 0; i < (int)Nspline; i++)
    {
        // for all B-splines (use symmetry)
        for (j = i; j <= i + (int)order and j < (int)Nspline; j++)
        {
            // determine relevant knots
            int ileft = j;
            int iright = i + order + 1;
            
            // "right" has to be smaller than "Rmax"
            if (iright > iknotmax)
                iright = iknotmax;
            
            // for all relevant knots
            for (iknot = ileft; iknot < iright; iknot++)
            {
                // get integration boundaries
                Complex xa = bspline_.t(iknot);
                Complex xb = bspline_.t(iknot+1);
                
                // throw away zero length intervals
                if (xa == xb)
                    continue;
                
                // results of the quadrature
                Complex integral;
                double logscale = 0.; // logarithm of the scale
                
                // use at least 2nd order
                int points = std::max(2, order + std::abs(a) + 1);
                
                // integrate
                integral = g_.quadMFP
                (
                    this, &RadialIntegrals::M_integrand,       // integrand pointer
                    points, iknot, xa, xb,                     // integration parameters
                    i, j, a, iknot, iknotmax, logscale         // data to pass to the integrator
                );
                
                // get the coordinates in m-matrix
                int x_1 = i;                // reference spline is i-th
                int y_1 = j - (i - order);
                int z_1 = iknot - i;
                
                // get the coordinates in m-matrix of the symmetric case
                int x_2 = j;                // reference spline is j-th
                int y_2 = i - (j - order);
                int z_2 = iknot - j;
                
                // save to m- and s-matrix
                Complex lg = ((integral == 0.) ? Complex(0.,special::constant::Inf) : std::log(integral) + logscale);
                m[(x_1 * (2 * order + 1) + y_1) * (order + 1) + z_1] = lg;
                m[(x_2 * (2 * order + 1) + y_2) * (order + 1) + z_2] = lg;
            }
        }
    }
    
    return m;
}


Complex RadialIntegrals::computeD_iknot (int i, int j, int iknot) const
{
    if (iknot < 0)
        iknot = bspline_.Nknot() - 1;
    
    // get interval boundaries
    Complex x1 = bspline_.t(iknot);
    Complex x2 = bspline_.t(iknot + 1);
    
    // throw away zero-length intervals
    if (x1 == x2)
        return 0;
    
    // get Gauss-Legendre nodes and weights for the interval [-1, 1]
    // - use at least 2nd order
    int points = std::max (2, bspline_.order());
    cArray xs = g_.p_points(points, x1, x2);
    cArray ws = g_.p_weights(points, x1, x2);
    
    // evaluate B-splines at Gauss-Legendre nodes
    Complex values_i[points], values_j[points];
    bspline_.dB(i, iknot, points, xs.data(), values_i);
    bspline_.dB(j, iknot, points, xs.data(), values_j);
    
    // result
    Complex res = 0;
    
    // accumulate the result
    for (int k = 0; k < points; k++)
        res += values_i[k] * values_j[k] * ws[k];
    
    return res;
}

Complex RadialIntegrals::computeD (int i, int j, int maxknot) const
{
    // get boundary iknots
    int left = std::max(i, j);
    int right = std::min(i, j) + bspline_.order();
    
    // cut at maxknot
    if (right > maxknot)
        right = maxknot;
    
    // the result
    Complex res = 0;
    
    // undergo integration on sub-intervals
    for (int iknot = left; iknot <= right; iknot++)
        res += computeD_iknot(i, j, iknot);
        
    return res;
}

Complex RadialIntegrals::computeM_iknot (int a, int i, int j, int iknot, Complex R) const
{
    // get interval boundaries
    Complex x1 = bspline_.t(iknot);
    Complex x2 = bspline_.t(iknot + 1);
    
    // throw away zero-length intervals
    if (x1 == x2)
        return 0;
    
    // get Gauss-Legendre nodes and weights for the interval [-1, 1]
    // - use at least 2nd order
    int points = std::max (2, bspline_.order() + std::abs(a) + 1);
    cArray xs = g_.p_points(points, x1, x2);
    cArray ws = g_.p_weights(points, x1, x2);
    
    // evaluate B-splines at Gauss-Legendre nodes
    Complex values_i[points], values_j[points];
    bspline_.B(i, iknot, points, xs.data(), values_i);
    bspline_.B(j, iknot, points, xs.data(), values_j);
    
    // result
    Complex res = 0;
    
    // accumulate the (damped) result
    if (R != 0.)
    {
        for (int k = 0; k < points; k++)
            res += values_i[k] * values_j[k] * pow(xs[k],a) * ws[k] * damp(xs[k],R);
    }
    else
    {
        for (int k = 0; k < points; k++)
            res += values_i[k] * values_j[k] * pow(xs[k],a) * ws[k];
    }
    
    return res;
}

Complex RadialIntegrals::computeM (int a, int i, int j, int maxknot) const
{
    // get boundary iknots
    int left = std::max(i, j);
    int right = std::min(i, j) + bspline_.order();
    
    // cut at maxknot
    if (maxknot != 0 and right > maxknot)
        right = maxknot;
    
    // the result
    Complex res = 0;
    
    // undergo integration on sub-intervals
    for (int iknot = left; iknot <= right; iknot++)
        res += computeM_iknot(a, i, j, iknot, bspline_.t(maxknot));
    
    return res;
}

void RadialIntegrals::setupOneElectronIntegrals ()
{
    // create file names for this radial integrals
    char D_name[20], S_name[20], Mm1_name[20], Mm1_tr_name[20], Mm2_name[20];
    snprintf (D_name,      sizeof(D_name),      "%d-D.hdf",      bspline_.order());
    snprintf (S_name,      sizeof(S_name),      "%d-S.hdf",      bspline_.order());
    snprintf (Mm1_name,    sizeof(Mm1_name),    "%d-Mm1.hdf",    bspline_.order());
    snprintf (Mm1_tr_name, sizeof(Mm1_tr_name), "%d-Mm1_tr.hdf", bspline_.order());
    snprintf (Mm2_name,    sizeof(Mm2_name),    "%d-Mm2.hdf",    bspline_.order());
    
    // load/compute derivative overlaps
    std::cout << "Loading/precomputing derivative overlaps... " << std::flush;
    D_.hdfload(D_name) or D_.populate (
        bspline_.order(), [=](int i, int j) -> Complex { return computeD(i, j, bspline_.Nknot() - 1); }
    ).hdfsave(D_name);
    
    // load/compute integral moments
    std::cout << "ok\n\nLoading/precomputing integral moments... " << std::flush;
    S_.hdfload(S_name) or S_.populate (
        bspline_.order(), [=](int m, int n) -> Complex { return computeM(0, m, n); }
    ).hdfsave(S_name);
    Mm1_.hdfload(Mm1_name) or Mm1_.populate (
        bspline_.order(), [=](int m, int n) -> Complex { return computeM(-1, m, n); }
    ).hdfsave(Mm1_name);
    Mm1_tr_.hdfload(Mm1_tr_name) or Mm1_tr_.populate (
        bspline_.order(),    [=](int m, int n) -> Complex { return computeM(-1, m, n, bspline_.Nreknot() - 1);}
    ).hdfsave(Mm1_tr_name);
    Mm2_.hdfload(Mm2_name) or Mm2_.populate (
        bspline_.order(), [=](int m, int n) -> Complex { return computeM(-2, m, n); }
    ).hdfsave(Mm2_name);
    std::cout << "ok\n\n";
}

void RadialIntegrals::setupTwoElectronIntegrals (Parallel const & par, CommandLine const & cmd, Array<bool> const & lambdas)
{
    // allocate storage
    R_tr_dia_.resize(lambdas.size());
    
/*#ifndef NO_OPENCL
    // prepare GPU kernels for computation of diagonal R-integrals
    setup_gpu_();
#endif*/
    
    // print information
#if defined(_OPENMP)
    #pragma omp parallel
    #pragma omp master
    std::cout << "Precomputing multipole integrals (λ = 0 .. " << lambdas.size() - 1 << ") using " 
              << omp_get_num_threads() << " threads." << std::endl;
#endif
    
    // shorthands
    int Nspline = bspline_.Nspline();
    int order = bspline_.order();
    
    // for all multipoles : compute / load
    for (int lambda = 0; lambda < (int)lambdas.size(); lambda++)
    {
        // this process will only compute a subset of radial integrals
        if (lambda % par.Nproc() != par.iproc())
            continue;
        
        // look for precomputed data on disk
        std::string R_name = format ("%d-R_tr_dia_%d.hdf", bspline_.order(), lambda);
        if (R_tr_dia_[lambda].hdfload(R_name))
        {
            std::cout << "\t- integrals for λ = " << lambda << " loaded from \"" << R_name << "\"\n";
            continue; // no need to compute
        }
        
        // logarithms of partial integral moments
        cArray Mtr_L, Mtr_mLm1;
        
        // compute partial moments
        Mtr_L    = std::move(computeMi( lambda,   bspline_.Nreknot() - 1));
        Mtr_mLm1 = std::move(computeMi(-lambda-1, bspline_.Nreknot() - 1));
        
/*#ifndef NO_OPENCL
        // set up outer quadrature rule for GPU
        std::size_t Nouter = bspline_.order() + lambda + 10;    // outer quadrature points count
        const double *pxOut, *pwOut;
        g_.GaussLegendreData::gauss_nodes_and_weights(Nouter, pxOut, pwOut);
        xOut0_.disconnect(); xOut0_ = rArrayView(Nouter, const_cast<double*>(pxOut)); xOut0_.connect(context_, CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY);
        wOut0_.disconnect(); wOut0_ = rArrayView(Nouter, const_cast<double*>(pwOut)); wOut0_.connect(context_, CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY);
        
        // set up inner quadrature rule for GPU
        std::size_t Ninner = bspline_.order() + lambda + 1;     // inner quadrature points count
        const double *pxIn, *pwIn;
        g_.GaussLegendreData::gauss_nodes_and_weights(Ninner, pxIn,  pwIn);
        xIn0_.disconnect(); xIn0_ = rArrayView(Nouter, const_cast<double*>(pxIn)); xIn0_.connect(context_, CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY);
        wIn0_.disconnect(); wIn0_ = rArrayView(Nouter, const_cast<double*>(pwIn)); wIn0_.connect(context_, CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY);
#endif*/
        
        // elements of R_tr
        lArray R_tr_i, R_tr_j, th_R_tr_i, th_R_tr_j;
        cArray R_tr_v, th_R_tr_v;
        
/*#ifdef NO_OPENCL*/
        # pragma omp parallel default(none) \
            private (th_R_tr_i, th_R_tr_j, th_R_tr_v) \
            firstprivate (Nspline, order, lambda, Mtr_L, Mtr_mLm1) \
            shared (R_tr_i, R_tr_j, R_tr_v)
        {
            // for all B-spline pairs
            # pragma omp for schedule(dynamic,1)
/*#endif*/
            for (int i = 0; i < Nspline; i++)
            for (int j = 0; j < Nspline; j++)
            {
                // for all nonzero, nonsymmetry R-integrals
                for (int k = i; k <= i + order and k < Nspline; k++) // enforce i ≤ k
                for (int l = j; l <= j + order and l < Nspline; l++) // enforce j ≤ l
                {
                    // skip symmetry ijkl <-> jilk (others are accounted for in the limits)
                    if (i > j and k > l)
                        continue;
                    
                    // evaluate B-spline integral
                    Complex Rijkl_tr = computeR (lambda, i, j, k, l, Mtr_L, Mtr_mLm1);
                    
                    // store all symmetries
                    allSymmetries(i, j, k, l, Rijkl_tr, th_R_tr_i, th_R_tr_j, th_R_tr_v);
                }
            }
            
            # pragma omp critical
            {
                // merge the thread local arrays
                R_tr_i.append(th_R_tr_i.begin(), th_R_tr_i.end());
                R_tr_j.append(th_R_tr_j.begin(), th_R_tr_j.end());
                R_tr_v.append(th_R_tr_v.begin(), th_R_tr_v.end());
            }
/*#ifdef NO_OPENCL*/
        }
/*#endif*/
        
        // create matrices and save them to disk; use only upper part of the matrix R as we haven't computed whole lower part anyway
        R_tr_dia_[lambda] = CooMatrix(Nspline*Nspline, Nspline*Nspline, R_tr_i, R_tr_j, R_tr_v).todia(upper);
        R_tr_dia_[lambda].hdfsave(R_name, true, 10);
        
        std::cout << "\t- integrals for λ = " << lambda << " computed\n";
    }
    
#ifndef NO_MPI
    // for all multipoles : synchronize
    if (par.active())
    {
        for (int lambda = 0; lambda < (int)lambdas.size(); lambda++)
        {
            // get owner process of this multipole
            int owner = lambda % par.Nproc();
            
            // get dimensions
            int diagsize = R_tr_dia_[lambda].diag().size();
            int datasize = R_tr_dia_[lambda].data().size();
            
            // owner will broadcast dimensions
            MPI_Bcast(&diagsize, 1, MPI_INT, owner, MPI_COMM_WORLD);
            MPI_Bcast(&datasize, 1, MPI_INT, owner, MPI_COMM_WORLD);
            
            // get arrays
            iArray diag = R_tr_dia_[lambda].diag();
            cArray data = R_tr_dia_[lambda].data();
            diag.resize(diagsize);
            data.resize(datasize);
            
            // master will broadcast arrays
            MPI_Bcast(&diag[0], diag.size(), MPI_INT, owner, MPI_COMM_WORLD);
            MPI_Bcast(&data[0], data.size(), MPI_DOUBLE_COMPLEX, owner, MPI_COMM_WORLD);
            
            // reconstruct objects
            R_tr_dia_[lambda] = SymDiaMatrix(Nspline * Nspline, diag, data);
            
            if (owner != par.iproc())
            {
                std::cout << "\t- integrals for λ = " << lambda << " acquired from process " << owner << "\n";
                
                // save to disk (if the file doesn't already exist)
                std::string R_name = format ("%d-R_tr_dia_%d.hdf", bspline_.order(), lambda);
                if (not HDFFile(R_name, HDFFile::readonly).valid())
                    R_tr_dia_[lambda].hdfsave(R_name, true, 10);
            }
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
    }
#endif
}

// kernels' source as byte array, generated by "xxd" from the CL source
char slater_cl [] = {
    #include "slater_cl.c"
    , 0x00 // terminate the string by zero
};

// pointer to the source; to be used in setup
char * slater_source = &slater_cl[0];

void RadialIntegrals::setup_gpu_ ()
{
    std::cout << "Setting up OpenCL environment" << std::endl;
    char text [1000];
    
    // use platform 0
    clGetPlatformIDs (1, &platform_, nullptr);
    clGetPlatformInfo (platform_, CL_PLATFORM_NAME, sizeof(text), text, nullptr);
    std::cout << "\tplatform: " << text << " ";
    clGetPlatformInfo (platform_, CL_PLATFORM_VENDOR, sizeof(text), text, nullptr);
    std::cout << "(" << text << ")" << std::endl;
    clGetPlatformInfo (platform_, CL_PLATFORM_VERSION, sizeof(text), text, nullptr);
    std::cout << "\tavailable version: " << text << std::endl;
    
    // use device 0
    clGetDeviceIDs (platform_, CL_DEVICE_TYPE_GPU, 1, &device_, nullptr);
//     clGetDeviceIDs (platform_, CL_DEVICE_TYPE_CPU, 1, &device_, nullptr);
    clGetDeviceInfo (device_, CL_DEVICE_NAME, sizeof(text), text, nullptr);
    std::cout << "\tdevice: " << text << " ";
    clGetDeviceInfo (device_, CL_DEVICE_VENDOR, sizeof(text), text, nullptr);
    std::cout << "(" << text << ")" << std::endl;
    cl_ulong size;
    clGetDeviceInfo (device_, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &size, 0);
    std::cout << "\tlocal memory size: " << size/1024 << " kiB" << std::endl;
    clGetDeviceInfo (device_, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &size, 0);
    std::cout << "\tglobal memory size: " << format("%.2f", size/pow(1024,3)) << " GiB " << std::endl;
    clGetDeviceInfo (device_, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_ulong), &size, 0);
    std::cout << "\tmax compute units: " << size << std::endl;
    clGetDeviceInfo (device_, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(cl_ulong), &size, 0);
    std::cout << "\tmax work group size: " << size << std::endl << std::endl;
    
    // create context and command queue
    context_ = clCreateContext (nullptr, 1, &device_, nullptr, nullptr, nullptr);
    queue_ = clCreateCommandQueue (context_, device_, 0, nullptr);
    
    // setup compile flags
    std::ostringstream flags;
    flags << " -cl-fast-relaxed-math ";
    flags << " -D ORDER="     << bspline_.order()   << " ";
    flags << " -D NSPLINE="   << bspline_.Nspline() << " ";
    flags << " -D NOUTMAX="   << bspline_.order()+maxlambda()+10 << " ";
    flags << " -D NINMAX="    << bspline_.order()+maxlambda()+1  << " ";
    flags << " -D IKNOTMAX="  << bspline_.Nreknot()-1 << " ";
    
    // build program
    program_ = clCreateProgramWithSource (context_, 1, const_cast<const char**>(&slater_source), nullptr, nullptr);
    
    std::cout << "Create program DONE.\n";
    std::cout << "Flags: " << flags.str() << std::endl;
    
    clBuildProgram (program_, 1, &device_, flags.str().c_str(), nullptr, nullptr);
    
    std::cout << "Build DONE.\n";
    
    cl_build_status status;
    clGetProgramBuildInfo(program_, device_, CL_PROGRAM_BUILD_STATUS, sizeof(status), &status, nullptr);
    if (status != CL_SUCCESS)
    {
        std::cout << std::endl << "Source:" << std::endl << slater_source << std::endl;
        std::cout << std::endl << "Command line:" << std::endl << flags.str() << std::endl << std::endl;
        
        char log [100000];
        clGetProgramBuildInfo(program_, device_, CL_PROGRAM_BUILD_LOG, sizeof(log), log, nullptr);
        std::cout << "clGetProgramBuildInfo: log" << std::endl << log << std::endl;
        
        throw exception ("Failed to initialize OpenCL.");
    }
    
    // set program entry point
    rint_ = clCreateKernel(program_, "R_integral", nullptr);
    
    std::cout << "Kernel DONE.\n";
    
    // create and connect B-spline knot array
    t_ = bspline_.t();
    t_.connect(context_, CL_MEM_COPY_HOST_PTR | CL_MEM_WRITE_ONLY);
    
    // create and connect the data array
    R_gpu_ = cArray(bspline_.Nreknot() - 1);
    R_gpu_.connect(context_, CL_MEM_COPY_HOST_PTR | CL_MEM_WRITE_ONLY);
    
    std::cout << "Setup DONE\n";
}
