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

#include <algorithm>
#include <complex>

#include <omp.h>

#include "arrays.h"
#include "bspline.h"
#include "gauss.h"
#include "moments.h"
#include "parallel.h"

/**
 * Potential suppressing factor. 
 * @param r Radius.
 * @param R Truncation radius.
 */
inline double damp(Complex r, Complex R)
{
    // if sufficiently far, return clean zero
    if (r.real() > R.real())
        return 0.;
    
    // else damp using tanh(x) distribution
    return tanh(0.125 * (R.real() - r.real()));
}

// ----------------------------------------------------------------------- //
//  Computation of B-spline moments                                        //
// ----------------------------------------------------------------------- //

void RadialIntegrals::M_integrand(int n, Complex * const restrict in, Complex * const restrict out, void *data) const
{
    // extract data
    int i = ((int*)data)[0];
    int j = ((int*)data)[1];
    int a = ((int*)data)[2];
    int iknot = ((int*)data)[3];
    int iknotmax = ((int*)data)[4];
    Complex R = bspline_.t(iknotmax);
    
    // evaluate B-splines
    Complex values_i[n], values_j[n];
    bspline_.B(i, iknot, n, in, values_i);
    bspline_.B(j, iknot, n, in, values_j);
    
    // fill output array
    int k;
    if (R != 0.)
    {
        for (k = 0; k < n; k++)
            out[k] = values_i[k] * values_j[k] * pow(in[k],a) * damp(in[k],R);
    }
    else
    {
        for (k = 0; k < n; k++)
            out[k] = values_i[k] * values_j[k] * pow(in[k],a);
    }
}

/** 
 * Compute two-spline integral moment of degree "a"
 * @param a Moment degree.
 * @param iknotmax Truncation knot.
 */
cArray RadialIntegrals::computeMi(int a, int iknotmax) const
{
    int Nspline = bspline_.Nspline();
    int order = bspline_.order();
    
    int i, j, iknot;
    size_t size = Nspline * (2 * order + 1) * (order + 1);
    cArray m(size);
    
    // for all B-splines
    for (i = 0; i < (int)Nspline; i++)
    {
        // for all B-splines (use symmetry)
        for (j = i; j <= i + (int)order and j < (int)Nspline; j++)
        {
            // determine relevant knots
            int ileft = (i < j) ? j : i;
            int iright = ((i < j) ? i : j) + order + 1;
            
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
                
                // compute the moment
                int points = order + abs(a) + 1;
                int data[5] = {i, j, a, iknot, iknotmax};
                Complex integral = g_.quadMFP (
                    this, &RadialIntegrals::M_integrand,
                    data, points, iknot, xa, xb
                );
                
                // get the coordinates in m-matrix
                int x_1 = i;                // reference spline is i-th
                int y_1 = j - (i - order);
                int z_1 = iknot - i;
                
                // get the coordinates in m-matrix of the symmetric case
                int x_2 = j;                // reference spline is j-th
                int y_2 = i - (j - order);
                int z_2 = iknot - j;
                
                // save to m-matrix
                m[(x_1 * (2 * order + 1) + y_1) * (order + 1) + z_1] = integral;
                m[(x_2 * (2 * order + 1) + y_2) * (order + 1) + z_2] = integral;
            }
        }
    }
    
    return m;
}


Complex RadialIntegrals::computeD_iknot(int i, int j, int iknot) const
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
    int points = bspline_.order();
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

Complex RadialIntegrals::computeD(int i, int j, int maxknot) const
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

Complex RadialIntegrals::computeM_iknot(int a, int i, int j, int iknot, Complex R) const
{
    // get interval boundaries
    Complex x1 = bspline_.t(iknot);
    Complex x2 = bspline_.t(iknot + 1);
    
    // throw away zero-length intervals
    if (x1 == x2)
        return 0;
    
    // get Gauss-Legendre nodes and weights for the interval [-1, 1]
    int points = bspline_.order() + std::abs(a) + 1;
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

Complex RadialIntegrals::computeM(int a, int i, int j, int maxknot) const
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

void RadialIntegrals::setupOneElectronIntegrals()
{
    std::cout << "Loading/precomputing derivative overlaps... ";
    
    D_.hdfload("D.hdf") or D_.populate (
        bspline_.order(), [=](int i, int j) -> Complex { return computeD(i, j, bspline_.Nknot() - 1); }
    ).hdfsave("D.hdf");
    
    std::cout << "ok\n\nLoading/precomputing integral moments... ";
    
    S_.hdfload("S.hdf") or S_.populate (
        bspline_.order(), [=](int m, int n) -> Complex { return computeM(0, m, n); }
    ).hdfsave("S.hdf");
    //
    Mm1_.hdfload("Mm1.hdf") or Mm1_.populate (
        bspline_.order(), [=](int m, int n) -> Complex { return computeM(-1, m, n); }
    ).hdfsave("Mm1.hdf");
    //
    Mm1_tr_.hdfload("Mm1_tr.hdf") or Mm1_tr_.populate (
        bspline_.order(),    [=](int m, int n) -> Complex { return computeM(-1, m, n, bspline_.Nreknot() - 1);}
    ).hdfsave("Mm1_tr.hdf");
    //
    Mm2_.hdfload("Mm2.hdf") or Mm2_.populate (
        bspline_.order(), [=](int m, int n) -> Complex { return computeM(-2, m, n); }
    ).hdfsave("Mm2.hdf");
    
    std::cout << "ok\n\n";
}

void RadialIntegrals::setupTwoElectronIntegrals(Parallel const & par, size_t maxlambda)
{
    // allocate storage
    R_tr_dia_.resize(maxlambda + 1);
    
    #pragma omp parallel
    {
        #pragma omp master
        {
            std::cout << "Precomputing multipole integrals (λ = 0 .. " 
                      << maxlambda
                      << ") using " 
                      << omp_get_num_threads() 
                      << " threads.\n";
        }
    }
    
    // shorthands
    int Nspline = bspline_.Nspline();
    int order = bspline_.order();
    
    // for all multipoles : compute / load
    for (int lambda = 0; lambda <= (int)maxlambda; lambda++)
    {
        // this process will only compute a subset of radial integrals
        if (lambda % par.Nproc() != par.iproc())
            continue;
        
        // look for precomputed data on disk
        std::ostringstream oss2;
        oss2 << "R_tr_dia_" << lambda << ".hdf";
        int R_integ_exists = R_tr_dia_[lambda].hdfload(oss2.str().c_str());
        
        if (R_integ_exists)
        {
            std::cout << "\t- integrals for λ = " << lambda << " loaded from \"" << oss2.str().c_str() << "\"\n";
            continue; // no need to compute
        }
        
        // precompute necessary partial integral moments
        // - truncated moments λ and -λ-1
        cArray Mtr_L    = computeMi( lambda,   bspline_.Nreknot()-1);
        cArray Mtr_mLm1 = computeMi(-lambda-1, bspline_.Nreknot()-1);
        
        // elements of R_tr
        lArray R_tr_i, R_tr_j, th_R_tr_i, th_R_tr_j;
        cArray R_tr_v, th_R_tr_v;
        
        # pragma omp parallel default(none) \
            private (th_R_tr_i, th_R_tr_j, th_R_tr_v) \
            firstprivate (Nspline, order, lambda, Mtr_L, Mtr_mLm1) \
            shared (R_tr_i, R_tr_j, R_tr_v)
        {
            // for all B-spline pairs
            # pragma omp for schedule(dynamic,1)
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
                    Complex Rijkl_tr = computeR(lambda, i, j, k, l, Mtr_L, Mtr_mLm1);
                    
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
        }
        
        // create matrices and save them to disk; use only upper part of the matrix R as we haven't computed whole lower part anyway
        R_tr_dia_[lambda] = CooMatrix(Nspline*Nspline, Nspline*Nspline, R_tr_i, R_tr_j, R_tr_v).todia(upper);
        R_tr_dia_[lambda].hdfsave(oss2.str().c_str(), true, 10);
        
        std::cout << "\t- integrals for λ = " << lambda << " computed\n";
    }
    
    // for all multipoles : synchronize
#ifndef NO_MPI
    if (par.active())
    {
        for (int lambda = 0; lambda <= (int)maxlambda; lambda++)
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
                std::cout << "\t- integrals for λ = " << lambda << " retrieved from process " << owner << "\n";
                
                // save to disk (if the file doesn't already exist)
                std::ostringstream oss3;
                oss3 << "R_tr_dia_" << lambda << ".hdf";
                if (not HDFFile(oss3.str().c_str(), HDFFile::readonly).valid())
                    R_tr_dia_[lambda].hdfsave(oss3.str().c_str(), true, 10);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
#endif
    
    std::cout << "\t- R_tr[λ] has " << R_tr_dia_[0].data().size() << " nonzero elements\n";
    // --------------------------------------------------------------------- //
}
