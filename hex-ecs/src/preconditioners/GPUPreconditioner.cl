//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2017, Jakub Benda, Charles University in Prague                    //
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

// Enable double precision (redundant in OpenCL 2.0).
#pragma OPENCL EXTENSION cl_khr_fp64: enable

// Necessary compile-time definitions:
// -D ORDER=... (implies local size in "mmul_2el")
// -D NSPLINE_ATOM=... (needed by "mmul_1el", "mmul_2el", "mul_ABt" and "kron_div")
// -D NSPLINE_PROJ=... (needed by "mmul_1el", "mmul_2el", "mul_ABt" and "kron_div")
// -D NREKNOT_ATOM=... (needed by "mmul_2el")
// -D NREKNOT_PROJ=... (needed by "mmul_2el")
// -D NLOCAL=... (local size in "scalar_product", "norm" and "mmul_1el")
// -D NBLOCK_SIZE=... (block size and local size in "mul_ABt")
// -D ANGULAR_BASIS_SIZE=... (number of coupled angular states, used in offset routines)
// -D PROJECTILE_BASIS_OFFSET=... (skipped projectile basis splines, used in "mmul_2el")
// -D NSRCSEG=...
// -D NDSTSEG=...
// -D ZA=... (atom nucleus charge, hydrogen = +1)
// -D ZP=... (projectile charge, electron = -1)

// The following section is only for the syntax highlighting purposes in KDevelop:
#ifndef ORDER

    #include <clc/clc.h>

    #ifndef Real
    #define Real float
    #endif

    #ifndef Complex
    #define Complex float2
    #endif

    #ifndef ORDER
    #define ORDER 4
    #endif

    #ifndef NLOCAL
    #define NLOCAL 1
    #endif

    #ifndef NSPLINE_ATOM
    #define NSPLINE_ATOM 1
    #endif

    #ifndef NSPLINE_PROJ
    #define NSPLINE_PROJ 1
    #endif

    #ifndef BLOCK_SIZE
    #define BLOCK_SIZE 16
    #endif

    #ifndef ZA
    #define ZA 1
    #endif

    #ifndef ZP
    #define ZP -1
    #endif

#endif

// Derived variables.
#define NROW (NSPLINE_ATOM * NSPLINE_PROJ)

// These are sometimes not defined for double.
// #define min(x,y) (x < y ? x : y)
// #define max(x,y) (x > y ? x : y)

/**
 * @brief Complex multiplication.
 * 
 * Multiplies two complex numbers and returns their product.
 */
Complex cmul (private Complex a, private Complex b)
{
    private Complex c;
    
    c.x = a.x * b.x - a.y * b.y;
    c.y = a.x * b.y + a.y * b.x;
    
    return c;
}

/**
 * @brief Complex division.
 * 
 * Divides two complex numbers and returns the fraction.
 * No overflow checking is done.
 */
Complex cdiv (private Complex a, private Complex b)
{
    private Complex c;
    private Real b2 = b.x * b.x + b.y * b.y;
    
    c.x = (a.x * b.x + a.y * b.y) / b2;
    c.y = (a.y * b.x - a.x * b.y) / b2;
    
    return c;
}

/**
 * @brief Integer power.
 * 
 * Fast integer power; uses only multiplications.
 * The number of multiplications is proportional to log(n).
 * 
 * This is faster than "pown" when all workitems have the same exponent.
 */
Real pow_int (private Real x, private int n)
{
    private Real value = 1; // = x^0
    
    do
    {
        if (n % 2 == 1)
            value *= x;
        
        n /= 2;
        x *= x;
    }
    while (n != 0);
    
    return value;
}

/**
 * @brief AXBY operation.
 * 
 * Computes the linear combination @f$ ax + by @f$ of vectors @f$ x @f$ and 
 * @f$ y @f$ and stores the result in @f$ x @f$.
 * \f[
 *         x_i = a x_i + b y_i
 * \f]
 * 
 * @param a Complex factor.
 * @param x Source and destination vector.
 * @param b Complex factor.
 * @param y Source vector.
 */
kernel void a_vec_b_vec (private Complex a, global Complex *x, private Complex b, global Complex *y)
{
    private int i = get_global_id(0);
    
    if (i < NROW)
        x[i] = cmul(a,x[i]) + cmul(b,y[i]);
}

/**
 * @brief Full scalar product.
 * 
 * Calculates element-wise product of the arrays and reduces them using local memory
 * to one number per work group. These intermediate numbers are then summed by CPU.
 * 
 * @param u Source vector.
 * @param v Source vector.
 * @param z Output vector for intermediate segment scalar products.
 */
kernel void scalar_product (global Complex *u, global Complex *v, global Complex *z)
{
    // position of this thread among other threads
    private int iglobal = get_global_id(0);
    private int ilocal = get_local_id(0);
    
    // calculate product of array elements
    local Complex uv[NLOCAL];
    uv[ilocal] = (iglobal < NROW ? cmul(u[iglobal],v[iglobal]) : (Complex)(0,0));
    barrier(CLK_LOCAL_MEM_FENCE);
    
    // reduce the per-element products
    int stride = 1;
    while (stride < NLOCAL)
    {
        if (ilocal % (2 * stride) == 0 && iglobal + stride < NROW)
            uv[ilocal] += uv[ilocal + stride];
        stride *= 2;
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    
    // write the reduced scalar product as computed by this group
    if (ilocal == 0)
        z[get_group_id(0)] = uv[0];
}

/**
 * @brief Full vector norm.
 * 
 * Calculates scalar product of the segments of the arrays belonging to different groups
 * as one number per work group. These intermediate numbers are then summed by CPU.
 * 
 * @param v Source vector.
 * @param z Output vector for intermediate segment norms.
 */
kernel void norm (global Complex *v, global Real *z)
{
    // position of this thread among other threads
    private int iglobal = get_global_id(0);
    private int ilocal = get_local_id(0);
    
    // get element to process by this thread
    private Complex vi = v[iglobal];
    
    // calculate squared modulus of an array element
    local Real vv[NLOCAL];
    vv[ilocal] = (iglobal < NROW ? vi.x * vi.x + vi.y * vi.y : 0);
    barrier(CLK_LOCAL_MEM_FENCE);
    
    // reduce the per-element products
    private int stride = 1;
    while (stride < NLOCAL)
    {
        if (ilocal % (2 * stride) == 0 && iglobal + stride < NROW)
            vv[ilocal] += vv[ilocal + stride];
        stride *= 2;
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    
    // write the reduced scalar product as computed by this group
    if (ilocal == 0)
        z[get_group_id(0)] = vv[0];
}

/**
 * @brief Multiplication by one-electron Hamiltonian matrix.
 * 
 * Multiplies given vector by one-electron part of the Hamiltonian matrix,
 * that can expressed as a Kronecker product of simple one-electron matrices.
 * \f[
 *     y_{ij} = \sum_{kl} \left(
 *         E S_{ik} S_{jl} - H_{ik}^{(1)} S_{jl} - S_{ik} H_{jl}^{(2)}
 *     \right) x_{kl}
 * \f]
 * @param E Total energy of the system.
 * @param Spa Row-padded upper overlap matrix for atomic electron.
 * @param Dpa Row-padded upper derivative overlap matrix for atomic electron.
 * @param M1pa Row-padded upper integral moment matrix (for r^1) for atomic electron.
 * @param M2pa Row-padded upper integral moment matrix (for r^2) for atomic electron.
 * @param Spp Row-padded upper overlap matrix for projectile electron.
 * @param Dpp Row-padded upper derivative overlap matrix for projectile electron.
 * @param M1pp Row-padded upper integral moment matrix (for r^1) for projectile electron.
 * @param M2pp Row-padded upper integral moment matrix (for r^2) for projectile electron.
 * @param l1 Angular momentum of atomic electron.
 * @param l2 Angular momentum of projectile electron.
 * @param x Source vector.
 * @param y Destination vector.
 */
kernel void mmul_1el
(
    // energy
    private Real E,
    // row-padded one-electron matrices (atom)
    global Complex const * const restrict Spa,
    global Complex const * const restrict Dpa,
    global Complex const * const restrict M1pa,
    global Complex const * const restrict M2pa,
    // row-padded one-electron matrices (projectile)
    global Complex const * const restrict Spp,
    global Complex const * const restrict Dpp,
    global Complex const * const restrict M1pp,
    global Complex const * const restrict M2pp,
    // angular momenta
    private int l1,
    private int l2,
    // source and target vector
    global Complex const * const restrict x,
    global Complex       * const restrict y
)
{
    // output vector element index
    private int i = get_global_id(0) / NSPLINE_PROJ;
    private int j = get_global_id(0) % NSPLINE_PROJ;
    
    // initialize the output element
    y[i * NSPLINE_PROJ + j] = 0;
    
    // for all source vector elements
    if (i < NSPLINE_ATOM)
    for (private int k = i - ORDER; k <= i + ORDER; k++) if (0 <= k && k < NSPLINE_ATOM)
    for (private int l = j - ORDER; l <= j + ORDER; l++) if (0 <= l && l < NSPLINE_PROJ)
    {
        // compute multi-indices
        private int ik = min(i,k) * (ORDER + 1) + abs(i - k);
        private int jl = min(j,l) * (ORDER + 1) + abs(l - j);
        
        // calculate the one-electron part of the hamiltonian matrix element Hijkl
        private Complex elem = E * cmul(Spa[ik],Spp[jl]);
        elem -= (Real)(0.5f) * (cmul(Dpa[ik],Spp[jl]) + cmul(Spa[ik],Dpp[jl]));
        elem -= (Real)(0.5f) * l1 * (l1 + 1) * cmul(M2pa[ik],Spp[jl]) + (Real)(0.5f) * l2 * (l2 + 1) * cmul(Spa[ik],M2pp[jl]);
        elem += ZA * cmul(M1pa[ik],Spp[jl]) - ZP * cmul(Spa[ik],M1pp[jl]);
        
        // multiply right-hand side by that matrix element
        y[i * NSPLINE_PROJ + j] += cmul(elem, x[k * NSPLINE_PROJ + l]);
    }
}

/**
 * @brief KPA preconditioner.
 * 
 * Divides the source vector by the KPA preconditioner term.
 * \f[
 *     y_{ij} = \frac{y_{ij}}{E - D_i - D_j}
 * \f]
 * @param E Total energy.
 * @param D1 Eigenvalues for first angular momentum.
 * @param D2 Eigenvalues for second angular momentum.
 * @param y Source/destination vector.
 */
kernel void kron_div
(
    private Complex E,
    global Complex const * const restrict D1,
    global Complex const * const restrict D2,
    global Complex       * const restrict y
)
{
    // get worker's offset index (0 <= j < NSPLINE_PROJ)
    private int j = get_global_id(0);
    
    // y = y / (E (I kron I) - (D1 kron I) - (I kron D2))
    for (private int i = 0; i < NSPLINE_ATOM; i++)
        y[i * NSPLINE_PROJ + j] = cdiv(y[i * NSPLINE_PROJ + j], E - D1[i] - D2[j]);
}

Real unrotate (Complex x)
{
    return 0; // FIXME
}

Real clamp (Real x, Real a, Real b)
{
    return max(a, min(x, b));
}

/**
 * @brief Multiply by two-electron integral matrix.
 * 
 * Multiplies the source vector by the modified two-electron matrix R,
 * where only the fully decoupled integrals are considered.
 */
kernel void mmul_2el_decoupled_splines
(
    // B-spline knots (atom and projectile)
    constant Complex const * const restrict ta,
    constant Complex const * const restrict tp,
    // multipole
    private int lambda,
    // angular integral
    private Real f,
    // one-electron moments (atomic electron)
    global Complex const * const restrict MLa,
    global Complex const * const restrict MmLm1a,
    // one-electron moments (projectile electron)
    global Complex const * const restrict MLp,
    global Complex const * const restrict MmLm1p,
    // source and target vector
    global Complex const * const restrict x,
    global Complex       * const restrict y
)
{
    // B-spline indices
    uint a = get_global_id(0) / NSPLINE_PROJ;
    uint b = get_global_id(0) % NSPLINE_PROJ;
    uint c = get_global_id(1) / NSPLINE_PROJ;
    uint d = get_global_id(1) % NSPLINE_PROJ;
    
    // B-spline bounding knots
    Real ta1 = unrotate(ta[a]), ta2 = unrotate(ta[a + ORDER + 1]);
    Real tb1 = unrotate(tp[b]), tb2 = unrotate(tp[b + ORDER + 1]);
    Real tc1 = unrotate(ta[c]), tc2 = unrotate(ta[c + ORDER + 1]);
    Real td1 = unrotate(tp[d]), td2 = unrotate(tp[d + ORDER + 1]);
    
    // Do nothing when there is no overlap.
    if (max(ta1,tc1) >= min(ta2,tc2) || max(tb1,td1) >= min(tb2,td2))
        return;
    
    Real t_xmin = clamp(max(ta1,tc1), rxmin_, rxmax_);
    Real t_xmax = clamp(min(ta2,tc2), rxmin_, rxmax_);
    Real t_ymin = clamp(max(tb1,td1), rymin_, rymax_);
    Real t_ymax = clamp(min(tb2,td2), rymin_, rymax_);
    
    if (t_ymax <= t_xmin)
    {
        Real scale = pow_int(t_ymax / t_xmax, lambda) / t_xmax;
        return scale * MmLm1a[lambda](a,c) * MLp[lambda](b,d);
    }
    
    if (t_xmax <= t_ymin)
    {
        Real scale = pow_int(t_xmax / t_ymax, lambda) / t_ymax;
        return scale * MLa[lambda](a,c) * MmLm1p[lambda](b,d);
    }
}

/**
 * @brief Multiply by two-electron integral matrix.
 * 
 * Multiplies the source vector by the modified two-electron matrix R,
 * where only the cell-decoupled part of non-decoupled integrals are considered.
 */
kernel void mmul_2el_decoupled_cells
(
    
)
{
    
}

/**
 * @brief Multiply by two-electron integral matrix.
 * 
 * Multiplies the source vector by the modified two-electron matrix R,
 * where only the cell-coupled part of non-decoupled integrals are considered.
 */
kernel void mmul_2el_coupled_cells
(
    
)
{
    
}
