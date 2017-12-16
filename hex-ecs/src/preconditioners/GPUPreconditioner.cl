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

// --------------------------------------------------------------------------------- //

// This section is here only for purposes of syntax highlighting in KDevelop.
#ifdef SYNTAX_CHECK

    #include <clc/clc.h>

    // Necessary compile-time definitions
    
    #define Long      long  // or int
    #define Real     float  // or double
    #define Complex float2  // or double2
    #define ORDER        4  // B-spline order
    #define NLOCAL       1  // local size in "scalar_product", "norm" and "mmul_1el"
    #define NSPLINE_ATOM 1  // needed by "mmul_1el", "mmul_2el", "mul_ABt" and "kron_div"
    #define NSPLINE_PROJ 1  // needed by "mmul_1el", "mmul_2el", "mul_ABt" and "kron_div"
    #define BLOCK_SIZE  16  // block size and local size in "mul_ABt"
    #define ANGULAR_BASIS_SIZE   1  // coupled number of angular states
    #define ZA           1  // nuclear charge, hydrogen = +1
    #define ZP          -1  // projectile charge, electron = -1
    #define RXMIN        0  // potential bound
    #define RXMAX        0  // potential bound
    #define RYMIN        0  // potential bound
    #define RYMAX        0  // potential bound
    #define ROTFACT      1  // rotation factor, cos(ecs_theta)
    
#endif

// --------------------------------------------------------------------------------- //

// Derived variables.
#define NROW (NSPLINE_ATOM * NSPLINE_PROJ)

// --------------------------------------------------------------------------------- //

// These functions are often not properly defined (at least for double):

//- Minimum.
//#define min(x,y) (x < y ? x : y)

//- Maximum.
//#define max(x,y) (x > y ? x : y)

//- Clamp function (left & right bound).
Real myclamp (Real x, Real a, Real b)
{
    return max(a, min(x, b));
}

// --------------------------------------------------------------------------------- //

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
 * @brief Multiplication by BlockSymBandMatrix.
 * 
 * Multiplies a vector by BlockSymBandMatrix,
 * \f[
 *     y = A \cdot x
 * \f]
 */
kernel void mmul_simple (global Complex *A, global Complex *x, global Complex *y)
{
    // output vector element index
    private int i = get_global_id(0) / NSPLINE_PROJ;
    private int j = get_global_id(0) % NSPLINE_PROJ;
    
    // initialize the output element
    private Complex result = 0;
    
    // for all source vector elements
    if (i < NSPLINE_ATOM)
    for (private int k = i - ORDER; k <= i + ORDER; k++) if (0 <= k && k < NSPLINE_ATOM)
    for (private int l = j - ORDER; l <= j + ORDER; l++) if (0 <= l && l < NSPLINE_PROJ)
    {
        // compute multi-indices
        private int ik = min(i,k) * (ORDER + 1) + abs(i - k);
        private int jl = min(j,l) * (ORDER + 1) + abs(l - j);
        
        // get the matrix element
        private Complex elem = A[ik * (ORDER + 1) * NSPLINE_PROJ + jl];
        
        // multiply right-hand side by that matrix element
        result += cmul(elem, x[k * NSPLINE_PROJ + l]);
    }
    
    // push result to global memory
    if (i < NSPLINE_ATOM)
    {
        y[i * NSPLINE_PROJ + j] = result;
    }
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
        elem -= ZA * (-1) * cmul(M1pa[ik],Spp[jl]) + ZA * ZP * cmul(Spa[ik],M1pp[jl]);
        
        // multiply right-hand side by that matrix element
        y[i * NSPLINE_PROJ + j] += cmul(elem, x[k * NSPLINE_PROJ + l]);
    }
}

/**
 * @brief Multiplication by one-electron Hamiltonian matrix.
 * 
 * This is a variant of @ref mmul_1el, which operates on segments.
 */
kernel void mmul_1el_seg
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
    global Complex       * const restrict y,
    // block and segment
    private int ill,
    private int k
)
{
    // output vector element index
    private int j = get_global_id(0);
    
    // for all relevant target vector segments
    if (j < NSPLINE_PROJ)
    for (private int i = k - ORDER; i <= k + ORDER; i++) if (0 <= i && i <= NSPLINE_ATOM)
    {
        // initialize the output element
        y[(ill * (2 * ORDER + 1) + k + ORDER - i) * NSPLINE_PROJ + j] = 0;
        
        // for all relevant source vector elements
        for (private int l = j - ORDER; l <= j + ORDER; l++) if (0 <= l && l < NSPLINE_PROJ)
        {
            // compute multi-indices
            private int ik = min(i,k) * (ORDER + 1) + abs(i - k);
            private int jl = min(j,l) * (ORDER + 1) + abs(l - j);
            
            // calculate the one-electron part of the hamiltonian matrix element Hijkl
            private Complex elem = E * cmul(Spa[ik],Spp[jl]);
            elem -= (Real)(0.5f) * (cmul(Dpa[ik],Spp[jl]) + cmul(Spa[ik],Dpp[jl]));
            elem -= (Real)(0.5f) * l1 * (l1 + 1) * cmul(M2pa[ik],Spp[jl]) + (Real)(0.5f) * l2 * (l2 + 1) * cmul(Spa[ik],M2pp[jl]);
            elem -= ZA * (-1) * cmul(M1pa[ik],Spp[jl]) + ZA * ZP * cmul(Spa[ik],M1pp[jl]);
            
            // multiply right-hand side by that matrix element
            y[(ill * (2 * ORDER + 1) + k + ORDER - i) * NSPLINE_PROJ + j] += cmul(elem, x[ill * NSPLINE_PROJ + l]);
        }
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

Real unrotateX (private Complex A)
{
    // left complex part
    if (A.x < RXMIN)   return RXMIN + (A.x - RXMIN) / ROTFACT;
    
    // right complex part
    if (A.x > RXMAX)   return RXMAX + (A.x - RXMAX) / ROTFACT;
    
    // only real part
    return A.x;
}

Real unrotateY (private Complex A)
{
    // left complex part
    if (A.x < RYMIN)   return RYMIN + (A.x - RYMIN) / ROTFACT;
    
    // right complex part
    if (A.x > RYMAX)   return RYMAX + (A.x - RYMAX) / ROTFACT;
    
    // only real part
    return A.x;
}

/**
 * @brief Multiply by two-electron integral matrix.
 * 
 * Multiplies the source vector by the modified two-electron matrix R,
 * where only the fully decoupled integrals are considered.
 */
kernel void mmul_2el_decoupled
(
    // B-spline knots (atom and projectile)
    constant Complex const * const restrict ta,
    constant Complex const * const restrict tp,
    // angular integrals
    constant Real  const * const restrict f,
    // multipole
    private int lambda,
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
    private int a = get_global_id(0) / NSPLINE_PROJ;
    private int b = get_global_id(0) % NSPLINE_PROJ;
    
    // B-spline bounding knots
    private Real ta1 = unrotateX(ta[a]), ta2 = unrotateX(ta[a + ORDER + 1]);
    private Real tb1 = unrotateY(tp[b]), tb2 = unrotateY(tp[b + ORDER + 1]);
    
    // loop over free B-spline indices
    for (private int c = a > ORDER ? a - ORDER : 0; c < NSPLINE_ATOM && c <= a + ORDER; c++)
    for (private int d = b > ORDER ? b - ORDER : 0; d < NSPLINE_PROJ && d <= b + ORDER; d++)
    {
        // B-spline bounding knots
        private Real tc1 = unrotateX(ta[c]), tc2 = unrotateX(ta[c + ORDER + 1]);
        private Real td1 = unrotateY(tp[d]), td2 = unrotateY(tp[d + ORDER + 1]);
        
        if (max(ta1,tc1) >= min(ta2,tc2) || max(tb1,td1) >= min(tb2,td2))
            continue;
        
        // restrict effective radius
        private Real t_xmin = myclamp(max(ta1,tc1), RXMIN, RXMAX);
        private Real t_xmax = myclamp(min(ta2,tc2), RXMIN, RXMAX);
        private Real t_ymin = myclamp(max(tb1,td1), RYMIN, RYMAX);
        private Real t_ymax = myclamp(min(tb2,td2), RYMIN, RYMAX);
        
        // decoupled y < x
        if (t_ymax <= t_xmin)
        {
            private Real scale = pow_int(t_ymax / t_xmax, lambda) / t_xmax;
            private Complex R = scale * cmul(MmLm1a[min(a,c) * (ORDER + 1) + abs_diff(a,c)], MLp[min(b,d) * (ORDER + 1) + abs_diff(b,d)]);
            
            y[a * NSPLINE_PROJ + b] += ZP * f[lambda] * cmul(R, x[c * NSPLINE_PROJ + d]);
        }
        
        // decoupled x < y
        if (t_xmax <= t_ymin)
        {
            private Real scale = pow_int(t_xmax / t_ymax, lambda) / t_ymax;
            private Complex R = scale * cmul(MLa[min(a,c) * (ORDER + 1) + abs_diff(a,c)], MmLm1p[min(b,d) * (ORDER + 1) + abs_diff(b,d)]);
            
            y[a * NSPLINE_PROJ + b] += ZP * f[lambda] * cmul(R, x[c * NSPLINE_PROJ + d]);
        }
    }
}

/**
 * @brief Multiply by two-electron integral matrix.
 * 
 * Similar to @ref mmul_2el_decoupled, but operating on segments.
 */
kernel void mmul_2el_decoupled_segment
(
    // B-spline knots (atom and projectile)
    constant Complex const * const restrict ta,
    constant Complex const * const restrict tp,
    // angular integrals
    constant Real  const * const restrict f,
    // multipole
    private int lambda,
    // one-electron moments (atomic electron)
    global Complex const * const restrict MLa,
    global Complex const * const restrict MmLm1a,
    // one-electron moments (projectile electron)
    global Complex const * const restrict MLp,
    global Complex const * const restrict MmLm1p,
    // source and target vector
    global Complex const * const restrict x,
    global Complex       * const restrict y,
    // segment
    private int k
)
{
    // B-spline indices
    private int j = get_global_id(0);
    
    if (j < NSPLINE_PROJ)
    for (private int i = k > ORDER ? k - ORDER : 0; i < NSPLINE_ATOM && i <= k + ORDER; i++)
    for (private int l = j > ORDER ? j - ORDER : 0; l < NSPLINE_PROJ && l <= j + ORDER; l++)
    {
        // B-spline bounding knots
        private Real ti1 = unrotateX(ta[i]), ti2 = unrotateX(ta[i + ORDER + 1]);
        private Real tj1 = unrotateY(tp[j]), tj2 = unrotateY(tp[j + ORDER + 1]);
        private Real tk1 = unrotateX(ta[k]), tk2 = unrotateX(ta[k + ORDER + 1]);
        private Real tl1 = unrotateY(tp[l]), tl2 = unrotateY(tp[l + ORDER + 1]);
        
        if (max(ti1,tk1) >= min(ti2,tk2) || max(tj1,tl1) >= min(tj2,tl2))
            continue;
        
        // restrict effective radius
        private Real t_xmin = myclamp(max(ti1,tk1), RXMIN, RXMAX);
        private Real t_xmax = myclamp(min(ti2,tk2), RXMIN, RXMAX);
        private Real t_ymin = myclamp(max(tj1,tl1), RYMIN, RYMAX);
        private Real t_ymax = myclamp(min(tj2,tl2), RYMIN, RYMAX);
        
        // multi-indices
        private int ik = min(i,k) * (ORDER + 1) + abs_diff(i,k);
        private int jl = min(j,l) * (ORDER + 1) + abs_diff(j,l);
        
        // radial integral
        private Complex R = 0;
        
        // decoupled y < x
        if (t_ymax <= t_xmin)
        {
            private Real scale = pow_int(t_ymax / t_xmax, lambda) / t_xmax;
            R = scale * cmul(MmLm1a[ik], MLp[jl]);
        }
        
        // decoupled x < y
        if (t_xmax <= t_ymin)
        {
            private Real scale = pow_int(t_xmax / t_ymax, lambda) / t_ymax;
            R = scale * cmul(MLa[ik], MmLm1p[jl]);
        }
        
        // update all angular blocks
        if (R.x != 0 || R.y != 0)
        for (private int ill  = 0; ill  < ANGULAR_BASIS_SIZE; ill++)
        for (private int illp = 0; illp < ANGULAR_BASIS_SIZE; illp++)
        {
            y[(ill * (2 * ORDER + 1) + i + ORDER - k) * NSPLINE_PROJ + j]
                += ZP * f[lambda] * cmul(R, x[illp * NSPLINE_PROJ + l]);
        }
    }
}

/**
 * @brief Multiply by two-electron integral matrix.
 * 
 * Multiplies the source vector by the modified two-electron matrix R,
 * where only the coupled integrals are considered.
 */
kernel void mmul_2el_coupled
(
    // angular integrals
    constant Real  const * const restrict f,
    // multipole
    private int lambda,
    // radial integrals (CSR storage)
    global Long    const * const restrict Rp,
    global Long    const * const restrict Ri,
    global Complex const * const restrict Rx,
    // source and target vector
    global Complex const * const restrict x,
    global Complex       * const restrict y
)
{
    private int a = get_global_id(0) / NSPLINE_PROJ;
    private int b = get_global_id(0) % NSPLINE_PROJ;
    
    for (private int c = a > ORDER ? a - ORDER : 0; c < NSPLINE_ATOM && c <= a + ORDER; c++)
    for (private int d = b > ORDER ? b - ORDER : 0; d < NSPLINE_PROJ && d <= b + ORDER; d++)
    {
        private int I = min(a,c) * (ORDER + 1) + abs_diff(a,c);
        private int J = min(b,d) * (ORDER + 1) + abs_diff(b,d);
        
        for (private int idx = Rp[I]; idx < Rp[I + 1]; idx++)
        {
            if (Ri[idx] == J)
            {
                y[a * NSPLINE_PROJ + b] += ZP * f[lambda] * cmul(Rx[idx], x[c * NSPLINE_PROJ + d]);
            }
        }
    }
}

/**
 * @brief Multiply by two-electron integral matrix.
 * 
 * Similar to @ref mmul_2el_coupled, only acting on segments.
 */
kernel void mmul_2el_coupled_segment
(
    // angular integrals
    constant Real  const * const restrict f,
    // multipole
    private int lambda,
    // radial integrals (CSR storage)
    global Long    const * const restrict Rp,
    global Long    const * const restrict Ri,
    global Complex const * const restrict Rx,
    // source and target vector
    global Complex const * const restrict x,
    global Complex       * const restrict y,
    // segment
    private int k
)
{
    private int j = get_global_id(0);
    
    if (j < NSPLINE_PROJ)
    for (private int i = k > ORDER ? k - ORDER : 0; i < NSPLINE_ATOM && i <= k + ORDER; i++)
    for (private int l = j > ORDER ? j - ORDER : 0; l < NSPLINE_PROJ && l <= j + ORDER; l++)
    {
        private int ik = min(i,k) * (ORDER + 1) + abs_diff(i,k);
        private int jl = min(j,l) * (ORDER + 1) + abs_diff(j,l);
        
        for (private int idx = Rp[ik]; idx < Rp[ik + 1]; idx++)
        {
            if (Ri[idx] == jl)
            {
                for (private int ill  = 0; ill  < ANGULAR_BASIS_SIZE; ill++)
                for (private int illp = 0; illp < ANGULAR_BASIS_SIZE; illp++)
                {
                    y[(ill * (2 * ORDER + 1) + i + ORDER - k) * NSPLINE_PROJ + j] += ZP *
                        f[(lambda * ANGULAR_BASIS_SIZE + ill) * ANGULAR_BASIS_SIZE + illp] *
                        cmul(Rx[idx], x[illp * NSPLINE_PROJ + l]);
                }
            }
        }
    }
}

/**
 * @brief Matrix-matrix multiplication.
 * 
 * General matrix-matrix multiplication.
 * \f[
 *     C_{ij} = \sum_{k} A_{ik} B_{kj}
 * \f]
 * @param m Row count of A.
 * @param n Column count of B.
 * @param k Column count of A.
 * @param A Input m-by-k matrix in row-major storage.
 * @param B Input k-by-n matrix in column-major storage (i.e., transposed).
 * @param C Output m-by-n matrix in row-major storage.
 */
kernel void mul_ABt
(
    private int m, private int n, private int k,
    global Complex const * const restrict A,
    global Complex const * const restrict B,
    global Complex       * const restrict C
)
{
    // destination block indices
    private int irow_block = get_group_id(0);
    private int icol_block = get_group_id(1);
    
    // global destination index
    private int irow_glob = get_global_id(0);
    private int icol_glob = get_global_id(1);
    
    // group worker threads; indices within the destination block
    private int irow_loc = get_local_id(0);
    private int icol_loc = get_local_id(1);
    
    // work arrays
    local Complex Aloc[BLOCK_SIZE][BLOCK_SIZE];
    local Complex Bloc[BLOCK_SIZE][BLOCK_SIZE];
    
    // aggregated scalar product of the destination element C[idy,idx]
    private Complex res = 0;
    
    // calculate number of blocks
    private int Nblock = ((k + BLOCK_SIZE - 1) / BLOCK_SIZE);
    
    // for all source blocks
    for (private int iblock = 0; iblock < Nblock; iblock++)
    {
        // get indices of the current elements in the matrices
        private int Arow = irow_glob, Acol = iblock * BLOCK_SIZE + icol_loc;
        private int Bcol = icol_glob, Brow = iblock * BLOCK_SIZE + irow_loc;
        
        // load elements of the source blocks into the local memory, pad by zeros
        barrier(CLK_LOCAL_MEM_FENCE);
        Aloc[icol_loc][irow_loc] = (Arow < m && Acol < k ? A[Arow * k + Acol] : (Complex)(0.,0.));
        Bloc[irow_loc][icol_loc] = (Brow < k && Bcol < n ? B[Brow + k * Bcol] : (Complex)(0.,0.));
        barrier(CLK_LOCAL_MEM_FENCE);
        
        // each group's thread will calculate one of the scalar products
        for (private int K = 0; K < BLOCK_SIZE; K++)
            if (iblock * BLOCK_SIZE + K < k)
                res += cmul(Aloc[K][irow_loc],Bloc[K][icol_loc]);
    }
    
    // store result to device memory
    if (irow_glob < m && icol_glob < n)
        C[irow_glob * n + icol_glob] = res;
}

/**
 * @brief Matrix-matrix multiplication.
 * 
 * General matrix-matrix multiplication.
 * \f[
 *     C_{ij} = \sum_{k} A_{ik} B_{kj}
 * \f]
 * @param m Row count of A.
 * @param n Column count of B.
 * @param k Column count of A.
 * @param A Input m-by-k matrix in column-major storage (i.e., transposed).
 * @param B Input k-by-n matrix in column-major storage (i.e., transposed).
 * @param C Output m-by-n matrix in row-major storage.
 */
kernel void mul_AtBt
(
    private int m, private int n, private int k,
    global Complex const * const restrict A,
    global Complex const * const restrict B,
    global Complex       * const restrict C
)
{
    // destination block indices
    private int irow_block = get_group_id(0);
    private int icol_block = get_group_id(1);
    
    // global destination index
    private int irow_glob = get_global_id(0);
    private int icol_glob = get_global_id(1);
    
    // group worker threads; indices within the destination block
    private int irow_loc = get_local_id(0);
    private int icol_loc = get_local_id(1);
    
    // work arrays
    local Complex Aloc[BLOCK_SIZE][BLOCK_SIZE];
    local Complex Bloc[BLOCK_SIZE][BLOCK_SIZE];
    
    // aggregated scalar product of the destination element C[idy,idx]
    private Complex res = 0;
    
    // calculate number of blocks
    private int Nblock = ((k + BLOCK_SIZE - 1) / BLOCK_SIZE);
    
    // for all source blocks
    for (private int iblock = 0; iblock < Nblock; iblock++)
    {
        // get indices of the current elements in the matrices
        private int Arow = irow_glob, Acol = iblock * BLOCK_SIZE + icol_loc;
        private int Bcol = icol_glob, Brow = iblock * BLOCK_SIZE + irow_loc;
        
        // load elements of the source blocks into the local memory, pad by zeros
        barrier(CLK_LOCAL_MEM_FENCE);
        Aloc[icol_loc][irow_loc] = (Arow < m && Acol < k ? A[Arow + k * Acol] : (Complex)(0.,0.));
        Bloc[irow_loc][icol_loc] = (Brow < k && Bcol < n ? B[Brow + k * Bcol] : (Complex)(0.,0.));
        barrier(CLK_LOCAL_MEM_FENCE);
        
        // each group's thread will calculate one of the scalar products
        for (private int K = 0; K < BLOCK_SIZE; K++)
            if (iblock * BLOCK_SIZE + K < k)
                res += cmul(Aloc[K][irow_loc],Bloc[K][icol_loc]);
    }
    
    // store result to device memory
    if (irow_glob < m && icol_glob < n)
        C[irow_glob * n + icol_glob] = res;
}
