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

// Necessary compile-time definitions:
// -D ORDER=... (implies local size in "mmul_2el")
// -D NSPLINE_ATOM=... (needed by "mmul_1el", "mmul_2el", "mul_ABt" and "kron_div")
// -D NSPLINE_PROJ=... (needed by "mmul_1el", "mmul_2el", "mul_ABt" and "kron_div")
// -D NREKNOT_ATOM=... (needed by "mmul_2el")
// -D NREKNOT_PROJ=... (needed by "mmul_2el")
// -D NLOCAL=... (local size in "scalar_product", "norm" and "mmul_1el")
// -D NBLOCK_SIZE=... (block size and local size in "mul_ABt")

// Enable double precision (redundant in OpenCL 2.0).
#pragma OPENCL EXTENSION cl_khr_fp64: enable

// Derived variables.
#define NROW            (NSPLINE_ATOM * NSPLINE_PROJ)
#define BLOCK_VOLUME    (BLOCK_SIZE * BLOCK_SIZE)
#define NUM_BLOCKS_X    ((NSPLINE_ATOM + BLOCK_SIZE - 1) / BLOCK_SIZE)
#define NUM_BLOCKS_Y    ((NSPLINE_PROJ + BLOCK_SIZE - 1) / BLOCK_SIZE)

/**
 * @brief Complex multiplication.
 * 
 * Multiplies two complex numbers and returns the product.
 */
double2 cmul (double2 a, double2 b)
{
    double2 c;
    c.x = a.x * b.x - a.y * b.y;
    c.y = a.x * b.y + a.y * b.x;
    return c;
}

/**
 * @brief Complex division.
 * 
 * Divides two complex numbers and returns the fraction. No overflow
 * checking is done.
 */
double2 cdiv (double2 a, double2 b)
{
    double2 c;
    double b2 = b.x * b.x + b.y * b.y;
    
    c.x = (a.x * b.x + a.y * b.y) / b2;
    c.y = (a.y * b.x - a.x * b.y) / b2;
    
    return c;
}

/**
 * @brief Integer power.
 * 
 * Fast integer power of arbitrary numerical type. Uses only
 * multiplications. The number of multiplications is proportional
 * to log(n).
 */
double pow_int (double x, unsigned n)
{
    double value = 1;
    
    do
    {
        if(n % 2 == 1)
            value *= x;
        
        n /= 2;
        x *= x;
    }
    while (n);
    
    return value;
}

/**
 * @brief AXBY operation.
 * 
 * Computes the linear combination @f$ ax + by @f$ of vectors @f$ x @f$ and 
 * @f$ y @f$ and stores the result in @f$ x @f$.
 * @param a Complex factor.
 * @param x Source and destination vector.
 * @param b Complex factor.
 * @param y Source vector.
 */
kernel void a_vec_b_vec (private double2 a, global double2 *x, private double2 b, global double2 *y)
{
    uint i = get_global_id(0);
    
    if (i < NROW)
        x[i] = cmul(a,x[i]) + cmul(b,y[i]);
}

/**
 * @brief Full scalar product.
 * 
 * Calculates element-wise product of the arrays and reduces them using local memory
 * to one number per work group. These intermediate numbers are the summed by CPU.
 * @param u Source vector.
 * @param v Source vector.
 * @param z Output vector for intermediate segment scalar products.
 */
kernel void scalar_product (global double2 *u, global double2 *v, global double2 *z)
{
    // position of this thread among other threads
    private int iglobal = get_global_id(0);
    private int ilocal = get_local_id(0);
    
    // calculate product of array elements
    local double2 uv[NLOCAL];
    uv[ilocal] = (iglobal < NROW ? cmul(u[iglobal],v[iglobal]) : (double2)(0.,0.));
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
 * as one number per work group. These intermediate numbers are the summed by CPU.
 * @param v Source vector.
 * @param z Output vector for intermediate segment norms.
 */
kernel void norm (global double2 *v, global double *z)
{
    // position of this thread among other threads
    private int iglobal = get_global_id(0);
    private int ilocal = get_local_id(0);
    
    // get element to process by this thread
    private double2 vi = v[iglobal];
    
    // calculate squared modulus of an array element
    local double vv[NLOCAL];
    vv[ilocal] = (iglobal < NROW ? vi.x * vi.x + vi.y * vi.y : 0.);
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
 * @param E Total energy of the system.
 * @param Sp Row-padded upper overlap matrix.
 * @param Dp Row-padded upper derivative overlap matrix.
 * @param M1p Row-padded upper integral moment matrix (for r^1).
 * @param M2p Row-padded upper integral moment matrix (for r^2).
 * @param l1 Angular momentum of first electron.
 * @param l2 Angular momentum of second electron.
 * @param x Source vector.
 * @param y Destination vector.
 */
kernel void mmul_1el
(
    // energy
    private double E,
    // row-padded one-electron matrices
    global double2 const * const restrict Sp,
    global double2 const * const restrict Dp,
    global double2 const * const restrict M1p,
    global double2 const * const restrict M2p,
    // angular momenta and nonzero multipoles
    private int l1,
    private int l2,
    // source and target vector
    global double2 const * const restrict x,
    global double2       * const restrict y
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
        private double2 elem = E * cmul(Sp[ik],Sp[jl]);
        elem -= 0.5 * (cmul(Dp[ik],Sp[jl]) + cmul(Sp[ik],Dp[jl]));
        elem -= 0.5 * l1 * (l1 + 1.) * cmul(M2p[ik],Sp[jl]) + 0.5 * l2 * (l2 + 1.) * cmul(Sp[ik],M2p[jl]);
        elem += cmul(M1p[ik],Sp[jl]) + cmul(Sp[ik],M1p[jl]);
        
        // multiply right-hand side by that matrix element
        y[i * NSPLINE_PROJ + j] += cmul(elem, x[k * NSPLINE_PROJ + l]);
    }
}

/**
 * @brief Multiplication by two-electron Hamiltonian matrix.
 * 
 * Multiplies vector by a two-electron Hamiltonian matrix @f$ R_{ijkl}^\lambda @f$.
 * Each multi-index @f$ (i,j) @f$ is assigned to a different thread.
 * This kernel is called for every multipole independently.
 * 
 * @param f Real prefactor (angular integral).
 * @param MiL Partial integral moments (r^lambda).
 * @param MimLm1 Partial integral moments (r^(-lambda-1)).
 * @param x Source vector.
 * @param y Destination vector.
 */
kernel void mmul_2el
(
    // B-spline knots
    constant double2 const * const restrict t,
    // multipole
    private int lambda,
    // angular integral
    private double f,
    // one-electron dull moments
    global double2 const * const restrict ML,
    global double2 const * const restrict MmLm1,
    // one-electron partial moments
    global double2 const * const restrict MiL,
    global double2 const * const restrict MimLm1,
    // source and target vector
    global double2 const * const restrict x,
    global double2       * const restrict y
)
{
    // output vector element index
    private int i = get_global_id(0) / NSPLINE_PROJ;
    private int j = get_global_id(0) % NSPLINE_PROJ;
    
    // auxiliary variables
    private double2 m_ik, m_jl;
    private double scale, tx, ty;
    
    // pointers to the needed partial integral moments
    global double2 const * restrict M_ik;
    global double2 const * restrict M_jl;
    
    // for all source vector elements
    if (i < NSPLINE_ATOM)
    for (private int k = i - ORDER; k <= i + ORDER; k++) if (0 <= k && k < NSPLINE_ATOM)
    for (private int l = j - ORDER; l <= j + ORDER; l++) if (0 <= l && l < NSPLINE_PROJ)
    {
        // matrix element
        private double2 elem = 0;
        
        // Are the integral moments completely decoupled, i.e. there is there no overlap between Bi, Bj, Bk and Bl?
        // In such cases we can compute the off-diagonal contribution just as a product of the two
        // (scaled) integral moments of order "lambda" and "-lambda-1", respectively.
        
        tx = t[min(i,k) + ORDER + 1].x;
        ty = t[min(j,l) + ORDER + 1].x;
        
        // (j,l) << (i,k)
        if (min(j,l) + ORDER < max(i,k))
        {
            scale = pow_int(ty / tx, lambda) / tx;
            elem += scale * MmLm1[min(i,k) * (ORDER + 1) + abs(i-k)] * ML[min(j,l) * (ORDER + 1) + abs(j-l)];
        }
        
        // (i,k) << (j,l)
        else if (min(i,k) + ORDER < max(j,l))
        {
            scale = pow_int(tx / ty, lambda) / ty;
            elem += scale * ML[min(i,k) * (ORDER + 1) + abs(i-k)] * MmLm1[min(j,l) * (ORDER + 1) + abs(j-l)];
        }
        
        // Further parts are a bit cryptical, because we are using precomputed
        // (partial, per knot) integral moments, which are quite compactly stored
        // in arrays M_L and M_mLm1 of shape [Nspline * (2*order+1) * (order+1)],
        // but the aim is straightforward: Just to sum the offdiagonal elements,
        // i.e. the products of two two-spline integrals, when ix != iy.
        
        // (i,k) ~ (j,l)
        else
        {
            // ix < iy
            M_ik = MiL    + (i * (2*ORDER+1) + k - (i-ORDER)) * (ORDER+1);
            M_jl = MimLm1 + (j * (2*ORDER+1) + l - (j-ORDER)) * (ORDER+1);
            for (private int ix = i; ix < min(i + ORDER + 1, NREKNOT_ATOM - 1); ix++) if (t[ix + 1].x > 0)
            {
                m_ik = M_ik[ix - i]; tx = t[ix + 1].x;
                
                for (private int iy = max(j, ix + 1); iy < min(j + ORDER + 1, NREKNOT_PROJ - 1); iy++)
                {
                    m_jl = M_jl[iy - j]; ty = t[iy + 1].x;
                    scale = pow_int(tx/ty,lambda)/ty;
                    elem += cmul(m_ik,m_jl) * scale;
                }
            }
            
            // ix > iy (by renaming the ix,iy indices)
            M_ik = MimLm1 + (i * (2*ORDER+1) + k - (i-ORDER)) * (ORDER+1);
            M_jl = MiL    + (j * (2*ORDER+1) + l - (j-ORDER)) * (ORDER+1);
            for (private int ix = j; ix < min(j + ORDER + 1, NREKNOT_PROJ - 1); ix++) if (t[ix + 1].x > 0)
            {
                m_jl = M_jl[ix - j]; tx = t[ix + 1].x;
                
                for (private int iy = max(i, ix + 1); iy < min(i + ORDER + 1, NREKNOT_ATOM - 1); iy++)
                {
                    m_ik = M_ik[iy - i]; ty = t[iy + 1].x;
                    scale = pow_int(tx/ty,lambda)/ty;
                    elem += cmul(m_ik,m_jl) * scale;
                }
            }
        }
        
        // multiply right-hand side by that matrix element
        y[i * NSPLINE_PROJ + j] -= f * cmul(elem, x[k * NSPLINE_PROJ + l]);
    }
}

/**
 * @brief Matrix-matrix multiplication.
 * 
 * General matrix-matrix multiplication.
 * @param A Input matrix (row-major storage).
 * @param B Input matrix (column-major storage).
 * @param C Outpu matrix (row-major storage).
 */
kernel void mul_ABt
(
    global double2 const * const restrict A, // input matrix A (row-major)
    global double2 const * const restrict B, // input matrix B (col-major)
    global double2       * const restrict C  // output matrix C (row-major)
)
{
    // work arrays
    local double2 Aloc[BLOCK_SIZE][BLOCK_SIZE];
    local double2 Bloc[BLOCK_SIZE][BLOCK_SIZE];
    
    // destination blocks
    private int idyblock = get_group_id(0); // row
    private int idxblock = get_group_id(1); // col
    
    // group worker threads
    private int iylocal = get_local_id(0); // row
    private int ixlocal = get_local_id(1); // col
    
    // aggregated scalar product of the destination element C[get_global_id(0),get_global_id(1)]
    private double2 res = 0;
    
    // for all source blocks
    for (private int iblock = 0; iblock < NUM_BLOCKS_X; iblock++)
    {
        // load source blocks into the local memory (WARNING: Should be padded by zeros: will segfault on CPU.)
        barrier(CLK_LOCAL_MEM_FENCE);
        Aloc[ixlocal][iylocal] = A[(idyblock * BLOCK_SIZE + iylocal) * NSPLINE_PROJ + (iblock * BLOCK_SIZE + ixlocal)];
        Bloc[iylocal][ixlocal] = B[(idxblock * BLOCK_SIZE + ixlocal) * NSPLINE_PROJ + (iblock * BLOCK_SIZE + iylocal)];
        barrier(CLK_LOCAL_MEM_FENCE);
        
        // each group's thread will calculate one of BLOCK_VOLUME scalar products
        for (private int k = 0; k < BLOCK_SIZE; k++)
            if (iblock * BLOCK_SIZE + k < NSPLINE_ATOM)
                res += cmul(Aloc[k][iylocal],Bloc[k][ixlocal]);
    }
    
    // store result to device memory
    if (get_global_id(0) < NSPLINE_ATOM && get_global_id(1) < NSPLINE_PROJ)
        C[get_global_id(0) * NSPLINE_PROJ + get_global_id(1)] = res;
}

/**
 * @brief KPA preconditioner.
 * 
 * Divides the source vector by the KPA preconditioner term.
 * @param E Total energy.
 * @param D1 Eigenvalues for first angular momentum.
 * @param D2 Eigenvalues for second angular momentum.
 * @param y Source/destination vector.
 */
kernel void kron_div
(
    private double2 E,
    global double2 const * const restrict D1,
    global double2 const * const restrict D2,
    global double2       * const restrict y
)
{
    // get worker's segment index (0 <= i < NSPLINE_ATOM)
    private int i = get_global_id(0);
    
    // y = y / (E (I kron I) - (D1 kron I) - (I kron D2))
    for (private int j = 0; j < NSPLINE_PROJ; j++)
        y[j * NSPLINE_ATOM + i] = cdiv(y[j * NSPLINE_ATOM + i], E - D1[j] - D2[i]);
}
