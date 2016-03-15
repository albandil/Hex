//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2016, Jakub Benda, Charles University in Prague                    //
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
// -D ANGULAR_BASIS_SIZE=... (number of coupled angular states, used in offset routines)
// -D PROJECTILE_BASIS_OFFSET=... (skipped projectile basis splines, used in "mmul_2el")
// -D NSRCSEG=...
// -D NDSTSEG=...

// Enable double precision (redundant in OpenCL 2.0).
#pragma OPENCL EXTENSION cl_khr_fp64: enable

// Derived variables.
#define NROW (NSPLINE_ATOM * NSPLINE_PROJ)

/**
 * @brief Complex multiplication.
 * 
 * Multiplies two complex numbers and returns the product.
 */
double2 cmul (private double2 a, private double2 b)
{
    private double2 c;
    
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
double2 cdiv (private double2 a, private double2 b)
{
    private double2 c;
    private double b2 = b.x * b.x + b.y * b.y;
    
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
double pow_int (private double x, private int n)
{
    private double value = 1; // = x^0
    
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
 * @param a Complex factor.
 * @param x Source and destination vector.
 * @param b Complex factor.
 * @param y Source vector.
 */
kernel void a_vec_b_vec (private double2 a, global double2 *x, private double2 b, global double2 *y)
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
 * as one number per work group. These intermediate numbers are then summed by CPU.
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
 * \f[
 *     y_{ij} = \sum_{kl} \left(
 *         E S_{ik} S_{jl} - H_{ik}^{(1)} S_{jl} - S_{ik} H_{jl}^{(2)}
 *     \right) x_{kl}
 * \f]
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
    // row-padded one-electron matrices (atom)
    global double2 const * const restrict Spa,
    global double2 const * const restrict Dpa,
    global double2 const * const restrict M1pa,
    global double2 const * const restrict M2pa,
    // row-padded one-electron matrices (projectile)
    global double2 const * const restrict Spp,
    global double2 const * const restrict Dpp,
    global double2 const * const restrict M1pp,
    global double2 const * const restrict M2pp,
    // angular momenta
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
        private double2 elem = E * cmul(Spa[ik],Spp[jl]);
        elem -= 0.5 * (cmul(Dpa[ik],Spp[jl]) + cmul(Spa[ik],Dpp[jl]));
        elem -= 0.5 * l1 * (l1 + 1.) * cmul(M2pa[ik],Spp[jl]) + 0.5 * l2 * (l2 + 1.) * cmul(Spa[ik],M2pp[jl]);
        elem += cmul(M1pa[ik],Spp[jl]) + cmul(Spa[ik],M1pp[jl]);
        
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
 * \f[
 *     y_{ij} = y_{ij} - \sum_{kl} f^{\lambda} R_{ijkl}^{\lambda} x_{kl}
 * \f]
 * @param f Real prefactor (angular integral).
 * @param MiL Partial integral moments (r^lambda).
 * @param MimLm1 Partial integral moments (r^(-lambda-1)).
 * @param x Source vector.
 * @param y Destination vector.
 */
kernel void mmul_2el
(
    // B-spline knots (atom and projectile)
    constant double2 const * const restrict ta,
    constant double2 const * const restrict tp,
    // multipole
    private int lambda,
    // angular integral
    private double f,
    // one-electron full and partial moments (atom)
    global double2 const * const restrict MLa,
    global double2 const * const restrict MmLm1a,
    global double2 const * const restrict MiLa,
    global double2 const * const restrict MimLm1a,
    // one-electron full and partial moments (projectile)
    global double2 const * const restrict MLp,
    global double2 const * const restrict MmLm1p,
    global double2 const * const restrict MiLp,
    global double2 const * const restrict MimLm1p,
    // two-electron diagonal contributions
    global double2 const * const restrict Rdia,
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
    for (private int k = max(0, i - ORDER); k <= min(i + ORDER, NSPLINE_ATOM - 1); k++)
    for (private int l = max(0, j - ORDER); l <= min(j + ORDER, NSPLINE_PROJ - 1); l++)
    {
        // matrix element
        private double2 elem = 0;
        
        // boundary knots of the participating B-splines
        private double ti1 = ta[i].x, ti2 = ta[i + ORDER + 1].x;
        private double tj1 = tp[j].x, tj2 = tp[j + ORDER + 1].x;
        private double tk1 = ta[k].x, tk2 = ta[k + ORDER + 1].x;
        private double tl1 = tp[l].x, tl2 = tp[l + ORDER + 1].x;
        
        // Are the integral moments completely decoupled, i.e. is there no overlap between Ba, Bb, Bc and Bd?
        // In such cases we can compute the off-diagonal contribution just as a product of the two
        // integral moments of order "lambda" and "-lambda-1", respectively. Moreover, in such
        // case there is no diagonal contribution, because there is no overlap between the _four_
        // participating B-splines.
        
        // (j,l) << (i,k)
        if (min(tj2,tl2) <= max(ti1,tk1))
        {
            tx = min(ti2,tk2);
            ty = min(tj2,tl2);
            scale = pow_int(ty / tx, lambda) / tx;
            elem += scale * cmul(MmLm1a[min(i,k) * (ORDER + 1) + abs(i-k)], MLp[min(j,l) * (ORDER + 1) + abs(j-l)]);
        }
        
        // (i,k) << (j,l)
        else if (min(ti2,tk2) <= max(tj1,tl1))
        {
            tx = min(ti2,tk2);
            ty = min(tj2,tl2);
            scale = pow_int(tx / ty, lambda) / ty;
            elem += scale * cmul(MLa[min(i,k) * (ORDER + 1) + abs(i-k)], MmLm1p[min(j,l) * (ORDER + 1) + abs(j-l)]);
        }
        
        // The rest allows overlap of the four B-splines and is used only by the first
        // (origin) panel.
        
        // Further parts are a bit cryptical, because we are using precomputed
        // (partial, per knot) integral moments, which are quite compactly stored
        // in arrays M_L and M_mLm1 of shape [Nspline * (2*order+1) * (order+1)],
        // but the aim is straightforward: Just to sum the offdiagonal elements,
        // i.e. the products of two two-spline integrals, when ix != iy.
        
        // (i,k) ~ (j,l)
        else
        {
            // retrieve diagonal contribution
            if (i <= j && i <= k && i <= l)
                elem += Rdia[((i * (ORDER+1) + (j-i)) * (ORDER+1) + (k-i)) * (ORDER+1) + (l-i)];
            else if (j <= i && j <= k && j <= l)
                elem += Rdia[((j * (ORDER+1) + (i-j)) * (ORDER+1) + (l-j)) * (ORDER+1) + (k-j)];
            else if (k <= i && k <= j && k <= l)
                elem += Rdia[((k * (ORDER+1) + (j-k)) * (ORDER+1) + (i-k)) * (ORDER+1) + (l-k)];
            else // (l <= i && l <= j && l <= k)
                elem += Rdia[((l * (ORDER+1) + (i-l)) * (ORDER+1) + (j-l)) * (ORDER+1) + (k-l)];
            
            // ix < iy
            M_ik = MiLa    + (i * (2*ORDER+1) + k - (i-ORDER)) * (ORDER+1);
            M_jl = MimLm1p + (j * (2*ORDER+1) + l - (j-ORDER)) * (ORDER+1);
            for (private int ix = i; ix < min(i + ORDER + 1, NREKNOT_ATOM - 1); ix++) if (ta[ix + 1].x > 0)
            {
                m_ik = M_ik[ix - i]; tx = ta[ix + 1].x;
                
                for (private int iy = max(j, ix + 1); iy < min(j + ORDER + 1, NREKNOT_PROJ - 1); iy++)
                {
                    m_jl = M_jl[iy - j]; ty = tp[iy + 1].x;
                    scale = pow_int(tx/ty,lambda)/ty;
                    elem += cmul(m_ik,m_jl) * scale;
                }
            }
            
            // ix > iy (by renaming the ix,iy indices)
            M_ik = MimLm1a + (i * (2*ORDER+1) + k - (i-ORDER)) * (ORDER+1);
            M_jl = MiLp    + (j * (2*ORDER+1) + l - (j-ORDER)) * (ORDER+1);
            for (private int ix = j; ix < min(j + ORDER + 1, NREKNOT_PROJ - 1); ix++) if (tp[ix + 1].x > 0)
            {
                m_jl = M_jl[ix - j]; tx = tp[ix + 1].x;
                
                for (private int iy = max(i, ix + 1); iy < min(i + ORDER + 1, NREKNOT_ATOM - 1); iy++)
                {
                    m_ik = M_ik[iy - i]; ty = ta[iy + 1].x;
                    scale = pow_int(tx/ty,lambda)/ty;
                    elem += cmul(m_ik,m_jl) * scale;
                }
            }
        }
        
        // multiply right-hand side by that matrix element
        y[i * (ulong)(NSPLINE_PROJ) + j] -= f * cmul(elem, x[k * (ulong)(NSPLINE_PROJ) + l]);
    }
}

kernel void mmul_2el_decoupled
(
    // B-spline knots (atom and projectile)
    constant double2 const * const restrict ta,
    constant double2 const * const restrict tp,
    // multipole
    private int lambda,
    // angular integral
    private double f,
    // one-electron moments
    global double2 const * const restrict MLa,
    global double2 const * const restrict MmLm1a,
    global double2 const * const restrict MLp,
    global double2 const * const restrict MmLm1p,
    // source and target vector
    global double2 const * const restrict x,
    global double2       * const restrict y
)
{
    // output vector element index
    private int i = get_global_id(0) / NSPLINE_PROJ;
    private int j = get_global_id(0) % NSPLINE_PROJ;
    
    // auxiliary variables
    private double scale, tx, ty;
    private double2 elem;
    
    // for all source vector elements
    if (i < NSPLINE_ATOM)
    for (private int k = max(0, i - ORDER); k <= min(i + ORDER, NSPLINE_ATOM - 1); k++)
    for (private int l = max(0, j - ORDER); l <= min(j + ORDER, NSPLINE_PROJ - 1); l++)
    {
        // matrix element
        elem = 0;
        
        // (j,l) << (i,k)
        if (min(j,l) + ORDER < max(i,k))
        {
            tx = min(ta[i + ORDER + 1].x,ta[k + ORDER + 1].x);
            ty = min(tp[j + ORDER + 1].x,tp[l + ORDER + 1].x);
            scale = pow_int(ty / tx, lambda) / tx;
            elem += scale * cmul(MmLm1a[min(i,k) * (ORDER + 1) + abs(i-k)], MLp[min(j,l) * (ORDER + 1) + abs(j-l)]);
        }
        
        // (i,k) << (j,l)
        if (min(i,k) + ORDER < max(j,l))
        {
            tx = min(ta[i + ORDER + 1].x,ta[k + ORDER + 1].x);
            ty = min(tp[j + ORDER + 1].x,tp[l + ORDER + 1].x);
            scale = pow_int(tx / ty, lambda) / ty;
            elem += scale * cmul(MLa[min(i,k) * (ORDER + 1) + abs(i-k)], MmLm1p[min(j,l) * (ORDER + 1) + abs(j-l)]);
        }
        
        // multiply right-hand side by the matrix element
        y[i * (ulong)(NSPLINE_PROJ) + j] -= f * cmul(elem, x[k * (ulong)(NSPLINE_PROJ) + l]);
    }
}

kernel void mmul_2el_coupled
(
    // B-spline knots (atom and projectile)
    constant double2 const * const restrict ta,
    constant double2 const * const restrict tp,
    // multipole
    private int lambda,
    // angular integral
    private double f,
    // one-electron partial moments (atom)
    global double2 const * const restrict MiLa,
    global double2 const * const restrict MimLm1a,
    // one-electron partial moments (projectile)
    global double2 const * const restrict MiLp,
    global double2 const * const restrict MimLm1p,
    // two-electron diagonal contributions
    global double2 const * const restrict Rdia,
    // source and target vector
    global double2 const * const restrict x,
    global double2       * const restrict y
)
{
    // output vector element index
    private int i = get_global_id(0) / (2 * ORDER + 1);
    private int j = get_global_id(0) % (2 * ORDER + 1) + i - ORDER;
    
    // auxiliary variables
    private double2 elem, m_ik, m_jl;
    private double scale, tx, ty;
    
    // pointers to the needed partial integral moments
    global double2 const * restrict M_ik;
    global double2 const * restrict M_jl;
    
    // for all source vector elements
    if (i < NSPLINE_ATOM && 0 <= j && j <= NSPLINE_PROJ)
    for (private int k = max(0, i - ORDER); k <= min(i + ORDER, NSPLINE_ATOM - 1); k++)
    for (private int l = max(0, j - ORDER); l <= min(j + ORDER, NSPLINE_PROJ - 1); l++)
    if (max(i,k) <= min(j,l) + ORDER && max(j,l) <= min(i,k) + ORDER)
    {
        // matrix element
        elem = 0;
        
        // (i,k) ~ (j,l)
        
        // retrieve diagonal contribution
        if (i <= j && i <= k && i <= l)
            elem += Rdia[((i * (ORDER+1) + (j-i)) * (ORDER+1) + (k-i)) * (ORDER+1) + (l-i)];
        else if (j <= i && j <= k && j <= l)
            elem += Rdia[((j * (ORDER+1) + (i-j)) * (ORDER+1) + (l-j)) * (ORDER+1) + (k-j)];
        else if (k <= i && k <= j && k <= l)
            elem += Rdia[((k * (ORDER+1) + (j-k)) * (ORDER+1) + (i-k)) * (ORDER+1) + (l-k)];
        else // (l <= i && l <= j && l <= k)
            elem += Rdia[((l * (ORDER+1) + (i-l)) * (ORDER+1) + (j-l)) * (ORDER+1) + (k-l)];
        
        // ix < iy
        M_ik = MiLa    + (i * (2*ORDER+1) + k - (i-ORDER)) * (ORDER+1);
        M_jl = MimLm1p + (j * (2*ORDER+1) + l - (j-ORDER)) * (ORDER+1);
        for (private int ix = i; ix < min(i + ORDER + 1, NREKNOT_ATOM - 1); ix++) if (ta[ix + 1].x > 0)
        {
            m_ik = M_ik[ix - i]; tx = ta[ix + 1].x;
            
            for (private int iy = max(j, ix + 1); iy < min(j + ORDER + 1, NREKNOT_PROJ - 1); iy++)
            {
                m_jl = M_jl[iy - j]; ty = tp[iy + 1].x;
                scale = pow_int(tx/ty,lambda)/ty;
                elem += cmul(m_ik,m_jl) * scale;
            }
        }
        
        // ix > iy (by renaming the ix,iy indices)
        M_ik = MimLm1a + (i * (2*ORDER+1) + k - (i-ORDER)) * (ORDER+1);
        M_jl = MiLp    + (j * (2*ORDER+1) + l - (j-ORDER)) * (ORDER+1);
        for (private int ix = j; ix < min(j + ORDER + 1, NREKNOT_PROJ - 1); ix++) if (tp[ix + 1].x > 0)
        {
            m_jl = M_jl[ix - j]; tx = tp[ix + 1].x;
            
            for (private int iy = max(i, ix + 1); iy < min(i + ORDER + 1, NREKNOT_ATOM - 1); iy++)
            {
                m_ik = M_ik[iy - i]; ty = ta[iy + 1].x;
                scale = pow_int(tx/ty,lambda)/ty;
                elem += cmul(m_ik,m_jl) * scale;
            }
        }
        
        // multiply right-hand side by that matrix element
        y[i * (ulong)(NSPLINE_PROJ) + j] -= f * cmul(elem, x[k * (ulong)(NSPLINE_PROJ) + l]);
    }
}

/**
 * @brief Multiplication by two-electron Hamiltonian matrix.
 * 
 * Multiplies vector by a two-electron Hamiltonian matrix @f$ R_{ijkl}^\lambda @f$.
 * Each multi-index @f$ (i,j) @f$ is assigned to a different thread.
 * This kernel is called for every multipole independently.
 * \f[
 *     y_{ij} = y_{ij} - \sum_{kl} f^{\lambda} R_{ijkl}^{\lambda} x_{kl}
 * \f]
 * @param f Real prefactor (angular integral).
 * @param MiL Partial integral moments (r^lambda).
 * @param MimLm1 Partial integral moments (r^(-lambda-1)).
 * @param x Source vector.
 * @param y Destination vector.
 */
kernel void mmul_2el_offset
(
    // B-spline knots (atom and projectile)
    constant double2 const * const restrict ta,
    constant double2 const * const restrict tp,
    // multipole
    private int lambda,
    // angular integral
    global double  const * const restrict f, private int foffset,
    // one-electron full and partial moments (atom)
    global double2 const * const restrict MLa,
    global double2 const * const restrict MmLm1a,
    global double2 const * const restrict MiLa,
    global double2 const * const restrict MimLm1a,
    // one-electron full and partial moments (projectile)
    global double2 const * const restrict MLp,
    global double2 const * const restrict MmLm1p,
    global double2 const * const restrict MiLp,
    global double2 const * const restrict MimLm1p,
    // two-electron diagonal contributions
    global double2 const * const restrict Rdia,
    // source and target vector
    global double2 const * const restrict x,
    global double2       * const restrict y,
    // angular block domain
    private short x_ang_begin,
    private short y_ang_begin
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
        
        // boundary knots of the participating B-splines
        private double ti1 = ta[i].x, ti2 = ta[i + ORDER + 1].x;
        private double tj1 = tp[j].x, tj2 = tp[j + ORDER + 1].x;
        private double tk1 = ta[k].x, tk2 = ta[k + ORDER + 1].x;
        private double tl1 = tp[l].x, tl2 = tp[l + ORDER + 1].x;
        
        // Are the integral moments completely decoupled, i.e. there is there no overlap between Ba, Bb, Bc and Bd?
        // In such cases we can compute the off-diagonal contribution just as a product of the two
        // integral moments of order "lambda" and "-lambda-1", respectively. Moreover, in such
        // case there is no diagonal contribution, because there is no overlap between the _four_
        // participating B-splines.
        
        // (j,l) << (i,k)
        if (min(tj2,tl2) <= max(ti1,tk1))
        {
            tx = min(ti2,tk2);
            ty = min(tj2,tl2);
            scale = pow_int(ty / tx, lambda) / tx;
            elem += scale * cmul(MmLm1a[min(i,k) * (ORDER + 1) + abs(i-k)], MLp[min(j,l) * (ORDER + 1) + abs(j-l)]);
        }
        
        // (i,k) << (j,l)
        else if (min(ti2,tk2) <= max(tj1,tl1))
        {
            tx = min(ti2,tk2);
            ty = min(tj2,tl2);
            scale = pow_int(tx / ty, lambda) / ty;
            elem += scale * cmul(MLa[min(i,k) * (ORDER + 1) + abs(i-k)], MmLm1p[min(j,l) * (ORDER + 1) + abs(j-l)]);
        }
        
        // The rest allows overlap of the four B-splines and is used only by the first
        // (origin) panel.
        
        // Further parts are a bit cryptical, because we are using precomputed
        // (partial, per knot) integral moments, which are quite compactly stored
        // in arrays M_L and M_mLm1 of shape [Nspline * (2*order+1) * (order+1)],
        // but the aim is straightforward: Just to sum the offdiagonal elements,
        // i.e. the products of two two-spline integrals, when ix != iy.
        
        // (i,k) ~ (j,l)
        else
        {
            // retrieve diagonal contribution
            if (i <= j && i <= k && i <= l)
                elem += Rdia[((i * (ORDER+1) + (j-i)) * (ORDER+1) + (k-i)) * (ORDER+1) + (l-i)];
            else if (j <= i && j <= k && j <= l)
                elem += Rdia[((j * (ORDER+1) + (i-j)) * (ORDER+1) + (l-j)) * (ORDER+1) + (k-j)];
            else if (k <= i && k <= j && k <= l)
                elem += Rdia[((k * (ORDER+1) + (j-k)) * (ORDER+1) + (i-k)) * (ORDER+1) + (l-k)];
            else // (l <= i && l <= j && l <= k)
                elem += Rdia[((l * (ORDER+1) + (i-l)) * (ORDER+1) + (j-l)) * (ORDER+1) + (k-l)];
            
            // ix < iy
            M_ik = MiLa    + (i * (2*ORDER+1) + k - (i-ORDER)) * (ORDER+1);
            M_jl = MimLm1p + (j * (2*ORDER+1) + l - (j-ORDER)) * (ORDER+1);
            for (private int ix = i; ix < min(i + ORDER + 1, NREKNOT_ATOM - 1); ix++) if (ta[ix + 1].x > 0)
            {
                m_ik = M_ik[ix - i]; tx = ta[ix + 1].x;
                
                for (private int iy = max(j, ix + 1); iy < min(j + ORDER + 1, NREKNOT_PROJ - 1); iy++)
                {
                    m_jl = M_jl[iy - j]; ty = tp[iy + 1].x;
                    scale = pow_int(tx/ty,lambda)/ty;
                    elem += cmul(m_ik,m_jl) * scale;
                }
            }
            
            // ix > iy (by renaming the ix,iy indices)
            M_ik = MimLm1a + (i * (2*ORDER+1) + k - (i-ORDER)) * (ORDER+1);
            M_jl = MiLp    + (j * (2*ORDER+1) + l - (j-ORDER)) * (ORDER+1);
            for (private int ix = j; ix < min(j + ORDER + 1, NREKNOT_PROJ - 1); ix++) if (tp[ix + 1].x > 0)
            {
                m_jl = M_jl[ix - j]; tx = tp[ix + 1].x;
                
                for (private int iy = max(i, ix + 1); iy < min(i + ORDER + 1, NREKNOT_ATOM - 1); iy++)
                {
                    m_ik = M_ik[iy - i]; ty = ta[iy + 1].x;
                    scale = pow_int(tx/ty,lambda)/ty;
                    elem += cmul(m_ik,m_jl) * scale;
                }
            }
        }
        
        // all-to-all-segments multiplication
        for (int ill  = y_ang_begin; ill  < min(y_ang_begin + NDSTSEG, ANGULAR_BASIS_SIZE); ill ++)
        for (int illp = x_ang_begin; illp < min(x_ang_begin + NSRCSEG, ANGULAR_BASIS_SIZE); illp++)
        {
            y[((ill - y_ang_begin) * (ulong)(NSPLINE_ATOM) + i) * NSPLINE_PROJ + j] -=
                f[foffset + ill * ANGULAR_BASIS_SIZE + illp]
                * cmul(elem, x[((illp - x_ang_begin) * (ulong)(NSPLINE_ATOM) + k) * NSPLINE_PROJ + l]);
        }
    }
}

kernel void mmul_2el_decoupled_offset
(
    // B-spline knots (atom and projectile)
    constant double2 const * const restrict ta,
    constant double2 const * const restrict tp,
    // multipole
    private int lambda,
    // angular integral
    global double  const * const restrict f, private int foffset,
    // one-electron moments (atom)
    global double2 const * const restrict MLa,
    global double2 const * const restrict MmLm1a,
    // one-electron moments (projectile)
    global double2 const * const restrict MLp,
    global double2 const * const restrict MmLm1p,
    // source and target vector
    global double2 const * const restrict x,
    global double2       * const restrict y,
    // angular block domain
    private short x_ang_begin,
    private short y_ang_begin
)
{
    // output vector element index
    private int i = get_global_id(0) / NSPLINE_PROJ;
    private int j = get_global_id(0) % NSPLINE_PROJ;
    
    // auxiliary variables
    private double scale, tx, ty;
    
    // for all source vector elements
    if (i < NSPLINE_ATOM)
    for (private int k = i - ORDER; k <= i + ORDER; k++) if (0 <= k && k < NSPLINE_ATOM)
    for (private int l = j - ORDER; l <= j + ORDER; l++) if (0 <= l && l < NSPLINE_PROJ)
    {
        // matrix element
        private double2 elem = 0;
        
        // Are the integral moments completely decoupled, i.e. there is there no overlap between Ba, Bb, Bc and Bd?
        // In such cases we can compute the off-diagonal contribution just as a product of the two
        // integral moments of order "lambda" and "-lambda-1", respectively. Moreover, in such
        // case there is no diagonal contribution, because there is no overlap between the _four_
        // participating B-splines.
        
        // (j,l) << (i,k)
        if (min(j,l) + ORDER < max(i,k))
        {
            tx = min(ta[i + ORDER + 1].x,ta[k + ORDER + 1].x);
            ty = min(tp[j + ORDER + 1].x,tp[l + ORDER + 1].x);
            scale = pow_int(ty / tx, lambda) / tx;
            elem += scale * cmul(MmLm1a[min(i,k) * (ORDER + 1) + abs(i-k)], MLp[min(j,l) * (ORDER + 1) + abs(j-l)]);
        }
        
        // (i,k) << (j,l)
        if (min(i,k) + ORDER < max(j,l))
        {
            tx = min(ta[i + ORDER + 1].x,ta[k + ORDER + 1].x);
            ty = min(tp[j + ORDER + 1].x,tp[l + ORDER + 1].x);
            scale = pow_int(tx / ty, lambda) / ty;
            elem += scale * cmul(MLa[min(i,k) * (ORDER + 1) + abs(i-k)], MmLm1p[min(j,l) * (ORDER + 1) + abs(j-l)]);
        }
        
        // all-to-all-segments multiplication
        for (int ill  = y_ang_begin; ill  < min(y_ang_begin + NDSTSEG, ANGULAR_BASIS_SIZE); ill ++)
        for (int illp = x_ang_begin; illp < min(x_ang_begin + NSRCSEG, ANGULAR_BASIS_SIZE); illp++)
        {
            y[((ill - y_ang_begin) * (ulong)(NSPLINE_ATOM) + i) * NSPLINE_PROJ + j] -=
                f[foffset + ill * ANGULAR_BASIS_SIZE + illp]
                * cmul(elem, x[((illp - x_ang_begin) * (ulong)(NSPLINE_ATOM) + k) * NSPLINE_PROJ + l]);
        }
    }
}

kernel void mmul_2el_coupled_offset
(
    // B-spline knots (atom and projectile)
    constant double2 const * const restrict ta,
    constant double2 const * const restrict tp,
    // multipole
    private int lambda,
    // angular integral
    global double  const * const restrict f, private int foffset,
    // one-electron partial moments (atom)
    global double2 const * const restrict MiLa,
    global double2 const * const restrict MimLm1a,
    // one-electron partial moments (projectile)
    global double2 const * const restrict MiLp,
    global double2 const * const restrict MimLm1p,
    // two-electron diagonal contributions
    global double2 const * const restrict Rdia,
    // source and target vector
    global double2 const * const restrict x,
    global double2       * const restrict y,
    // angular block domain
    private short x_ang_begin,
    private short y_ang_begin
)
{
    // output vector element index
    private int i = get_global_id(0) / (2 * ORDER + 1);
    private int j = get_global_id(0) % (2 * ORDER + 1) + i - ORDER;
    
    // auxiliary variables
    private double2 m_ik, m_jl;
    private double scale, tx, ty;
    
    // pointers to the needed partial integral moments
    global double2 const * restrict M_ik;
    global double2 const * restrict M_jl;
    
    // for all source vector elements
    if (i < NSPLINE_ATOM && 0 <= j && j <= NSPLINE_PROJ)
    for (private int k = i - ORDER; k <= i + ORDER; k++) if (0 <= k && k < NSPLINE_ATOM)
    for (private int l = j - ORDER; l <= j + ORDER; l++) if (0 <= l && l < NSPLINE_PROJ)
    if (max(i,k) <= min(j,l) + ORDER && max(j,l) <= min(i,k) + ORDER)
    {
        // matrix element
        private double2 elem = 0;
        
        // (i,k) ~ (j,l)
        
        // retrieve diagonal contribution
        if (i <= j && i <= k && i <= l)
            elem += Rdia[((i * (ORDER+1) + (j-i)) * (ORDER+1) + (k-i)) * (ORDER+1) + (l-i)];
        else if (j <= i && j <= k && j <= l)
            elem += Rdia[((j * (ORDER+1) + (i-j)) * (ORDER+1) + (l-j)) * (ORDER+1) + (k-j)];
        else if (k <= i && k <= j && k <= l)
            elem += Rdia[((k * (ORDER+1) + (j-k)) * (ORDER+1) + (i-k)) * (ORDER+1) + (l-k)];
        else // (l <= i && l <= j && l <= k)
            elem += Rdia[((l * (ORDER+1) + (i-l)) * (ORDER+1) + (j-l)) * (ORDER+1) + (k-l)];
        
        // ix < iy
        M_ik = MiLa    + (i * (2*ORDER+1) + k - (i-ORDER)) * (ORDER+1);
        M_jl = MimLm1p + (j * (2*ORDER+1) + l - (j-ORDER)) * (ORDER+1);
        for (private int ix = i; ix < min(i + ORDER + 1, NREKNOT_ATOM - 1); ix++) if (ta[ix + 1].x > 0)
        {
            m_ik = M_ik[ix - i]; tx = ta[ix + 1].x;
            
            for (private int iy = max(j, ix + 1); iy < min(j + ORDER + 1, NREKNOT_PROJ - 1); iy++)
            {
                m_jl = M_jl[iy - j]; ty = tp[iy + 1].x;
                scale = pow_int(tx/ty,lambda)/ty;
                elem += cmul(m_ik,m_jl) * scale;
            }
        }
        
        // ix > iy (by renaming the ix,iy indices)
        M_ik = MimLm1a + (i * (2*ORDER+1) + k - (i-ORDER)) * (ORDER+1);
        M_jl = MiLp    + (j * (2*ORDER+1) + l - (j-ORDER)) * (ORDER+1);
        for (private int ix = j; ix < min(j + ORDER + 1, NREKNOT_PROJ - 1); ix++) if (tp[ix + 1].x > 0)
        {
            m_jl = M_jl[ix - j]; tx = tp[ix + 1].x;
            
            for (private int iy = max(i, ix + 1); iy < min(i + ORDER + 1, NREKNOT_ATOM - 1); iy++)
            {
                m_ik = M_ik[iy - i]; ty = ta[iy + 1].x;
                scale = pow_int(tx/ty,lambda)/ty;
                elem += cmul(m_ik,m_jl) * scale;
            }
        }
        
        // all-to-all-segments multiplication
        for (int ill  = y_ang_begin; ill  < min(y_ang_begin + NDSTSEG, ANGULAR_BASIS_SIZE); ill ++)
        for (int illp = x_ang_begin; illp < min(x_ang_begin + NSRCSEG, ANGULAR_BASIS_SIZE); illp++)
        {
            y[((ill - y_ang_begin) * (ulong)(NSPLINE_ATOM) + i) * NSPLINE_PROJ + j] -=
                f[foffset + ill * ANGULAR_BASIS_SIZE + illp]
                * cmul(elem, x[((illp - x_ang_begin) * (ulong)(NSPLINE_ATOM) + k) * NSPLINE_PROJ + l]);
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
 * @param B Input k-by-n matrix in column-major(!) storage (i.e., transposed).
 * @param C Output m-by-n matrix in row-major storage.
 */
kernel void mul_ABt
(
    private int m, private int n, private int k,
    global double2 const * const restrict A,
    global double2 const * const restrict B,
    global double2       * const restrict C
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
    local double2 Aloc[BLOCK_SIZE][BLOCK_SIZE];
    local double2 Bloc[BLOCK_SIZE][BLOCK_SIZE];
    
    // aggregated scalar product of the destination element C[idy,idx]
    private double2 res = 0;
    
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
        Aloc[icol_loc][irow_loc] = (Arow < m && Acol < k ? A[Arow * k + Acol] : (double2)(0.,0.));
        Bloc[irow_loc][icol_loc] = (Brow < k && Bcol < n ? B[Brow + k * Bcol] : (double2)(0.,0.));
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
    private double2 E,
    global double2 const * const restrict D1,
    global double2 const * const restrict D2,
    global double2       * const restrict y
)
{
    // get worker's offset index (0 <= j < NSPLINE_PROJ)
    private int j = get_global_id(0);
    
    // y = y / (E (I kron I) - (D1 kron I) - (I kron D2))
    for (private int i = 0; i < NSPLINE_ATOM; i++)
        y[i * NSPLINE_PROJ + j] = cdiv(y[i * NSPLINE_PROJ + j], E - D1[i] - D2[j]);
}
