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
// -D ORDER=...
// -D NSPLINE=...
// -D NREKNOT=...
// -D NLOCAL=...

// Enable double precision.
#pragma OPENCL EXTENSION cl_khr_fp64: enable

// Derived variables.
#define NDIAG   ((2*ORDER+1)*(2*ORDER+1))
#define NROW    (NSPLINE * NSPLINE)
#define NCOL    (NSPLINE * NSPLINE)

/**
 * @brief Complex multiplication.
 * 
 * Multiplies two complex numbers and returns the product.
 */
inline double2 cmul (double2 a, double2 b)
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
inline double2 cdiv (double2 a, double2 b)
{
    double2 c;
    double b2 = b.x * b.x + b.y * b.y;
    
    c.x = (a.x * b.x + a.y * b.y) / b2;
    c.y = (a.y * b.x - a.x * b.y) / b2;
    
    return c;
}

/**
 * @brief AXBY operation.
 * 
 * Computes the linear combination @f$ ax + by @f$ of vectors @f$ x @f$ and 
 * @f$ y @f$ and stores the result in @f$ x @f$.
 */
kernel void a_vec_b_vec (private double2 a, global double2 *x, private double2 b, global double2 *y)
{
    uint i = get_global_id(0);
    
    if (i < NROW)
        x[i] = cmul(a,x[i]) + cmul(b,y[i]);
}

/**
 * @brief Vector-vector multiplication.
 * 
 * Computes per-element vector-vector multiplication @f$ a * b @f$ and stores
 * the resulting vector in @f$ c @f$,
 */
kernel void vec_mul_vec (global double2 *a, global double2 *b, global double2 *c)
{
    uint i = get_global_id(0);
    
    if (i < NROW)
        c[i] = cmul(a[i],b[i]);
}

/**
 * @brief Per-element square of a vector.
 * 
 * Similarly to @ref vec_mul_vec, computes a per-element square of
 * a given vector @f$ v @f$ and stores the result in the vector 
 * given as the second argument.
 */
kernel void vec_norm (global double2 *v, global double *n)
{
    uint i = get_global_id(0);
    double2 vi = v[i];
    n[i] = vi.x * vi.x + vi.y * vi.y;
}

/**
 * @brief Full scalar product.
 */
kernel void scalar_product (global double2 *u, global double2 *v, global double2 *z)
{
    uint iglobal = get_global_id(0);
    uint ilocal = get_local_id(0);
    
    local double2 uv[NLOCAL];
    
    if (iglobal < NROW)
        uv[ilocal] = cmul(u[iglobal],v[iglobal]);
    else
        uv[ilocal] = 0;
    
    // wait on completition of local write instructions
    barrier(CLK_LOCAL_MEM_FENCE);
    
    // reduce the per-element products
    int stride = 1;
    while (stride < NLOCAL)
    {
        if (ilocal % (2 * stride) == 0 && iglobal + stride < NROW)
            uv[ilocal] += uv[ilocal + stride];
        
        stride *= 2;
        
        // wait on completition of local write instructions
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    
    // write the reduced scalar product as computed by this group
    if (ilocal == 0)
        z[get_group_id(0)] = uv[0];
}

/**
 * @brief Full vector norm.
 */
kernel void norm (global double2 *v, global double *z)
{
    uint iglobal = get_global_id(0);
    uint ilocal = get_local_id(0);
    
    local double vv[NLOCAL];
    
    if (iglobal < NROW)
    {
        double2 vi = v[iglobal];
        vv[ilocal] = vi.x * vi.x + vi.y * vi.y;
    }
    else
        vv[ilocal] = 0;
    
    // wait for completition of local write instructions
    barrier(CLK_LOCAL_MEM_FENCE);
    
    // reduce the per-element products
    int stride = 1;
    while (stride < NLOCAL)
    {
        if (ilocal % (2 * stride) == 0 && iglobal + stride < NROW)
            vv[ilocal] += vv[ilocal + stride];
        
        stride *= 2;
        
        // wait for completition of local write instructions
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    
    // write the reduced scalar product as computed by this group
    if (ilocal == 0)
        z[get_group_id(0)] = vv[0];
}

/**
 * @brief Multiplication by a superblock.
 * 
 * This routine will multiply super-segment of the right-hand side by a super-block
 * of the matrix that corresponds to given angular momenta.
 */
kernel void mmul
(
    // energy
    double E,
    // row-padded one-electron matrices
    global double2 const * const restrict Sp,
    global double2 const * const restrict Dp,
    global double2 const * const restrict M1p,
    global double2 const * const restrict M2p,
    // angular momenta and nonzero multipoles
    int l1, int l2, int Nlambdas,
    global int     const * const restrict lambdas,
    // precomputed angular integrals
    global double  const * const restrict fs,
    // one-electron partial moments
    global double2 const * const restrict MiL,
    global double2 const * const restrict MimLm1,
    // source and target vector
    global double2 const * const restrict x,
    global double2       * const restrict y
)
{
    // block row index
    int i = get_group_id(0);
    
    // clear output vector
    for (int first_j = 0; first_j < NSPLINE; first_j += NLOCAL)
    {
        // get line index of the current worker
        int j = first_j + get_local_id(0);
        
        // erase the element
        if (j < NSPLINE)
            y[i * NSPLINE + j] = 0;
    }
    
    // for all block column indices
    for (int k = i - ORDER; k <= i + ORDER; k++) if (0 <= k && k < NSPLINE)
    {
        // loop over line groups
        for (int first_j = 0; first_j < NSPLINE; first_j += NLOCAL) if (first_j + get_local_id(0) < NSPLINE)
        {
            // get line index of the current worker
            int j = first_j + get_local_id(0);
            
            // scalar product of this block's line and the source vector segment
            double2 prod = 0;
            
            // for all elements of this block's line
            for (int l = j - ORDER; l <= j + ORDER; l++) if (0 <= l && l < NSPLINE)
            {
                // compute multi-indices
                int ik = min(i,k) * (ORDER + 1) + abs(i - k);
                int jl = min(j,l) * (ORDER + 1) + abs(l - j);
                
                // calculate the one-electron part of the hamiltonian matrix element Hijkl
                double2 elem = 0;
                elem += E * cmul(Sp[ik],Sp[jl]);
                elem -= 0.5 * (cmul(Dp[ik],Sp[jl]) + cmul(Sp[ik],Dp[jl]));
                elem -= 0.5 * l1 * (l1 + 1.) * cmul(M2p[ik],Sp[jl]) + 0.5 * l2 * (l2 + 1.) * cmul(Sp[ik],M2p[jl]);
                elem += cmul(M1p[ik],Sp[jl]) + cmul(Sp[ik],M1p[jl]);
                
                // calculate the simplified two-electron part of the hamiltonian matrix element Hijkl
                for (int ilambda = 0; ilambda < Nlambdas; ilambda++)
                {
                    // get pointers to the needed partial integral moments
                    global double2 const * const MiL_ik    = MiL    + ((lambdas[ilambda] * NSPLINE + i) * (2*ORDER+1) + k - (i-ORDER)) * (ORDER+1);
                    global double2 const * const MimLm1_ik = MimLm1 + ((lambdas[ilambda] * NSPLINE + i) * (2*ORDER+1) + k - (i-ORDER)) * (ORDER+1);
                    global double2 const * const MiL_jl    = MiL    + ((lambdas[ilambda] * NSPLINE + j) * (2*ORDER+1) + l - (j-ORDER)) * (ORDER+1);
                    global double2 const * const MimLm1_jl = MimLm1 + ((lambdas[ilambda] * NSPLINE + j) * (2*ORDER+1) + l - (j-ORDER)) * (ORDER+1);
                    
                    // ix < iy
                    for (int ix = i; ix <= i + ORDER && ix < NREKNOT - 1; ix++)
                    for (int iy = max(j,ix+1); iy <= j + ORDER && iy < NREKNOT - 1; iy++)
                    {
                        double2 m_ik = MiL_ik[ix - i], m_jl = MimLm1_jl[iy - j];
                        
                        // multiply real x real (merge exponents)
                        if (m_ik.y == 0 && m_jl.y == 0)
                        {
                            elem.x -= fs[ilambda] * exp(m_ik.x + m_jl.x);
                        }
                        
                        // multiply other cases
                        else
                        {
                            double2 M_ik = 0; if (m_ik.y == 0) M_ik.x = exp(m_ik.x); else M_ik = m_ik;
                            double2 M_jl = 0; if (m_jl.y == 0) M_jl.x = exp(m_jl.x); else M_jl = m_jl;
                            elem -= fs[ilambda] * cmul(M_ik,M_jl);
                        }
                    }
                    
                    // ix > iy (by renaming the ix,iy indices)
                    for (int ix = j; ix <= j + ORDER && ix < NREKNOT - 1; ix++)
                    for (int iy = max(i,ix+1); iy <= i + ORDER && iy < NREKNOT - 1; iy++)
                    {
                        double2 m_jl = MiL_jl[ix - j], m_ik = MimLm1_ik[iy - i];
                        
                        // multiply real x real (merge exponents)
                        if (m_ik.y == 0 && m_jl.y == 0)
                        {
                            elem.x -= fs[ilambda] * exp(m_ik.x + m_jl.x);
                        }
                        
                        // multiply other cases
                        else
                        {
                            double2 M_ik = 0; if (m_ik.y == 0) M_ik.x = exp(m_ik.x); else M_ik = m_ik;
                            double2 M_jl = 0; if (m_jl.y == 0) M_jl.x = exp(m_jl.x); else M_jl = m_jl;
                            elem -= fs[ilambda] * cmul(M_ik,M_jl);
                        }
                    }
                } // end for (ilambda)
                
                // multiply right-hand side by that matrix element
                prod += cmul(elem, x[k * NSPLINE + l]);
                
            } // end for (l)
            
            // update result vector
            y[i * NSPLINE + j] += prod;
            
        } // end for (turn) ... group of lines to process by the work-group
        
        // wait for completition of this block
        barrier(CLK_LOCAL_MEM_FENCE);
        
    } // end for (k) ... block index in row
}

/**
 * KPA preconditioner
 * 
 * 1) kron_dot1(CS1,CS2,r,C)
 * 2) kron_dot2(CS1,CS2,t,C)
 * 3) kron_div (E,D1,D2,t)
 * 4) kron_dot1(SC1,SC2,t,C)
 * 5) kron_dot2(SC1,SC2,z,C)
 */
kernel void kron_dot1
(
    global double2 const * const restrict A,
    global double2 const * const restrict B,
    global double2 const * const restrict v,
    global double2       * const restrict C
)
{
    // get worker's segment index (0 <= i < NSPLINE)
    int seg = get_global_id(0);
    
    // for all rows of B
    for (int irow = 0; irow < NSPLINE; irow++)
    {
        // scalar product of the current row of B and the worker's segment
        double2 prod = 0;
        for (int icol = 0; icol < NSPLINE; icol++)
            prod = prod + cmul(B[irow * NSPLINE + icol],v[seg * NSPLINE + icol]);
        C[irow * NSPLINE + seg] = prod;
    }
}
kernel void kron_dot2
(
    global double2 const * const restrict A,
    global double2 const * const restrict B,
    global double2       * const restrict w,
    global double2 const * const restrict C
)
{
    // get worker's segment index (0 <= i < NSPLINE)
    int seg = get_global_id(0);
    
    // for all rows of C
    for (int irow = 0; irow < NSPLINE; irow++)
    {
        // scalar product of the current row of A and C
        double2 prod = 0;
        for (int icol = 0; icol < NSPLINE; icol++)
            prod = prod + cmul(A[seg * NSPLINE + icol],C[irow * NSPLINE + icol]);
        w[seg * NSPLINE + irow] = prod;
    }
}
kernel void kron_div
(
    double2 E,
    global double2 const * const restrict D1,
    global double2 const * const restrict D2,
    global double2       * const restrict y
)
{
    // get worker's segment index (0 <= i < NSPLINE)
    int i = get_global_id(0);
    
    // y = y / (E (I kron I) - (D1 kron I) - (I kron D2))
    for (int j = 0; j < NSPLINE; j++)
        y[j * NSPLINE + i] = cdiv(y[j * NSPLINE + i], E - D1[j] - D2[i]);
}
