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

#pragma OPENCL EXTENSION cl_khr_fp64: enable

// -D ORDER=...
// #define ORDER 5

// -D NSPLINE=...
// #define NSPLINE 117

// -D DIAGONALS=...
// #define DIAGONALS \
//     -590, -589, -588, -587, -586, -585, -584, -583, -582, -581, -580, \
//     -473, -472, -471, -470, -469, -468, -467, -466, -465, -464, -463, \
//     -356, -355, -354, -353, -352, -351, -350, -349, -348, -347, -346, \
//     -239, -238, -237, -236, -235, -234, -233, -232, -231, -230, -229, \
//     -122, -121, -120, -119, -118, -117, -116, -115, -114, -113, -112, \
//       -5,   -4,   -3,   -2,   -1,    0,    1,    2,    3,    4,    5, \
//      112,  113,  114,  115,  116,  117,  118,  119,  120,  121,  122, \
//      229,  230,  231,  232,  233,  234,  235,  236,  237,  238,  239, \
//      346,  347,  348,  349,  350,  351,  352,  353,  354,  355,  356, \
//      463,  464,  465,  466,  467,  468,  469,  470,  471,  472,  473, \
//      580,  581,  582,  583,  584,  585,  586,  587,  588,  589,  590
//
// -D NLOCAL=...
// #define NLOCAL 64

// Derived variables.
#define NDIAG   ((2*ORDER+1)*(2*ORDER+1))
#define NROW    (NSPLINE * NSPLINE)
#define NCOL    (NSPLINE * NSPLINE)

/**
 * @brief Complex multiplication.
 * 
 * Multiplies two complex numbers and returns the product.
 */
inline double2 cpxm (double2 a, double2 b)
{
    double2 c;
    c.x = a.x * b.x - a.y * b.y;
    c.y = a.x * b.y + a.y * b.x;
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
    x[i] = cpxm(a,x[i]) + cpxm(b,y[i]);
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
    c[i] = cpxm(a[i],b[i]);
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
    double2 ui = u[iglobal];
    double2 vi = v[iglobal];
    
    uint ilocal = get_local_id(0);
    local double2 uv[NLOCAL];
    uv[ilocal].x = ui.x * vi.x - ui.y * vi.y;
    uv[ilocal].y = ui.x * vi.y + ui.y * vi.x;
    
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
    double2 vi = v[iglobal];
    
    uint ilocal = get_local_id(0);
    local double vv[NLOCAL];
    vv[ilocal] = vi.x * vi.x + vi.y * vi.y;
    
    // wait on completition of local write instructions
    barrier(CLK_LOCAL_MEM_FENCE);
    
    // reduce the per-element products
    int stride = 1;
    while (stride < NLOCAL)
    {
        if (ilocal % (2 * stride) == 0 && iglobal + stride < NROW)
            vv[ilocal] += vv[ilocal + stride];
        
        stride *= 2;
        
        // wait on completition of local write instructions
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    
    // write the reduced scalar product as computed by this group
    if (ilocal == 0)
        z[get_group_id(0)] = vv[0];
}

/**
 * @brief DIA-matrix-vector multiplication.
 * 
 * The input matrix is supplied in the form of zero-padded concatenated diagonals.
 * For example, the matrix
 * @f[
 *     A = \pmatrix {
 *         a_{11} & a_{12} &        &        &        \cr
 *         a_{21} & a_{22} & a_{23} &        &        \cr
 *                & a_{32} & a_{33} & a_{34} &        \cr
 *                &        & a_{43} & a_{44} & a_{45} \cr
 *                &        &        & a_{54} & a_{55} \cr
 *     }
 * @f]
 * would be column-padded to
 * @f[
 *     A' = \pmatrix {
 *              0 & a_{11} & a_{12} \cr
 *         a_{21} & a_{22} & a_{23} \cr
 *         a_{32} & a_{33} & a_{34} \cr
 *         a_{43} & a_{44} & a_{45} \cr
 *         a_{54} & a_{55} &      0 \cr
 *     }
 * @f]
 * and sent in as a 1D array
 * @f[
 *     A'' = [0, a_{21} , a_{32}, a_{43}, a_{54}; a_{11}, a_{22}, a_{33}, a_{44}, a_{55}; a_{12}, a_{23}, a_{34}, a_{45}, 0] \ .
 * @f]
 */
kernel void DIA_dot_vec (global double2 *A, global double2 *x, global double2 *y)
{
    // diagonal labels
    const int diagonals [] = { DIAGONALS };
    
    // matrix row for this local thread
    int irow = get_global_id(0);
    
    // only do something if this matrix row exists
    if (irow < NROW)
    {
        // scalar product
        double2 yloc = 0;
        
        // for all diagonals
        for (int idiag = 0; idiag < NDIAG; idiag++)
        {
            // get column index for all threads
            int icol = irow + diagonals[idiag];
            
            // multiply the elements
            if (0 <= icol && icol < NCOL)
            {
                yloc += cpxm(A[idiag * NROW + irow], x[icol]);
            }
        }
        
        // copy results to global memory
        y[irow] = yloc;
    }
}
