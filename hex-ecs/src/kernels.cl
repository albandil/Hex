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
 * @brief CSR-matrix-vector multiplication.
 * 
 * Multiplies the given vector @f$ x @f$ by a compressed sparse row-major
 * matrix given as three arrays (row pointers, column indices and matrix
 * elements). The resulting vector is stored in @f$ y @f$.
 */
kernel void CSR_dot_vec (global long *Ap, global long *Ai, global double2 *Ax, global double2 *x, global double2 *y)
{
    uint i = get_global_id(0);
    double2 sprod = 0.;
    
    for (int idx = Ap[i]; idx < Ap[i + 1]; idx++)
        sprod += cpxm(Ax[idx],x[Ai[idx]]);
    
    y[i] = sprod;
}

/*kernel void DIA_dot_vec (global double2 *A, global double2 *x, global double2 *y)
{
    // some compile-time constants ------------------------------------- //
    //
    
    // NOTE : Need to be modified by the program if
    //                 order != 5        or
    //               Nspline != 114
    
    #define Ndiag   121    // = (2*order + 1)^2
    #define Nlocal  128    // = Ndiag rounded to multiple of 64
    #define Nrow    13689  // = Nspline^2
    #define Ncol    13689  // = Nspline^2
    
    // diagonal labels
    const int diagonals [] = {
        -590, -589, -588, -587, -586, -585, -584, -583, -582, -581, -580,
        -473, -472, -471, -470, -469, -468, -467, -466, -465, -464, -463,
        -356, -355, -354, -353, -352, -351, -350, -349, -348, -347, -346,
        -239, -238, -237, -236, -235, -234, -233, -232, -231, -230, -229,
        -122, -121, -120, -119, -118, -117, -116, -115, -114, -113, -112,
          -5,   -4,   -3,   -2,   -1,    0,    1,    2,    3,    4,    5,
         112,  113,  114,  115,  116,  117,  118,  119,  120,  121,  122,
         229,  230,  231,  232,  233,  234,  235,  236,  237,  238,  239,
         346,  347,  348,  349,  350,  351,  352,  353,  354,  355,  356,
         463,  464,  465,  466,  467,  468,  469,  470,  471,  472,  473,
         580,  581,  582,  583,  584,  585,  586,  587,  588,  589,  590
    };    
                                                                         //
    // ----------------------------------------------------------------- //
    
    // get group indices (numbering the rows)
    int irow = get_group_id(0);
    
    // get global indices (numbering the consecutive elements of matrix)
    int ielem = get_global_id(0) % get_local_size(0) + irow * Ndiag;
    
    // get local indices (numbering the diagonals)
    int ilocal = get_local_id(0);
        
    // get column indices
    int icol = irow + diagonals[ilocal];
    
    // define and clear the temporary result array
    local double2 tmp[Nlocal];
    tmp[ilocal] = 0.;
    
    // multiply matrix row by the vector (per component, mask by column range)
    if (ilocal < Ndiag && 0 <= icol && icol < Ncol)
        tmp[ilocal] = cpxm(A[ielem],x[icol]);
    
    // wait for finish of local memory writes
    barrier(CLK_LOCAL_MEM_FENCE);
    
    // reduce the dot product by bisection
    uint stride = 1, size = Nlocal;
    while (size > 1)
    {
        // sum elements
        tmp[ilocal] += tmp[ilocal + stride];
        
        // update stride and size
        stride *= 2;
        size /= 2;
        
        // wait for finish of local memory writes
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    
    // group master will store the result of the dot product
    if (ilocal == 0)
        y[irow] = tmp[0];
}*/
