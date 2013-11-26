/**
 * @brief Complex multiplication.
 */
inline double2 complex_multiply (double2 a, double2 b)
{
    double2 c;
    c.x = a.x * b.x - a.y * b.y;
    c.y = a.x * b.y + a.y * b.x;
    return c;
}

/**
 * @brief AXBY operation.
 */
kernel void a_vec_b_vec (private double2 a, global double2 *x, private double2 b, global double2 *y)
{
    uint i = get_global_id(0);
    x[i] = complex_multiply(a,x[i]) + complex_multiply(b,y[i]);
}

/**
 * @brief Vector-vector multiplication.
 */
kernel void vec_mul_vec (global double2 *a, global double2 *b, global double2 *c)
{
    uint i = get_global_id(0);
    c[i] = complex_multiply(a[i],b[i]);
}

/**
 * @brief Per-element square of a vector.
 */
kernel void vec_norm (global double2 *v, global double *n)
{
    uint i = get_global_id(0);
    double2 vi = v[i];
    n[i] = vi.x * vi.x + vi.y * vi.y;
}

/**
 * @brief CSR-matrix-vector multiplication.
 */
kernel void CSR_dot_vec (global long *Ap, global long *Ai, global double2 *Ax, global double2 *x, global double2 *y)
{
    uint i = get_global_id(0);
    double2 sprod = 0.;
    
    for (int idx = Ap[i]; idx < Ap[i + 1]; idx++)
        sprod += complex_multiply(Ax[idx],x[Ai[idx]]);
    
    y[i] = sprod;
}

/**
 * @brief DIA-matrix-vector multiplication.
 */
kernel void DIA_dot_vec (global double2 *A, global double2 *x, global double2 *y)
{
    // B-spline order and derived numbers to be supplied before compilation
    const uint order = 5;
    const uint Ndiag = 121;
    const uint Nlocal = 128;
    
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
    
    // get group indices (numbering the rows)
    int irow = get_group_id(0);
    int Nrow = get_num_groups(0);
    
    // get local indices (numbering the diagonals)
    int ilocal = get_local_id(0);
    
    // get column index
    int icol = irow + diagonals[ilocal];
    
    // --------------------------------------------------------------------- //
    
    // define and clear the temporary result array
    local double2 tmp[Nlocal];
    tmp[ilocal] = 0.;
    
    // multiply matrix row by the vector (per component)
    if (ilocal < Ndiag)
    {
        tmp[ilocal] = A[ilocal] * x[icol];
    }
    
    // filter out NaN-s
    if (!(tmp[ilocal] == tmp[ilocal]))
    {
        tmp[ilocal] = 0.;
    }
    
    // wait for finish of local memory writes
    barrier(CLK_LOCAL_MEM_FENCE);
    
    // reduce the dot product by bisection
    uint stride = 1, size = Nlocal;
    while (size > 1)
    {
        if (ilocal % (2 * stride) == 0)
            tmp[ilocal] += tmp[ilocal + stride];
        
        stride *= 2;
        size /= 2;
    }
    
    // --------------------------------------------------------------------- //
    
    // group master will store the result of the dot product
    if (ilocal == 0)
    {
        y[irow] += tmp[0];
    }
}
