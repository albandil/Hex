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

// Define the following macros:
//   -D ORDER=...
//   -D NSPLINE=...

// ------------------------------------------------------------------------- //
// Complex arithmetic.                                                       //
// ------------------------------------------------------------------------- //

double2 cmul (const double2 a, const double2 b)
{
    // c = a * b
    double2 c;
    c.x = a.x * b.x - a.y * b.y;
    c.y = a.x * b.y + a.y * b.x;
    return c;
}

double2 cdiv (const double2 a, const double2 b)
{
    // c = a / b
    double2 c;
    c.x = (a.x * b.x + a.y * b.y) / (b.x * b.x + b.y * b.y);
    c.y = (a.y * b.x - a.x * b.y) / (b.x * b.x + b.y * b.y);
    return c;
}

double2 cpow (const double2 a, const double x)
{
    // b = a ^ x
    double2 b;
    const double phi = x * atan2(a.y,a.x);
    const double r = pow(a.x * a.x + a.y * a.y, 0.5 * x);
    b.x = r * cos(phi);
    b.y = r * sin(phi);
    return b;
}

// ------------------------------------------------------------------------- //
// B-spline evaluation.                                                      //
// ------------------------------------------------------------------------- //

/**
 * @brief Evaluate single B-spline.
 * 
 * This function evaluates a single B-spline on a given knot.
 */
double2 B_single
(
    constant double2 const * const restrict t,
    const int i,
    const int iknot,
    const double2 x
)
{
    // FIXME Assuming ORDER = 4
    
    // compute zero-order B-splines (k = 0)
    double2 b0, b1, b2, b3, b4;
    b0.x = (i + 0 == iknot ? 1.0 : 0.0);
    b1.x = (i + 1 == iknot ? 1.0 : 0.0);
    b2.x = (i + 2 == iknot ? 1.0 : 0.0);
    b3.x = (i + 3 == iknot ? 1.0 : 0.0);
    b4.x = (i + 4 == iknot ? 1.0 : 0.0);
    
    // compute first-order B-splines (k = 1)
    b0 = ((t[i+1] == t[i+0]) ? 0.0 : cmul(b0,cdiv(x-t[i+0],t[i+1]-t[i+0])))
       + ((t[i+2] == t[i+1]) ? 0.0 : cmul(b1,cdiv(t[i+2]-x,t[i+2]-t[i+1])));
    b1 = ((t[i+2] == t[i+1]) ? 0.0 : cmul(b1,cdiv(x-t[i+1],t[i+2]-t[i+1])))
       + ((t[i+3] == t[i+2]) ? 0.0 : cmul(b2,cdiv(t[i+3]-x,t[i+3]-t[i+2])));
    b2 = ((t[i+3] == t[i+2]) ? 0.0 : cmul(b2,cdiv(x-t[i+2],t[i+3]-t[i+2])))
       + ((t[i+4] == t[i+3]) ? 0.0 : cmul(b3,cdiv(t[i+4]-x,t[i+4]-t[i+3])));
    b3 = ((t[i+4] == t[i+3]) ? 0.0 : cmul(b3,cdiv(x-t[i+3],t[i+4]-t[i+3])))
       + ((t[i+5] == t[i+4]) ? 0.0 : cmul(b4,cdiv(t[i+5]-x,t[i+5]-t[i+4])));
    
    // compute second-order B-splines (k = 2)
    b0 = ((t[i+2] == t[i+0]) ? 0.0 : cmul(b0,cdiv(x-t[i+0],t[i+2]-t[i+0])))
       + ((t[i+3] == t[i+1]) ? 0.0 : cmul(b1,cdiv(t[i+3]-x,t[i+3]-t[i+1])));
    b1 = ((t[i+3] == t[i+1]) ? 0.0 : cmul(b1,cdiv(x-t[i+1],t[i+3]-t[i+1])))
       + ((t[i+4] == t[i+2]) ? 0.0 : cmul(b2,cdiv(t[i+4]-x,t[i+4]-t[i+2])));
    b2 = ((t[i+4] == t[i+2]) ? 0.0 : cmul(b2,cdiv(x-t[i+2],t[i+4]-t[i+2])))
       + ((t[i+5] == t[i+3]) ? 0.0 : cmul(b3,cdiv(t[i+5]-x,t[i+5]-t[i+3])));
    
    // compute third-order B-splines (k = 3)
    b0 = ((t[i+3] == t[i+0]) ? 0.0 : cmul(b0,cdiv(x-t[i+0],t[i+3]-t[i+0])))
       + ((t[i+4] == t[i+1]) ? 0.0 : cmul(b1,cdiv(t[i+4]-x,t[i+4]-t[i+1])));
    b1 = ((t[i+4] == t[i+1]) ? 0.0 : cmul(b1,cdiv(x-t[i+1],t[i+4]-t[i+1])))
       + ((t[i+5] == t[i+2]) ? 0.0 : cmul(b2,cdiv(t[i+5]-x,t[i+5]-t[i+2])));
    
    // compute fourth-order B-spline (k = 4)
    b0 = ((t[i+4] == t[i+0]) ? 0.0 : cmul(b0,cdiv(x-t[i+0],t[i+4]-t[i+0])))
       + ((t[i+5] == t[i+1]) ? 0.0 : cmul(b1,cdiv(t[i+5]-x,t[i+5]-t[i+1])));
    
    return b0;
}

// ------------------------------------------------------------------------- //
// Radial integral functions.                                                //
// ------------------------------------------------------------------------- //

/**
 * @brief Damping factor.
 * 
 * Factor that will damp the potential.
 */
double damp (const double2 b, const double2 a, const double R)
{
    // compute hyperradius
    const double r = hypot(b.x, a.x);
    
    // if sufficiently far, return clean zero
    // else damp using tanh(x) distribution
    return return (r > R) ? 0.0 : tanh(0.125 * (R - r));
}

/**
 * @brief Triangular part of the two-electron ("Slater") integral.
 * 
 * This kernel should be called in the following way:
 * - One work-group per knot.
 * - Work-group size equal to (or larger than) the outer quadrature node count.
 * 
 * @param idx    R-integral four-index (i,j,k,l).
 * @param lambda Multipole moment (angular momentum transfer).
 * @param t      Knot sequence.
 * @param xIn0   Inner Gauss-Legendre quadrature nodes.
 * @param wIn0   Inner Gauss-Legendre quadrature weights.
 * @param xOut0  Outer Gauss-Legendre quadrature nodes.
 * @param wOut0  Outer Gauss-Legendre quadrature weights.
 * @param R_out  Output array for knot contributions.
 * 
 * The length of the array R_out should be equal to (or larger than) the
 * work-group count.
 */
kernel void R_integral
(
    const int4 idx,
    const int lambda,
    constant double2 const * const restrict t,
    constant double  const * const restrict xIn0,
    constant double  const * const restrict wIn0,
    constant double  const * const restrict xOut0,
    constant double  const * const restrict wOut0,
    global   double2       * const restrict R_out
)
{
    // get knot and outer quadrature node index
    const int iknot = get_group_id(0);
    const int ioutp = get_local_id(0);
    
    // get evaluation points count
    const int nIn  = ORDER + lambda + 1;
    const int nOut = ORDER + lambda + 10;
    
    // resulting contribution to the integral
    // - one knot per work-group
    // - one quadrature point per in-group (local) worker
    local double2 R_local[NOUTMAX];
    
    // compute only if the interval is non-empty
    if (!(t[iknot+1] == t[iknot]) && ioutp < NOUTMAX)
    {
        // evaluation points and weights for (outer) interval
        //     x in (t[iknot],t[iknot+1])
        const double2 xOut = (t[iknot+1] - t[iknot]) * xOut0[u] + t[iknot];
        const double2 wOut = (t[iknot+1] - t[iknot]) * wOut0[u];
        
        // evaluate outer B-splines
        const double2 values_i = B_single(t, idx.x, iknot, xOut[ioutp]);
        const double2 values_j = B_single(t, idx.y, iknot, xOut[ioutp]);
        
        // for all inner points
        double2 evalIn = 0.;
        for (int v = 0; v < nIn; v++)
        {
            // evaluation points and weights for (inner) interval
            //     y in (t[iknot],x)
            const double2 xIn = (xOut[ioutp] - t[iknot]) * xIn0[v] + t[iknot];
            const double2 wIn = (xOut[ioutp] - t[iknot]) * wIn0[v];
            
            // evaluate inner B-splines
            const double values_k = B_single(t, idx.z, iknot, xIn[v]);
            const double values_l = B_single(t, idx.w, iknot, xIn[v]);
            
            // update inner integral
            evalIn += cmul(wIn,cmul(values_k,cmul(values_l,cpow(cdiv(xIn,xOut),lambda)))) * damp(xIn, 0, t[IKNOTMAX].x);
        }
        
        // evaluate outer integral
        R_local[ioutp] = cmul(wOut,cmul(values_i,cmul(values_j,cdiv(evalIn,xOut)))) * damp(0, xOut, t[IKNOTMAX].x);
    }
    
    // wait for local memory flush
    barrier(CLK_LOCAL_MEM_FENCE);
    
    // master worker will reduce the contributions from this work-group
    if (ioutp == 0)
    {
        double2 R = 0.;
        
        for (int ipt = 0; ipt < nOut; ipt++)
            R = R + R_local[ipt];
        
        R_out[iknot] = R;
    }
}
