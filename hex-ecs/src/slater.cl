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

#pragma OPENCL EXTENSION cl_khr_fp64: enable

// Define the following macros:
//   -D ORDER=...
//   -D IKNOTMAX=...  (Nreknot-1)

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

double2 cexp (const double2 z)
{
    double2 w;
    w.x = cos(z.y);
    w.y = sin(z.y);
    return w * exp(z.x);
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
 * This function evaluates a single B-spline on a given knot and returns
 * its value as a complex number.
 * 
 * @note This function is implemented for fourth-order B-spline only.
 * 
 * @param t Complex array containing B-spline knots.
 * @param i B-spline index.
 * @param iknot Which interval is it that we evaluate the B-spline on.
 * @param x Complex coordinate from the interval given by "iknot".
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
    double2 b0 = 0.0, b1 = 0.0, b2 = 0.0, b3 = 0.0, b4 = 0.0;
    b0.x = (i + 0 == iknot ? 1.0 : 0.0);
    b1.x = (i + 1 == iknot ? 1.0 : 0.0);
    b2.x = (i + 2 == iknot ? 1.0 : 0.0);
    b3.x = (i + 3 == iknot ? 1.0 : 0.0);
    b4.x = (i + 4 == iknot ? 1.0 : 0.0);
    
    // compute first-order B-splines (k = 1)
    b0 = ((t[i+1].x == t[i+0].x) ? 0.0 : cmul(b0,cdiv(x-t[i+0],t[i+1]-t[i+0])))
       + ((t[i+2].x == t[i+1].x) ? 0.0 : cmul(b1,cdiv(t[i+2]-x,t[i+2]-t[i+1])));
    b1 = ((t[i+2].x == t[i+1].x) ? 0.0 : cmul(b1,cdiv(x-t[i+1],t[i+2]-t[i+1])))
       + ((t[i+3].x == t[i+2].x) ? 0.0 : cmul(b2,cdiv(t[i+3]-x,t[i+3]-t[i+2])));
    b2 = ((t[i+3].x == t[i+2].x) ? 0.0 : cmul(b2,cdiv(x-t[i+2],t[i+3]-t[i+2])))
       + ((t[i+4].x == t[i+3].x) ? 0.0 : cmul(b3,cdiv(t[i+4]-x,t[i+4]-t[i+3])));
    b3 = ((t[i+4].x == t[i+3].x) ? 0.0 : cmul(b3,cdiv(x-t[i+3],t[i+4]-t[i+3])))
       + ((t[i+5].x == t[i+4].x) ? 0.0 : cmul(b4,cdiv(t[i+5]-x,t[i+5]-t[i+4])));
    
    // compute second-order B-splines (k = 2)
    b0 = ((t[i+2].x == t[i+0].x) ? 0.0 : cmul(b0,cdiv(x-t[i+0],t[i+2]-t[i+0])))
       + ((t[i+3].x == t[i+1].x) ? 0.0 : cmul(b1,cdiv(t[i+3]-x,t[i+3]-t[i+1])));
    b1 = ((t[i+3].x == t[i+1].x) ? 0.0 : cmul(b1,cdiv(x-t[i+1],t[i+3]-t[i+1])))
       + ((t[i+4].x == t[i+2].x) ? 0.0 : cmul(b2,cdiv(t[i+4]-x,t[i+4]-t[i+2])));
    b2 = ((t[i+4].x == t[i+2].x) ? 0.0 : cmul(b2,cdiv(x-t[i+2],t[i+4]-t[i+2])))
       + ((t[i+5].x == t[i+3].x) ? 0.0 : cmul(b3,cdiv(t[i+5]-x,t[i+5]-t[i+3])));
    
    // compute third-order B-splines (k = 3)
    b0 = ((t[i+3].x == t[i+0].x) ? 0.0 : cmul(b0,cdiv(x-t[i+0],t[i+3]-t[i+0])))
       + ((t[i+4].x == t[i+1].x) ? 0.0 : cmul(b1,cdiv(t[i+4]-x,t[i+4]-t[i+1])));
    b1 = ((t[i+4].x == t[i+1].x) ? 0.0 : cmul(b1,cdiv(x-t[i+1],t[i+4]-t[i+1])))
       + ((t[i+5].x == t[i+2].x) ? 0.0 : cmul(b2,cdiv(t[i+5]-x,t[i+5]-t[i+2])));
    
    // compute fourth-order B-spline (k = 4)
    b0 = ((t[i+4].x == t[i+0].x) ? 0.0 : cmul(b0,cdiv(x-t[i+0],t[i+4]-t[i+0])))
       + ((t[i+5].x == t[i+1].x) ? 0.0 : cmul(b1,cdiv(t[i+5]-x,t[i+5]-t[i+1])));
    
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
    return (r > R) ? 0.0 : tanh(0.125 * (R - r));
}

/**
 * @brief Two-electron ("Slater") integrals.
 * 
 * This kernel should be called once to compute all integrals.
 * Alternatively, it can be called as a series to apply consecutively
 * on given chunks of the index array "idx_R". The chunk size
 * is always equal to the global work item count and the chunk
 * origin is given by the argument "offset".
 * 
 * @param lambda    Multipole moment (angular momentum transfer).
 * @param offset    Which work-chunk of integral to compute.
 * @param t         Knot sequence.
 * @param xIn0      Inner Gauss-Legendre quadrature nodes.
 * @param wIn0      Inner Gauss-Legendre quadrature weights.
 * @param xOut0     Outer Gauss-Legendre quadrature nodes.
 * @param wOut0     Outer Gauss-Legendre quadrature weights.
 * @param Mtr_L     Per-spline and -knot integral moments (power: lambda)
 * @param Mtr_mLm1  Per-spline and -knot integral moments (power: -lambda-1)
 * @param idx_R     In-out array. On entry every 128 bits contain 4 integer indices
 *                  of the integral to compute. On return, the same indices contain
 *                  one double precision complex number, which is the value of the
 *                  integral.
 */
kernel void R_integral
(
    const int lambda,
    const ulong offset,
    constant double2 const * const restrict t,
    constant double  const * const restrict xIn0,
    constant double  const * const restrict wIn0,
    constant double  const * const restrict xOut0,
    constant double  const * const restrict wOut0,
    global double2 const * const restrict Mtr_L,
    global double2 const * const restrict Mtr_mLm1,
    global double2 * const restrict idx_R
)
{
    //
    // preparations
    //
    
    // get worker's ID (= which two-electron integral to compute)
    const ulong iwork = offset + get_global_id(0);
    
    // get integral index 4-tuple
    const int4 idx = as_int4(idx_R[iwork]);
    
    //
    // the off-diagonal contribution to the integral
    //
    
    double2 Rtr_Labcd_offdiag = 0.0;
    
    global double2 const * const restrict Mtr_L_ac    = Mtr_L    + (idx.x * (2*ORDER+1) + idx.z - (idx.x-ORDER)) * (ORDER+1);
    global double2 const * const restrict Mtr_mLm1_ac = Mtr_mLm1 + (idx.x * (2*ORDER+1) + idx.z - (idx.x-ORDER)) * (ORDER+1);
    global double2 const * const restrict Mtr_L_bd    = Mtr_L    + (idx.y * (2*ORDER+1) + idx.w - (idx.y-ORDER)) * (ORDER+1);
    global double2 const * const restrict Mtr_mLm1_bd = Mtr_mLm1 + (idx.y * (2*ORDER+1) + idx.w - (idx.y-ORDER)) * (ORDER+1);
    
    for (int ix = 0; ix < IKNOTMAX; ix++)
    {
        for (int iy = ix + 1; iy < IKNOTMAX; iy++)
        {
            // ix < iy
            if (idx.x <= ix && ix <= idx.x + ORDER && idx.y <= iy && iy <= idx.y + ORDER)
            {
                const double2 lg = Mtr_L_ac[ix - idx.x] + Mtr_mLm1_bd[iy - idx.y];
                if (isfinite(lg.y))
                    Rtr_Labcd_offdiag += cexp(lg);
            }
            
            // ix > iy (by renaming the ix,iy indices)
            if (idx.y <= ix && ix <= idx.y + ORDER && idx.x <= iy && iy <= idx.x + ORDER)
            {
                const double2 lg = Mtr_L_bd[ix - idx.y] + Mtr_mLm1_ac[iy - idx.x];
                if (isfinite(lg.y))
                    Rtr_Labcd_offdiag += cexp(lg);
            }
        }
    }
    
    //
    // the diagonal contribution to the integral
    // 
    
    double2 Rtr_Labcd_diag = 0.0;
    
    // get evaluation points count
    const int nIn  = ORDER + lambda + 1;
    const int nOut = ORDER + lambda + 10;
    
    // for all knots
    for (int iknot = 0; iknot < IKNOTMAX; iknot++)
    {
        if (t[iknot+1].x != t[iknot].x)
        {
            for (int u = 0; u < nOut; u++)
            {
                // evaluation points and weights for (outer) interval
                //     x in (t[iknot],t[iknot+1])
                const double2 xOut = 0.5 * ((1.0 - xOut0[u]) * t[iknot] + (1.0 + xOut0[u]) * t[iknot+1]);
                const double2 wOut = 0.5 * (t[iknot+1] - t[iknot]) * wOut0[u];
                
                // evaluate outer B-splines
                const double2 values_i1 = B_single(t, idx.x, iknot, xOut);
                const double2 values_j1 = B_single(t, idx.z, iknot, xOut);
                const double2 values_i2 = B_single(t, idx.y, iknot, xOut);
                const double2 values_j2 = B_single(t, idx.w, iknot, xOut);
                
                // for all inner points
                double2 evalIn1 = 0.0, evalIn2 = 0.0;
                for (int v = 0; v < nIn; v++)
                {
                    // evaluation points and weights for (inner) interval
                    //     y in (t[iknot],x)
                    const double2 xIn = 0.5 * ((1.0 - xIn0[v]) * t[iknot] + (1.0 + xIn0[v]) * xOut);
                    const double2 wIn = 0.5 * (xOut - t[iknot]) * wIn0[v];
                    
                    // evaluate inner B-splines
                    const double2 values_k1 = B_single(t, idx.y, iknot, xIn);
                    const double2 values_l1 = B_single(t, idx.w, iknot, xIn);
                    const double2 values_k2 = B_single(t, idx.x, iknot, xIn);
                    const double2 values_l2 = B_single(t, idx.z, iknot, xIn);
                    
                    // update inner integral
                    evalIn1 += cmul(wIn,cmul(values_k1,cmul(values_l1,cpow(cdiv(xIn,xOut),lambda)))) * damp(xIn, 0, t[IKNOTMAX].x);
                    evalIn2 += cmul(wIn,cmul(values_k2,cmul(values_l2,cpow(cdiv(xIn,xOut),lambda)))) * damp(xIn, 0, t[IKNOTMAX].x);
                }
                
                // evaluate outer integral
                Rtr_Labcd_diag += (
                    cmul(wOut,cmul(values_i1,cmul(values_j1,cdiv(evalIn1,xOut))))
                  + cmul(wOut,cmul(values_i2,cmul(values_j2,cdiv(evalIn2,xOut))))
                ) * damp(0, xOut, t[IKNOTMAX].x);
            }
        }
    }
    
    //
    // store the result to the "in-out" array
    //
    
    idx_R[iwork] = Rtr_Labcd_offdiag + Rtr_Labcd_diag;
}
