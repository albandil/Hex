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

inline double2 cmul (double2 a, double2 b)
{
    // c = a * b
    double2 c;
    c.x = a.x * b.x - a.y * b.y;
    c.y = a.x * b.y + a.y * b.x;
    return c;
}

double2 cdiv (double2 a, double2 b)
{
    // c = a / b
    double2 c;
    c.x = (a.x * b.x + a.y * b.y) / (b.x * b.x + b.y * b.y);
    c.y = (a.y * b.x - a.x * b.y) / (b.x * b.x + b.y * b.y);
    return c;
}

double2 cpow (double2 a, double x)
{
    // b = a ^ x
    double2 b;
    double phi = x * atan2(a.y,a.x), r = pow(a.x * a.x + a.y * a.y, 0.5 * x);
    b.x = r * cos(phi);
    b.y = r * sin(phi);
    return b;
}

// ------------------------------------------------------------------------- //
// B-spline evaluation.                                                      //
// ------------------------------------------------------------------------- //

void B
(
    constant double2 const * const restrict t,
    int i, int iknot, int n,
    double2 const * const restrict x,
    double2       * const restrict y
)
{
    double2 z[ORDER + 1];
    
    for (int u = 0; u < n; u++)
    {
        // initialize temporary array "z" with zero-order B-splines
        for (int j = 0; j <= ORDER; j++)
        {
            z[j] = ((i + j == iknot) ? (1) : (0));
        }
        
        // compose higher-order B-splines
        for (int k = 1; k <= ORDER; k++)
        {
            for (int j = 0; j <= ORDER - k; j++)
            {
                double2 res = 0;
                
                if (t[i+j+k].x != t[i+j].x || t[i+j+k].y != t[i+j].y)
                    res = res + cmul(z[j], cdiv(x[u] - t[i+j], t[i+j+k] - t[i+j]));
                
                if (t[i+j+k+1].x != t[i+j+1].x || t[i+j+k+1].y != t[i+j+1].y)
                    res = res + cmul(z[j+1], cdiv(t[i+j+k+1] - x[u], t[i+j+k+1] - t[i+j+1]));
                
                z[j] = res;
            }
        }
        
        // save result
        y[u] = z[0];
    }
}

// ------------------------------------------------------------------------- //
// Radial integral functions.                                                //
// ------------------------------------------------------------------------- //

double damp (double2 b, double2 a, double R)
{
    // compute hyperradius
    double r = hypot(b.x, a.x);
    
    // if sufficiently far, return clean zero
    if (r > R)
        return 0.;
    
    // else damp using tanh(x) distribution
    return tanh(0.125 * (R - r));
}

kernel void R_integral
(
// R-integral four-index
    const int4 idx,
// multipole moment (angular momentum transfer)
    const int lambda,
// knot sequence
    constant double2 const * const t,
// inner Gauss-Legendre quadrature rule
    constant double const * const xIn0,
    constant double const * const wIn0,
// outer Gauss-Legendre quadrature rule
    constant double const * const xOut0,
    constant double const * const wOut0,
// output array for knot contributions
    global double2 * const restrict R_out
)
{
    // get knot
    const int iknot = get_global_id(0);
    
    // resulting contribution to the integral
    double2 R = 0;
    
    // compute only if the interval is non-empty
    if (t[iknot+1].x != t[iknot].x || t[iknot+1].y != t[iknot].y)
    {
        // get evaluation points count
        const int nIn  = ORDER + lambda + 1;
        const int nOut = ORDER + lambda + 10;
        
        // evaluation points and weights for interval
        //     x in (t[iknot],t[iknot+1])
        double2 xOut[NOUTMAX], wOut[NOUTMAX];
        for (int u = 0; u < nOut; u++)
        {
            xOut[u] = t[iknot] + (t[iknot+1] - t[iknot]) * xOut0[u];
            wOut[u] = (t[iknot+1] - t[iknot]) * wOut0[u];
        }
        
        // evaluate outer B-splines
        double2 values_i[NOUTMAX], values_j[NOUTMAX];
        B(t, idx.x, iknot, nOut, xOut, values_i);
        B(t, idx.y, iknot, nOut, xOut, values_j);
        
        // evaluate outer integral, fill output array
        for (int u = 0; u < nOut; u++)
        {
            // evaluation points and weights for interval
            //     y in (t[iknot],x)
            double2 xIn[NINMAX], wIn[NINMAX];
            for (int v = 0; v < nIn; v++)
            {
                xIn[v] = t[iknot] + (xOut[u] - t[iknot]) * xIn0[v];
                wIn[v] = (xOut[u] - t[iknot]) * wIn0[v];
            }
            
            // evaluate B-splines
            double2 values_k[NINMAX], values_l[NINMAX];
            B(t, idx.z, iknot, nIn, xIn, values_k);
            B(t, idx.w, iknot, nIn, xIn, values_l);
            
            // evaluate inner integral
            double2 evalIn = 0;
            for (int v = 0; v < nIn; v++)
                evalIn += cmul(wIn[v], cmul(values_k[v], cmul(values_l[v], cpow(cdiv(xIn[v],xOut[u]),lambda)))) * damp(xIn[v], 0, t[IKNOTMAX].x);
            
            // evaluate outer integral
            R += cmul(wOut[u], cmul(values_i[u], cmul(values_j[u], cdiv(evalIn,xOut[u])))) * damp(0, xOut[u], t[IKNOTMAX].x);
        }
    }
    
    // store the result
    R_out[iknot] = R;
}
