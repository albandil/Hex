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

#ifndef _MOMTS_H_
#define _MOMTS_H_

#include <cmath>
#include <vector>
#include <complex>

#include "arrays.h"
#include "bspline.h"
#include "gauss.h"
#include "specf.h"

class weightEdgeDamp
{
public:
    double operator() (Complex z) const
    {
        double R0 = Bspline::ECS().R0();
        
        // this will suppress function value from R0+1 onwards
        // which is useful for expanding (divergent) Ricatti-Bessel function
        return (z.imag() == 0.) ? (1+tanh(R0 - 5 - z.real()))/2 : 0.;
    }
};

class weightEndDamp
{
public:
    double operator() (Complex z) const
    {
        double Rmax = Bspline::ECS().Rmax();
        
        // whis will suppress function value at Rmax
        // which is useful for expanding everywhere-nonzero hydrogenic function
        return tanh(Rmax - z.real());
    }
};

/**
 * Compute derivative overlap of B-splines @f$ B_i @f$ and @f$ B_j @f$
 * over the knot "iknot", using Gauss-Legendre integration.
 * @param i      B-spline index.
 * @param j      B-spline index.
 * @param iknot  Interval index.
 */
Complex computeD_iknot(int i, int j, int iknot);

/**
 * Compute derivative overlap for B-splines @f$ B_i @f$ and @f$ B_j @f$.
 * @param i B-spline index.
 * @param j B-spline index.
 * @param maxknot Right-most knot of any integration.
 */
Complex computeD(int i, int j, int maxknot = Bspline::ECS().Nknot() - 1);

/**
 * Compute integral moment of coordinate power between the B-splines
 * @f$ B_i @f$ and @f$ B_j @f$
 * over the knot "iknot", using Gauss-Legendre integration.
 * @param a      Exponent.
 * @param i      B-spline index.
 * @param j      B-spline index.
 * @param iknot  Interval index.
 * @param R      Potential truncation point.
 */
Complex computeM_iknot(int a, int i, int j, int iknot, Complex R);

/**
 * Compute integral moment of coordinate power between the B-splines
 * @f$ B_i @f$ and @f$ B_j @f$
 * @param a Exponent.
 * @param i B-spline index.
 * @param j B-spline index.
 * @param maxknot Right-most knot of any integration.
 */
Complex computeM(int a, int i, int j, int maxknot = 0);

/**
 * Compute integral moment of degree "a" for every B-spline pair and every
 * interknot sub-interval. Store in 1-D array of shape
 * @code
 * [ Nspline × Nspline × Nintval ]
 * @endcode
 * Zero entries are stores as well, to allow for better caching.
 * 
 * @param a Moment degree.
 */
cArray computeMi(int a, int iknotmax = 0);

/** Compute P-overlaps
 * Compute overlap vector of B-splines vs. hydrogen Pnl function.
 * @param n Principal quantum number.
 * @param l Orbital quantum number.
 * @param weightf Weight function to multiply every value of the hydrogenic function.
 *                It is expected to have the "double operator() (Complex z)" interface,
 *                where the sent value is the complex coordinate.
 */
template <class Functor> cArray overlapP(int n, int l, Functor weightf)
{
    cArray res(Bspline::ECS().Nspline());
    
    // per interval
    int points = 20;
    
    // evaluated B-spline and hydrogenic functions (auxiliary variables)
    cArray evalB(points);
    cArray evalP(points);
    
    // for all knots
    for (int iknot = 0; iknot < Bspline::ECS().Nknot() - 1; iknot++)
    {
        // skip zero length intervals
        if (Bspline::ECS().t(iknot) == Bspline::ECS().t(iknot+1))
            continue;
        
        // which points are to be used here?
        cArray xs = p_points(points, Bspline::ECS().t(iknot), Bspline::ECS().t(iknot+1));
        cArray ws = p_weights(points, Bspline::ECS().t(iknot), Bspline::ECS().t(iknot+1));
        
        // evaluate the hydrogenic function
        std::transform(
            xs.begin(), xs.end(), evalP.begin(),
            [ = ](Complex x) -> Complex {
                gsl_sf_result R;
                if (gsl_sf_hydrogenicR_e(n, l, 1., x.real(), &R) == GSL_EUNDRFLW)
                    return 0.;
                else
                    return weightf(x) * x * R.val;
            }
        );
        
        // for all relevant B-splines
        for (int ispline = std::max(iknot-Bspline::ECS().order(),0); ispline < Bspline::ECS().Nspline() and ispline <= iknot; ispline++)
        {
            // evaluate the B-spline
            Bspline::ECS().B(ispline, iknot, points, xs.data(), evalB.data());
            
            // sum with weights
            Complex sum = 0.;
            for (int ipoint = 0; ipoint < points; ipoint++)
                sum += ws[ipoint] * evalP[ipoint] * evalB[ipoint];
            
            // store the overlap
            res[ispline] += sum;
        }
    }
    
    return res;
}

/** Compute j-overlaps
 * Compute overlap integrals for Riccati-Bessel function.
 * @param maxL2 Maximal degree of the Riccati-Bessel function.
 * @param vk Array containing linear momenta.
 * @param weightf Weight function to multiply every value of the Bessel function.
 *                It is expected to have the "double operator() (Complex z)" interface,
 *                where the sent value is the complex coordinate.
 * @return Array of shape [vk.size() × (maxell + 1) × Nspline] in column-major format.
 */
template <class Functor> cArray overlapj(int maxell, std::vector<double> vk, Functor weightf)
{
    // shorthands
    int Nenergy = vk.size();
    int Nspline = Bspline::ECS().Nspline();
    int Nknot = Bspline::ECS().Nknot();
    int order = Bspline::ECS().order();
    
    // reserve space for the output array
    size_t size = Nspline * Nenergy * (maxell + 1);
    cArray res(size);
    
    // per interval
    int points = 20;
        
    // for all knots
    # pragma omp parallel for
    for (int iknot = 0; iknot < Nknot - 1; iknot++)
    {
        // skip zero length intervals
        if (Bspline::ECS().t(iknot) == Bspline::ECS().t(iknot+1))
            continue;
        
        // which points are to be used here?
        cArray xs = p_points(points, Bspline::ECS().t(iknot), Bspline::ECS().t(iknot+1));
        cArray ws = p_weights(points, Bspline::ECS().t(iknot), Bspline::ECS().t(iknot+1));
        
        // evaluate relevant B-splines on this knot
        cArrays evalB(Nspline);
        for (int ispline = std::max(iknot-order,0); ispline < Nspline and ispline <= iknot; ispline++)
        {
            evalB[ispline] = cArray(points);
            Bspline::ECS().B(ispline, iknot, points, xs.data(), evalB[ispline].data());
        }
        
        // for all linear momenta (= energies)
        for (int ie = 0; ie < Nenergy; ie++)
        {
            // evaluate the Riccati-Bessel function for this knot and energy and for all angular momenta
            cArrays evalj(points);
            for (int ipoint = 0; ipoint < points; ipoint++)
                evalj[ipoint] = weightf(xs[ipoint]) * ric_jv(maxell, vk[ie] * xs[ipoint]);
            
            // for all angular momenta
            for (int l = 0; l <= maxell; l++)
            {
                // for all relevant B-splines
                for (int ispline = std::max(iknot-order,0); ispline < Nspline and ispline <= iknot; ispline++)
                {
                    // sum with weights
                    Complex sum = 0.;
                    for (int ipoint = 0; ipoint < points; ipoint++)
                        sum += ws[ipoint] * evalj[ipoint][l] * evalB[ispline][ipoint];
                    
                    // store the overlap; keep the shape Nmomenta × Nspline × (maxl+1)
                    res[(ie * (maxell + 1) + l) * Nspline + ispline] += sum;
                }
            }
        }
    }
    
    return res;
}

#endif
