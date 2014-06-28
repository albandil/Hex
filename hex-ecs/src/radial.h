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

#ifndef HEX_MOMENTS
#define HEX_MOMENTS

#include <cmath>
#include <vector>

#include "arrays.h"
#include "bspline.h"
#include "complex.h"
#include "gauss.h"
#include "input.h"
#include "parallel.h"
#include "special.h"
#include "matrix.h"

class weightEdgeDamp
{
    public:
        
        // constructor
        weightEdgeDamp(Bspline const & bspline)
            : bspline_(bspline) {}
        
        // weight function
        double operator() (Complex z) const
        {
            double R0 = bspline_.R0();
            
            // this will suppress function value from R0+1 onwards
            // which is useful for expanding (divergent) Ricatti-Bessel function
            return (z.imag() == 0.) ? (1+tanh(R0 - 5 - z.real()))/2 : 0.;
        }
    
    private:
        
        /// B-spline environment
        Bspline const & bspline_;
};

class weightEndDamp
{
    public:
        
        // constructor
        weightEndDamp(Bspline const & bspline)
            : bspline_(bspline) {}
        
        // weight function
        double operator() (Complex z) const
        {
            double Rmax = bspline_.Rmax();
            
            // whis will suppress function value at Rmax
            // which is useful for expanding everywhere-nonzero hydrogenic function
            return tanh(Rmax - z.real());
        }
        
    private:
        
        /// B-spline environment
        Bspline const & bspline_;
};

class RadialIntegrals
{
    public:
        
        // constructor
        RadialIntegrals(Bspline const & bspline)
            : bspline_(bspline), g_(bspline), D_(bspline.Nspline()), S_(bspline.Nspline()),
              Mm1_(bspline.Nspline()), Mm1_tr_(bspline.Nspline()), Mm2_(bspline.Nspline()) {}
        
        // public callable members
        void setupOneElectronIntegrals();
        void setupTwoElectronIntegrals(Parallel const & par, CommandLine const & cmd, Array<bool> const & lambdas);
        
        /**
         * Compute derivative overlap of B-splines @f$ B_i @f$ and @f$ B_j @f$
         * over the knot "iknot", using Gauss-Legendre integration.
         * @param i      B-spline index.
         * @param j      B-spline index.
         * @param iknot  Interval index.
         */
        Complex computeD_iknot(int i, int j, int iknot) const;
        
        /**
         * Compute derivative overlap for B-splines @f$ B_i @f$ and @f$ B_j @f$.
         * @param i B-spline index.
         * @param j B-spline index.
         * @param maxknot Right-most knot of any integration.
         */
        Complex computeD(int i, int j, int maxknot = -1) const;
        
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
        Complex computeM_iknot(int a, int i, int j, int iknot, Complex R) const;
        
        /**
         * Compute integral moment of coordinate power between the B-splines
         * @f$ B_i @f$ and @f$ B_j @f$
         * @param a Exponent.
         * @param i B-spline index.
         * @param j B-spline index.
         * @param maxknot Right-most knot of any integration.
         */
        Complex computeM(int a, int i, int j, int maxknot = 0) const;
        
        /**
         * Compute logarithms of integral moment of degree "a" for every B-spline pair and every
         * interknot sub-interval. Store in 1-D array of shape
         * @code
         * [ Nspline × Nspline × Nintval ]
         * @endcode
         *  
         * @param a Moment degree.
         * @param iknotmax Index of knot that terminates the integration range.
         */
        cArray computeMi(int a, int iknotmax = 0) const;
        
        rArray computeScale (int a, int iknotmax = 0) const;
        
        void M_integrand (int n, Complex *in, Complex *out, int i, int j, int a, int iknot, int iknotmax, double& logscale) const;
        
        /**
         * Compute the two-electron (Slater-type) four-B-spline integral.
         * 
         * @param lambda Multipole degree.
         * @param i First (x-dependent) B-spline index.
         * @param j Second (y-dependent) B-spline index.
         * @param k Third (x-dependent) B-spline index.
         * @param l Fourth (y-dependent) B-spline index.
         * @param Mtr_L Logarithms of R₀-truncated partial moments.
         * @param Mtr_mLm1 Logarithms of R₀-truncated partial moments.
         * 
         * Given the R-type integral symmetry, following calls will produce identical results:
         * 
         * @code
         * computeR(lambda, i, j, k, l);
         * computeR(lambda, j, i, l, k);
         * computeR(lambda, k, j, i, l);
         * computeR(lambda, i, l, k, j);
         * computeR(lambda, k, l, i, j);
         * @endcode
         * 
         * @return Value of Rtr.
         */
        Complex computeR (
            int lambda,
            int i, int j, int k, int l,
            cArray const & Mtr_L,
            cArray const & Mtr_mLm1
        ) const;
        
        Complex computeRdiag (int L, int a, int b, int c, int d, int iknot, int iknotmax) const;
        Complex computeRtri (int L, int k, int l, int m, int n, int iknot, int iknotmax) const;
        void R_inner_integrand (int n, Complex* in, Complex* out, int i, int j, int L, int iknot, int iknotmax, Complex x) const;
        void R_outer_integrand (int n, Complex* in, Complex* out, int i, int j, int k, int l, int L, int iknot, int iknotmax) const;
        
        void allSymmetries (
            int i, int j, int k, int l,
            Complex Rijkl_tr,
            NumberArray<long> & R_tr_i,
            NumberArray<long> & R_tr_j,
            NumberArray<Complex> & R_tr_v
        ) const;
        
        /** Compute P-overlaps
         * Compute overlap vector of B-splines vs. hydrogen Pnl function.
         * @param n Principal quantum number.
         * @param l Orbital quantum number.
         * @param weightf Weight function to multiply every value of the hydrogenic function.
         *                It is expected to have the "double operator() (Complex z)" interface,
         *                where the sent value is the complex coordinate.
         */
        template <class Functor> cArray overlapP (int n, int l, Functor weightf) const
        {
            cArray res(bspline_.Nspline());
            
            // per interval
            int points = 20;
            
            // evaluated B-spline and hydrogenic functions (auxiliary variables)
            cArray evalB(points);
            cArray evalP(points);
            
            // for all knots
            for (int iknot = 0; iknot < bspline_.Nknot() - 1; iknot++)
            {
                // skip zero length intervals
                if (bspline_.t(iknot) == bspline_.t(iknot+1))
                    continue;
                
                // which points are to be used here?
                cArray xs = g_.p_points(points, bspline_.t(iknot), bspline_.t(iknot+1));
                cArray ws = g_.p_weights(points, bspline_.t(iknot), bspline_.t(iknot+1));
                
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
                for (int ispline = std::max(iknot-bspline_.order(),0); ispline < bspline_.Nspline() and ispline <= iknot; ispline++)
                {
                    // evaluate the B-spline
                    bspline_.B(ispline, iknot, points, xs.data(), evalB.data());
                    
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

        /**
         * @brief Compute j-overlaps
         * 
         * Compute B-spline overlap integrals for Riccati-Bessel function.
         * 
         * @param maxell Maximal degree of the Riccati-Bessel function.
         * @param vk Array containing linear momenta.
         * @param weightf Weight function to multiply every value of the Bessel function.
         *                It is expected to have the "double operator() (Complex z)" interface,
         *                where the sent value is the complex coordinate.
         * @return Array of shape [vk.size() × (maxell + 1) × Nspline] in column-major format.
         */
        template <class Functor> cArray overlapj (int maxell, const rArrayView vk, Functor weightf) const
        {
            // shorthands
            int Nenergy = vk.size();
            int Nspline = bspline_.Nspline();
            int Nknot = bspline_.Nknot();
            int order = bspline_.order();
            
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
                if (bspline_.t(iknot) == bspline_.t(iknot+1))
                    continue;
                
                // which points are to be used here?
                cArray xs = g_.p_points(points, bspline_.t(iknot), bspline_.t(iknot+1));
                cArray ws = g_.p_weights(points, bspline_.t(iknot), bspline_.t(iknot+1));
                
                // evaluate relevant B-splines on this knot
                cArrays evalB(Nspline);
                for (int ispline = std::max(iknot-order,0); ispline < Nspline and ispline <= iknot; ispline++)
                {
                    evalB[ispline] = cArray(points);
                    bspline_.B(ispline, iknot, points, xs.data(), evalB[ispline].data());
                }
                
                // for all linear momenta (= energies)
                for (int ie = 0; ie < Nenergy; ie++)
                {
                    // evaluate the Riccati-Bessel function for this knot and energy and for all angular momenta
                    cArrays evalj(points);
                    for (int ipoint = 0; ipoint < points; ipoint++)
                        evalj[ipoint] = weightf(xs[ipoint]) * special::ric_jv(maxell, vk[ie] * xs[ipoint]);
                    
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
        
        Bspline const & bspline() const { return bspline_; }
        
        SymDiaMatrix const & D() const { return D_; }
        SymDiaMatrix const & S() const { return S_; }
        SymDiaMatrix const & Mm1() const { return Mm1_; }
        SymDiaMatrix const & Mm1_tr() const { return Mm1_tr_; }
        SymDiaMatrix const & Mm2() const { return Mm2_; }
        
        SymDiaMatrix const & R_tr_dia(unsigned i) const
        {
            assert(i < R_tr_dia_.size());
            return R_tr_dia_[i];
        }
        
        size_t maxlambda() const { return R_tr_dia_.size() - 1; }
        
    private:
        
        // B-spline environment
        Bspline const & bspline_;
        
        // Gauss-Legendre integrator
        GaussLegendre g_;
        
        //
        // matrices
        //
        
        SymDiaMatrix D_, S_, Mm1_, Mm1_tr_, Mm2_;
        Array<SymDiaMatrix> R_tr_dia_;
};

#endif
