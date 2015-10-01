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

#ifndef HEX_MOMENTS
#define HEX_MOMENTS

#include <cmath>
#include <vector>

#include "arrays.h"
#include "bspline.h"
#include "gauss.h"
#include "io.h"
#include "parallel.h"
#include "special.h"
#include "matrix.h"

#define EXPANSION_QUADRATURE_POINTS 20

/**
 * Potential suppressing factor. 
 * @param y Radial coordinate of some electron.
 * @param x Radial coordinate of the other electron.
 * @param R Truncation radius.
 */
inline double damp (Complex y, Complex x, Complex R)
{
    // compute hyperradius
    double r = std::hypot(x.real(), y.real());
    
    // if sufficiently far, return clean zero
    if (r > R.real())
        return 0.;
    
    // else damp using tanh(x) distribution
    return std::tanh(0.125 * (R.real() - r));
}

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
            return (z.imag() == 0.) ? (1 + std::tanh(R0 - 5 - z.real()))/2 : 0.;
        }
    
    private:
        
        /// B-spline environment
        Bspline const & bspline_;
};

class weightEndDamp
{
    public:
        
        // constructor
        weightEndDamp (Bspline const & bspline)
            : bspline_(bspline) {}
        
        // weight function
        double operator() (Complex z) const
        {
            double Rmax = bspline_.Rmax();
            
            // whis will suppress function value at Rmax
            // which is useful for expanding everywhere-nonzero hydrogenic function
            return std::tanh(Rmax - z.real());
        }
        
    private:
        
        /// B-spline environment
        Bspline const & bspline_;
};

class RadialIntegrals
{
    public:
        
        // constructor
        RadialIntegrals (const Bspline& bspline_atom, const Bspline& bspline_proj, int Nlambdas);
        
        // public callable members
        void setupOneElectronIntegrals (Parallel const & par, CommandLine const & cmd);
        void setupTwoElectronIntegrals (Parallel const & par, CommandLine const & cmd);
        
        // verbosity control
        void verbose (bool v) { verbose_ = v; }
        
        /**
         * @brief Partial derivative overlap.
         * 
         * Compute derivative overlap of B-splines @f$ B_i @f$ and @f$ B_j @f$
         * over the knot "iknot", using Gauss-Legendre integration.
         * @param i      B-spline index.
         * @param j      B-spline index.
         * @param iknot  Interval index.
         */
        Complex computeD_iknot
        (
            Bspline const & bspline, GaussLegendre const & g,
            int i, int j, int iknot
        ) const;
        
        /**
         * @brief Derivative overlap.
         * 
         * Compute derivative overlap for B-splines @f$ B_i @f$ and @f$ B_j @f$.
         * @param i B-spline index.
         * @param j B-spline index.
         * @param maxknot Right-most knot of any integration.
         */
        Complex computeD
        (
            Bspline const & bspline, GaussLegendre const & g,
            int i, int j, int maxknot = -1
        ) const;
        
        /**
         * @brief Partial integral moment.
         * 
         * Compute integral moment of coordinate power between the B-splines
         * @f$ B_i @f$ and @f$ B_j @f$
         * over the knot "iknot", using Gauss-Legendre integration.
         * @param a      Exponent.
         * @param i      B-spline index.
         * @param j      B-spline index.
         * @param iknot  Interval index.
         * @param R      Potential truncation point.
         */
        Complex computeM_iknot
        (
            Bspline const & bspline, GaussLegendre const & g,
            int a, int i, int j,
            int iknot, Complex R, double scale
        ) const;
        
        /**
         * @brief Integral moments.
         * 
         * Compute integral moment of coordinate power between the B-splines
         * @f$ B_i @f$ and @f$ B_j @f$
         * @param a Exponent.
         * @param i B-spline index.
         * @param j B-spline index.
         * @param maxknot Right-most knot of any integration.
         */
        Complex computeM
        (
            Bspline const & bspline, GaussLegendre const & g,
            int a, int i, int j,
            int maxknot = -1, bool scale = false
        ) const;
        
        /**
         * @brief Partial integral moments.
         * 
         * Compute logarithms of integral moment of degree "a" for every B-spline pair and every
         * interknot sub-interval. Store in 1-D array of shape
         * @code
         * [ Nspline × Nspline × Nintval ]
         * @endcode
         *  
         * @param a Moment degree.
         * @param iknotmax Index of knot that terminates the integration range.
         */
        cArray computeMi
        (
            Bspline const & bspline, GaussLegendre const & g,
            int a, int iknotmax = 0
        ) const;
        
        /**
         * @brief Integrand used in calculation of integral moments.
         * 
         * This function is used by @ref computeMi to calculate partial integral moments
         * of a coordinate between two B-splines.
         */
        void Mi_integrand
        (
            int n, Complex *in, Complex *out,
            Bspline const & bspline_ij, int i, int j,
            int a, int iknot, int iknotmax
        ) const;
        
        cArray diagonalR (int lambda) const;
        
        /**
         * @brief Two-electron integral.
         * 
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
        Complex computeR
        (
            int lambda,
            int i, int j, int k, int l,
            bool simple = false
        ) const;
        
        /**
         * @brief Diagonal contribution to R-integral.
         * 
         * Calculates integral
         * @f[
         *     R_{ijkl}^\lambda = \int\limits_{t_i}^{t_{i+1}} \int\limits_{t_i}^{t_{i+1}}
         *     B_a(r_1) B_b(r_2) \frac{r_<^\lambda}{r_>^{\lambda+1}} B_c(r_1) B_d(r_2)
         *     \mathrm{d}r_1 \mathrm{d}r_2 \,.
         * @f]
         * 
         * The integral is computed as a sum of two triangular integrals, see @ref computeRtri.
         * 
         * @param L The multipole (@f$ \lambda @f$).
         * @param a First B-spline index.
         * @param b Second B-spline index.
         * @param c Third B-spline index.
         * @param d Fourth B-spline index.
         * @param iknot Integration interval (t[iknot] ... t[iknot+1]).
         * @param iknotmax Truncation knot for use in damping factor.
         */
        Complex computeRdiag (int L, int a, int b, int c, int d, int iknot, int iknotmax) const;
        
        /**
         * @brief Triangular R-integral.
         * 
         * Calculates integral
         * @f[
         *     R_{klmn}^{\lambda\triangle} = \int\limits_{t_i}^{t_{i+1}}
         *     B_k(r_1) B_l(r_1) \frac{1}{r_1} \left( \int\limits_{t_i}^{r_1}
         *     \left(\frac{r_2}{r_1}\right)^\lambda B_m(r_2) B_n(r_2)
         *     \mathrm{d}r_2 \right) \mathrm{d}r_1 \,,
         * @f]
         * 
         * which is used to calculate a diagonal contribution to the full R-integral,
         * see @ref computeRdiag. Inner integrand is a polynomial, so
         * a simple Gauss-Legendre quadrature is used (with a sufficient order). The outer
         * integral is polynomial only when the integration starts from zero (@f$ t_i = 0 @f$).
         * In such case it has degree equal to the combined order of the B-splines.
         * It is assumed that also elsewhere the integrand can be well approximated
         * by a polynomial of the same degree.
         * 
         * Uses functions @ref R_outer_integrand and @ref R_inner_integrand.
         * 
         * @param L The multipole (@f$ \lambda @f$).
         * @param a First B-spline index.
         * @param b Second B-spline index.
         * @param c Third B-spline index.
         * @param d Fourth B-spline index.
         * @param iknot Integration interval (t[iknot] ... t[iknot+1]).
         * @param iknotmax Truncation knot for use in damping factor.
         */
        Complex computeRtri
        (
            int L,
            int k, int l,
            int m, int n,
            int iknot, int iknotmax
        ) const;
        
        /**
         * @brief Evaluates inner integrand of the triangular R-integral.
         * 
         * Calculates for a set of coordinates @f$ r_2 @f$ (and given coordinate @f$ r_1 @f$) the expression
         * @f[
         *     g_{ij}^\lambda(r_1;r_2) = B_i(r_2) B_j(r_2) \left(\frac{r_2}{r_1}\right)^\lambda \xi(r_2) \,,
         * @f]
         * where @f$ \xi(r) @f$ is the damping function; see @ref damp.
         * 
         * @param n Number of evaluation points.
         * @param in Pointer to input array of (complex) evaluation points.
         * @param out Pointer to output array of (complex) evaluations.
         * @param i First B-spline index.
         * @param j Second B-spline index.
         * @param L Multipole (@f$ \lambda @f$).
         * @param iknot Integration interval.
         * @param iknotmax Damping parameter.
         * @param x Fixed value of @f$ r_1 @f$.
         */
        void R_inner_integrand
        (
            int n, Complex* in, Complex* out,
            int i, int j,
            int L, int iknot, int iknotmax, Complex x
        ) const;
        
        /**
         * @brief Evaluates outer integrand of the triangular R-integral.
         * 
         * Calculates for a set of coordinates @f$ r_1 @f$ the expression
         * @f[
         *     f_{ijkl}^\lambda(r_1) = B_i(r_1) B_j(r_1) \frac{1}{r_1} \xi(r_1) \int g_{kl}^\lambda(r_1;r_2) \mathrm{d}r_2 \,,
         * @f]
         * where @f$ \xi(r) @f$ is the damping function (see @ref damp) and @f$ g_{kl}^\lambda(x,y) @f$
         * is the inner integrand (see @ref R_inner_integrand).
         * 
         * @param n Number of evaluation points.
         * @param in Pointer to input array of (complex) evaluation points.
         * @param out Pointer to output array of (complex) evaluations.
         * @param i First B-spline index.
         * @param j Second B-spline index.
         * @param L Multipole (@f$ \lambda @f$).
         * @param iknot Integration interval.
         * @param iknotmax Damping parameter.
         * @param x Fixed value of @f$ r_1 @f$.
         */
        void R_outer_integrand
        (
            int n, Complex* in, Complex* out,
            int i, int j,
            int k, int l,
            int L, int iknot, int iknotmax
        ) const;
        
        /** 
         * @brief Compute P-overlaps
         * 
         * Compute overlap vector of B-splines vs. hydrogen Pnl function.
         * 
         * @param bspline B-spline basis.
         * @param g Gauss-Legendre integrator adapted to the B-spline basis.
         * @param n Principal quantum number.
         * @param l Orbital quantum number.
         * @param weightf Weight function to multiply every value of the hydrogenic function.
         *                It is expected to have the "double operator() (Complex z)" interface,
         *                where the sent value is the complex coordinate.
         */
        template <class Functor> cArray overlapP (Bspline const & bspline, GaussLegendre const & g, int n, int l, Functor weightf) const
        {
            // result
            cArray res (bspline.Nspline());
            
            // per interval
            int points = EXPANSION_QUADRATURE_POINTS;
            
            // evaluated B-spline and hydrogenic functions
            cArray evalB (points), evalP (points);
            
            // quadrature nodes and weights
            cArray xs (points), ws (points);
            
            // for all knots
            for (int iknot = 0; iknot < bspline.Nknot() - 1; iknot++)
            {
                // skip zero length intervals
                if (bspline.t(iknot) == bspline.t(iknot+1))
                    continue;
                
                // which points are to be used here?
                g.scaled_nodes_and_weights(points, bspline.t(iknot), bspline.t(iknot+1), &xs[0], &ws[0]);
                
                // evaluate the hydrogenic function
                std::transform
                (
                    xs.begin(), xs.end(), evalP.begin(),
                    [ = ](Complex x) -> Complex
                    {
                        gsl_sf_result R;
                        if (gsl_sf_hydrogenicR_e(n, l, 1., x.real(), &R) == GSL_EUNDRFLW)
                            return 0.;
                        else
                            return weightf(x) * x * R.val;
                    }
                );
                
                // for all relevant B-splines
                for (int ispline = std::max(iknot-bspline.order(),0); ispline < bspline.Nspline() and ispline <= iknot; ispline++)
                {
                    // evaluate the B-spline
                    bspline.B(ispline, iknot, points, xs.data(), evalB.data());
                    
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
        template <class Functor> cArray overlapj (Bspline const & bspline, GaussLegendre const & g, int maxell, const rArrayView vk, Functor weightf) const
        {
            // shorthands
            int Nenergy = vk.size();
            int Nspline = bspline.Nspline();
            int Nknot = bspline.Nknot();
            int order = bspline.order();
            
            // reserve space for the output array
            std::size_t size = Nspline * Nenergy * (maxell + 1);
            cArray res (size);
            
            // per interval
            int points = EXPANSION_QUADRATURE_POINTS;
            
            // quadrature weights and nodes
            cArray xs (points), ws (points);
            
            // for all knots
            # pragma omp parallel for firstprivate (xs,ws)
            for (int iknot = 0; iknot < Nknot - 1; iknot++)
            {
                // skip zero length intervals
                if (bspline.t(iknot) == bspline.t(iknot+1))
                    continue;
                
                // which points are to be used here?
                g.scaled_nodes_and_weights(points, bspline.t(iknot), bspline.t(iknot+1), &xs[0], &ws[0]);
                
                // evaluate relevant B-splines on this knot
                cArrays evalB(Nspline);
                for (int ispline = std::max(iknot-order,0); ispline < Nspline and ispline <= iknot; ispline++)
                {
                    evalB[ispline] = cArray(points);
                    bspline.B(ispline, iknot, points, xs.data(), evalB[ispline].data());
                }
                
                // for all linear momenta (= energies)
                for (int ie = 0; ie < Nenergy; ie++)
                {
                    // evaluate the Riccati-Bessel function for this knot and energy and for all angular momenta
                    cArrays evalj(points);
                    for (int ipoint = 0; ipoint < points; ipoint++)
                    {
                        // compute the damping factor
                        double damp = weightf(xs[ipoint]);
                        
                        // if the factor is numerical zero, do not evaluate the function, just allocate zeros
                        if (damp == 0)
                        {
                            evalj[ipoint] = cArray(maxell + 1);
                            continue;
                        }
                        
                        // evaluate all Riccati-Bessel functions in point
                        evalj[ipoint] = damp * special::ric_jv(maxell, vk[ie] * xs[ipoint]);
                        
                        // clear all possible NaN entries (these may occur for far radii, where should be zero)
                        for (int l = 0; l <= maxell; l++) if (not Complex_finite(evalj[ipoint][l]))
                            evalj[ipoint][l] = 0.;
                    }
                    
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
        
        /// Return reference to the B-spline object.
        Bspline const & bspline_atom () const { return bspline_atom_; }
        Bspline const & bspline_proj () const { return bspline_proj_; }
        
        /// Return the Gauss-Legendre integrator object.
        GaussLegendre const & gaussleg_atom () const { return g_atom_; }
        GaussLegendre const & gaussleg_proj () const { return g_proj_; }
        
        /// Return reference to the precomputed derivative overlap matrix.
        SymBandMatrix<Complex> const & D_atom () const { return D_atom_; }
        SymBandMatrix<Complex> const & D_proj () const { return D_proj_; }
        
        /// Return reference to the precomputed overlap matrix.
        SymBandMatrix<Complex> const & S_atom () const { return S_atom_; }
        SymBandMatrix<Complex> const & S_proj () const { return S_proj_; }
        
        /// Return reference to the precomputed integral moment matrix of order -1.
        SymBandMatrix<Complex> const & Mm1_atom () const { return Mm1_atom_; }
        SymBandMatrix<Complex> const & Mm1_proj () const { return Mm1_proj_; }
        
        /// Return reference to the precomputed integral moment matrix of order -1, truncated at the end of the real grid.
        SymBandMatrix<Complex> const & Mm1_tr_atom () const { return Mm1_tr_atom_; }
        SymBandMatrix<Complex> const & Mm1_tr_proj () const { return Mm1_tr_proj_; }
        
        /// Return reference to the precomputed integral moment matrix of order -2.
        SymBandMatrix<Complex> const & Mm2_atom () const { return Mm2_atom_; }
        SymBandMatrix<Complex> const & Mm2_proj () const { return Mm2_proj_; }
        
        /// Return reference to the precomputed array of diagonal contributions to two-electron integrals.
        cArray const & R_tr_dia_diag (unsigned i) const
        {
            assert(i < R_tr_dia_diag_.size());
            return R_tr_dia_diag_[i];
        }
        
        /// Return reference to the precomputed matrix of two-electron integrals for given multipole.
        BlockSymBandMatrix<Complex> const & R_tr_dia (unsigned i) const
        {
            assert(i < R_tr_dia_.size());
            return R_tr_dia_[i];
        }
        
        /// Return reference to precomputed full (scaled) integral moments of order L.
        SymBandMatrix<Complex> const & Mtr_L_atom (int L) const { return Mtr_L_atom_[L]; }
        SymBandMatrix<Complex> const & Mtr_L_proj (int L) const { return Mtr_L_proj_[L]; }
        
        /// Return reference to precomputed full (scaled) integral moments of order -L-1.
        SymBandMatrix<Complex> const & Mtr_mLm1_atom (int L) const { return Mtr_mLm1_atom_[L]; }
        SymBandMatrix<Complex> const & Mtr_mLm1_proj (int L) const { return Mtr_mLm1_proj_[L]; }
        
        /// Return view of precomputed partial integral moments of order L.
        cArrayView Mitr_L_atom (int L) const
        {
            if (L < 0)
                return Mitr_L_atom_;
            
            std::size_t mi_size = bspline_atom_.Nspline() * (2 * bspline_atom_.order() + 1) * (bspline_atom_.order() + 1);
            
            return cArrayView
            (
                Mitr_L_atom_,
                L * mi_size,
                mi_size
            );
        }
        cArrayView Mitr_L_proj (int L) const
        {
            if (L < 0)
                return Mitr_L_proj_;
            
            std::size_t mi_size = bspline_proj_.Nspline() * (2 * bspline_proj_.order() + 1) * (bspline_proj_.order() + 1);
            
            return cArrayView
            (
                Mitr_L_proj_,
                L * mi_size,
                mi_size
            );
        }
        
        /// Return view of precomputed partial integral moments of order -L-1.
        cArrayView Mitr_mLm1_atom (int L) const
        {
            if (L < 0)
                return Mitr_mLm1_atom_;
            
            std::size_t mi_size = bspline_atom_.Nspline() * (2 * bspline_atom_.order() + 1) * (bspline_atom_.order() + 1);
            
            return cArrayView
            (
                Mitr_mLm1_atom_,
                L * mi_size,
                mi_size
            );
        }
        cArrayView Mitr_mLm1_proj (int L) const
        {
            if (L < 0)
                return Mitr_mLm1_proj_;
            
            std::size_t mi_size = bspline_proj_.Nspline() * (2 * bspline_proj_.order() + 1) * (bspline_proj_.order() + 1);
            
            return cArrayView
            (
                Mitr_mLm1_proj_,
                L * mi_size,
                mi_size
            );
        }
        
        /**
         * @brief Calculate particular sub-matrix of the radial integrals matrix.
         * 
         * Calculate particular sub-matrix of the radial integrals matrix (with block indices "i" and "k")
         * and return it in a form of a dense array (copying structure of the overlap matrix).
         */
        SymBandMatrix<Complex> calc_R_tr_dia_block (unsigned lambda, int i, int k, bool simple = false) const;
        
        /**
         * @brief Multiply vector by matrix of two-electron integrals.
         * 
         * This routine will multiply several source vectors by the matrix of two-electron
         * integrals for given multipole 'lambda'. The matrix elements are
         * computed anew and applied directly to the vectors to minimize memory
         * requirements. This method - instead of caching the whole integral
         * matrix in memory or on disk - is used in the 'lightweight' mode,
         * which can be requested by the command line option --lightweight.
         */
        void apply_R_matrix
        (
            unsigned lambda,
            Complex a, const cArrayView src,
            Complex b,       cArrayView dst,
            bool simple = false
        ) const;
        
        /// Return maximal multipole, for which there are precomputed two-electron integrals.
        int maxlambda () const { return Nlambdas_ - 1; }
        
    private:
        
        // B-spline environment
        Bspline const & bspline_atom_;
        Bspline const & bspline_proj_;
        
        // Gauss-Legendre integrator
        GaussLegendre g_atom_;
        GaussLegendre g_proj_;
        
        // knot that terminates multipole scaled region
        int lastscaled_;
        
        // one-electron moment and overlap matrices
        SymBandMatrix<Complex> D_atom_, S_atom_, Mm1_atom_, Mm1_tr_atom_, Mm2_atom_;
        SymBandMatrix<Complex> D_proj_, S_proj_, Mm1_proj_, Mm1_tr_proj_, Mm2_proj_;
        
        // one-electron full integral moments for various orders (used to calculate R-integrals)
        std::vector<SymBandMatrix<Complex>> Mtr_L_atom_, Mtr_mLm1_atom_;
        std::vector<SymBandMatrix<Complex>> Mtr_L_proj_, Mtr_mLm1_proj_;
        
        // partial one-electron integral moments for various orders (used to calculate R-integrals)
        cArray Mitr_L_atom_, Mitr_mLm1_atom_;
        cArray Mitr_L_proj_, Mitr_mLm1_proj_;
        
        // diagonal contributions to R_tr_dia
        cArrays R_tr_dia_diag_;
        
        // two-electron integral matrices
        Array<BlockSymBandMatrix<Complex>> R_tr_dia_;
        
        // verbose output
        bool verbose_;
        
        // number of multipole matrices
        int Nlambdas_;
        
        // projectile basis shift
        int proj_basis_shift_;
};

#endif
