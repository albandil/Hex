//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2016, Jakub Benda, Charles University in Prague                    //
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
#include <functional>
#include <vector>

#include "hex-arrays.h"
#include "hex-symbandmatrix.h"
#include "hex-special.h"

#include "bspline.h"
#include "gauss.h"
#include "io.h"
#include "parallel.h"

#define EXPANSION_QUADRATURE_POINTS 20

/**
 * Potential suppressing factor. 
 * @param y Radial coordinate of some electron.
 * @param x Radial coordinate of the other electron.
 * @param R Truncation radius.
 */
inline Real damp (Complex y, Complex x, Complex R)
{
    // compute hyperradius
    Real r = std::max(x.real(), y.real());
    
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
        Real operator() (Complex z) const
        {
            Real R0 = bspline_.R0();
            
            // this will suppress function value from R0-5 onwards
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
        Real operator() (Complex z) const
        {
            Real Rmax = bspline_.Rmax();
            
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
        RadialIntegrals
        (
            Bspline const & bspline_inner,
            Bspline const & bspline_full,
            int Nlambdas
        );
        
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
        
        /// Calculate inter-basis overlap between two B-splines.
        Complex computeS12
        (
            GaussLegendre const & g,
            Bspline const & bspline1,
            Bspline const & bspline2,
            int i, int j
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
            int iknot, Complex R, Real scale
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
         */
        cArray computeMi
        (
            Bspline const & bspline,
            GaussLegendre const & g,
            int a
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
         * @brief Compute B-spline overlaps of arbitrary function.
         * 
         * @param bspline B-spline basis.
         * @param g Gauss-Legendre integrator adapted to the B-spline basis.
         * @param funct Some one-dimensional function (Complex -> Complex).
         * @param weightf Weight function to multiply every value of the hydrogenic function (Complex -> Real).
         */
        cArray overlap (Bspline const & bspline, GaussLegendre const & g, std::function<Complex(Complex)> funct, std::function<Real(Complex)> weightf) const;
        
        /** 
         * @brief Compute P-overlaps
         * 
         * Compute overlap vector of B-splines vs. hydrogen Pnl function.
         * 
         * @param bspline B-spline basis.
         * @param g Gauss-Legendre integrator adapted to the B-spline basis.
         * @param n Principal quantum number.
         * @param l Orbital quantum number.
         * @param weightf Weight function to multiply every value of the hydrogenic function (Complex -> Real).
         */
        cArray overlapP (Bspline const & bspline, GaussLegendre const & g, int n, int l, std::function<Real(Complex)> weightf) const;

        /**
         * @brief Compute j-overlaps
         * 
         * Compute B-spline overlap integrals for Riccati-Bessel function.
         * 
         * @param maxell Maximal degree of the Riccati-Bessel function.
         * @param vk Array containing linear momenta.
         * @param weightf Weight function to multiply every value of the Bessel function (Complex -> Real).
         * @return Array of shape [vk.size() × (maxell + 1) × Nspline] in column-major format.
         */
        cArray overlapj (Bspline const & bspline, GaussLegendre const & g, int maxell, const rArrayView vk, std::function<Real(Complex)> weightf) const;
        
        /// Return reference to the B-spline object.
        Bspline const & bspline_inner () const { return bspline_inner_; }
        Bspline const & bspline_full  () const { return bspline_full_ ; }
        
        /// Return the Gauss-Legendre integrator object.
        GaussLegendre const & gaussleg_inner () const { return g_inner_; }
        GaussLegendre const & gaussleg_full  () const { return g_full_ ; }
        
        /// Return reference to the precomputed derivative overlap matrix.
        SymBandMatrix<Complex> const & D_inner () const { return D_inner_; }
        SymBandMatrix<Complex> const & D_full  () const { return D_full_; }
        
        /// Return reference to the precomputed overlap matrix.
        SymBandMatrix<Complex> const & S_inner () const { return S_inner_; }
        SymBandMatrix<Complex> const & S_full  () const { return S_full_ ; }
        
        /// Return reference to the inter-basis overlaps.
//         RowMatrix<Complex> const & S12 () const { return S12_; }
//         RowMatrix<Complex> const & S21 () const { return S21_; }
        
        /// Return reference to the precomputed integral moment matrix of order -1.
        SymBandMatrix<Complex> const & Mm1_inner () const { return Mm1_inner_; }
        SymBandMatrix<Complex> const & Mm1_full  () const { return Mm1_full_ ; }
        
        /// Return reference to the precomputed integral moment matrix of order -1, truncated at the end of the real grid.
        SymBandMatrix<Complex> const & Mm1_tr_inner () const { return Mm1_tr_inner_; }
        SymBandMatrix<Complex> const & Mm1_tr_full  () const { return Mm1_tr_full_ ; }
        
        /// Return reference to the precomputed integral moment matrix of order -2.
        SymBandMatrix<Complex> const & Mm2_inner () const { return Mm2_inner_; }
        SymBandMatrix<Complex> const & Mm2_full  () const { return Mm2_full_ ; }
        
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
        SymBandMatrix<Complex> const & Mtr_L_inner (int L) const { return Mtr_L_inner_[L]; }
        SymBandMatrix<Complex> const & Mtr_L_full  (int L) const { return Mtr_L_full_ [L]; }
        
        /// Return reference to precomputed full (scaled) integral moments of order -L-1.
        SymBandMatrix<Complex> const & Mtr_mLm1_inner (int L) const { return Mtr_mLm1_inner_[L]; }
        SymBandMatrix<Complex> const & Mtr_mLm1_full  (int L) const { return Mtr_mLm1_full_ [L]; }
        
        /// Return view of precomputed partial integral moments of order L.
        cArrayView Mitr_L_inner (int L) const
        {
            if (L < 0)
                return Mitr_L_inner_;
            
            std::size_t mi_size = bspline_inner_.Nspline() * (2 * bspline_inner_.order() + 1) * (bspline_inner_.order() + 1);
            
            return cArrayView
            (
                Mitr_L_inner_,
                L * mi_size,
                mi_size
            );
        }
        
        cArrayView Mitr_L_full (int L) const
        {
            if (L < 0)
                return Mitr_L_full_;
            
            std::size_t mi_size = bspline_full_.Nspline() * (2 * bspline_full_.order() + 1) * (bspline_full_.order() + 1);
            
            return cArrayView
            (
                Mitr_L_full_,
                L * mi_size,
                mi_size
            );
        }
        
        /// Return view of precomputed partial integral moments of order -L-1.
        cArrayView Mitr_mLm1_inner (int L) const
        {
            if (L < 0)
                return Mitr_mLm1_inner_;
            
            std::size_t mi_size = bspline_inner_.Nspline() * (2 * bspline_inner_.order() + 1) * (bspline_inner_.order() + 1);
            
            return cArrayView
            (
                Mitr_mLm1_inner_,
                L * mi_size,
                mi_size
            );
        }
        
        cArrayView Mitr_mLm1_full (int L) const
        {
            if (L < 0)
                return Mitr_mLm1_full_;
            
            std::size_t mi_size = bspline_full_.Nspline() * (2 * bspline_full_.order() + 1) * (bspline_full_.order() + 1);
            
            return cArrayView
            (
                Mitr_mLm1_full_,
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
        Bspline const & bspline_inner_;
        Bspline const & bspline_full_;
        
        // Gauss-Legendre integrator
        GaussLegendre g_inner_;
        GaussLegendre g_full_;
        
        // one-electron moment and overlap matrices
        SymBandMatrix<Complex> D_inner_, S_inner_, Mm1_inner_, Mm1_tr_inner_, Mm2_inner_;
        SymBandMatrix<Complex> D_full_, S_full_, Mm1_full_, Mm1_tr_full_, Mm2_full_;
        
        // inter-basis overlaps
        RowMatrix<Complex> S12_, S21_;
        
        // one-electron full integral moments for various orders (used to calculate R-integrals)
        std::vector<SymBandMatrix<Complex>> Mtr_L_inner_, Mtr_mLm1_inner_;
        std::vector<SymBandMatrix<Complex>> Mtr_L_full_, Mtr_mLm1_full_;
        
        // partial one-electron integral moments for various orders (used to calculate R-integrals)
        cArray Mitr_L_inner_, Mitr_mLm1_inner_;
        cArray Mitr_L_full_, Mitr_mLm1_full_;
        
        // diagonal contributions to R_tr_dia
        cArrays R_tr_dia_diag_;
        
        // two-electron integral matrices
        Array<BlockSymBandMatrix<Complex>> R_tr_dia_;
        
        // verbose output
        bool verbose_;
        
        // number of multipole matrices
        int Nlambdas_;
};

#endif
