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

#ifndef HEX_ECS_RADIAL_H
#define HEX_ECS_RADIAL_H

// --------------------------------------------------------------------------------- //

#include <cmath>
#include <functional>
#include <vector>

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"
#include "hex-symbandmatrix.h"
#include "hex-special.h"

// --------------------------------------------------------------------------------- //

#include "bspline.h"
#include "gauss.h"
#include "io.h"
#include "parallel.h"

// --------------------------------------------------------------------------------- //

// This is a hard-coded number of points for Gauss-Legendre quadrature of the radial functions.
#define EXPANSION_QUADRATURE_POINTS 20

// --------------------------------------------------------------------------------- //

// Complex -> Complex function.
typedef std::function<Complex(Complex)> CCFunction;

// Complex x Complex -> Complex function.
typedef std::function<Complex(Complex,Complex)> C2CFunction;

// --------------------------------------------------------------------------------- //

/**
 * @brief Class that calculates and manages the radial integrals for given B-spline bases.
 * 
 * This class contains important routines for calculation of overlap matrices and matrix elements of various
 * functions for given B-spline bases. The calculation of the overlaps uses Gauss-Legendre quadrature,
 * that is to be specified in a form of a preallocated @ref GaussLegendre object.
 * 
 * The overlaps considered are one-electron integrals of various one-dimensional functions
 * (or given arbitrary function) and two-electron overlaps of the multipole elements (or given
 * adbitrary function).
 * 
 * There are two variants of the B-spline basis specified for every axis: the inner basis and the
 * full basis. These bases can be identical, when no channel reduction is used in the solution
 * of the scattering equations (typically for above-ionization energies). However, the inner and
 * outer bases can be different, when the channel reduction will take place.
 * The inner basis then contains real B-spline knots only and will be used for the full-channel
 * calculation of the inner problem. The full basis contains the inner basis knots but is generally much longer
 * and terminated by the ECS complex grid section. This basis is used for solution of energy-allowed channels.
 * The energy-forbidden channels are thus solved only in the inner basis, which is all right,
 * because the energy-forbidden channel functions exponentially decrease with the radial distance.
 */
class RadialIntegrals
{
    public:
        
        /**
         * @brief Constructor for idential bases.
         * 
         * This constructor will initialize the RadialIntegrals object assuming that the radial
         * B-spline bases for both axes are identical.
         */
        RadialIntegrals
        (
            Bspline const & bspline_inner,
            Bspline const & bspline_full,
            int Nlambdas
        );
        
        /**
         * @brief Constructor for distinct bases.
         * 
         * This constructor will initialize the RadialIntegrals object with different B-spline
         * bases for every axis.
         */
        RadialIntegrals
        (
            Bspline const & bspline_x_inner,
            Bspline const & bspline_x_full,
            Bspline const & bspline_y_inner,
            Bspline const & bspline_y_full,
            int Nlambdas
        );
        
        /**
         * @brief Calculate one-electron integral matrices.
         * 
         * This function will calculate the common one-electron overlap matrices
         * S, D, Mm1, Mm2. The parameters @c sharedscratch and @c iammaster are used
         * when writing the results to disk to avoid multiple processes writing
         * the same file.
         */
        void setupOneElectronIntegrals (bool sharedscratch = false, bool iammaster = true);
        
        /**
         * @brief Calculate one-electron integral matrices.
         * 
         * This function calls the other @ref setupOneElectronIntegrals with the
         * appropriate parameters.
         */
        void setupOneElectronIntegrals (Parallel const & par, CommandLine const & cmd);
        
        /**
         * @brief Calculate the two-electron integral matrix.
         * 
         * This function will calculate the two-electron multipole overlap
         * matrices.
         */
        void setupTwoElectronIntegrals (Parallel const & par, CommandLine const & cmd);
        
        /**
         * @brief Verbosity control.
         * 
         * Setting this to false will inhibit standard output messages from this class.
         */
        void verbose (bool v) { verbose_ = v; }
        
        /**
         * @brief Maximal multipole moment.
         *
         * Returns maximal multipole (also known as "angular momentum transfer"),
         * for which there are precomputed two-electron integrals.
         */
        int maxlambda () const { return Nlambdas_ - 1; }
        
        // compute overlap matrix of radial function
        Complex computeOverlapMatrixElement_iknot (Bspline const & bspline, GaussLegendre const & g, int i, int j, CCFunction func, int iknot) const;
        Complex computeOverlapMatrixElement (Bspline const & bspline, GaussLegendre const & g, int i, int j, CCFunction func) const;
        SymBandMatrix<Complex> computeOverlapMatrix (Bspline const & bspline, CCFunction func) const;
        
        /**
         * @brief Compute overlap matrix of two B-spline bases.
         * 
         * This function returns a matrix of B-spline overlaps across two bases.
         * The overlap integrals are calculated only between the given bounds
         * @c R_left and @c R_right (these are the real parts of the potentially
         * complex coordinates). Both bounds must be aligned to some B-spline knot.
         */
        static CsrMatrix<LU_int_t,Complex> computeOverlapMatrix
        (
            Bspline const & bspline1,
            Bspline const & bspline2,
            Real R_left,
            Real R_right
        );
        
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
            Bspline const & bspline,
            GaussLegendre const & g,
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
            Bspline const & bspline,
            GaussLegendre const & g,
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
            Bspline const & bspline,
            GaussLegendre const & g,
            int a, int i, int j,
            int iknot,
            CCFunction weight,
            Real scale
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
            Bspline const & bspline,
            GaussLegendre const & g,
            int a, int i, int j,
            bool truncate = false,
            bool scale = false
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
        
        /**
         * @brief Two-electron integral for multipole @f$ lambda @f$ - diagonal knots contribution.
         * 
         * The two-electron integral calculated using @ref computeR is calculated
         * as a sum of contributions from individual 2D parcels bounded by Cartesian product
         * of the relevant knots. The two-electron integral on a chosen parcel often separates
         * to a product of one-dimensional integrals. This function returns the contribution from
         * all 2D knot parcels, where it does NOT separate.
         */
        cArray diagonalR (int lambda) const;
        
        /**
         * @brief Two-electron integral for multipole @f$ lambda @f$.
         * 
         * Compute the two-electron (Slater-type) four-B-spline integral.
         * Uses the "diagonal" contribution from @ref diagonalR and the off-diagonal
         * contributions calculated from products of the partial moments.
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
        Complex computeRdiag
        (
            int L,
            int a, int b, int c, int d,
            int iknot, int iknotmax
        ) const;
        
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
         * @brief Calculate particular sub-matrix of the radial integrals matrix.
         * 
         * Calculate particular sub-matrix of the radial integrals matrix (with block indices "i" and "k")
         * and return it in a form of a dense array (copying structure of the overlap matrix).
         */
        SymBandMatrix<Complex> calc_R_tr_dia_block
        (
            unsigned lambda,
            int i, int k,
            bool simple = false
        ) const;
        
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
        
        /** 
         * @brief Compute B-spline overlaps of arbitrary one-dimensional function.
         * 
         * @param bspline B-spline basis.
         * @param g Gauss-Legendre integrator adapted to the B-spline basis.
         * @param funct Some one-dimensional function (Complex -> Complex).
         */
        cArray overlap
        (
            Bspline const & bspline,
            GaussLegendre const & g,
            CCFunction funct
        ) const;
        
        /** 
         * @brief Compute B-spline overlaps of arbitrary two-dimensional function.
         * 
         * @param bsplinex B-spline basis for x-axis.
         * @param gx Gauss-Legendre integrator adapted to the B-spline basis for x-axis.
         * @param bspliney B-spline basis for y-axis.
         * @param gy Gauss-Legendre integrator adapted to the B-spline basis for y-axis.
         * @param funct Some one-dimensional function (Complex x Complex -> Complex).
         */
        cArray overlap
        (
            Bspline const & bsplinex,
            GaussLegendre const & gx,
            Bspline const & bspliney,
            GaussLegendre const & gy,
            C2CFunction funct
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
         * @param weightf Weight function to multiply every value of the hydrogenic function (Complex -> Real).
         */
        cArray overlapP
        (
            Bspline const & bspline,
            GaussLegendre const & g,
            int n, int l
        ) const;
        
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
        cArray overlapj
        (
            Bspline const & bspline,
            GaussLegendre const & g,
            int maxell,
            const rArrayView vk,
            bool fast_bessel = false
        ) const;
        
        // Return reference to the B-spline object.
        Bspline const & bspline_inner_x () const { return bspline_inner_x_; }
        Bspline const & bspline_full_x  () const { return bspline_full_x_ ; }
        Bspline const & bspline_inner_y () const { return bspline_inner_y_; }
        Bspline const & bspline_full_y  () const { return bspline_full_y_ ; }
        
        // Return the Gauss-Legendre integrator object.
        GaussLegendre const & gaussleg_inner_x () const { return g_inner_x_; }
        GaussLegendre const & gaussleg_full_x  () const { return g_full_x_ ; }
        GaussLegendre const & gaussleg_inner_y () const { return g_inner_y_; }
        GaussLegendre const & gaussleg_full_y  () const { return g_full_y_ ; }
        
        #define OneElectronMatrixAccessors(M) \
            SymBandMatrix<Complex> const & M##_inner_x () const { return M##_inner_x_; } \
            SymBandMatrix<Complex> const & M##_full_x  () const { return M##_full_x_; } \
            SymBandMatrix<Complex> const & M##_inner_y () const { return M##_inner_y_; } \
            SymBandMatrix<Complex> const & M##_full_y  () const { return M##_full_y_; } \
            Complex M##_inner_x (std::size_t i, std::size_t j) const { return M##_inner_x_(i,j); } \
            Complex M##_full_x  (std::size_t i, std::size_t j) const { return M##_full_x_(i,j);  } \
            Complex M##_inner_y (std::size_t i, std::size_t j) const { return M##_inner_y_(i,j); } \
            Complex M##_full_y  (std::size_t i, std::size_t j) const { return M##_full_y_(i,j);  }
        
        // Access the precomputed one-electron overlap matrices.
        OneElectronMatrixAccessors(D)
        OneElectronMatrixAccessors(S)
        OneElectronMatrixAccessors(Mm1_tr)
        OneElectronMatrixAccessors(Mm2)
        
        #define OneElectronMatrixArrayAccessors(M) \
            SymBandMatrix<Complex> const & M##_inner_x (int L) const { return M##_inner_x_[L]; } \
            SymBandMatrix<Complex> const & M##_full_x  (int L) const { return M##_full_x_ [L]; } \
            SymBandMatrix<Complex> const & M##_inner_y (int L) const { return M##_inner_y_[L]; } \
            SymBandMatrix<Complex> const & M##_full_y  (int L) const { return M##_full_y_ [L]; } \
            Complex M##_inner_x (int L, std::size_t i, std::size_t j) const { return M##_inner_x_[L](i,j); } \
            Complex M##_full_x  (int L, std::size_t i, std::size_t j) const { return M##_full_x_ [L](i,j); } \
            Complex M##_inner_y (int L, std::size_t i, std::size_t j) const { return M##_inner_y_[L](i,j); } \
            Complex M##_full_y  (int L, std::size_t i, std::size_t j) const { return M##_full_y_ [L](i,j); }
        
        // Access the precomputed scaled full integral moments of order L / -L-1.
        OneElectronMatrixArrayAccessors(Mtr_L)
        OneElectronMatrixArrayAccessors(Mtr_mLm1)
        
        #define OneElectronPartialMatrixAccessors(M) \
            cArrayView M##_inner_x (int L = -1) const \
            { \
                if (L < 0) return M##_inner_x_; \
                std::size_t mi_size = bspline_inner_x_.Nspline() * (2 * bspline_inner_x_.order() + 1) * (bspline_inner_x_.order() + 1); \
                return cArrayView (M##_inner_x_, L * mi_size, mi_size); \
            } \
            cArrayView M##_inner_y (int L = -1) const \
            { \
                if (L < 0) return M##_inner_y_; \
                std::size_t mi_size = bspline_inner_y_.Nspline() * (2 * bspline_inner_y_.order() + 1) * (bspline_inner_y_.order() + 1); \
                return cArrayView (M##_inner_y_, L * mi_size, mi_size); \
            } \
            cArrayView M##_full_x (int L = -1) const \
            { \
                if (L < 0) return M##_full_x_; \
                std::size_t mi_size = bspline_full_x_.Nspline() * (2 * bspline_full_x_.order() + 1) * (bspline_full_x_.order() + 1); \
                return cArrayView (M##_full_x_, L * mi_size, mi_size); \
            } \
            cArrayView M##_full_y (int L = -1) const \
            { \
                if (L < 0) return M##_full_y_; \
                std::size_t mi_size = bspline_full_y_.Nspline() * (2 * bspline_full_y_.order() + 1) * (bspline_full_y_.order() + 1); \
                return cArrayView (M##_full_y_, L * mi_size, mi_size); \
            }
        
        // Access the precomputed scaled partial overlap matrices of order L / -L-1.
        OneElectronPartialMatrixAccessors(Mitr_L)
        OneElectronPartialMatrixAccessors(Mitr_mLm1)
        
        // Return reference to the precomputed array of diagonal contributions to two-electron integrals.
        cArray const & R_tr_dia_diag (unsigned i) const { return R_tr_dia_diag_[i]; }
        
        // Return reference to the precomputed matrix of two-electron integrals for given multipole.
        BlockSymBandMatrix<Complex> const & R_tr_dia (unsigned i) const { return R_tr_dia_[i]; }
        
    private:
        
        // B-spline environment
        Bspline bspline_inner_x_, bspline_full_x_;
        Bspline bspline_inner_y_, bspline_full_y_;
        
        // Gauss-Legendre integrator
        GaussLegendre g_inner_x_, g_inner_y_;
        GaussLegendre g_full_x_, g_full_y_;
        
        // one-electron moment and overlap matrices
        SymBandMatrix<Complex> D_inner_x_, S_inner_x_, Mm1_inner_x_, Mm1_tr_inner_x_, Mm2_inner_x_;
        SymBandMatrix<Complex> D_inner_y_, S_inner_y_, Mm1_inner_y_, Mm1_tr_inner_y_, Mm2_inner_y_;
        SymBandMatrix<Complex> D_full_x_, S_full_x_, Mm1_full_x_, Mm1_tr_full_x_, Mm2_full_x_;
        SymBandMatrix<Complex> D_full_y_, S_full_y_, Mm1_full_y_, Mm1_tr_full_y_, Mm2_full_y_;
        
        // one-electron full integral moments for various orders (used to calculate R-integrals)
        std::vector<SymBandMatrix<Complex>> Mtr_L_inner_x_, Mtr_mLm1_inner_x_;
        std::vector<SymBandMatrix<Complex>> Mtr_L_inner_y_, Mtr_mLm1_inner_y_;
        std::vector<SymBandMatrix<Complex>> Mtr_L_full_x_, Mtr_mLm1_full_x_;
        std::vector<SymBandMatrix<Complex>> Mtr_L_full_y_, Mtr_mLm1_full_y_;
        
        // partial one-electron integral moments for various orders (used to calculate R-integrals)
        cArray Mitr_L_inner_x_, Mitr_mLm1_inner_x_;
        cArray Mitr_L_inner_y_, Mitr_mLm1_inner_y_;
        cArray Mitr_L_full_x_, Mitr_mLm1_full_x_;
        cArray Mitr_L_full_y_, Mitr_mLm1_full_y_;
        
        // diagonal contributions to R_tr_dia
        cArrays R_tr_dia_diag_;
        
        // two-electron integral matrices
        std::vector<BlockSymBandMatrix<Complex>> R_tr_dia_;
        
        // verbose output
        bool verbose_;
        
        // number of multipole matrices
        int Nlambdas_;
};

#endif // HEX_ECS_RADIAL_H
