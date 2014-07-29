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

#ifndef HEX_AMPLITUDES
#define HEX_AMPLITUDES

#include <vector>

#include "angular.h"
#include "arrays.h"
#include "bspline.h"
#include "radial.h"

class Amplitudes
{
    public:
        
        Amplitudes
        (
            Bspline const & bspline, InputFile const & inp, Parallel const & par,
            AngularBasis const & ang
        ) : bspline_(bspline), rad_(bspline_), inp_(inp), par_(par), ang_(ang)
        {
            // nothing to do
        }
        
        void extract ();
        
        void writeSQL_files ();
        
        void writeICS_files ();
        
    private:
        
        // Λ[ie] indexed by (ni,li,2*ji,2*mi,nf,lf,2*jf,2*mf) and (L,S,l')
        std::map <
            std::tuple<int,int,int,int,int,int,int,int>,
            std::map<std::tuple<int,int,int>,cArray>
        > Lambda_LSlp;
        
        // T[ie] indexed by (ni,li,2*ji,2*mi,nf,lf,2*jf,2*mf) and (l',j')
        std::map <
            std::tuple<int,int,int,int,int,int,int,int>,
            std::map<std::tuple<int,int>,cArray>
        > Tmat_jplp;
        
        // σ[ie] indexed by (ni,li,2*ji,2*mi,nf,lf,2*jf,2*mf)
        std::map <
            std::tuple<int,int,int,int,int,int,int,int>,
            rArray
        > sigma;
        
        /**
         * @brief Extract radial part of scattering amplitude.
         * 
         * Compute radial integrals for evaluation of the discrete T-matrices,
         * @f[
         *    \Lambda_l^{(1)LMS} = 
         *      \int_0^{R_0}
         *        P_{n_f l_f}(r_1)
         *        \mathcal{W}\left[
         *           \psi_{\ell_1 \ell_2}^{LMS}(r_1,\bullet),
         *            \hat{j}_l(k_f\bullet)
         *         \right]_{R_0}
         *        \mathrm{d}r_1 \ .
         *  @f]
         */
        std::map<std::tuple<int,int,int>,cArray> computeLambda_
        (
            std::tuple<int,int,int,int,int,int,int,int> transition,
            rArray const & kf
        );
        
        std::map<std::tuple<int,int>,cArray> computeTmat_
        (
            std::tuple<int,int,int,int,int,int,int,int> transition,
            rArray const & kf
        );
        
        rArray computeSigma_
        (
            std::tuple<int,int,int,int,int,int,int,int> transition,
            rArray const & kf
        );
        
        /**
         * @brief Extract radial part of ionization amplitude.
         * 
         * Compute hyperangular integrals for evaluation of the ionization T-matrices,
         * @f[
         *     \Xi_{\ell_1 \ell_2}^{LS}(k_1,k_2) =
         *       \int_0^{\pi/2}
         *         \left(
         *           F_{\ell_1}(k_1,r_1) F_{\ell_2}(k_2,r_2)
         *           \frac{\partial}{\partial\rho}
         *           \psi^{LS}_{\ell_1\ell_2}(r_1,r_2)
         *           -
         *           \psi^{LS}_{\ell_1\ell_2}(r_1,r_2)
         *           \frac{\partial}{\partial\rho}
         *           F_{\ell_1}(k_1,r_1) F_{\ell_2}(k_2,r_2)
         *         \right)
         *         \rho\mathrm{d}\alpha \ ,
         * @f]
         * where @f$ r_1 = \rho \cos\alpha @f$ and @f$ r_2 = \rho \sin\alpha @f$.
         */
        /*cArrays computeXi
        (
            // TODO
        );*/
        
        // B-spline environment
        Bspline const & bspline_;
        
        // radial integrals
        RadialIntegrals const & rad_;
        
        // input data
        InputFile const & inp_;
        
        // parallel environment
        Parallel const & par_;
        
        // angular basis
        AngularBasis const & ang_;
};

#endif
