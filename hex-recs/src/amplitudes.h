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
        
        /**
         * @brief Transition description.
         * 
         * The structure contains ll discrete quantum numbers identificating
         * the transition
         * \f[
         *     | n_i l_i j_i m_i \rangle | \mathbf{k}_i \sigma_i \rangle
         *     \rightarrow
         *     | n_f l_f j_f m_f \rangle | \mathbf{k}_f \sigma_f \rangle \,.
         * \f]
         * 
         * The structure also defined an ordering operator so that it can be used
         * as an index in std::map container.
         */
        typedef struct tTransition
        {
            // initial state
            int ni, li, two_ji, two_mi, two_si;
            
            // final state
            int nf, lf, two_jf, two_mf, two_sf;
            
            // necessary ordering operator so that we can use this structure as an index in std::map
            bool operator < (tTransition const & t) const
            {
                #define Compare(x) if (x < t.x) return true; if (x > t.x) return false;
                
                // perform standard lexicographic comparison (by components)
                Compare(ni); Compare(li); Compare(two_ji); Compare(two_mi); Compare(two_si);
                Compare(nf); Compare(lf); Compare(two_jf); Compare(two_mf); Compare(two_sf);
                
                // they are equal
                return false;
            }
        }
        Transition;
        
        // shorteded data types
        typedef std::map<std::tuple<int,int,int>,cArray> ThreeIntComplexMap;
        typedef std::map<std::tuple<int,int>,cArray> TwoIntComplexMap;
        
        // Λ[ie] indexed by (ni,li,2*ji,2*mi,nf,lf,2*jf,2*mf) and (L,S,l')
        std::map<Transition,ThreeIntComplexMap> Lambda_LSlp;
        
        // T[ie] indexed by (ni,li,2*ji,2*mi,nf,lf,2*jf,2*mf) and (l',j')
        std::map<Transition,TwoIntComplexMap> Tmat_jplp;
        
        // σ[ie] indexed by (ni,li,2*ji,2*mi,nf,lf,2*jf,2*mf)
        std::map<Transition,rArray> sigma;
        
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
        ThreeIntComplexMap computeLambda_ (Transition, rArray const &);
        
        TwoIntComplexMap computeTmat_ (Transition, rArray const &);
        
        rArray computeSigma_ (Transition, rArray const &);
        
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
