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

/**
 * @brief Extraction of amplitudes from solutions.
 * 
 * This class is responsible for the last stage of the computation:
 * the extraction of T-matrices from the calculated solutions. Given
 * the input file (and some other environment classes) it will load
 * solution one after one and extract all requested amplitudes.
 * Unavailable solutions are skipped with a warning message, so it is
 * possible to run the extraction even when some other process is yet
 * doing the main computation.
 * 
 * The class also offers the possibility to write out the data: as an
 * SQL batch file for use in hex-db, or as a text file containing
 * cross section column data.
 */
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
        
        /**
         * @brief Extract the amplitudes.
         * 
         * This is the main routine of the class. It will do the extraction of
         * T-matrices, which involves some linear algebra (~ radial integrals in
         * B-spline basis) and computations of various coupling coefficients. In the
         * end the class will contain T-matrices and partial cross sections for
         * transitions requested in the input file.
         */
        void extract ();
        
        /**
         * @brief Write SQL batch file.
         * 
         * This routine will write out T-matrices in the form of SQL batch file
         * usable in hex-db. The file's name will be "tmat-J-M.sql", where J and M
         * will be substituted by actual values of these total quantum numbers.
         * See the documentation of hex-db for details on the output format.
         */
        void writeSQL_files ();
        
        /**
         * @brief Write integral cross sections to file.
         * 
         * This routine will write the partial integral cross sections (i.e. for
         * the given total quantum numbers) for every transition and energy into
         * a text column-formatted file "ics-J-M.dat". (J and M will be replaced
         * by the actual values of these quantum numbers.) The format is as follows:
         * The first column lists the energies, the next columns are cross sections
         * for every possible transition allowed by input file. The first row is
         * the header, which is commented out by hash symbol ('#') for easy use in
         * gnuplot. An example of a name of a column is
         @verbatim
              2p1/2(-1/2)⁺-2p3/2(1/2)⁻
         @endverbatim
         * which is a transition from
         * @f$ |nljm\rangle |\mathbf{k}\sigma\rangle = |2,1,1/2,-1/2\rangle |\mathbf{k}_i +\rangle @f$
         * to
         * @f$ |nljm\rangle |\mathbf{k}\sigma\rangle = |2,1,1/2,-1/2\rangle |\mathbf{k}_f -\rangle @f$.
         * The signs (+/-) are the projectile spin orientation.
         */
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
         * Compute radial integrals for evaluation of the discrete T-matrices.
         */
        ThreeIntComplexMap computeLambda_ (Transition, rArray const &);
        
        /**
         * @brief Evaluate T-matrices.
         * 
         * This function will used precomputed Lambda-s to compute the
         * discrete scattering T-matrices. It it necessary to call
         * @ref computeLambda_ for the given trnasition first and store that
         * precomputed Lambda dataset to 'Lambda_LSlp' with the proper
         * key (= Transition object).
         * 
         * The T-matrices are stored in TwoIntComplexMap as complex array
         * (an item per impact energy) indexed by @f$ j' @f$ and @f$ l' @f$
         * (the angular momenta of the projectile).
         */
        TwoIntComplexMap computeTmat_ (Transition, rArray const &);
        
        /**
         * @brief Evaluate cross sections.
         * 
         * This routine will reuse the T-matrices calculated by @ref computeTmat_
         * and compute the cross sections.
         */
        rArray computeSigma_ (Transition, rArray const &);
        
        // B-spline environment
        Bspline const & bspline_;
        
        // radial integrals
        RadialIntegrals rad_;
        
        // input data
        InputFile const & inp_;
        
        // parallel environment
        Parallel const & par_;
        
        // angular basis
        AngularBasis const & ang_;
};

#endif
