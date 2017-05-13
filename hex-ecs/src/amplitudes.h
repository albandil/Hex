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

#ifndef HEX_AMPLITUDES
#define HEX_AMPLITUDES

#include <vector>

#include "hex-arrays.h"
#include "hex-chebyshev.h"

#include "bspline.h"
#include "io.h"
#include "parallel.h"
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
            Bspline const & bspline_inner,
            Bspline const & bspline_full,
            InputFile const & inp,
            Parallel const & par,
            CommandLine const & cmd,
            std::vector<std::pair<int,int>> const & ang
        );
        
        /**
         * @brief Extract the amplitudes.
         * 
         * This is the main routine of the class. It will do the extraction of
         * T-matrices, which involves some linear algebra (~ radial integrals in
         * B-spline basis) and computations of various coupling coefficients. In the
         * end the class will contain T-matrices and partial cross sections for
         * transitions requested in the input file.
         */
        void extract (std::string directory = ".");
        
        /**
         * @brief Write SQL batch file.
         * 
         * This routine will write out T-matrices in the form of SQL batch file
         * usable in hex-db. The file's name will be "tmat-n-L-S-Pi.sql", where the symbols
         * will be substituted by actual values of these total quantum numbers.
         * See the documentation of hex-db for details on the output format.
         */
        void writeSQL_files (std::string directory = ".");
        
        /**
         * @brief Write integral cross sections to file.
         * 
         * This routine will write the partial integral cross sections (i.e. for
         * the given total quantum numbers) for every transition and energy into
         * a text column-formatted file "ics-n-L-S-Pi.dat". (Symbols will be replaced
         * by the actual values of these quantum numbers.) The format is as follows:
         * The first column lists the energies, the next columns are cross sections
         * for every possible transition allowed by input file. The first row is
         * the header, which is commented out by hash symbol ('#') for easy use in
         * gnuplot. An example of a name of a column is
         @verbatim
              2p(-1)-2p(1)
         @endverbatim
         * which is a transition from
         * @f$ |nlm\rangle |\mathbf{k}\rangle = |2,1,-1\rangle |\mathbf{k}_i\rangle @f$
         * to
         * @f$ |nlm\rangle |\mathbf{k}\rangle = |2,1,1\rangle |\mathbf{k}_f\rangle @f$.
         */
        void writeICS_files (std::string directory = ".");
        
        /**
         * @brief Set verbosity level.
         * 
         * This function can be used to mute the text output of this class.
         */
        void verbose (bool b) { verbose_ = b; }
        
    private:
        
        /**
         * @brief Transition description.
         * 
         * The structure contains ll discrete quantum numbers identificating
         * the transition
         * \f[
         *     | n_i l_i m_i \rangle | \mathbf{k}_i \rangle
         *     \rightarrow
         *     | n_f l_f m_f \rangle | \mathbf{k}_f \rangle \,.
         * \f]
         * 
         * The structure also defined an ordering operator so that it can be used
         * as an index in std::map container.
         */
        typedef struct tTransition
        {
            // initial state
            int ni, li, mi;
            
            // final state
            int nf, lf, mf;
            
            // necessary ordering operator so that we can use this structure as an index in std::map
            bool operator < (tTransition const & t) const
            {
                #define Compare(x) if (x < t.x) return true; if (x > t.x) return false;
                
                // perform standard lexicographic comparison (by components)
                Compare(ni); Compare(li); Compare(mi);
                Compare(nf); Compare(lf); Compare(mf);
                
                // they are equal
                return false;
            }
        }
        Transition;
        
        // Λ[ie] for both spins indexed by (ni,li,mi,nf,lf,mf) and l'
        std::map<Transition,std::vector<std::pair<cArray,cArray>>> Lambda_Slp;
        
        // T[ie] for both spins indexed by (ni,li,mi,nf,lf,mf) and l'
        std::map<Transition,std::vector<std::pair<cArray,cArray>>> Tmat_Slp;
        
        // σ[ie] for both spins indexed by (ni,li,mi,nf,lf,mf)
        std::map<Transition,std::pair<rArray,rArray>> sigma_S;
        
        // Ξ[ie] (lists of Chebyshev coefficients) indexed by (ni,li,mi) and (l1,l2)
        std::map<Transition,std::vector<std::pair<cArray,cArray>>> Xi_Sl1l2;
        
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
         *           \hat{j}_l(k_f\bullet)
         *        \right]_{R_0}
         *        \mathrm{d}r_1 \ .
         * @f]
         */
        void computeLambda_ (Transition T, BlockArray<Complex> & solution, int ie, int Spin);
        
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
        void computeTmat_ (Transition T);
        
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
         * 
         * @param bspline The B-spline environment to use when evaluating the solutions.
         * @param maxell Maximal angular momentum of the free electrons.
         * @param L Total angular momentum (partial wave).
         * @param Spin Conserved total spin (partial wave).
         * @param Pi Total parity (partial wave).
         * @param ni Initial atomic state - principal quantum number.
         * @param li Initial atomic state - orbital quantum number.
         * @param mi Initial atomic state - magnetic quantum number.
         * @param Ei Initial projectile energies.
         * @param ics Ionization cross section (on return).
         * @param coupled_states List of all coupled-state angular momenta pairs.
         * @return Vector of radial integrals.
         */
        void computeXi_ (Transition T, BlockArray<Complex> & solution, int ie, int Spin);
        
        /**
         * @brief Evaluate cross sections.
         * 
         * This routine will reuse the T-matrices calculated by @ref computeTmat_
         * and compute the cross sections.
         */
        void computeSigma_ (Transition T);
        
        Chebyshev<double,Complex> fcheb (cArrayView const & PsiSc, Real kmax, int l1, int l2);
        void computeSigmaIon_ (Transition T);
        
        cArray readAtomPseudoState (int lf, int ichan) const;
        
        // B-spline environment
        Bspline const & bspline_inner_;
        Bspline const & bspline_full_;
        
        // radial integrals
        RadialIntegrals rad_;
        
        // input data
        InputFile const & inp_;
        
        // parallel environment
        Parallel const & par_;
        
        // command line switches
        CommandLine const & cmd_;
        
        // angular basis
        std::vector<std::pair<int,int>> const & ang_;
        
        // solution reader
        SolutionIO reader_;
        
        // standard output verbosity
        bool verbose_;
};

#endif
