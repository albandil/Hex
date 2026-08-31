//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2017, Jakub Benda, Charles University in Prague                    //
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

#ifndef HEX_ANG_H
#define HEX_ANG_H

#include "hex-arrays.h"

#include "inout.h"

/**
 * @brief Angular basis.
 * 
 * List of coupled angular states included in the calculation.
 *
 * The class keeps two lists. @ref states_full is every coupled state (l1,l2) that
 * the total quantum numbers admit; @ref states is the subset that is actually solved
 * for. The two are the same unless the exchange symmetry is folded in (see
 * @ref folded), in which case the solved list holds only the states with l1 <= l2
 * and the remaining states are reconstructed from their mirror counterparts.
 *
 * The full list is always ordered so that the solved states come first. An index into
 * the solved list is therefore also a valid index into the full list and refers to the
 * same angular state, which lets the two be used interchangeably wherever only solved
 * blocks are addressed.
 */
class AngularBasis
{
    public:

        /// Constructor.
        AngularBasis (InputFile const & inp, bool fold = false);
        AngularBasis (AngularBasis const & ang)
            : L_(ang.L_), S_(ang.S_), Pi_(ang.Pi_), maxell_(ang.maxell_),
            maxlambda_(ang.maxlambda_), folded_(ang.folded_), states_(ang.states_),
            states_full_(ang.states_full_), mirror_(ang.mirror_), f_(ang.f_) {}

        /// List of coupled angular states that are solved for.
        std::vector<std::pair<int,int>> const & states () const { return states_; }

        /// List of all coupled angular states, those from @ref states first.
        std::vector<std::pair<int,int>> const & states_full () const { return states_full_; }

        /**
         * @brief Index of the mirror block.
         *
         * Returns the position of the state (l2,l1) in @ref states_full, given the
         * position of the state (l1,l2) there. States with l1 = l2 are their own
         * mirrors. When the mirror state is not a part of the basis at all (which
         * happens for a basis restricted by 'exchange = 0' in the input file), -1
         * is returned.
         */
        int mirror (int ill) const { return mirror_[ill]; }

        /**
         * @brief Whether the exchange symmetry is folded into the equations.
         *
         * When true, @ref states is only the l1 <= l2 half of @ref states_full and the
         * solution segments of the other half are given by @ref exchange_sign times the
         * transpose of the segments of their mirror blocks.
         */
        bool folded () const { return folded_; }

        /**
         * @brief Sign relating a block to its mirror.
         *
         * The (anti)symmetry of the wave function with respect to the exchange of the
         * two electrons gives
         * @f[
         *     \psi_{l_2 l_1}(r_1,r_2) = (-1)^{l_1+l_2+L+S} \psi_{l_1 l_2}(r_2,r_1) \ ,
         * @f]
         * and because every state of the basis satisfies l1 + l2 = 2n + L + Pi, the sign
         * reduces to (-1)^(S+Pi), which is the same for all blocks of the basis.
         */
        Real exchange_sign () const { return (S_ + Pi_) % 2 == 0 ? 1.0_r : -1.0_r; }

        /// Highest orbital number.
        int maxell () const { return maxell_; }

        /// Highest multipole.
        int maxlambda () const { return maxlambda_; }

        /// Angular integrals. Both indices refer to @ref states_full.
        Real f (int ill, int illp, int lambda) const
        {
            return f_[(lambda * states_full_.size() + ill) * states_full_.size() + illp];
        }

        /// Angular integrals.
        rArray const & f () const { return f_; }

        /// Total angular momentum.
        int L () const { return L_; }

        /// Total spin.
        int S () const { return S_; }
        int & S () { return S_; }

        /// Total parity.
        int Pi () const { return Pi_; }

    private:

        // Quantum numbers.
        int L_, S_, Pi_, nL_;

        // Highest orbital number.
        int maxell_;

        // Highest multipole.
        int maxlambda_;

        // Whether the l1 > l2 states have been excluded from 'states_'.
        bool folded_;

        // List of the coupled angular states that are solved for.
        std::vector<std::pair<int,int>> states_;

        // List of all coupled angular states, the solved ones first.
        std::vector<std::pair<int,int>> states_full_;

        // Position of the mirror state in 'states_full_' for every state there.
        std::vector<int> mirror_;

        // Angular integrals.
        rArray f_;
};

/**
 * @brief Mirror a solution segment.
 *
 * Rewrites the solution segment @c src of some angular block (l1,l2) into the segment
 * @c dst of the mirror block (l2,l1), which amounts to interchanging the two electrons.
 * The inner region is a square matrix of B-spline coefficients and is transposed; the
 * two groups of asymptotic channels describe the escape of either of the electrons, so
 * they change places. The sign that also relates the two blocks is
 * @ref AngularBasis::exchange_sign and is *not* applied here.
 *
 * Both segments have the layout
 * <pre>
 *   [ Nspline_inner x Nspline_inner ][ Nchan1 x Nspline_outer ][ Nchan2 x Nspline_outer ]
 * </pre>
 * except that the mirror block has the channel counts swapped. The two electrons must
 * share the same B-spline basis, otherwise the two layouts are not compatible at all.
 */
void mirror_segment
(
    const cArrayView src,
    cArrayView dst,
    int Nspline_inner,
    int Nspline_outer,
    int Nchan1,
    int Nchan2
);

#endif
