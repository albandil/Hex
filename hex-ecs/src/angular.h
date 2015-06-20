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

#ifndef HEX_ANGULAR_H
#define HEX_ANGULAR_H

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

/**
 * @brief Angular basis class.
 * 
 * The set of equations created by discretization of the ECS Schrödinger equation
 * in the basis of B-splines has a block structure determined by the angular momenta
 * of the two electrons. Every solution vector and right-hand side can be, similarly,
 * split into segments indexed by angular momenta.
 * 
 * In the case of the vectors (the solution or the right-hand side) the index is a
 * pair (ℓ₁,ℓ₂). In the matrix every blocks in consistently indexed by row index
 * (ℓ₁,ℓ₂) and column index (ℓ₁',ℓ₂'). The indices (ℓ₁,ℓ₂) and (ℓ₂,ℓ₁) which differ
 * only by swapping of the electron angular momenta, number linearly dependent states
 * due to the overall spatial symmetry or antisymmetry.
 * 
 * This class contains the full list of all angular momentum pairs that are considered
 * in the calculation. Also it offers some interface functions that relate linearly
 * dependent angular states to so called "basic symmetries", where by convention
 * ℓ₁ <= ℓ₂.
 * 
 * The class offers a limited STL-container-like interface through the functions
 * 'size', 'begin', 'end' and 'operator[]'.
 */
class AngularBasis
{
    public:
        
        /// Constructor.
        AngularBasis (int L, int S, int Pi, int nL) : L_(L), S_(S), Pi_(Pi), nL_(nL)
        {
            std::cout << "Setting up the coupled angular states..." << std::endl;
            
            // for given L, Π and levels list all available (ℓ₁ℓ₂) pairs
            for (int ell = 0; ell <= nL_; ell++)
            {
                std::cout << "\t-> [" << ell << "] ";
                
                // get sum of the angular momenta for this angular level
                int sum = 2 * ell + L_ + Pi_;
                
                // for all angular momentum pairs that do compose L
                for (int l1 = ell + Pi_; l1 <= sum - ell - Pi_; l1++)
                {
                    std::cout << "(" << l1 << "," << sum - l1 << ") ";
                    ll_.push_back(std::make_pair(l1, sum - l1));
                }
                std::cout << std::endl;
            }
            
            // sanity check
            assert((int)ll_.size() == (nL_ + 1) * (L_ + 1 - Pi_));
            
            std::cout << "\t-> The matrix of the set contains " << (nL_ + 1) * (L_ + 1 - Pi_)
                    << " diagonal blocks, of which " << (nL_ + 1) * ((L_ + 2 - Pi_) / 2)
                    << " are independent and will be solved for." << std::endl << std::endl;
        }
        
        /// Check that chosen angular pair is a basic symmetry.
        bool is_basic_symmetry (int ill) const
        {
            return ll_[ill].first <= ll_[ill].second;
        }
        
        /// Get index of basic symmetry in the full list of all angular states.
        int basic_symmetry (int ill) const
        {
            int ill0 = ill;
            if (not is_basic_symmetry(ill))
            {
                int igroup = ill / (L_ + 1 - Pi_);
                int istate = ill % (L_ + 1 - Pi_);
                ill0 = igroup * (L_ + 1 - Pi_) + (L_ + 1 - Pi_ - istate - 1);
            }
            return ill0;
        }
        
        /// Get index of the basic symmetry in hypothetic list of basic symmetries.
        int basic_symmetry_index (int ill) const
        {
            int igroup = ill / (L_ + 1 - Pi_);
            int istate = ill % (L_ + 1 - Pi_);
            return igroup * ((L_ + 2 - Pi_) / 2) + (L_ + 1 - Pi_ - istate - 1);
        }
        
        /// Get basic symmetry relative sign.
        int basic_symmetry_sign (int ill) const
        {
            return is_basic_symmetry(ill) ? +1 : ((S_ + Pi_) % 2 == 0 ? +1 : -1);
        }
        
        /// Get spin.
        int spin () const { return S_; }
        
        /// Set spin.
        void spin (int S) { S_ = S; }
        
        /// Get number of basic symmetries.
        int basic_size () const { return (nL_ + 1) * ((L_ + 2 - Pi_) / 2); }
        
    private:
        
        /// Total angular momentum.
        int L_;
        
        /// Total spin.
        int S_;
        
        /// Total parity.
        int Pi_;
        
        /// Number of groups of angular momentum pairs with given sum.
        int nL_;
        
        /// List of coupled angular states.
        std::vector<std::pair<int,int>> ll_;
    
    public:
        
        // STL container-like interface.
            
            /// Beginning of the list of angular momentum pairs.
            std::vector<std::pair<int,int>>::const_iterator begin () const { return ll_.begin(); }
            
            /// (One element past the) End of the list of angular momentum pairs.
            std::vector<std::pair<int,int>>::const_iterator end () const { return ll_.end(); }
            
            /// Size of the list of angular momentum pairs.
            int size () const { return ll_.size(); }
            
            /// Element access.
            std::pair<int,int> const & operator[] (int i) const { return ll_[i]; }
};

#endif // HEX_ANGULAR_H
