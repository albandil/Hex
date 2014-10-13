//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2014, Jakub Benda, Charles University in Prague                    //
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

#ifndef HEX_CCC_POTENTIAL
#define HEX_CCC_POTENTIAL

#include "arrays.h"
#include "basis.h"
// #include "gauss.h"
#include "quadrature.h"

class PotentialMatrix
{
    public:
        
        /**
         * @brief Construct the matrix.
         * 
         * This function will precompute entries of the potential matrix based
         * on the chosen basis and quadrature rule. This is one of the most
         * time expensive operations because for every element of the matrix
         * a double integral has to be evaluated.
         * @param basis LaguerreBasis object containing the bases.
         * @param quadrature QuadratureRule object containing the quadrature settings.
         * @param J Total angular momentum.
         * @param S Total spin.
         * @param Pi Total parity.
         */
        PotentialMatrix
        (
            LaguerreBasis const & basis,
            QuadratureRule const & quadrature,
            int J, int S, int Pi
        );
        
        /**
         * @brief Return reference to the matrix of the potential.
         */
        RowMatrix<double> const & matrix() const;
        
    private:
        
        /// Basis.
        LaguerreBasis const & basis_;
        
        /// Quadrature rule.
        QuadratureRule const & quadrature_;
        
        /// The (symmetrical) matrix of the potential.
        RowMatrix<double> matrix_;
        
        
//         GaussLegendre<double> g_;
        
        /**
         * @brief Double integral @f$ I_{\mathrm{dir}} @f$.
         */
        //@{
        double ComputeIdir_nested
        (
            int lambda,
            int L, int i, int l, double k,
            int Lp, int ip, int lp, double kp
        ) const;
        double ComputeIdir_Romberg
        (
            int lambda,
            int L, int i, int l, double k,
            int Lp, int ip, int lp, double kp
        ) const;
        double ComputeIdir_semisymbolic
        (
            int lambda,
            symbolic::poly const & xi, int l, double k,
            symbolic::poly const & xip, int lp, double kp
        ) const;
        //@}
};

class MatrixEquation
{
    public:
        
        // constructor
        MatrixEquation
        (
            QuadratureRule const & quadrature,
            PotentialMatrix const & potential
        );
        
        // solve the equation
        rArray solve () const;
};

#endif /* HEX_CCC_POTENTIAL */

