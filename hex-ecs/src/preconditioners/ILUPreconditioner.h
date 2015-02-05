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

#ifndef HEX_ILUPRECONDITIONER_H
#define HEX_ILUPRECONDITIONER_H

#include "../preconditioners.h"

/**
 * @brief ILU-preconditioned CG-based preconditioner.
 * 
 * Enhances CGPreconditioner conjugate gradients solver by incomplete LU factorization
 * preconditioning. This is done by redefining virtual function CG_prec. The factorization
 * is drop tolerance based and is computed by UMFPACK.
 */
class ILUCGPreconditioner : public CGPreconditioner
{
    public:
        
        static const std::string name;
        static const std::string description;
        
        ILUCGPreconditioner
        (
            Parallel const & par,
            InputFile const & inp,
            std::vector<std::pair<int,int>> const & ll,
            Bspline const & bspline,
            CommandLine const & cmd
        ) : CGPreconditioner(par, inp, ll, bspline, cmd), droptol_(cmd.droptol),
            csr_blocks_(ll.size()), lu_(ll.size())
        {
            // nothing more to do
        }
        
        // reuse parent definitions
        virtual RadialIntegrals const & rad () const { return CGPreconditioner::rad(); }
        virtual void multiply (const cArrayView p, cArrayView q) const { CGPreconditioner::multiply(p,q); }
        virtual void rhs (cArray & chi, int ienergy, int instate, int Spin) const { CGPreconditioner::rhs(chi,ienergy,instate,Spin); }
        virtual void precondition (const cArrayView r, cArrayView z) const { CGPreconditioner::precondition(r,z); }
        virtual void setup () { CGPreconditioner::setup(); }
        
        // declare own definitions
        virtual void update (double E);
        
        // inner CG callback (needed by parent)
        virtual void CG_prec (int iblock, const cArrayView r, cArrayView z) const;
        
    protected:
        
        // drop tolarance for the factorizations
        double droptol_;
        
        // diagonal CSR block for every coupled state
        mutable std::vector<CsrMatrix> csr_blocks_;
        
        // LU decompositions of the CSR blocks
        mutable std::vector<CsrMatrix::LUft> lu_;
};

#endif