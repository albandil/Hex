//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2025, Jakub Benda, Charles University in Prague                    //
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

#include <iostream>

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"
#include "hex-blas.h"
#include "hex-hydrogen.h"
#include "hex-misc.h"

// --------------------------------------------------------------------------------- //

#include "gauss.h"
#include "radial.h"

// --------------------------------------------------------------------------------- //

#include "DiagPreconditioner.h"

// --------------------------------------------------------------------------------- //

std::string DiagPreconditioner::description () const
{
    return "Jacobi preconditioner.";
}

void DiagPreconditioner::setup ()
{
    NoPreconditioner::setup();

    if (verbose_)
        std::cout << "Set up diag preconditioner" << std::endl << std::endl;

    RadialIntegrals const& rad = rad_inner();

    SymBandMatrix<Complex> const& S = rad.S();
    SymBandMatrix<Complex> const& D = rad.D();
    SymBandMatrix<Complex> const& Mm1 = rad.Mm1();
    SymBandMatrix<Complex> const& Mm2 = rad.Mm2();

    size_t Nspline = rad.bspline().Nspline();

    diagonal = cBlockArray(ang_->states().size());

    // calculate diagonal
    for (int ill = 0; ill < ang_->states().size(); ill++)
    {
        int l1 = ang_->states()[ill].first;
        int l2 = ang_->states()[ill].second;

        diagonal[ill].resize(Nspline*Nspline);

        for (int i = 0; i < Nspline; i++)
        {
            cArrayView segment = diagonal[ill].slice(i*Nspline, (i + 1)*Nspline);

            for (int j = 0; j < Nspline; j++)
            {
                segment[j] = 0.5_r * D(i,i) * S(j,j) - inp_->Za            * Mm1(i,i) * S(j,j) + 0.5_r * l1 * (l1 + 1) * Mm2(i,i) * S(j,j)
                           + 0.5_r * D(j,j) * S(i,i) + inp_->Za * inp_->Zp * Mm1(j,j) * S(i,i) + 0.5_r * l2 * (l2 + 1) * Mm2(j,j) * S(i,i);
            }

            for (int lambda = 0; lambda <= ang_->maxlambda(); lambda++)
            {
                double f = ang_->f(ill,ill,lambda);

                if (f != 0)
                {
                    //segment -= inp_->Zp * f * rad.calc_R_tr_dia_block(lambda, i, i).data().slice(0, Nspline);

                    SymBandMatrix<Complex> R = rad.calc_R_tr_dia_block(lambda, i, i);

                    for (int j = 0; j < Nspline; j++)
                        segment[j] -= inp_->Zp * f * R(j, j);
                }
            }
        }
    }
}

void DiagPreconditioner::precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const
{
    SymBandMatrix<Complex> const& S = rad_inner().S();

    size_t Nspline = rad_inner().bspline().Nspline();

    for (int ill = 0; ill < ang_->states().size(); ill++)
    {
        for (int i = 0; i < Nspline; i++)
        {
            for (int j = 0; j < Nspline; j++)
            {
                z[ill][i*Nspline + j] = r[ill][i*Nspline + j]/(E_*S(i,i)*S(j,j) - diagonal[ill][i*Nspline + j]);

                if (not (z[ill][i*Nspline + j] == z[ill][i*Nspline + j]))
                    z[ill][i*Nspline + j] = r[ill][i*Nspline + j];
            }
        }
    }
}

// --------------------------------------------------------------------------------- //

addClassToParentRunTimeSelectionTable(PreconditionerBase, DiagPreconditioner)

// --------------------------------------------------------------------------------- //

