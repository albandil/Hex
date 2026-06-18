//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2026, Jakub Benda, Charles University in Prague                    //
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

#include "hex-hydrogen.h"
#include "hex-special.h"

#include "eigchans.h"

void Eigchans::calc(InputFile const& inp, std::vector<std::pair<int,int>> const & ang)
{
    for (auto [ni, li, mi] : inp.instates)
    {
        // skip already processed shell
        if (eigmoms.find(ni) != eigmoms.end())
            continue;

        // calculate the 1/2r² potential matrix
        ColMatrix<Real> M(ang.size(), ang.size());

        for (std::size_t ill = 0; ill < ang.size(); ill++)
        {
            auto [l1, l2] = ang[ill];

            M(ill, ill) = l2*(l2 + 1);

            for (std::size_t illp = 0; illp < ang.size(); illp++)
            {
                auto [l1p, l2p] = ang[illp];

                auto f = special::computef(1, l1, l2, l1p, l2p, inp.L);
                auto d = Hydrogen::radialDipole(inp.Za, ni, l1, l1p);

                M(ill, illp) += 2*f*d;
            }
        }

        // diagonalize the matrix
        M.diagonalize(eigvals[ni], &eigvecs[ni]);

        // evaluate the generalized angular momenta
        eigmoms[ni] = -0.5_r + sqrt(0.25_z + eigvals[ni]);
    }
}
