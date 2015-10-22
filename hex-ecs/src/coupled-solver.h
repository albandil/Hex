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

#ifndef HEX_ECS_COUPLED_SOLVER_H
#define HEX_ECS_COUPLED_SOLVER_H

#include <vector>

#include "hex-itersolve.h"

#include "ang.h"
#include "bspline.h"
#include "io.h"
#include "parallel.h"

#include "radial.h"

class CoupledSolver
{
    public:
        
        /// Constructor.
        CoupledSolver
        (
            CommandLine & cmd,
            InputFile const & inp,
            Parallel const & par,
            AngularBasis const & ang,
            std::vector<Bspline> const & bspline,
            std::vector<Bspline> const & bspline_full
        );
        
        /// Pick preconditioner for a given panel.
        void choose_preconditioner (int ipanel);
        
        /// Precompute (or load) the radial data.
        void setup_preconditioner ();
        
        /// Find solution of the scattering equations on chosen panel.
        void solve ();
        
        /// Release resources.
        void finish ();
        
        /// Construct right-hand side.
        void rhs (cArray & rhs, int ie, int instate) const;
        RadialIntegrals rad_;
        
    private:
        
        /// Command line parameters.
        CommandLine & cmd_;
        
        /// Input file settings.
        InputFile const & inp_;
        
        /// Parallel environment.
        Parallel const & par_;
        
        /// Angular basis.
        AngularBasis ang_;
        
        /// Radial bases.
        std::vector<Bspline> const & bspline_;
};

#endif // HEX_ECS_COUPLED_SOLVER_H
