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

#include "radial.h"

// --------------------------------------------------------------------------------- //

#include "DOMPreconditioner.h"

// --------------------------------------------------------------------------------- //

DOMPreconditioner::DOMPreconditioner
(
    Parallel const & par,
    InputFile const & inp,
    AngularBasis const & ll,
    Bspline const & bspline_inner,
    Bspline const & bspline_full,
    CommandLine const & cmd
) : NoPreconditioner(par, inp, ll, bspline_inner, bspline_full, cmd)
{
    
}

std::string DOMPreconditioner::description () const
{
    return "Domain decompositon preconditioner.";
}

void DOMPreconditioner::setup ()
{
    
}

void DOMPreconditioner::update (Real E)
{
    
}

void DOMPreconditioner::precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const
{
    int xpanels = 2, ypanels = 2;
    int npanels = xpanels * ypanels;
    
    // solutions on individual sub-domains
    std::vector<PanelSolution> p (npanels);
    
    // split the residual among sub-domains
    for (int ixpanel = 0; ixpanel < xpanels; ixpanel++)
    for (int iypanel = 0; iypanel < ypanels; iypanel++)
        /* TODO */;
    
    // find the solution on sub-domains
    for (int ixpanel = 0; ixpanel < xpanels; ixpanel++)
    for (int iypanel = 0; iypanel < ypanels; iypanel++)
        solvePanel(xpanels, p, ixpanel, iypanel);
    
    // collect the solution from sub-domains
    for (int ixpanel = 0; ixpanel < xpanels; ixpanel++)
    for (int iypanel = 0; iypanel < ypanels; iypanel++)
        /* TODO */;
}

void DOMPreconditioner::finish ()
{
    
}

void DOMPreconditioner::solvePanel (int n, std::vector<PanelSolution> & p, int i, int j) const
{
    // create B-spline bases for this panel
    Bspline xspline_inner
    (
        rad().bspline_inner_x().order(),
        rad().bspline_inner_x().ECStheta(),
        rad().bspline_inner_x().cknots1(),
        rad().bspline_inner_x().rknots(),
        rad().bspline_inner_x().cknots2()
    );
    Bspline xspline_full
    (
        rad().bspline_full_x().order(),
        rad().bspline_full_x().ECStheta(),
        rad().bspline_full_x().cknots1(),
        rad().bspline_full_x().rknots(),
        rad().bspline_full_x().cknots2()
    );
    Bspline yspline_inner
    (
        rad().bspline_inner_y().order(),
        rad().bspline_inner_y().ECStheta(),
        rad().bspline_inner_y().cknots1(),
        rad().bspline_inner_y().rknots(),
        rad().bspline_inner_y().cknots2()
    );
    Bspline yspline_full
    (
        rad().bspline_full_y().order(),
        rad().bspline_full_y().ECStheta(),
        rad().bspline_full_y().cknots1(),
        rad().bspline_full_y().rknots(),
        rad().bspline_full_y().cknots2()
    );
    
    // calculate radial integrals
    RadialIntegrals rint
    (
        xspline_inner,
        xspline_full,
        yspline_inner,
        yspline_full,
        rad().maxlambda()
    );
    rint.setupOneElectronIntegrals(cmd_, par_);
    rint.setupTwoElectronIntegrals(cmd_, par_);
    
    // contruct the surrogate source
    // TODO
    
    // construct the matrix of the equations
    // TODO
    
    // solve the system with the combined right-hand side
    // TODO
    
    // evaluate the outgoing field by subtraction
    // TODO
}

// --------------------------------------------------------------------------------- //

addClassToParentRunTimeSelectionTable(PreconditionerBase, DOMPreconditioner)

// --------------------------------------------------------------------------------- //
