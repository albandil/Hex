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

#ifndef HEX_DOM_PRECONDITIONER_H
#define HEX_DOM_PRECONDITIONER_H

// --------------------------------------------------------------------------------- //

#include <array>
#include <vector>

// --------------------------------------------------------------------------------- //

#include "hex-csrmatrix.h"
#include "hex-densematrix.h"

// --------------------------------------------------------------------------------- //

#include "NoPreconditioner.h"

// --------------------------------------------------------------------------------- //

/**
 * @brief Domain decomposition preconditioner.
 * 
 * This preconditioner splits the simulated domain into several non-overlapping subdomains
 * and solves the smaller per-subdomain problems in a sequence. The numbering of
 * the domains is as follows:
 * 
 * @verbatim
 * 
 * ┌──┬──┬──┐
 * │2 │5 │8 │
 * ├──┼──┼──┤
 * │1 │4 │7 │
 * ├──┼──┼──┤
 * │0 │3 │6 │
 * └──┴──┴──┘
 * 
 * @endverbatim
 */
class DOMPreconditioner : public NoPreconditioner
{
    public:
        
        // run-time selection mechanism
        preconditionerRunTimeSelectionDefinitions(DOMPreconditioner, "DOM")
        
        // default constructor needed by the RTS mechanism
        DOMPreconditioner () {}
        
        // constructor
        DOMPreconditioner
        (
            CommandLine  const & cmd,
            InputFile    const & inp,
            Parallel     const & par,
            AngularBasis const & ang,
            Bspline const & bspline_inner,
            Bspline const & bspline_full,
            Bspline const & bspline_panel_x,
            Bspline const & bspline_panel_y
        );
        
        // description of the preconditioner
        virtual std::string description () const;
        
        // reuse parent definitions
        using NoPreconditioner::rhs;
        using NoPreconditioner::multiply;
        
        // declare own definitions
        virtual void setup ();
        virtual void update (Real E);
        virtual void precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const;
        virtual void finish ();
    
    protected:
        
        // neighbour panels
        enum Neighbours
        {
            Left   = 0,
            Right  = 1,
            Down   = 2,
            Up     = 3,
            nNbrs  = 4
        };
        
        // get reverse direction for a given direction
        int rev (int dir) const
        {
            switch (dir)
            {
                case Left  : return Right;
                case Right : return Left;
                case Down  : return Up;
                case Up    : return Down;
            };
            
            return nNbrs;
        }
        
        // solutions on the sub-domains
        class PanelSolution
        {
            public:
                
                PanelSolution
                (
                    int ixpanel, int iypanel,
                    int xpanels, int ypanels,
                    int order,
                    Real theta,
                    Bspline const & xspline, Bspline const & yspline,
                    rArray cxspline1_inner, rArray rxspline_inner, rArray cxspline2_inner,
                    rArray cyspline1_inner, rArray ryspline_inner, rArray cyspline2_inner,
                    rArray cxspline1_full,  rArray rxspline_full,  rArray cxspline2_full,
                    rArray cyspline1_full,  rArray ryspline_full,  rArray cyspline2_full,
                    int Nang
                );
                
                Bspline xspline_inner;  // inner x-axis B-spline basis
                Bspline yspline_inner;  // inner y-axis B-spline basis
                
                Bspline xspline_full;   // full x-axis B-spline basis
                Bspline yspline_full;   // full y-axis B-spline basis
                
                cBlockArray r;  // original source
                cBlockArray z;  // solution
                
                int ixpanel;   // which panel (x-dir)
                int iypanel;   // which panel (y-dir)
                
                int xoffset;   // x-offset of the real basis of panel
                int yoffset;   // y-offset of the real basis of panel
                
                int minpxspline, maxpxspline; // B-splines that have a counterpart in the global basis (x-dir)
                int minpyspline, maxpyspline; // B-splines that have a counterpart in the global basis (y-dir)
        };
        
        // find solution on a sub-domain
        void solvePanel
        (
            int cycle, int cycles,
            std::vector<PanelSolution> & p,
            int i, int j
        ) const;
        
        // add neighbour field interfaces
        void correctSource
        (
            cBlockArray & chi,
            std::vector<CooMatrix<LU_int_t,Complex>> & G,
            std::vector<PanelSolution> const & panels,
            int ipanel, int jpanel
        ) const;
        
        // evaluate matrix element
        Complex couplingMatrixElement
        (
            int ill, int illp,
            int i, int j, int k, int l
        ) const;
        
        // get knot sub-sequences
        void knotSubsequence
        (
            int ipanel,
            int npanels,
            Bspline const & bspline,
            rArray & rknots,
            rArray & cknots1,
            rArray & cknots2
        ) const;
        
        // interpolate residual to sub-domains
        void splitResidual (cBlockArray const & r, std::vector<PanelSolution> & p) const;
        
        // interpolate solution from sub-domains
        void collectSolution (cBlockArray & z, std::vector<PanelSolution> & p) const;
        
        // gap of real knots between the panel seam and the complex absorption layer
        int gap_;
};

// --------------------------------------------------------------------------------- //

#endif
