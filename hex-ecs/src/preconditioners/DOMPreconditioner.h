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

// --------------------------------------------------------------------------------- //

#include "NoPreconditioner.h"

// --------------------------------------------------------------------------------- //

/**
 * @brief Multiplicative Schwarz domain decomposition preconditioner.
 * 
 * This preconditioner splits the simulated domain into several overlapping subdomains
 * and solves the smaller per-subdomain problems in a sequence.
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
                    int order,
                    Real theta,
                    Bspline const & xspline, Bspline const & yspline,
                    rArray cxspline1_inner, rArray rxspline_inner, rArray cxspline2_inner,
                    rArray cyspline1_inner, rArray ryspline_inner, rArray cyspline2_inner,
                    rArray cxspline1_full,  rArray rxspline_full,  rArray cxspline2_full,
                    rArray cyspline1_full,  rArray ryspline_full,  rArray cyspline2_full,
                    int Nang
                );
                
                bool mapToPanel
                (
                    unsigned   ixspline, unsigned   iyspline,
                    unsigned & pxspline, unsigned & pyspline
                ) const;
                
                bool mapFromPanel
                (
                    unsigned & ixspline, unsigned & iyspline,
                    unsigned   pxspline, unsigned   pyspline
                ) const;
                
                Bspline xspline_inner;  // inner x-axis B-spline basis
                Bspline yspline_inner;  // inner y-axis B-spline basis
                
                Bspline xspline_full;   // full x-axis B-spline basis
                Bspline yspline_full;   // full y-axis B-spline basis
                
                CsrMatrix<LU_int_t,Complex> SaF, SbF;   // overlaps of panel and full basis
                CsrMatrix<LU_int_t,Complex> Saa, Sbb;   // panel B-spline self-overlap matrices
                
                std::shared_ptr<LUft> lu_Saa, lu_Sbb;   // LU decomposition of the panel overlaps
                
                cBlockArray r;  // original source
                cBlockArray z;  // solution
                
                std::array<cBlockArray,nNbrs> ssrc;  // surrogate sources from neighbour panels
                std::array<cBlockArray,nNbrs> outf;  // outgoing field to neighbour panels
                
                unsigned xoffset;   // x-offset of the real basis of panel
                unsigned yoffset;   // y-offset of the real basis of panel
        };
        
        // find solution on a sub-domain
        void solvePanel (int n, std::vector<PanelSolution> & p, int i, int j) const;
        
        // construct the surrogate source for panel's boundary
        void surrogateSource (PanelSolution * panel, int direction, PanelSolution * neighbour) const;
        
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
};

// --------------------------------------------------------------------------------- //

#endif
