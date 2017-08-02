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

#include "hex-csrmatrix.h"
#include "hex-itersolve.h"
#include "hex-misc.h"
#include "hex-vtkfile.h"

// --------------------------------------------------------------------------------- //

#include "radial.h"

// --------------------------------------------------------------------------------- //

#include "DOMPreconditioner.h"

// --------------------------------------------------------------------------------- //

//#define DOM_DEBUG

// --------------------------------------------------------------------------------- //

DOMPreconditioner::DOMPreconditioner
(
    CommandLine  const & cmd,
    InputFile    const & inp,
    Parallel     const & par,
    AngularBasis const & ang,
    Bspline const & bspline_inner,
    Bspline const & bspline_full,
    Bspline const & bspline_panel_x,
    Bspline const & bspline_panel_y
) : NoPreconditioner
    (
        cmd, inp, par, ang,
        bspline_inner, bspline_full,
        bspline_panel_x, bspline_panel_y
    ),
    gap_(0)
{
    // nothing more to do
}

std::string DOMPreconditioner::description () const
{
    return "Domain decompositon preconditioner.";
}

void DOMPreconditioner::setup ()
{
    NoPreconditioner::setup();
}

void DOMPreconditioner::update (Real E)
{
    NoPreconditioner::update(E);
}

void DOMPreconditioner::precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const
{
    // number of initial states (right-hand sides)
    Nini_ = r[0].size() / block_rank_[0];
    
    // B-spline parameters
    int order = rad_inner().bspline().order();
    Real theta = rad_inner().bspline().ECStheta();
    
    // construct sub-domain bases
    std::vector<PanelSolution> p;
    for (int ixpanel = 0; ixpanel < cmd_->dom_x_panels; ixpanel++)
    for (int iypanel = 0; iypanel < cmd_->dom_y_panels; iypanel++)
    {
        // knot sub-sequences
        rArray rxknots, cxknots1, cxknots2;
        rArray ryknots, cyknots1, cyknots2;
        
        // calculate the knot sub-sequences
        knotSubsequence(ixpanel, cmd_->dom_x_panels, rad_inner().bspline(), rxknots, cxknots1, cxknots2);
        knotSubsequence(iypanel, cmd_->dom_y_panels, rad_inner().bspline(), ryknots, cyknots1, cyknots2);
        
        if (ixpanel == 0) cxknots1.drop();
        if (iypanel == 0) cyknots1.drop();
        
        // create the B-spline bases for the panel
        p.emplace_back
        (
            ixpanel, iypanel,
            cmd_->dom_x_panels,
            cmd_->dom_y_panels,
            order,
            theta,
            rad_inner().bspline_x(),
            rad_inner().bspline_y(),
            cxknots1, rxknots, cxknots2,
            cyknots1, ryknots, cyknots2,
            cxknots1, rxknots, cxknots2,
            cyknots1, ryknots, cyknots2,
            ang_->states().size(),
            Nini_
        );
        
        // get the panel data structure
        PanelSolution & panel = p.back();
        
        // get the current B-spline objects
        Bspline const & xspline = panel.xspline_inner;
        Bspline const & yspline = panel.yspline_inner;
        
        std::cout << std::endl << std::endl;
        std::cout << "\tPanel (" << ixpanel << "," << iypanel << ")" << std::endl;
        std::cout << "\t  x basis : " << xspline.Rmin() << " " << xspline.R1() << " " << xspline.R2() << " " << xspline.Rmax() << std::endl;
        std::cout << "\t  y basis : " << yspline.Rmin() << " " << yspline.R1() << " " << yspline.R2() << " " << yspline.Rmax() << std::endl;
        std::cout << "\t  x reknot : " << xspline.rknots() << std::endl;
        std::cout << "\t  y reknot : " << yspline.rknots() << std::endl;
        std::cout << "\t  x exclusive B-splines : " << panel.minpxspline << " ... " << panel.maxpxspline
                  << " corresponding to " << panel.minpxspline + panel.xoffset << " ... " << panel.maxpxspline + panel.xoffset << std::endl;
        std::cout << "\t  y exclusive B-splines : " << panel.minpyspline << " ... " << panel.maxpyspline
                  << " corresponding to " << panel.minpyspline + panel.yoffset << " ... " << panel.maxpyspline + panel.yoffset << std::endl;
    }
    
    // interpolate the residual into sub-domains
    splitResidual(r, p);
    
    // reset the solution
    for (cArray & Z : z)
        Z.fill(0.);
    
    // find the solution on sub-domains
    std::cout << std::endl;
    int cycles = cmd_->dom_sweeps > 0 ? cmd_->dom_sweeps : std::max(cmd_->dom_x_panels,cmd_->dom_y_panels);
    for (int cycle = 0; cycle < cycles; cycle++)
    {
        for (int ixpanel = 0; ixpanel < cmd_->dom_x_panels; ixpanel++)
        for (int iypanel = 0; iypanel < cmd_->dom_y_panels; iypanel++)
        {
            solvePanel(cycle, cycles, p, ixpanel, iypanel);
            
#ifdef DOM_DEBUG
            collectSolution(z, p);
            cArray res (z[0].size());
            A_blocks_[0].dot(1.0, z[0], 0.0, res);
            res -= r[0];
            
            {
                VTKRectGridFile vtk;
                rArray gridx = linspace(0., rad_inner().bspline_x().Rmax(), 1001);
                rArray gridy = linspace(0., rad_inner().bspline_y().Rmax(), 1001);
                cArray eval = Bspline::zip(rad_inner().bspline_x(), rad_inner().bspline_y(), res, gridx, gridy);
                vtk.setGridX(gridx);
                vtk.setGridY(gridy);
                vtk.setGridZ(rArray{0});
                vtk.appendVector2DAttribute("residual", realpart(eval), imagpart(eval));
                vtk.writePoints(format("dom-%d-%d-%d-res.vtk", cycle, ixpanel, iypanel));
            }
            {
                VTKRectGridFile vtk;
                rArray gridx = linspace<Real>(0, rad_inner().bspline_x().Nspline(), rad_inner().bspline_x().Nspline() + 1);
                rArray gridy = linspace<Real>(0, rad_inner().bspline_y().Nspline(), rad_inner().bspline_y().Nspline() + 1);
                vtk.setGridX(gridx);
                vtk.setGridY(gridy);
                vtk.setGridZ(rArray{0});
                vtk.appendVector2DAttribute("residual", realpart(res), imagpart(res));
                vtk.writeCells(format("dom-components-%d-%d-%d-res.vtk", cycle, ixpanel, iypanel));
            }
#endif
        }
#ifdef DOM_DEBUG
        collectSolution(z, p); // <-- DEBUG
#endif
    }
    
    // interpolate the solution from sub-domains
    collectSolution(z, p);
    
    std::cout << std::endl << std::endl;
    std::cout << "\t    -> DOM preconditioning done!    ";
}

void DOMPreconditioner::finish ()
{
    NoPreconditioner::finish();
}

void DOMPreconditioner::solvePanel (int cycle, int cycles, std::vector<PanelSolution> & p, int ipanel, int jpanel) const
{
    // get reference to the current panel
    PanelSolution * pCentre = &p[ipanel * cmd_->dom_y_panels + jpanel];
    
    std::cout << std::endl;
    std::cout << "\t-----------------------------------------------" << std::endl;
    std::cout << "\tSolve panel " << ipanel << " " << jpanel
              << " (bases " << format("%04x %04x", pCentre->xspline_inner.hash(), pCentre->yspline_inner.hash())
              << ", sweep " << cycle + 1 << " of " << cycles << ")"
              << std::endl;
    std::cout << "\t-----------------------------------------------" << std::endl;
    std::cout << std::endl;
    
    // Sometimes all we need to do is to mirror panel solutions that have been already found.
    // This is possible only when calculating with exchange enabled, when we have
    // access to all angular symmetries. Of course, both particles must be electrons.
    /*if (inp_->exchange and inp_->Zp == -1 and cmd_->dom_x_panels == cmd_->dom_y_panels and ipanel > jpanel)
    {
        PanelSolution * pMirror = &p[jpanel * cmd_->dom_y_panels + ipanel];
        
        std::cout << "\tMirroring (" << jpanel << "," << ipanel << ") ..." << std::endl;
        
        // for all angular components of the central and mirror panel
        for (unsigned ill = 0; ill < ang_->states().size(); ill++)
        for (unsigned illp = 0; illp < ang_->states().size(); illp++)
        {
            // get particles' angular momenta
            int l1 = ang_->states()[ill].first;
            int l2 = ang_->states()[ill].second;
            int l1p = ang_->states()[illp].first;
            int l2p = ang_->states()[illp].second;
            
            // if this is the correct symmetry, use the mirror panel
            if (l1 == l2p and l2 == l1p)
            {
                assert(pCentre->z[ill].size() == pMirror->z[illp].size());
                
                // copy data from the mirror panel and change sign
                pCentre->z[ill] = pMirror->z[illp];
                pCentre->z[ill] *= (l1 + l2 + inp_->L + ang_->S()) % 2 == 0 ? 1.0_r : -1.0_r;
                
                // mirror the elements of the solution (i.e. swap B-spline indices)
                transpose(pCentre->z[ill], pCentre->xspline_inner.Nspline(), pCentre->yspline_inner.Nspline());
            }
        }
        
        return;
    }*/
    
    // create the preconditioner object
    PreconditionerBase * prec = PreconditionerBase::Choose
    (
        cmd_->dom_preconditioner,
        *cmd_, *inp_, *par_, *ang_,
        rad_inner().bspline(),  // inner region basis
        rad_full().bspline(),   // full domain basis
        pCentre->xspline_full,  // panel x basis
        pCentre->yspline_full   // panel y basis
    );
    
    std::cout.imbue(std::locale(std::locale::classic(), new MyNumPunct));
    std::cout << "\tPanel hamiltonian size: "
              << ang_->states().size() * pCentre->xspline_full.Nspline() * pCentre->yspline_full.Nspline()
              << std::endl;
    std::cout.imbue(std::locale::classic());
    
    // construct the matrix of the equations etc.
    prec->verbose(false);
    prec->setup();
    prec->update(E_);
    
    // number of B-splines in both directions of this panel
    int Nxspline = pCentre->xspline_inner.Nspline();
    int Nyspline = pCentre->yspline_inner.Nspline();
    
    // get right-hand side and solution arrays
    cBlockArray & psi = pCentre->z;
    cBlockArray   chi = pCentre->r;
    
    // correct right-hand side from the neighbour domains
    correctSource(chi, p, ipanel, jpanel);
    
#ifdef DOM_DEBUG
    cBlockArray col (chi.size()), mul (chi.size());
    for (unsigned ill = 0; ill < chi.size(); ill++)
    {
        {
            VTKRectGridFile vtk;
            rArray gridx = linspace(pCentre->xspline_inner.Rmin(), pCentre->xspline_inner.Rmax(), 2001);
            rArray gridy = linspace(pCentre->yspline_inner.Rmin(), pCentre->yspline_inner.Rmax(), 2001);
            cArray eval = Bspline::zip(pCentre->xspline_inner, pCentre->yspline_inner, chi[ill], gridx, gridy);
            vtk.setGridX(gridx);
            vtk.setGridY(gridy);
            vtk.setGridZ(rArray{0});
            vtk.appendVector2DAttribute("rhs", realpart(eval), imagpart(eval));
            vtk.writePoints(format("dom-%d-%d-%d-chi-%d.vtk", cycle, ipanel, jpanel, ill));
        }
        {
            VTKRectGridFile vtk;
            rArray gridx = linspace<Real>(0, pCentre->xspline_inner.Nspline(), pCentre->xspline_inner.Nspline() + 1);
            rArray gridy = linspace<Real>(0, pCentre->yspline_inner.Nspline(), pCentre->yspline_inner.Nspline() + 1);
            vtk.setGridX(gridx);
            vtk.setGridY(gridy);
            vtk.setGridZ(rArray{0});
            vtk.appendVector2DAttribute("rhs", realpart(chi[ill]), imagpart(chi[ill]));
            vtk.writeCells(format("dom-components-%d-%d-%d-chi-%d.vtk", cycle, ipanel, jpanel, ill));
        }
    }
#endif
    
    // reset the solution
    for (cArray & segment : psi)
        segment.fill(0.0);
    
    // conjugate gradients callbacks
    auto apply_preconditioner = [&prec](cBlockArray const & r, cBlockArray & z) -> void
    {
        prec->precondition(r, z);
    };
    auto matrix_multiply = [&prec](cBlockArray const & p, cBlockArray & q) -> void
    {
        prec->multiply(p, q);
    };
    auto compute_norm = [](cBlockArray const & A) -> Real
    {
        Real X = 0;
        for (cArray const & a : A)
            X += a.sqrnorm();
        return std::sqrt(X);
    };
    auto scalar_product = [](cBlockArray const & A, cBlockArray const & B) -> Complex
    {
        assert(A.size() == B.size());
        Complex X = 0;
        for (std::size_t i = 0; i < A.size(); i++)
            X += (A[i] | B[i]);
        return X;
    };
    auto axby_operation = [](Complex a, cBlockArray & A, Complex b, cBlockArray const & B) -> void
    {
        assert(A.size() == B.size());
        for (std::size_t i = 0; i < A.size(); i++)
        {
            assert(A[i].size() == B[i].size());
            for (std::size_t j = 0; j < A[i].size(); j++)
            {
                A[i][j] = a * A[i][j] + b * B[i][j];
            }
        }
    };
    auto new_array = [Nxspline,Nyspline,this](std::size_t N, std::string name) -> cBlockArray
    {
        cBlockArray A (N);
        for (cArray & a : A)
            a.resize(Nxspline * Nyspline * Nini_);
        return A;
    };
    auto process_solution = [](unsigned iteration, cBlockArray const & x) -> void
    {
        // nothing
    };
    
    // solve the system
    ConjugateGradients<Complex, cBlockArray, cBlockArray&> CG;
    CG.apply_preconditioner = apply_preconditioner;
    CG.matrix_multiply      = matrix_multiply;
    CG.verbose              = true;
    CG.compute_norm         = compute_norm;
    CG.scalar_product       = scalar_product;
    CG.axby                 = axby_operation;
    CG.new_array            = new_array;
    CG.process_solution     = process_solution;
    CG.reset();
    
    std::cout << "\tPanel solution" << std::endl;
    std::cout << "\t   i | time        | residual        | min  max  avg  block precond. iter." << std::endl;
    CG.solve(chi, psi, cmd_->prec_itertol, 0, 1000);

#ifdef DOM_DEBUG
    for (unsigned ill = 0; ill < psi.size(); ill++)
    {
        {
            VTKRectGridFile vtk;
            rArray gridx = linspace(pCentre->xspline_inner.Rmin(), pCentre->xspline_inner.Rmax(), 2001);
            rArray gridy = linspace(pCentre->yspline_inner.Rmin(), pCentre->yspline_inner.Rmax(), 2001);
            cArray eval = Bspline::zip(pCentre->xspline_inner, pCentre->yspline_inner, psi[ill], gridx, gridy);
            vtk.setGridX(gridx);
            vtk.setGridY(gridy);
            vtk.setGridZ(rArray{0});
            vtk.appendVector2DAttribute("solution", realpart(eval), imagpart(eval));
            vtk.writePoints(format("dom-%d-%d-%d-psi-%d.vtk", cycle, ipanel, jpanel, ill));
        }
        {
            VTKRectGridFile vtk;
            rArray gridx = linspace<Real>(0, pCentre->xspline_inner.Nspline(), pCentre->xspline_inner.Nspline() + 1);
            rArray gridy = linspace<Real>(0, pCentre->yspline_inner.Nspline(), pCentre->yspline_inner.Nspline() + 1);
            vtk.setGridX(gridx);
            vtk.setGridY(gridy);
            vtk.setGridZ(rArray{0});
            vtk.appendVector2DAttribute("solution", realpart(psi[ill]), imagpart(psi[ill]));
            vtk.writeCells(format("dom-components-%d-%d-%d-psi-%d.vtk", cycle, ipanel, jpanel, ill));
        }
    }
#endif
    
    prec->finish();
    delete prec;
}

void DOMPreconditioner::correctSource
(
    cBlockArray & chi,
    std::vector<PanelSolution> const & panels,
    int ipanel, int jpanel
) const
{
    PanelSolution const & p = panels[ipanel * cmd_->dom_y_panels + jpanel];
    
    int order = rad_inner().bspline().order();
    
    // loop over all B-spline elements of the target panel
    for (int pi = 0; pi < p.xspline_inner.Nspline(); pi++)
    for (int pj = 0; pj < p.yspline_inner.Nspline(); pj++)
    {
        // transform position from local index space to global index space
        int i = pi + p.xoffset;
        int j = pj + p.yoffset;
        
        // determine position of the central (target) B-spline index
        bool tgtInP0 = p.minpxspline + order < pi and pi < p.maxpxspline - order and
                       p.minpyspline + order < pj and pj < p.maxpyspline - order;
        bool tgtInP  = p.minpxspline <= pi and pi <= p.maxpxspline and
                       p.minpyspline <= pj and pj <= p.maxpyspline;
        bool tgtInQ  = tgtInP and not tgtInP0;
        bool tgtInR  = p.minpxspline - order <= pi and pi <= p.maxpxspline + order and
                       p.minpyspline - order <= pj and pj <= p.maxpyspline + order and not tgtInP;
        assert(!tgtInQ || !tgtInR);
        
        // loop over all surrounding (source) B-spline elements
        if (tgtInQ or tgtInR)
        for (int pk = std::max(pi - order, 0); pk <= std::min(pi + order, p.xspline_inner.Nspline() - 1); pk++)
        for (int pl = std::max(pj - order, 0); pl <= std::min(pj + order, p.yspline_inner.Nspline() - 1); pl++)
        {
            // transform position from local index space to global index space
            int k = pk + p.xoffset;
            int l = pl + p.yoffset;
            
            // determine position of the source B-spline index
            bool srcInP0 = p.minpxspline + order < pk and pk < p.maxpxspline - order and
                           p.minpyspline + order < pl and pl < p.maxpyspline - order;
            bool srcInP  = p.minpxspline <= pk and pk <= p.maxpxspline and
                           p.minpyspline <= pl and pl <= p.maxpyspline;
            bool srcInQ  = srcInP and not srcInP0;
            bool srcInR  = p.minpxspline - order <= pk and pk <= p.maxpxspline + order and
                           p.minpyspline - order <= pl and pl <= p.maxpyspline + order and not srcInP;
            assert(!srcInQ || !srcInR);
            
            // loop over all neighbour panels
            for (int kpanel = std::max(0, ipanel - 1); kpanel <= std::min(cmd_->dom_x_panels - 1, ipanel + 1); kpanel++)
            for (int lpanel = std::max(0, jpanel - 1); lpanel <= std::min(cmd_->dom_y_panels - 1, jpanel + 1); lpanel++)
            {
                // get neighbour panel info structure
                PanelSolution const & n = panels[kpanel * cmd_->dom_y_panels + lpanel];
                
                // transform source B-spline index to neighbour panel index space
                int nk = k - n.xoffset;
                int nl = l - n.yoffset;
                
                // locate the source B-spline in the neighbour panel
                bool nP = n.minpxspline <= nk and nk <= n.maxpxspline and
                          n.minpyspline <= nl and nl <= n.maxpyspline;
                bool nR = n.minpxspline - order <= nk and nk <= n.maxpxspline + order and
                          n.minpyspline - order <= nl and nl <= n.maxpyspline + order and not nP;
                
                // skip this panel if the source B-spline is not included in this panel
                if (not nP and not nR)
                    continue;
                
                // update the field source elements
                for (unsigned ill  = 0; ill  < chi.size(); ill ++)
                for (unsigned illp = 0; illp < chi.size(); illp++)
                {
                    Complex Aijkl = couplingMatrixElement(ill, illp, i, j, k, l);
                    
                    // for all initial states (right-hand sides)
                    for (int ini = 0; ini < Nini_; ini++)
                    {
                        std::size_t offset  = chi[ill ].size() * ini / Nini_;
                        std::size_t poffset = p.z[illp].size() * ini / Nini_;
                        std::size_t noffset = n.z[illp].size() * ini / Nini_;
                        
                        Complex & rhs = chi[ill ][pi * p.yspline_inner.Nspline() + pj + offset];
                        Complex   own = p.z[illp][pk * p.yspline_inner.Nspline() + pl + poffset];
                        Complex   nbr = n.z[illp][nk * n.yspline_inner.Nspline() + nl + noffset];
                        
                        // update the field source with the neighbour field (use only incoming component)
                        if (tgtInQ and srcInR and nP)
                            rhs -= Aijkl * (nbr - own);
                        
                        // update the field source with the neighbour field (already just the incoming component)
                        if (tgtInR and srcInQ and nR)
                            rhs += Aijkl * nbr;
                    }
                }
            }
        }
    }
}

Complex DOMPreconditioner::couplingMatrixElement
(
    int ill, int illp,
    int i, int j, int k, int l
) const
{
    // matrix element coupling B-splines across the panels boundary
    Complex Aijkl = 0;
    
    // angular diagonal contribution
    if (ill == illp)
    {
        int l1 = ang_->states()[ill].first;
        int l2 = ang_->states()[ill].second;
        
        Aijkl += E_ * rad_inner().S_x()(i,k) * rad_inner().S_y()(j,l);
        
        Aijkl -= 0.5_r * rad_inner().D_x()(i,k) * rad_inner().S_y()(j,l);
        Aijkl -= 0.5_r * rad_inner().S_x()(i,k) * rad_inner().D_y()(j,l);
        
        Aijkl -= (0.5_r * l1 * (l1 + 1)) * rad_inner().Mm2_x()(i,k) * rad_inner().S_y()(j,l);
        Aijkl -= (0.5_r * l2 * (l2 + 1)) * rad_inner().S_x()(i,k) * rad_inner().Mm2_y()(j,l);
        
        Aijkl -= rad_inner().Mm1_x()(i,k) * rad_inner().S_y()(j,l) * (inp_->Za * -1.0_r);
        Aijkl -= rad_inner().S_x()(i,k) * rad_inner().Mm1_y()(j,l) * (inp_->Za * inp_->Zp);
    }
    
    // angular off-diagonal contribution
    for (int lambda = 0; lambda <= rad_inner().maxlambda(); lambda++)
    {
        if (ang_->f(ill, illp, lambda) != 0)
        {
            Aijkl += inp_->Zp * ang_->f(ill, illp, lambda) * rad_inner().computeR(lambda,i,j,k,l);
        }
    }
    
    return Aijkl;
}

void DOMPreconditioner::knotSubsequence (int i, int n, Bspline const & bspline, rArray & rknots, rArray & cknots1, rArray & cknots2) const
{
    // boundaries of the original real grid
    int iR1 = bspline.iR1();
    int iR2 = bspline.iR2();
    
    // boundaries of the sub-grid
    int iRa = iR1 + (iR2 - iR1) * i / n;
    int iRb = iR1 + (iR2 - iR1) * (i + 1) / n;
    
    // Add some more real knots to both sides. The reason for this is to conserve
    // the total number of B-splines. For order = 2 we need at least this situation
    // at the end of a grid
    //                   |___
    //                ___/   \___
    //    ...+...+...+...+...+...+...+...+~~~
    //                   |
    // and this situation at the beginning of the next panel's grid
    //                   |    ___
    //                   |___/   \___
    //        ~~~+...+...+...+...+...+...+...+...
    //                   |
    // where the '+' symbols denote the knots, the vertical line is the panel separation
    // knot, the dots are real intervals and '~' marks the complex grid part. Displayed
    // are the last B-spline owned by first grid and the first B-spline owned by the
    // following grid. In the global basis the latter immediately follows the former.
    // To allow proper continuity and coupling between the grids it was necessary to
    // add '2*order' knots after the separation knot in the first grid, and 'order' knots
    // before the separation knot in the connected grid. This results in 'order'
    // additional B-spline elements both in the former and latter grid. Of course,
    // then there are also arbitrary number of complex B-splines in the damping area.
    // For obscure numerical reasons there can be also any number of additional real knots
    // apart from those pictured above.
    //
    // In brief:
    //    - Unless this is to be the first grid, add 1*order knots (+ arbitrary gap) at the beginning of the grid.
    //    - Unless this is to be the last grid, add 2*order knots (+ arbitrary gap) at the end of the grid.
    
    iRa = std::max(iR1, iRa - gap_ - 1 * bspline.order());
    iRb = std::min(iR2, iRb + gap_ + 2 * bspline.order());
    
    // new real knot sub-sequence
    rknots = inp_->rknots.slice(iRa, iRb + 1);
    
    // trailing complex knots
    cknots2 = inp_->cknots + rknots.back();
    
    // leading complex knots (need to reverse order to preserve interval stretching)
    cknots1 = -inp_->cknots;
    std::reverse(cknots1.begin(), cknots1.end());
    cknots1 += rknots.front();
}

void DOMPreconditioner::splitResidual (cBlockArray const & r, std::vector<PanelSolution> & panels) const
{
    for (int ixpanel = 0; ixpanel < cmd_->dom_x_panels; ixpanel++)
    for (int iypanel = 0; iypanel < cmd_->dom_y_panels; iypanel++)
    for (unsigned ill = 0; ill < r.size(); ill++)
    {
        PanelSolution & p = panels[ixpanel * cmd_->dom_y_panels + iypanel];
        
        std::size_t Nfxspline = rad_inner().bspline_x().Nspline();
        std::size_t Nfyspline = rad_inner().bspline_y().Nspline();
        std::size_t Npxspline = p.xspline_inner.Nspline();
        std::size_t Npyspline = p.yspline_inner.Nspline();
        
        p.r[ill].resize(Npxspline * Npyspline * Nini_);
        p.r[ill].fill(0.);
        
        // for all elements of the residual
        for (int ixspline = 0; ixspline < (int)Nfxspline; ixspline++)
        for (int iyspline = 0; iyspline < (int)Nfyspline; iyspline++)
        {
            // get corresponding indices in panel
            int pxspline = ixspline - p.xoffset;
            int pyspline = iyspline - p.yoffset;
            
            // skip elements that do not belong to this panel
            if (pxspline < p.minpxspline or pxspline > p.maxpxspline or pyspline < p.minpyspline or pyspline > p.maxpyspline)
                continue;
            
            // copy the element to the panel residual
            for (int ini = 0; ini < Nini_; ini++)
            {
                p.r[ill][(ini * Npxspline + pxspline) * Npyspline + pyspline]
                    = r[ill][(ini * Nfxspline + ixspline) * Nfyspline + iyspline];
            }
        }
        
#ifdef DOM_DEBUG
        {
            VTKRectGridFile vtk;
            rArray gridx = linspace(p.xspline_inner.Rmin(), p.xspline_inner.Rmax(), 2001);
            rArray gridy = linspace(p.yspline_inner.Rmin(), p.yspline_inner.Rmax(), 2001);
            cArray eval = Bspline::zip(p.xspline_inner, p.yspline_inner, p.r[ill], gridx, gridy);
            vtk.setGridX(gridx);
            vtk.setGridY(gridy);
            vtk.setGridZ(rArray{0});
            vtk.appendVector2DAttribute("split", realpart(eval), imagpart(eval));
            vtk.writePoints(format("dom-%d-%d-%d-split-%d.vtk", 0, ixpanel, iypanel, ill));
        }
        {
            VTKRectGridFile vtk;
            rArray gridx = linspace<Real>(0, p.xspline_inner.Nspline(), p.xspline_inner.Nspline() + 1) + p.xoffset;
            rArray gridy = linspace<Real>(0, p.yspline_inner.Nspline(), p.yspline_inner.Nspline() + 1) + p.yoffset;
            vtk.setGridX(gridx);
            vtk.setGridY(gridy);
            vtk.setGridZ(rArray{0});
            vtk.appendVector2DAttribute("split", realpart(p.r[ill]), imagpart(p.r[ill]));
            vtk.writeCells(format("dom-components-%d-%d-%d-split-%d.vtk", 0, ixpanel, iypanel, ill));
        }
#endif
    }
}

void DOMPreconditioner::collectSolution (cBlockArray & z, std::vector<PanelSolution> & panels) const
{
    for (unsigned ill = 0; ill < z.size(); ill++)
    {
        z[ill].fill(0.0);
        
        for (int ixpanel = 0; ixpanel < cmd_->dom_x_panels; ixpanel++)
        for (int iypanel = 0; iypanel < cmd_->dom_y_panels; iypanel++)
        {
            PanelSolution & p = panels[ixpanel * cmd_->dom_y_panels + iypanel];
            
            std::size_t Nfxspline = rad_inner().bspline_x().Nspline();
            std::size_t Nfyspline = rad_inner().bspline_y().Nspline();
            std::size_t Npxspline = p.xspline_inner.Nspline();
            std::size_t Npyspline = p.yspline_inner.Nspline();
            
            // for all elements of the panel solution
            for (int pxspline = p.minpxspline; pxspline <= p.maxpxspline; pxspline++)
            for (int pyspline = p.minpyspline; pyspline <= p.maxpyspline; pyspline++)
            {
                // get corresponding indices in global basis
                int ixspline = pxspline + p.xoffset;
                int iyspline = pyspline + p.yoffset;
                
                // update collected solution
                for (int ini = 0; ini < Nini_; ini++)
                {
                    z[ill][(ini * Nfxspline + ixspline) * Nfyspline + iyspline]
                        = p.z[ill][(ini * Npxspline + pxspline) * Npyspline + pyspline];
                }
            }
        }
        
        
#ifdef DOM_DEBUG
        static int n = 0;
        {
            VTKRectGridFile vtk;
            rArray grid = linspace(0., rad_inner().bspline().Rmax(), 2001);
            cArray eval = Bspline::zip(rad_inner().bspline(), rad_inner().bspline(), z[ill], grid, grid);
            vtk.setGridX(grid);
            vtk.setGridY(grid);
            vtk.setGridZ(rArray{0});
            vtk.appendVector2DAttribute("combined_wavefunction", realpart(eval), imagpart(eval));
            vtk.writePoints(format("combined-%d-%d.vtk", ill, n++));
        }
        {
            VTKRectGridFile vtk;
            rArray grid = linspace<Real>(0, rad_inner().bspline().Nspline(), rad_inner().bspline().Nspline() + 1);
            vtk.setGridX(grid);
            vtk.setGridY(grid);
            vtk.setGridZ(rArray{0});
            vtk.appendVector2DAttribute("combined_wavefunction", realpart(z[ill]), imagpart(z[ill]));
            vtk.writeCells(format("components-%d-%d.vtk", ill, n++));
        }
#endif
    }
}

DOMPreconditioner::PanelSolution::PanelSolution
(
    int ix, int iy, int xpanels, int ypanels,
    int order,
    Real theta,
    Bspline const & xspline, Bspline const & yspline,
    rArray cxspline1_inner, rArray rxspline_inner, rArray cxspline2_inner,
    rArray cyspline1_inner, rArray ryspline_inner, rArray cyspline2_inner,
    rArray cxspline1_full,  rArray rxspline_full,  rArray cxspline2_full,
    rArray cyspline1_full,  rArray ryspline_full,  rArray cyspline2_full,
    int Nang, int Nini
) : xspline_inner (order, theta, cxspline1_inner, rxspline_inner, cxspline2_inner),
    yspline_inner (order, theta, cyspline1_inner, ryspline_inner, cyspline2_inner),
    xspline_full (order, theta, cxspline1_full, rxspline_full, cxspline2_full),
    yspline_full (order, theta, cyspline1_full, ryspline_full, cyspline2_full),
    r (Nang), z (Nang), ixpanel(ix), iypanel(iy)
{
    // allocate memory for residual and solution
    for (int i = 0; i < Nang; i++) r[i].resize(xspline_inner.Nspline() * yspline_inner.Nspline() * Nini);
    for (int i = 0; i < Nang; i++) z[i].resize(xspline_inner.Nspline() * yspline_inner.Nspline() * Nini);
    
    // calculate panel-to-full real basis offset
    xoffset = xspline.knot(xspline_inner.R1()) - xspline_inner.iR1();
    yoffset = yspline.knot(yspline_inner.R1()) - yspline_inner.iR1();
    
    int gap = 0;
    
    // this panel's B-splines that have a counterpart in the global basis; the rest is specific to this panel
    minpxspline = (ixpanel == 0 ? 0 : xspline_inner.iR1() + gap + order);
    minpyspline = (iypanel == 0 ? 0 : yspline_inner.iR1() + gap + order);
    maxpxspline = (ixpanel + 1 == xpanels ? xspline_inner.Nspline() - 1 : xspline_inner.iR2() - 2*order - 1 - gap);
    maxpyspline = (iypanel + 1 == ypanels ? yspline_inner.Nspline() - 1 : yspline_inner.iR2() - 2*order - 1 - gap);
}

// --------------------------------------------------------------------------------- //

addClassToParentRunTimeSelectionTable(PreconditionerBase, DOMPreconditioner)

// --------------------------------------------------------------------------------- //
