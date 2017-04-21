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

#include <csignal>

// --------------------------------------------------------------------------------- //

#include "radial.h"

// --------------------------------------------------------------------------------- //

#include "DOMPreconditioner.h"

// --------------------------------------------------------------------------------- //

#define DOM_DEBUG

/// DEBUG : fixed division into sub-panels
int xpanels = 2, ypanels = 2;

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
    gap_(2*bspline_inner.order())
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
    // B-spline parameters
    int order = rad_inner().bspline().order();
    Real theta = rad_inner().bspline().ECStheta();
    
    // construct sub-domain bases
    std::vector<PanelSolution> p;
    for (int ixpanel = 0; ixpanel < xpanels; ixpanel++)
    for (int iypanel = 0; iypanel < ypanels; iypanel++)
    {
        // knot sub-sequences
        rArray rxknots, cxknots1, cxknots2;
        rArray ryknots, cyknots1, cyknots2;
        
        // calculate the knot sub-sequences
        knotSubsequence(ixpanel, xpanels, rad_inner().bspline(), rxknots, cxknots1, cxknots2);
        knotSubsequence(iypanel, ypanels, rad_inner().bspline(), ryknots, cyknots1, cyknots2);
        
        if (ixpanel == 0) cxknots1.drop();
        if (iypanel == 0) cyknots1.drop();
        
        // create the B-spline bases for the panel
        p.emplace_back
        (
            ixpanel, iypanel,
            order,
            theta,
            rad_inner().bspline_x(),
            rad_inner().bspline_y(),
            cxknots1, rxknots, cxknots2,
            cyknots1, ryknots, cyknots2,
            cxknots1, rxknots, cxknots2,
            cyknots1, ryknots, cyknots2,
            ang_->states().size()
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
    
    // find the solution on sub-domains
    std::cout << std::endl;
    int cycles = 2;
    for (int cycle = 0; cycle < cycles; cycle++)
    for (int ixpanel = 0; ixpanel < xpanels; ixpanel++)
    for (int iypanel = 0; iypanel < ypanels; iypanel++)
    {
        solvePanel(cycle, cycles, p, ixpanel, iypanel);
        collectSolution(z, p);
#ifdef DOM_DEBUG
        
        cArray res (z[0].size());
        A_blocks_[0].dot(1.0, z[0], 0.0, res);
        res -= r[0];
        
        std::ofstream ofs (format("dom-%d-%d-%d-res.vtk", cycle, ixpanel, iypanel));
        writeVTK_points
        (
            ofs,
            Bspline::zip
            (
                rad_inner().bspline_x(),
                rad_inner().bspline_y(),
                res,
                linspace(0., rad_inner().bspline_x().Rmax(), 1001),
                linspace(0., rad_inner().bspline_y().Rmax(), 1001)
            ),
            linspace(0., rad_inner().bspline_x().Rmax(), 1001),
            linspace(0., rad_inner().bspline_y().Rmax(), 1001),
            rArray{ 0. }
        );
        ofs.close();
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
    std::cout << std::endl;
    std::cout << "\t-------------------------------------------" << std::endl;
    std::cout << "\tSolve panel " << ipanel << " " << jpanel << " (cycle " << cycle + 1 << " of " << cycles << ")" << std::endl;
    std::cout << "\t-------------------------------------------" << std::endl;
    std::cout << std::endl;
    
    // get reference to the current panel
    PanelSolution * pCentre = &p[ipanel * ypanels + jpanel];
    
    // create the preconditioner object
    PreconditionerBase * prec = PreconditionerBase::Choose
    (
        "ILU", // TODO : Make (wisely) runtime selectable.
        *cmd_, *inp_, *par_, *ang_,
        rad_inner().bspline(),  // inner region basis
        rad_full().bspline(),   // full domain basis
        pCentre->xspline_full,  // panel x basis
        pCentre->yspline_full   // panel y basis
    );
    
    std::cout << "\tPanel hamiltonian size: "
              << ang_->states().size() * pCentre->xspline_full.Nspline() * pCentre->yspline_full.Nspline()
              << std::endl;
    
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
        rArray grid_x = linspace(pCentre->xspline_inner.Rmin(), pCentre->xspline_inner.Rmax(), 1001);
        rArray grid_y = linspace(pCentre->yspline_inner.Rmin(), pCentre->yspline_inner.Rmax(), 1001);
        std::ofstream ofs (format("dom-%d-%d-%d-chi-%d.vtk", cycle, ipanel, jpanel, ill));
        writeVTK_points
        (
            ofs,
            Bspline::zip
            (
                pCentre->xspline_inner,
                pCentre->yspline_inner,
                chi[ill],
                grid_x,
                grid_y
            ),
            grid_x,
            grid_y,
            rArray{ 0. }
        );
        ofs.close();
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
    auto new_array = [Nxspline,Nyspline](std::size_t N, std::string name) -> cBlockArray
    {
        cBlockArray A (N);
        for (cArray & a : A)
            a.resize(Nxspline * Nyspline);
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
    CG.solve(chi, psi, cmd_->itertol, 0, 1000);
    
#ifdef DOM_DEBUG
    int order = pCentre->xspline_inner.order();
    std::ofstream ofs;
    ofs.open(format("psi-int-%d-%d-%d-%d.txt", cycle, ipanel, jpanel, 0));
    for (int iy = 0; iy < pCentre->yspline_inner.Nspline(); iy++)
    {
        for (int ix = 0; ix < pCentre->xspline_inner.Nspline(); ix++)
        {
            if (ix >= pCentre->minpxspline and ix <= pCentre->maxpxspline and iy >= pCentre->minpyspline and iy <= pCentre->maxpyspline)
                ofs << psi[0][ix * pCentre->yspline_inner.Nspline() + iy].real() << "\t";
            else
                ofs << 0 << "\t";
        }
        ofs << std::endl;
    }
    ofs.close();
    ofs.open(format("psi-cpx-%d-%d-%d-%d.txt", cycle, ipanel, jpanel, 0));
    for (int iy = 0; iy < pCentre->yspline_inner.Nspline(); iy++)
    {
        for (int ix = 0; ix < pCentre->xspline_inner.Nspline(); ix++)
        {
            if (ix < pCentre->xspline_inner.iR1() or ix >= pCentre->xspline_inner.iR2() - order or
                iy < pCentre->yspline_inner.iR1() or iy >= pCentre->yspline_inner.iR2() - order)
                ofs << psi[0][ix * pCentre->yspline_inner.Nspline() + iy].real() << "\t";
            else
                ofs << 0 << "\t";
        }
        ofs << std::endl;
    }
    ofs.close();
    ofs.open(format("psi-gap-%d-%d-%d-%d.txt", cycle, ipanel, jpanel, 0));
    for (int iy = 0; iy < pCentre->yspline_inner.Nspline(); iy++)
    {
        for (int ix = 0; ix < pCentre->xspline_inner.Nspline(); ix++)
        {
            if (ix >= pCentre->minpxspline and ix <= pCentre->maxpxspline and iy >= pCentre->minpyspline and iy <= pCentre->maxpyspline)
                ofs << 0 << "\t";
            else if (ix < pCentre->xspline_inner.iR1() or ix >= pCentre->xspline_inner.iR2() - order or
                iy < pCentre->yspline_inner.iR1() or iy >= pCentre->yspline_inner.iR2() - order)
                ofs << 0 << "\t";
            else
                ofs << psi[0][ix * pCentre->yspline_inner.Nspline() + iy].real() << "\t";
        }
        ofs << std::endl;
    }
    ofs.close();
    
    rArray grid_x = linspace(pCentre->xspline_inner.Rmin(), pCentre->xspline_inner.Rmax(), 1001);
    rArray grid_y = linspace(pCentre->yspline_inner.Rmin(), pCentre->yspline_inner.Rmax(), 1001);
    for (unsigned ill = 0; ill < psi.size(); ill++)
    {
        std::ofstream ofs (format("dom-%d-%d-%d-psi-%d.vtk", cycle, ipanel, jpanel, ill));
        writeVTK_points
        (
            ofs,
            Bspline::zip
            (
                pCentre->xspline_inner,
                pCentre->yspline_inner,
                psi[ill],
                grid_x,
                grid_y
            ),
            grid_x,
            grid_y,
            rArray{ 0. }
        );
        ofs.close();
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
    PanelSolution const & p = panels[ipanel * ypanels + jpanel];
    
    int order = rad_inner().bspline().order();
    
    // loop over all potential neighbour panels (skip self and edge-adjacent neighbours)
    for (int kpanel = std::max(0, ipanel - 1); kpanel <= std::min(xpanels - 1, ipanel + 1); kpanel++)
    for (int lpanel = std::max(0, jpanel - 1); lpanel <= std::min(ypanels - 1, jpanel + 1); lpanel++)
    if ((kpanel != ipanel) or (lpanel != jpanel))
    {
        PanelSolution const & n = panels[kpanel * ypanels + lpanel];
        
        // shared B-spline bounds
        int bxmin = kpanel == ipanel ? p.xoffset + 0
            : std::max(p.xspline_inner.iR1() + p.xoffset, n.xspline_inner.iR1() + n.xoffset);
        int bymin = lpanel == jpanel ? p.yoffset + 0
            : std::max(p.yspline_inner.iR1() + p.yoffset, n.yspline_inner.iR1() + n.yoffset);
        int bxmax = kpanel == ipanel ? p.xoffset + p.xspline_inner.Nspline() - 1
            : std::min(p.xspline_inner.iR2() - order - 1 + p.xoffset, n.xspline_inner.iR2() - order - 1 + n.xoffset);
        int bymax = lpanel == jpanel ? p.yoffset + p.yspline_inner.Nspline() - 1
            : std::min(p.yspline_inner.iR2() - order - 1 + p.yoffset, n.yspline_inner.iR2() - order - 1 + n.yoffset);
        
        // for all (target) B-splines shared by these two panels
        for (int i = bxmin; i <= bxmax; i++)
        for (int j = bymin; j <= bymax; j++)
        {
            // central panel local indices
            int pi = i - p.xoffset;
            int pj = j - p.yoffset;
            
            // is this B-spline owned by the central panel?
            bool centraltgt =
                (lpanel == jpanel and p.minpxspline <= pi and pi <= p.maxpxspline) or
                (kpanel == ipanel and p.minpyspline <= pj and pj <= p.maxpyspline);
            
            // for all (source) B-splines shared by these two panels
            for (int k = std::max(bxmin, i - order); k <= std::min(bxmax, i + order); k++)
            for (int l = std::max(bymin, j - order); l <= std::min(bymax, j + order); l++)
            {
                // central panel local indices
                int pk = k - p.xoffset;
                int pl = l - p.yoffset;
                
                // neighbour panel local indices
                int nk = k - n.xoffset;
                int nl = l - n.yoffset;
                
                // is this B-spline owned by the central panel?
                bool centralsrc = 
                    (lpanel == jpanel and p.minpxspline <= pk and pk <= p.maxpxspline) or
                    (kpanel == ipanel and p.minpyspline <= pl and pl <= p.maxpyspline);
                
                // update the source
                if (centralsrc xor centraltgt)
                for (unsigned ill  = 0; ill  < chi.size(); ill ++)
                for (unsigned illp = 0; illp < chi.size(); illp++)
                {
                    Complex coupl = couplingMatrixElement(ill, illp, i, j, k, l);
                    Complex & rhs = chi[ill ][pi * p.yspline_inner.Nspline() + pj];
                    Complex  lsrc = p.z[illp][pk * p.yspline_inner.Nspline() + pl];
                    Complex   src = n.z[illp][nk * n.yspline_inner.Nspline() + nl];
                    
                    if (centraltgt) rhs -= coupl * (src - lsrc);
                    if (centralsrc) rhs += coupl * src;
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
        
        Aijkl -= rad_inner().Mm1_x()(i,k) * rad_inner().S_y()(j,l) * (-1.0_r);
        Aijkl -= rad_inner().S_x()(i,k) * rad_inner().Mm1_y()(j,l) * inp_->Zp;
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
    
    // add some more real knots to both sides, if available
    // - 'order' knots to tail to enable coupling to the next basis (if any)
    // - 'order' knots to front to enable coupling to the previous basis (if any)
    // - 'order' knots to front to add cut B-splines to the basis (if any)
    iRa = std::max(iR1, iRa - gap_ - bspline.order());
    iRb = std::min(iR2, iRb + gap_);
    
    // new knot sub-sequence
    rknots = inp_->rknots.slice(iRa, iRb + 1);
    cknots1 = inp_->cknots - inp_->cknots.back() + rknots.front();
    cknots2 = inp_->cknots + rknots.back();
}

void DOMPreconditioner::splitResidual (cBlockArray const & r, std::vector<PanelSolution> & panels) const
{
    for (int ixpanel = 0; ixpanel < xpanels; ixpanel++)
    for (int iypanel = 0; iypanel < ypanels; iypanel++)
    for (unsigned ill = 0; ill < r.size(); ill++)
    {
        PanelSolution & p = panels[ixpanel * ypanels + iypanel];
        
        std::size_t Nfxspline = rad_inner().bspline_x().Nspline();
        std::size_t Nfyspline = rad_inner().bspline_y().Nspline();
        std::size_t Npxspline = p.xspline_inner.Nspline();
        std::size_t Npyspline = p.yspline_inner.Nspline();
        
        p.r[ill].resize(Npxspline * Npyspline);
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
            p.r[ill][pxspline * Npyspline + pyspline] = r[ill][ixspline * Nfyspline + iyspline];
        }
    }
}

void DOMPreconditioner::collectSolution (cBlockArray & z, std::vector<PanelSolution> & panels) const
{
    for (unsigned ill = 0; ill < z.size(); ill++)
    {
        z[ill].fill(0.0);
        
        for (int ixpanel = 0; ixpanel < xpanels; ixpanel++)
        for (int iypanel = 0; iypanel < ypanels; iypanel++)
        {
            PanelSolution & p = panels[ixpanel * ypanels + iypanel];
            
            std::size_t Nfyspline = rad_inner().bspline_y().Nspline();
            std::size_t Npyspline = p.yspline_inner.Nspline();
            
            // for all elements of the panel solution
            for (int pxspline = p.minpxspline; pxspline <= p.maxpxspline; pxspline++)
            for (int pyspline = p.minpyspline; pyspline <= p.maxpyspline; pyspline++)
            {
                // get corresponding indices in global basis
                int ixspline = pxspline + p.xoffset;
                int iyspline = pyspline + p.yoffset;
                
                // update collected solution
                z[ill][ixspline * Nfyspline + iyspline] = p.z[ill][pxspline * Npyspline + pyspline];
            }
        }
        
#ifdef DOM_DEBUG
        static int n = 0;
        std::ofstream out (format("combined-%d-%d.vtk", ill, n++));
        writeVTK_points
        (
            out,
            Bspline::zip
            (
                rad_inner().bspline(),
                rad_inner().bspline(),
                z[ill],
                linspace(rad_inner().bspline().Rmin(), rad_inner().bspline().Rmax(), 2001),
                linspace(rad_inner().bspline().Rmin(), rad_inner().bspline().Rmax(), 2001)
            ),
            linspace(rad_inner().bspline().Rmin(), rad_inner().bspline().Rmax(), 2001),
            linspace(rad_inner().bspline().Rmin(), rad_inner().bspline().Rmax(), 2001),
            rArray{ 0. }
        );
#endif
    }
}

DOMPreconditioner::PanelSolution::PanelSolution
(
    int ix, int iy,
    int order,
    Real theta,
    Bspline const & xspline, Bspline const & yspline,
    rArray cxspline1_inner, rArray rxspline_inner, rArray cxspline2_inner,
    rArray cyspline1_inner, rArray ryspline_inner, rArray cyspline2_inner,
    rArray cxspline1_full,  rArray rxspline_full,  rArray cxspline2_full,
    rArray cyspline1_full,  rArray ryspline_full,  rArray cyspline2_full,
    int Nang
) : xspline_inner (order, theta, cxspline1_inner, rxspline_inner, cxspline2_inner),
    yspline_inner (order, theta, cyspline1_inner, ryspline_inner, cyspline2_inner),
    xspline_full (order, theta, cxspline1_full, rxspline_full, cxspline2_full),
    yspline_full (order, theta, cyspline1_full, ryspline_full, cyspline2_full),
    r (Nang), z (Nang), ixpanel(ix), iypanel(iy)
{
    // allocate memory for residual and solution
    for (int i = 0; i < Nang; i++) r[i].resize(xspline_inner.Nspline() * yspline_inner.Nspline());
    for (int i = 0; i < Nang; i++) z[i].resize(xspline_inner.Nspline() * yspline_inner.Nspline());
    
    // calculate panel-to-full real basis offset
    xoffset = xspline.knot(xspline_inner.R1()) - xspline_inner.iR1();
    yoffset = yspline.knot(yspline_inner.R1()) - yspline_inner.iR1();
    
    int gap = 2*order;
    
    // B-splines that have a counterpart in the global basis
    minpxspline = (ixpanel == 0 ? 0 : xspline_inner.iR1() + gap);
    minpyspline = (iypanel == 0 ? 0 : yspline_inner.iR1() + gap);
    maxpxspline = (ixpanel + 1 == xpanels ? xspline_inner.Nspline() - 1 : xspline_inner.iR2() - order - 1 - gap);
    maxpyspline = (iypanel + 1 == ypanels ? yspline_inner.Nspline() - 1 : yspline_inner.iR2() - order - 1 - gap);
}

// --------------------------------------------------------------------------------- //

addClassToParentRunTimeSelectionTable(PreconditionerBase, DOMPreconditioner)

// --------------------------------------------------------------------------------- //
