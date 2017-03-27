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
    gap_(bspline_inner.order() + 1)
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
    int Nspline = rad_inner().bspline().Nspline();
    int order = rad_inner().bspline().order();
    Real theta = rad_inner().bspline().ECStheta();
    
    // one-electron overlap matrix will be useful
    CsrMatrix<LU_int_t,Complex> SFF_csr = rad_inner_->S().tocoo<LU_int_t>().tocsr();
    std::unique_ptr<LUft> SFF_lu;
    SFF_lu.reset(LUft::Choose("lapack"));
    SFF_lu->factorize(SFF_csr);
    
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
        
        std::cout << std::endl;
        std::cout << "\tPanel (" << ixpanel << "," << iypanel << ")" << std::endl;
        std::cout << "\t  x basis : " << xspline.Rmin() << " " << xspline.R1() << " " << xspline.R2() << " " << xspline.Rmax() << std::endl;
        std::cout << "\t  y basis : " << yspline.Rmin() << " " << yspline.R1() << " " << yspline.R2() << " " << yspline.Rmax() << std::endl;
        
        // full panel basis overlap matrices
        panel.Sxx = RadialIntegrals::computeOverlapMatrix
        (
            xspline,
            xspline,
            xspline.Rmin(),
            xspline.Rmax()
        );
        panel.Syy = RadialIntegrals::computeOverlapMatrix
        (
            yspline,
            yspline,
            yspline.Rmin(),
            yspline.Rmax()
        );
        
        // full basis <-> panel basis overlaps with non-zero support only on real grids out of gap
        panel.SFx = RadialIntegrals::computeOverlapMatrix
        (
            rad_inner().bspline_x(),
            xspline,
            ixpanel     == 0       ? 0              : xspline.t(xspline.iR1() + order + 1).real(),
            ixpanel + 1 == xpanels ? xspline.Rmax() : xspline.t(xspline.iR2() - order - 1).real()
        );
        panel.SFy = RadialIntegrals::computeOverlapMatrix
        (
            rad_inner().bspline_y(),
            yspline,
            iypanel     == 0       ? 0              : yspline.t(yspline.iR1() + order + 1).real(),
            iypanel + 1 == ypanels ? yspline.Rmax() : yspline.t(yspline.iR2() - order - 1).real()
        );
        
        // factorize the square overlap matrices
        panel.lu_Sxx.reset(LUft::Choose("lapack"));
        panel.lu_Syy.reset(LUft::Choose("lapack"));
        panel.lu_Sxx->factorize(panel.Sxx);
        panel.lu_Syy->factorize(panel.Syy);
        
        // calculate reconstruction matrices for this panel
        panel.SFFm1SFx = ColMatrix<Complex>(Nspline, xspline.Nspline());
        panel.SFFm1SFy = ColMatrix<Complex>(Nspline, yspline.Nspline());
        SFF_lu->solve(panel.SFx.tocoo().transpose().torow().data(), panel.SFFm1SFx.data(), xspline.Nspline());
        SFF_lu->solve(panel.SFy.tocoo().transpose().torow().data(), panel.SFFm1SFy.data(), yspline.Nspline());
    }
    
    // interpolate the residual into sub-domains
    splitResidual(r, p);
    
    // find the solution on sub-domains
    std::cout << std::endl;
    int cycles = 2;
    for (int cycle = 0; cycle < cycles; cycle++)
    for (int ixpanel = 0; ixpanel < xpanels; ixpanel++)
    for (int iypanel = 0; iypanel < ypanels; iypanel++)
        solvePanel(cycle, cycles, xpanels, p, ixpanel, iypanel);
    
    // interpolate the solution from sub-domains
    collectSolution(z, p);
    
    std::cout << std::endl << std::endl;
    std::cout << "\t    -> DOM preconditioning done!    ";
}

void DOMPreconditioner::finish ()
{
    NoPreconditioner::finish();
}

void DOMPreconditioner::solvePanel (int cycle, int cycles, int n, std::vector<PanelSolution> & p, int i, int j) const
{
    std::cout << std::endl;
    std::cout << "\t-------------------------------------------" << std::endl;
    std::cout << "\tSolve panel " << i << " " << j << " (cycle " << cycle + 1 << " of " << cycles << ")" << std::endl;
    std::cout << "\t-------------------------------------------" << std::endl;
    std::cout << std::endl;
    
    // get reference to the current panel and its neighbours
    PanelSolution * pCentre = &p[i * n + j];
    PanelSolution * pNbr [nNbrs];
    pNbr[Left]  = (i     == 0 ? nullptr : &p[(i - 1) * n + j]);
    pNbr[Right] = (i + 1 == n ? nullptr : &p[(i + 1) * n + j]);
    pNbr[Down]  = (j     == 0 ? nullptr : &p[i * n + (j - 1)]);
    pNbr[Up]    = (j + 1 == n ? nullptr : &p[i * n + (j + 1)]);
    
    // B-spline order
    int order = pCentre->xspline_inner.order();
    
    // construct the surrogate sources for all boundaries shared with other panels
    for (int nbr = 0; nbr < nNbrs; nbr++) if (pNbr[nbr] != nullptr)
        surrogateSource(cycle, cycles, pCentre, nbr, pNbr[nbr]);
    
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
    
    // get radial integrals
    RadialIntegrals const & rad = dynamic_cast<NoPreconditioner*>(prec)->rad_panel();
    
    // number of B-splines in both directions of this panel
    int Nxspline = pCentre->xspline_inner.Nspline();
    int Nyspline = pCentre->yspline_inner.Nspline();
    
    // get right-hand side and solution arrays
    cBlockArray & psi = pCentre->z;
    cBlockArray   chi = pCentre->r;
    
    // reset the solution
    for (cArray & segment : psi)
        segment.fill(0.0);
    
    // add surrogate sources
    if (pNbr[Left] != nullptr)
    {
        int iR1 = pCentre->xspline_inner.iR1();
        std::cout << "\tAdd surrogate source from LEFT" << std::endl;
        for (unsigned ill = 0; ill < chi.size(); ill++)
        for (int i = iR1 - order; i <= iR1 + order; i++)
        for (int j = 0; j < Nyspline; j++)
            chi[ill][i * Nyspline + j] += rad.S_x(iR1, i) * pCentre->ssrc[Left][ill][j];
    }
    if (pNbr[Right] != nullptr)
    {
        int iR2 = pCentre->xspline_inner.iR2();
        std::cout << "\tAdd surrogate source from RIGHT" << std::endl;
        for (unsigned ill = 0; ill < chi.size(); ill++)
        for (int i = iR2 - 2*order - 1; i < iR2; i++)
        for (int j = 0; j < Nyspline; j++)
            chi[ill][i * Nyspline + j] += rad.S_x(iR2 - order - 1, i) * pCentre->ssrc[Right][ill][j];
    }
    if (pNbr[Down] != nullptr)
    {
        int iR1 = pCentre->yspline_inner.iR1();
        std::cout << "\tAdd surrogate source from DOWN" << std::endl;
        for (unsigned ill = 0; ill < chi.size(); ill++)
        for (int i = 0; i < Nxspline; i++)
        for (int j = iR1 - order; j <= iR1 + order; j++)
            chi[ill][i * Nyspline + j] += rad.S_y(iR1, j) * pCentre->ssrc[Down][ill][i];
    }
    if (pNbr[Up] != nullptr)
    {
        int iR2 = pCentre->yspline_inner.iR2();
        std::cout << "\tAdd surrogate source from UP" << std::endl;
        for (unsigned ill = 0; ill < chi.size(); ill++)
        for (int i = 0; i < Nxspline; i++)
        for (int j = iR2 - 2*order - 1; j < iR2; j++)
            chi[ill][i * Nyspline + j] += rad.S_y(iR2 - order - 1, j) * pCentre->ssrc[Up][ill][i];
    }
    
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
    
#ifdef DOM_DEBUG
    for (unsigned ill = 0; ill < chi.size(); ill++)
    {
        std::ofstream ofs0 (format("dom-%d-%d-%d-chi0-%d.vtk", cycle, i, j, ill));
        writeVTK_points
        (
            ofs0,
            Bspline::zip
            (
                pCentre->xspline_inner,
                pCentre->yspline_inner,
                chi[ill],
                linspace(pCentre->xspline_inner.Rmin(), pCentre->xspline_inner.Rmax(), 1001),
                linspace(pCentre->yspline_inner.Rmin(), pCentre->yspline_inner.Rmax(), 1001)
            ),
            linspace(pCentre->xspline_inner.Rmin(), pCentre->xspline_inner.Rmax(), 1001),
            linspace(pCentre->yspline_inner.Rmin(), pCentre->yspline_inner.Rmax(), 1001),
            rArray{ 0. }
        );
        ofs0.close();
        
        cArray c (Nxspline * Nyspline);
        cArray d (Nxspline * Nyspline);
        pCentre->lu_Syy->solve(chi[ill], c, Nxspline);
        transpose(c, Nyspline, Nxspline);
        pCentre->lu_Sxx->solve(c, d, Nyspline);
        transpose(d, Nxspline, Nyspline);
        
        std::ofstream ofs (format("dom-%d-%d-%d-chi-%d.vtk", cycle, i, j, ill));
        writeVTK_points
        (
            ofs,
            Bspline::zip
            (
                pCentre->xspline_inner,
                pCentre->yspline_inner,
                d,
                linspace(pCentre->xspline_inner.Rmin(), pCentre->xspline_inner.Rmax(), 1001),
                linspace(pCentre->yspline_inner.Rmin(), pCentre->yspline_inner.Rmax(), 1001)
            ),
            linspace(pCentre->xspline_inner.Rmin(), pCentre->xspline_inner.Rmax(), 1001),
            linspace(pCentre->yspline_inner.Rmin(), pCentre->yspline_inner.Rmax(), 1001),
            rArray{ 0. }
        );
        ofs.close();
    }
#endif
    
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
    rArray grid_x = linspace(pCentre->xspline_inner.Rmin(), pCentre->xspline_inner.Rmax(), 1001);
    if (grid_x.front() < 50. and 50. < grid_x.back())
        grid_x.insert(std::lower_bound(grid_x.begin(), grid_x.end(), 50.), 50.);
    rArray grid_y = linspace(pCentre->yspline_inner.Rmin(), pCentre->yspline_inner.Rmax(), 1001);
    if (grid_y.front() < 50. and 50. < grid_y.back())
        grid_y.insert(std::lower_bound(grid_y.begin(), grid_y.end(), 50.), 50.);
    for (unsigned ill = 0; ill < psi.size(); ill++)
    {
        std::ofstream ofs (format("dom-%d-%d-%d-psi-%d.vtk", cycle, i, j, ill));
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
    
    // evaluate the field at the edge of the domain
    if (pNbr[Left] != nullptr)
    {
        int iR1 = pCentre->xspline_inner.iR1();
        Real R1a = pCentre->xspline_inner.t(iR1 + order + 1).real();
        
        for (unsigned ill = 0; ill < chi.size(); ill++)
        {
            pCentre->outf[Left][ill].fill(0);
            for (int i = iR1 + 1; i < iR1 + order + 1; i++)
            for (int j = 0; j < Nyspline; j++)
                pCentre->outf[Left][ill][j] += psi[ill][i * Nyspline + j] * pCentre->xspline_inner.bspline(i, iR1 + order + 1, order, R1a);
            
#ifdef DOM_DEBUG
            write_array
            (
                linspace(pCentre->yspline_inner.Rmin(), pCentre->yspline_inner.Rmax(), 10001),
                pCentre->yspline_inner.zip
                (
                    pCentre->outf[Left][ill],
                    linspace(pCentre->yspline_inner.Rmin(), pCentre->yspline_inner.Rmax(), 10001)
                ),
                format("dom-%d-%d-%d-outf-left-%d.dat", cycle, i, j, ill)
            );
#endif
        }
    }
    if (pNbr[Right] != nullptr)
    {
        int iR2 = pCentre->xspline_inner.iR2();
        Real R2a = pCentre->xspline_inner.t(iR2 - order - 1).real();
        
        for (unsigned ill = 0; ill < chi.size(); ill++)
        {
            pCentre->outf[Right][ill].fill(0);
            for (int i = iR2 - 2*order - 1; i < iR2 - order - 1; i++)
            for (int j = 0; j < Nyspline; j++)
                pCentre->outf[Right][ill][j] += psi[ill][i * Nyspline + j] * pCentre->xspline_inner.bspline(i, iR2 - order - 1, order, R2a);
            
#ifdef DOM_DEBUG
            write_array
            (
                linspace(pCentre->yspline_inner.Rmin(), pCentre->yspline_inner.Rmax(), 10001),
                pCentre->yspline_inner.zip
                (
                    pCentre->outf[Right][ill],
                    linspace(pCentre->yspline_inner.Rmin(), pCentre->yspline_inner.Rmax(), 10001)
                ),
                format("dom-%d-%d-%d-outf-right-%d.dat", cycle, i, j, ill)
            );
#endif
        }
    }
    if (pNbr[Down] != nullptr)
    {
        int iR1 = pCentre->yspline_inner.iR1();
        Real R1a = pCentre->yspline_inner.t(iR1 + order + 1).real();
        
        for (unsigned ill = 0; ill < chi.size(); ill++)
        {
            pCentre->outf[Down][ill].fill(0);
            for (int i = 0; i < Nxspline; i++)
            for (int j = iR1 + 1; j < iR1 + order + 1; j++)
                pCentre->outf[Down][ill][i] += psi[ill][i * Nyspline + j] * pCentre->yspline_inner.bspline(j, iR1 + order + 1, order, R1a);
            
#ifdef DOM_DEBUG
            write_array
            (
                linspace(pCentre->xspline_inner.Rmin(), pCentre->xspline_inner.Rmax(), 10001),
                pCentre->xspline_inner.zip
                (
                    pCentre->outf[Down][ill],
                    linspace(pCentre->xspline_inner.Rmin(), pCentre->xspline_inner.Rmax(), 10001)
                ),
                format("dom-%d-%d-%d-outf-down-%d.dat", cycle, i, j, ill)
            );
#endif
        }
    }
    if (pNbr[Up] != nullptr)
    {
        int iR2 = pCentre->yspline_inner.iR2();
        Real R2a = pCentre->yspline_inner.t(iR2 - order - 1).real();
        
        for (unsigned ill = 0; ill < chi.size(); ill++)
        {
            pCentre->outf[Up][ill].fill(0);
            for (int i = 0; i < Nxspline; i++)
            for (int j = iR2 - 2*order - 1; j < iR2 - order - 1; j++)
                pCentre->outf[Up][ill][i] += psi[ill][i * Nyspline + j] * pCentre->yspline_inner.bspline(j, iR2 - order - 1, order, R2a);
            
#ifdef DOM_DEBUG
            write_array
            (
                linspace(pCentre->xspline_inner.Rmin(), pCentre->xspline_inner.Rmax(), 10001),
                pCentre->xspline_inner.zip
                (
                    pCentre->outf[Up][ill],
                    linspace(pCentre->xspline_inner.Rmin(), pCentre->xspline_inner.Rmax(), 10001)
                ),
                format("dom-%d-%d-%d-outf-up-%d.dat", cycle, i, j, ill)
            );
#endif
        }
    }
    
    prec->finish();
}

void DOMPreconditioner::surrogateSource (int cycle, int cycles, PanelSolution * panel, int direction, PanelSolution * neighbour) const
{
    //
    // What needs to be yet taken care of here:
    //
    //  1. The angular states are not coupled by the multi-pole potential.
    //     This may corrupt the accuracy.
    //
    //  2. The surrogate sources are calculated even in the complex parts of the grid of the other dimension.
    //     These tails may contribute to other edges, which is not desirable.
    //
    
    if (direction == Left)
    {
        std::cout << "\tConstruct surrogate source from LEFT" << std::endl;
        
        // original B-spline basis
        Bspline const & org_xspline = panel->xspline_inner;
        
        // short B-spline basis
        Bspline xspline
        (
            org_xspline.order(),
            org_xspline.ECStheta(),
            org_xspline.cknots1(),
            org_xspline.rknots(),
            org_xspline.cknots2()
        );
        Bspline const & yspline = panel->yspline_inner;
        
        // basis sizes
        int Nxspline = xspline.Nspline();
        int Nyspline = yspline.Nspline();
        
        // other quantities
        int order = xspline.order();
        int iR1x = xspline.iR1();
        Real R1xa = xspline.t(iR1x + order + 1).real();
        
        // calculate radial integrals
        RadialIntegrals rint (xspline, yspline, ang_->maxlambda() + 1);
        rint.verbose(false);
        rint.setupOneElectronIntegrals(*par_, *cmd_);
        rint.setupTwoElectronIntegrals(*par_, *cmd_);
        //rint.subsetTwoElectronIntegrals(*par_, *cmd_, rad_inner());
        
        // for all angular states
        for (unsigned ill = 0; ill < ang_->states().size(); ill++)
        {
            // angular momenta
            int l1 = ang_->states()[ill].first;
            int l2 = ang_->states()[ill].second;
            
            // construct matrix of the free equations
            BlockSymBandMatrix<Complex> A =
                kron(Complex(E_) * rint.S_x(), rint.S_y())
                - kron(0.5_z * rint.D_x(), rint.S_y())
                - kron(0.5_z * rint.S_x(), rint.D_y())
                - kron(Complex(0.5 * l1 * (l1 + 1)) * rint.Mm2_x(), rint.S_y())
                - kron(Complex(0.5 * l2 * (l2 + 1)) * rint.S_x(), rint.Mm2_y())
                + kron(rint.Mm1_x(), rint.S_y())
                + kron(rint.S_x(), rint.Mm1_y());
            
            for (int lambda = 0; lambda <= ang_->maxlambda(); lambda++)
                A -= Complex(ang_->f(ill, ill, lambda)) * rint.R_tr_dia(lambda);
            
            // resize matrix for additional elements
            CooMatrix<LU_int_t,Complex> A_coo = A.tocoo<LU_int_t>();
            CooMatrix<LU_int_t,Complex> A_coo0 = A_coo;
            A_coo.resize
            (
                A_coo.rows() + Nyspline,
                A_coo.cols() + Nyspline
            );
            
            // add coupling to surrogate source
            for (int i = iR1x - order; i <= iR1x + order; i++)
            for (int j = 0; j < Nyspline; j++)
            {
                A_coo.add
                (
                    i * Nyspline + j,
                    Nxspline * Nyspline + j,
                    -rint.S_x(i, iR1x)
                );
            }
            
            // add coupling to outgoing field
            for (int j = 0; j < Nyspline; j++)
            for (int k = iR1x + 1; k < iR1x + order + 1; k++)
            {
                A_coo.add
                (
                    Nxspline * Nyspline + j,
                    k * Nyspline + j,
                    xspline.bspline(k, iR1x + order + 1, order, R1xa)
                );
            }
            
            // set up the right-hand side
            cArray rhs (Nxspline * Nyspline + Nyspline);
            cArrayView (rhs, Nxspline * Nyspline, Nyspline) = neighbour->outf[Right][ill] - panel->outf[Left][ill];
            
#ifdef DOM_DEBUG
            write_array
            (
                linspace(yspline.Rmin(), yspline.Rmax(), 10001),
                yspline.zip
                (
                    neighbour->outf[Right][ill] - panel->outf[Left][ill],
                    linspace(yspline.Rmin(), yspline.Rmax(), 10001)
                ),
                format("dom-%d-%d-%d-left-nbrf-%d.txt", cycle, panel->ixpanel, panel->iypanel, ill)
            );
#endif
            
            // factorize the block matrix and solve
            cArray solution (rhs.size());
            CsrMatrix<LU_int_t,Complex> A_csr = A_coo.tocsr();
            std::unique_ptr<LUft> A_lu;
            A_lu.reset(LUft::Choose("umfpack")); // TODO : Make (wisely) runtime selectable.
            A_lu->factorize(A_csr);
            A_lu->solve(rhs, solution, 1);
            
            // update the surrogate source
            panel->ssrc[Left][ill] += cArrayView (solution, Nxspline * Nyspline, Nyspline);
            
#ifdef DOM_DEBUG
            std::ofstream ofsp (format("dom-%d-%d-%d-left-ispsi-%d.vtk", cycle, panel->ixpanel, panel->iypanel, ill));
            writeVTK_points
            (
                ofsp,
                Bspline::zip
                (
                    xspline,
                    yspline,
                    solution.slice(0, Nxspline * Nyspline),
                    linspace(xspline.Rmin(), xspline.Rmax(), 1001),
                    linspace(yspline.Rmin(), yspline.Rmax(), 1001)
                ),
                linspace(xspline.Rmin(), xspline.Rmax(), 1001),
                linspace(yspline.Rmin(), yspline.Rmax(), 1001),
                rArray{ 0 }
            );
            ofsp.close();
#endif
        }
    }
    
    if (direction == Right)
    {
        std::cout << "\tConstruct surrogate source from RIGHT" << std::endl;
        
        // original B-spline basis
        Bspline const & org_xspline = panel->xspline_inner;
        
        // short B-spline basis
        Bspline xspline
        (
            org_xspline.order(),
            org_xspline.ECStheta(),
            org_xspline.cknots1(),
            org_xspline.rknots(),
            org_xspline.cknots2()
        );
        Bspline const & yspline = panel->yspline_inner;
        
        // basis sizes
        int Nxspline = xspline.Nspline();
        int Nyspline = yspline.Nspline();
        
        // other quantities
        int order = xspline.order();
        int iR2x = xspline.iR2();
        Real R2xa = xspline.t(iR2x - order - 1).real();
        
        // calculate radial integrals
        RadialIntegrals rint (xspline, yspline, ang_->maxlambda() + 1);
        rint.verbose(false);
        rint.setupOneElectronIntegrals(*par_, *cmd_);
        rint.setupTwoElectronIntegrals(*par_, *cmd_);
        
        // for all angular states
        for (unsigned ill = 0; ill < ang_->states().size(); ill++)
        {
            // angular momenta
            int l1 = ang_->states()[ill].first;
            int l2 = ang_->states()[ill].second;
            
            // construct matrix of the free equations
            BlockSymBandMatrix<Complex> A =
                kron(Complex(E_) * rint.S_x(), rint.S_y())
                - kron(0.5_z * rint.D_x(), rint.S_y())
                - kron(0.5_z * rint.S_x(), rint.D_y())
                - kron(Complex(0.5 * l1 * (l1 + 1)) * rint.Mm2_x(), rint.S_y())
                - kron(Complex(0.5 * l2 * (l2 + 1)) * rint.S_x(), rint.Mm2_y())
                + kron(rint.Mm1_x(), rint.S_y())
                + kron(rint.S_x(), rint.Mm1_y());
            
            for (int lambda = 0; lambda <= ang_->maxlambda(); lambda++)
                A -= Complex(ang_->f(ill, ill, lambda)) * rint.R_tr_dia(lambda);
            
            // resize matrix for additional elements
            CooMatrix<LU_int_t,Complex> A_coo = A.tocoo<LU_int_t>();
            CooMatrix<LU_int_t,Complex> A_coo0 = A_coo;
            A_coo.resize
            (
                A_coo.rows() + Nyspline,
                A_coo.cols() + Nyspline
            );
            
            // add coupling to surrogate source
            for (int i = iR2x - 2*order - 1; i < iR2x; i++)
            for (int j = 0; j < Nyspline; j++)
            {
                A_coo.add
                (
                    i * Nyspline + j,
                    Nxspline * Nyspline + j,
                    -rint.S_x(i, iR2x - order - 1)
                );
            }
            
            // add coupling to outgoing field
            for (int j = 0; j < Nyspline; j++)
            for (int k = iR2x - 2*order - 1; k < iR2x - order - 1; k++)
            {
                A_coo.add
                (
                    Nxspline * Nyspline + j,
                    k * Nyspline + j,
                    xspline.bspline(k, iR2x - order - 1, order, R2xa)
                );
            }
            
            // set up the right-hand side
            cArray rhs (Nxspline * Nyspline + Nyspline);
            cArrayView (rhs, Nxspline * Nyspline, Nyspline) = neighbour->outf[Left][ill] - panel->outf[Right][ill];
            
#ifdef DOM_DEBUG
            write_array
            (
                linspace(yspline.Rmin(), yspline.Rmax(), 10001),
                yspline.zip
                (
                    neighbour->outf[Left][ill] - panel->outf[Right][ill],
                    linspace(yspline.Rmin(), yspline.Rmax(), 10001)
                ),
                format("dom-%d-%d-%d-right-nbrf-%d.txt", cycle, panel->ixpanel, panel->iypanel, ill)
            );
#endif
            
            // factorize the block matrix and solve
            cArray solution (rhs.size());
            CsrMatrix<LU_int_t,Complex> A_csr = A_coo.tocsr();
            std::unique_ptr<LUft> A_lu;
            A_lu.reset(LUft::Choose("umfpack")); // TODO : Make (wisely) runtime selectable.
            A_lu->factorize(A_csr);
            A_lu->solve(rhs, solution, 1);
            
            // update the surrogate source
            panel->ssrc[Right][ill] += cArrayView (solution, Nxspline * Nyspline, Nyspline);
            
#ifdef DOM_DEBUG
            std::ofstream ofsp (format("dom-%d-%d-%d-right-ispsi-%d.vtk", cycle, panel->ixpanel, panel->iypanel, ill));
            writeVTK_points
            (
                ofsp,
                Bspline::zip
                (
                    xspline,
                    yspline,
                    solution.slice(0, Nxspline * Nyspline),
                    linspace(xspline.Rmin(), xspline.Rmax(), 1001),
                    linspace(yspline.Rmin(), yspline.Rmax(), 1001)
                ),
                linspace(xspline.Rmin(), xspline.Rmax(), 1001),
                linspace(yspline.Rmin(), yspline.Rmax(), 1001),
                rArray{ 0 }
            );
            ofsp.close();
#endif
        }
    }
    
    if (direction == Down)
    {
        std::cout << "\tConstruct surrogate source from DOWN" << std::endl;
        
        // original B-spline basis
        Bspline const & org_yspline = panel->yspline_inner;
        
        // short B-spline basis
        Bspline yspline
        (
            org_yspline.order(),
            org_yspline.ECStheta(),
            org_yspline.cknots1(),
            org_yspline.rknots(),
            org_yspline.cknots2()
        );
        Bspline const & xspline = panel->xspline_inner;
        
        // basis sizes
        int Nxspline = xspline.Nspline();
        int Nyspline = yspline.Nspline();
        
        // other quantities
        int order = yspline.order();
        int iR1y = yspline.iR1();
        Real R1ya = yspline.t(iR1y + order + 1).real();
        
        // calculate radial integrals
        RadialIntegrals rint (xspline, yspline, ang_->maxlambda() + 1);
        rint.verbose(false);
        rint.setupOneElectronIntegrals(*par_, *cmd_);
        rint.setupTwoElectronIntegrals(*par_, *cmd_);
        
        // for all angular states
        for (unsigned ill = 0; ill < ang_->states().size(); ill++)
        {
            // angular momenta
            int l1 = ang_->states()[ill].first;
            int l2 = ang_->states()[ill].second;
            
            // construct matrix of the free equations
            BlockSymBandMatrix<Complex> A =
                kron(Complex(E_) * rint.S_x(), rint.S_y())
                - kron(0.5_z * rint.D_x(), rint.S_y())
                - kron(0.5_z * rint.S_x(), rint.D_y())
                - kron(Complex(0.5 * l1 * (l1 + 1)) * rint.Mm2_x(), rint.S_y())
                - kron(Complex(0.5 * l2 * (l2 + 1)) * rint.S_x(), rint.Mm2_y())
                + kron(rint.Mm1_x(), rint.S_y())
                + kron(rint.S_x(), rint.Mm1_y());
            
            for (int lambda = 0; lambda <= ang_->maxlambda(); lambda++)
                A -= Complex(ang_->f(ill, ill, lambda)) * rint.R_tr_dia(lambda);
            
            // resize matrix for additional elements
            CooMatrix<LU_int_t,Complex> A_coo = A.tocoo<LU_int_t>();
            CooMatrix<LU_int_t,Complex> A_coo0 = A_coo;
            A_coo.resize
            (
                A_coo.rows() + Nxspline,
                A_coo.cols() + Nxspline
            );
            
            // add coupling to surrogate source
            for (int i = 0; i < Nxspline; i++)
            for (int j = iR1y - order; j < iR1y + order + 1; j++)
            {
                A_coo.add
                (
                    i * Nyspline + j,
                    Nxspline * Nyspline + i,
                    -rint.S_y(j, iR1y)
                );
            }
            
            // add coupling to outgoing field
            for (int i = 0; i < Nxspline; i++)
            for (int l = iR1y + 1; l < iR1y + order + 1; l++)
            {
                A_coo.add
                (
                    Nxspline * Nyspline + i,
                    i * Nyspline + l,
                    yspline.bspline(l, iR1y + order + 1, order, R1ya)
                );
            }
            
            // set up the right-hand side
            cArray rhs (Nxspline * Nyspline + Nxspline);
            cArrayView (rhs, Nxspline * Nyspline, Nxspline) = neighbour->outf[Up][ill] - panel->outf[Down][ill];
            
#ifdef DOM_DEBUG
            write_array
            (
                linspace(xspline.Rmin(), xspline.Rmax(), 10001),
                xspline.zip
                (
                    neighbour->outf[Up][ill] - panel->outf[Down][ill],
                    linspace(xspline.Rmin(), xspline.Rmax(), 10001)
                ),
                format("dom-%d-%d-%d-down-nbrf-%d.txt", cycle, panel->ixpanel, panel->iypanel, ill)
            );
#endif
            
            // factorize the block matrix and solve
            cArray solution (rhs.size());
            CsrMatrix<LU_int_t,Complex> A_csr = A_coo.tocsr();
            std::unique_ptr<LUft> A_lu;
            A_lu.reset(LUft::Choose("umfpack")); // TODO : Make (wisely) runtime selectable.
            A_lu->factorize(A_csr);
            A_lu->solve(rhs, solution, 1);
            
            // update the surrogate source
            panel->ssrc[Down][ill] += cArrayView (solution, Nxspline * Nyspline, Nxspline);
            
#ifdef DOM_DEBUG
            std::ofstream ofsp (format("dom-%d-%d-%d-down-ispsi-%d.vtk", cycle, panel->ixpanel, panel->iypanel, ill));
            writeVTK_points
            (
                ofsp,
                Bspline::zip
                (
                    xspline,
                    yspline,
                    solution.slice(0, Nxspline * Nyspline),
                    linspace(xspline.Rmin(), xspline.Rmax(), 1001),
                    linspace(yspline.Rmin(), yspline.Rmax(), 1001)
                ),
                linspace(xspline.Rmin(), xspline.Rmax(), 1001),
                linspace(yspline.Rmin(), yspline.Rmax(), 1001),
                rArray{ 0 }
            );
            ofsp.close();
#endif
        }
    }
    
    if (direction == Up)
    {
        std::cout << "\tConstruct surrogate source from UP" << std::endl;
        
        // original B-spline basis
        Bspline const & org_yspline = panel->yspline_inner;
        
        // short B-spline basis
        Bspline yspline
        (
            org_yspline.order(),
            org_yspline.ECStheta(),
            org_yspline.cknots1(),
            org_yspline.rknots(),
            org_yspline.cknots2()
        );
        Bspline const & xspline = panel->xspline_inner;
        
        // basis sizes
        int Nxspline = xspline.Nspline();
        int Nyspline = yspline.Nspline();
        
        // other quantities
        int order = yspline.order();
        int iR2y = yspline.iR2();
        Real R2ya = yspline.t(iR2y - order - 1).real();
        
        // calculate radial integrals
        RadialIntegrals rint (xspline, yspline, ang_->maxlambda() + 1);
        rint.verbose(false);
        rint.setupOneElectronIntegrals(*par_, *cmd_);
        rint.setupTwoElectronIntegrals(*par_, *cmd_);
        
        // for all angular states
        for (unsigned ill = 0; ill < ang_->states().size(); ill++)
        {
            // angular momenta
            int l1 = ang_->states()[ill].first;
            int l2 = ang_->states()[ill].second;
            
            // construct matrix of the free equations
            BlockSymBandMatrix<Complex> A =
                kron(Complex(E_) * rint.S_x(), rint.S_y())
                - kron(0.5_z * rint.D_x(), rint.S_y())
                - kron(0.5_z * rint.S_x(), rint.D_y())
                - kron(Complex(0.5 * l1 * (l1 + 1)) * rint.Mm2_x(), rint.S_y())
                - kron(Complex(0.5 * l2 * (l2 + 1)) * rint.S_x(), rint.Mm2_y())
                + kron(rint.Mm1_x(), rint.S_y())
                + kron(rint.S_x(), rint.Mm1_y());
            
            for (int lambda = 0; lambda <= ang_->maxlambda(); lambda++)
                A -= Complex(ang_->f(ill, ill, lambda)) * rint.R_tr_dia(lambda);
            
            // resize matrix for additional elements
            CooMatrix<LU_int_t,Complex> A_coo = A.tocoo<LU_int_t>();
            CooMatrix<LU_int_t,Complex> A_coo0 = A_coo;
            A_coo.resize
            (
                A_coo.rows() + Nxspline,
                A_coo.cols() + Nxspline
            );
            
            // add coupling to surrogate source
            for (int i = 0; i < Nxspline; i++)
            for (int j = iR2y - 2*order - 1; j < iR2y; j++)
            {
                A_coo.add
                (
                    i * Nyspline + j,
                    Nxspline * Nyspline + i,
                    -rint.S_y(j, iR2y - order - 1)
                );
            }
            
            // add coupling to outgoing field
            for (int i = 0; i < Nxspline; i++)
            for (int l = iR2y - 2*order - 1; l < iR2y - order - 1; l++)
            {
                A_coo.add
                (
                    Nxspline * Nyspline + i,
                    i * Nyspline + l,
                    yspline.bspline(l, iR2y - order - 1, order, R2ya)
                );
            }
            
            // set up the right-hand side
            cArray rhs (Nxspline * Nyspline + Nxspline);
            cArrayView (rhs, Nxspline * Nyspline, Nxspline) = neighbour->outf[Down][ill] - panel->outf[Up][ill];
            
#ifdef DOM_DEBUG
            write_array
            (
                linspace(xspline.Rmin(), xspline.Rmax(), 10001),
                xspline.zip
                (
                    neighbour->outf[Down][ill] - panel->outf[Up][ill],
                    linspace(xspline.Rmin(), xspline.Rmax(), 10001)
                ),
                format("dom-%d-%d-%d-up-nbrf-%d.txt", cycle, panel->ixpanel, panel->iypanel, ill)
            );
#endif
            
            // factorize the block matrix and solve
            cArray solution (rhs.size());
            CsrMatrix<LU_int_t,Complex> A_csr = A_coo.tocsr();
            std::unique_ptr<LUft> A_lu;
            A_lu.reset(LUft::Choose("umfpack")); // TODO : Make (wisely) runtime selectable.
            A_lu->factorize(A_csr);
            A_lu->solve(rhs, solution, 1);
            
            // update the surrogate source
            panel->ssrc[Up][ill] += cArrayView (solution, Nxspline * Nyspline, Nxspline);
            
#ifdef DOM_DEBUG
            std::ofstream ofsp (format("dom-%d-%d-%d-up-ispsi-%d.vtk", cycle, panel->ixpanel, panel->iypanel, ill));
            writeVTK_points
            (
                ofsp,
                Bspline::zip
                (
                    xspline,
                    yspline,
                    solution.slice(0, Nxspline * Nyspline),
                    linspace(xspline.Rmin(), xspline.Rmax(), 1001),
                    linspace(yspline.Rmin(), yspline.Rmax(), 1001)
                ),
                linspace(xspline.Rmin(), xspline.Rmax(), 1001),
                linspace(yspline.Rmin(), yspline.Rmax(), 1001),
                rArray{ 0 }
            );
            ofsp.close();
#endif
        }
    }
}

void DOMPreconditioner::knotSubsequence (int i, int n, Bspline const & bspline, rArray & rknots, rArray & cknots1, rArray & cknots2) const
{
    // boundaries of the original real grid
    Real R1 = bspline.R1();
    Real R2 = bspline.R2();
    
    // boundaries of the sub-grid
    Real Ra = R1 + (R2 - R1) * i / n;
    Real Rb = R1 + (R2 - R1) * (i + 1) / n;
    
    // indices of the knots at coordinates Ra and Rb
    int iRa = bspline.knot(Ra);
    int iRb = bspline.knot(Rb);
    
    // add some more knots to both sides for surrogate source (if available)
    iRa = std::max(bspline.iR1(), iRa - bspline.order() - 1);
    iRb = std::min(bspline.iR2(), iRb + bspline.order() + 1);
    
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
        
        int order = rad_full().bspline().order(), pxspline, pyspline;
        std::size_t Nfyspline = rad_inner().bspline_y().Nspline();
        std::size_t Npxspline = p.xspline_inner.Nspline();
        std::size_t Npyspline = p.yspline_inner.Nspline();
        
        // allocate memory or erase the old residual
        p.r[ill].resize(Npxspline * Npyspline);
        p.r[ill].fill(0.);
        
        // copy only elements corresponding to the real B-splines (and avoid gap)
        for (int ixspline = 0; ixspline + order < rad_inner().bspline_x().iR2(); ixspline++)
        for (int iyspline = 0; iyspline + order < rad_inner().bspline_y().iR2(); iyspline++)
        if (p.mapToPanel(ixspline, iyspline, pxspline, pyspline))
        {
            if
            (
                // left & down : the B-spline must start after gap
                    (ixpanel == 0 or p.xspline_full.iR1() + gap_ + order + 1 <= pxspline) and
                    (iypanel == 0 or p.yspline_full.iR1() + gap_ + order + 1 <= pyspline) and
                // right & up : the B-spline must end before gap
                    (ixpanel + 1 == xpanels or pxspline + gap_ + 2 * (order + 1) <= p.xspline_full.iR2()) and
                    (iypanel + 1 == ypanels or pyspline + gap_ + 2 * (order + 1) <= p.yspline_full.iR2())
            )
            {
                // copy the element to the panel residual
                p.r[ill][pxspline * Npyspline + pyspline] = r[ill][ixspline * Nfyspline + iyspline];
            }
        }
    }
}

void DOMPreconditioner::collectSolution (cBlockArray & z, std::vector<PanelSolution> & panels) const
{
    for (unsigned ill = 0; ill < z.size(); ill++)
    {
        // reset solution
        z[ill].fill(0.0);
        
        // combine panel solutions
        for (PanelSolution const & panel : panels)
            z[ill] += kron_dot(panel.SFFm1SFx, panel.SFFm1SFy, panel.z[ill]);
        
#ifdef DOM_DEBUG
        std::ofstream out (format("combined-%d.vtk", ill));
        writeVTK_points
        (
            out,
            Bspline::zip
            (
                rad_inner().bspline(),
                rad_inner().bspline(),
                z[ill],
                linspace(rad_inner().bspline().Rmin(), rad_inner().bspline().Rmax(), 1001),
                linspace(rad_inner().bspline().Rmin(), rad_inner().bspline().Rmax(), 1001)
            ),
            linspace(rad_inner().bspline().Rmin(), rad_inner().bspline().Rmax(), 1001),
            linspace(rad_inner().bspline().Rmin(), rad_inner().bspline().Rmax(), 1001),
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
    // allocate memory for surrogate sources and outgoing fields
    ssrc[Left] .resize(Nang); for (cArray & arr : ssrc[Left] ) arr.resize(yspline_inner.Nspline());
    outf[Left] .resize(Nang); for (cArray & arr : outf[Left] ) arr.resize(yspline_inner.Nspline());
    ssrc[Right].resize(Nang); for (cArray & arr : ssrc[Right]) arr.resize(yspline_inner.Nspline());
    outf[Right].resize(Nang); for (cArray & arr : outf[Right]) arr.resize(yspline_inner.Nspline());
    ssrc[Down] .resize(Nang); for (cArray & arr : ssrc[Down] ) arr.resize(xspline_inner.Nspline());
    outf[Down] .resize(Nang); for (cArray & arr : outf[Down] ) arr.resize(xspline_inner.Nspline());
    ssrc[Up]   .resize(Nang); for (cArray & arr : ssrc[Up]   ) arr.resize(xspline_inner.Nspline());
    outf[Up]   .resize(Nang); for (cArray & arr : outf[Up]   ) arr.resize(xspline_inner.Nspline());
    
    // allocate memory for residual and solution
    for (int i = 0; i < Nang; i++) r[i].resize(xspline_inner.Nspline() * yspline_inner.Nspline());
    for (int i = 0; i < Nang; i++) z[i].resize(xspline_inner.Nspline() * yspline_inner.Nspline());
    
    // calculate panel-to-full real basis offset
    xoffset = xspline.knot(xspline_inner.R1()) - xspline_inner.iR1();
    yoffset = yspline.knot(yspline_inner.R1()) - yspline_inner.iR1();
}

bool DOMPreconditioner::PanelSolution::mapToPanel
(
    int   ixspline, int   iyspline,
    int & pxspline, int & pyspline
) const
{
    pxspline = ixspline - xoffset;
    pyspline = iyspline - yoffset;
    
    return pxspline >= xspline_inner.iR1() and
           pyspline >= yspline_inner.iR1() and
           pxspline + xspline_inner.order() + 1 <= xspline_inner.iR2() and
           pyspline + yspline_inner.order() + 1 <= yspline_inner.iR2();
}

bool DOMPreconditioner::PanelSolution::mapFromPanel
(
    int & ixspline, int & iyspline,
    int   pxspline, int   pyspline
) const
{
    pxspline = ixspline + xoffset;
    pyspline = iyspline + yoffset;
    
    return true;
}

// --------------------------------------------------------------------------------- //

addClassToParentRunTimeSelectionTable(PreconditionerBase, DOMPreconditioner)

// --------------------------------------------------------------------------------- //
