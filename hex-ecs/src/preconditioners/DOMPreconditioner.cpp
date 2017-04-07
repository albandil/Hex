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
    std::array<PanelSolution*,nNbrs> pNbr;
    pNbr[Left]  = (i     == 0 ? nullptr : &p[(i - 1) * n + j]);
    pNbr[Right] = (i + 1 == n ? nullptr : &p[(i + 1) * n + j]);
    pNbr[Down]  = (j     == 0 ? nullptr : &p[i * n + (j - 1)]);
    pNbr[Up]    = (j + 1 == n ? nullptr : &p[i * n + (j + 1)]);
    
    // B-spline order
    int order = pCentre->xspline_inner.order();
    
    // construct the surrogate sources for all boundaries shared with other panels
    surrogateSource(pCentre, pNbr);
    
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
        
        for (unsigned ill = 0; ill < chi.size(); ill++)
        for (int i = iR1 - order; i <= iR1 + order; i++)
        for (int j = 0; j < Nyspline; j++)
        for (int l = (j < order ? 0 : j - order); l <= j + order and l < Nyspline; l++)
            chi[ill][i * Nyspline + j] += rad.S_x(i, iR1) * rad.S_y(j, l) * pCentre->ssrc[Left][ill][l];
    }
    if (pNbr[Right] != nullptr)
    {
        int iR2 = pCentre->xspline_inner.iR2();
        
        for (unsigned ill = 0; ill < chi.size(); ill++)
        for (int i = iR2 - 2*order - 1; i < iR2; i++)
        for (int j = 0; j < Nyspline; j++)
        for (int l = (j < order ? 0 : j - order); l <= j + order and l < Nyspline; l++)
            chi[ill][i * Nyspline + j] += rad.S_x(i, iR2 - order - 1) * rad.S_y(j, l) * pCentre->ssrc[Right][ill][l];
    }
    if (pNbr[Down] != nullptr)
    {
        int iR1 = pCentre->yspline_inner.iR1();
        
        for (unsigned ill = 0; ill < chi.size(); ill++)
        for (int i = 0; i < Nxspline; i++)
        for (int j = iR1 - order; j <= iR1 + order; j++)
        for (int k = (i < order ? 0 : i - order); k <= i + order and k < Nxspline; k++)
            chi[ill][i * Nyspline + j] += rad.S_x(i, k) * rad.S_y(iR1, j) * pCentre->ssrc[Down][ill][k];
    }
    if (pNbr[Up] != nullptr)
    {
        int iR2 = pCentre->yspline_inner.iR2();
        
        for (unsigned ill = 0; ill < chi.size(); ill++)
        for (int i = 0; i < Nxspline; i++)
        for (int j = iR2 - 2*order - 1; j < iR2; j++)
        for (int k = (i < order ? 0 : i - order); k <= i + order and k < Nxspline; k++)
            chi[ill][i * Nyspline + j] += rad.S_x(i, k) * rad.S_y(iR2 - order - 1, j) * pCentre->ssrc[Up][ill][k];
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

void DOMPreconditioner::surrogateSource
(
    PanelSolution * panel,
    std::array<PanelSolution*,nNbrs> & pNbr
) const
{
    std::cout << "\tConstruct surrogate source" << std::endl;
    
    // B-spline parameters
    unsigned order = panel->xspline_inner.order();
    Real ecstheta = panel->xspline_inner.ECStheta();
    
    // empty deleter used below to prevent deallocation of non-owned pointers
    auto nodelete = [](Bspline*){};
    
    // B-spline bases for all four edges
    std::array<std::shared_ptr<Bspline>,nNbrs> xspline, yspline;
    if (pNbr[Left] != nullptr)
    {
        rArray const & cknots1 = panel->xspline_inner.cknots1();
        rArray const & rknots  = panel->xspline_inner.rknots ();
        rArray const & cknots2 = panel->xspline_inner.cknots2();
        
        /*xspline[Left].reset
        (
            new Bspline
            (
                order,
                ecstheta,
                cknots1,
                rknots.slice(0, 2*order + 3),
                cknots2 - cknots2.front() + rknots.front(2*order + 2)
            )
        );*/
        
        xspline[Left].reset(&panel->xspline_inner, nodelete);
        yspline[Left].reset(&panel->yspline_inner, nodelete);
    }
    if (pNbr[Right] != nullptr)
    {
        rArray const & cknots1 = panel->xspline_inner.cknots1();
        rArray const & rknots  = panel->xspline_inner.rknots ();
        rArray const & cknots2 = panel->xspline_inner.cknots2();
        
        /*xspline[Right].reset
        (
            new Bspline
            (
                order,
                ecstheta,
                cknots1 - cknots1.back() + rknots.back(2*order + 2),
                rknots.slice(rknots.size() - 2*order - 3, rknots.size()),
                cknots2
            )
        );*/
        
        xspline[Right].reset(&panel->xspline_inner, nodelete);
        yspline[Right].reset(&panel->yspline_inner, nodelete);
    }
    if (pNbr[Down] != nullptr)
    {
        rArray const & cknots1 = panel->yspline_inner.cknots1();
        rArray const & rknots  = panel->yspline_inner.rknots ();
        rArray const & cknots2 = panel->yspline_inner.cknots2();
        
        xspline[Down].reset(&panel->xspline_inner, nodelete);
        yspline[Down].reset(&panel->yspline_inner, nodelete);
        
        /*yspline[Down].reset
        (
            new Bspline
            (
                order,
                ecstheta,
                cknots1,
                rknots.slice(0, 2*order + 3),
                cknots2 - cknots2.front() + rknots.front(2*order + 2)
            )
        );*/
    }
    if (pNbr[Up] != nullptr)
    {
        rArray const & cknots1 = panel->yspline_inner.cknots1();
        rArray const & rknots  = panel->yspline_inner.rknots ();
        rArray const & cknots2 = panel->yspline_inner.cknots2();
        
        xspline[Up].reset(&panel->xspline_inner, nodelete);
        yspline[Up].reset(&panel->yspline_inner, nodelete);
        
        /*yspline[Up].reset
        (
            new Bspline
            (
                order,
                ecstheta,
                cknots1 - cknots1.back() + rknots.back(2*order + 2),
                rknots.slice(rknots.size() - 2*order - 3, rknots.size()),
                cknots2
            )
        );*/
    }
    
    // calculate radial integrals
    std::array<std::shared_ptr<RadialIntegrals>,nNbrs> rint;
    for (int nbr = 0; nbr < nNbrs; nbr++) if (pNbr[nbr] != nullptr)
    {
        rint[nbr].reset(new RadialIntegrals(*xspline[nbr], *yspline[nbr], ang_->maxlambda() + 1));
        rint[nbr]->verbose(false);
        rint[nbr]->setupOneElectronIntegrals(*par_, *cmd_);
        rint[nbr]->setupTwoElectronIntegrals(*par_, *cmd_);
    }
    
    // B-spline basis sizes and other dimensions
    std::array<std::size_t,nNbrs> Nxspline, Nyspline;
    std::array<std::size_t,nNbrs> iR1x, iR2x, iR1y, iR2y;
    for (int nbr = 0; nbr < nNbrs; nbr++) if (pNbr[nbr] != nullptr)
    {
        // number of B-splines in the basis
        Nxspline[nbr] = xspline[nbr]->Nspline();
        Nyspline[nbr] = yspline[nbr]->Nspline();
        
        // B-spline knots that bound the real grids
        iR1x[nbr] = xspline[nbr]->iR1();
        iR2x[nbr] = xspline[nbr]->iR2();
        iR1y[nbr] = yspline[nbr]->iR1();
        iR2y[nbr] = yspline[nbr]->iR2();
    }
    
    // indices of the first and (one after) the last B-splines that contribute to the matching line
    std::array<std::size_t,nNbrs> iBeg, iEnd, ssrcsize;
    for (int nbr = 0; nbr < nNbrs; nbr++) if (pNbr[nbr] != nullptr)
    {
        // vertical edges
        if (nbr == Left or nbr == Right)
        {
            iBeg[nbr] = (pNbr[Down] ? iR1y[nbr] + 1 : 0);
            iEnd[nbr] = (pNbr[Up]   ? iR2y[nbr] - order - 1 : iR2y[nbr] - 1);
            ssrcsize[nbr] = iEnd[nbr] - iBeg[nbr];
        }
        
        // horizontal edges
        if (nbr == Down or nbr == Up)
        {
            iBeg[nbr] = (pNbr[Left]  ? iR1x[nbr] + 1 : 0);
            iEnd[nbr] = (pNbr[Right] ? iR2x[nbr] - order - 1 : iR2x[nbr] - 1);
            ssrcsize[nbr] = iEnd[nbr] - iBeg[nbr];
        }
    }
    
    // knot indices of the matching line
    std::array<std::size_t,nNbrs> iMline;
    iMline[Left]  = (pNbr[Left]  ? iR1x[Left]  + order + 1 : iR1x[Left] );
    iMline[Right] = (pNbr[Right] ? iR2x[Right] - order - 1 : iR2x[Right]);
    iMline[Down]  = (pNbr[Down]  ? iR1y[Down]  + order + 1 : iR1y[Down] );
    iMline[Up]    = (pNbr[Up]    ? iR2y[Up]    - order - 1 : iR2y[Up]   );
    
    // coordinate positions of the matching line
    std::array<Real,nNbrs> Mline;
    Mline[Left]  = (pNbr[Left]  ? xspline[Left] ->t(iMline[Left]) .real() : 0);
    Mline[Right] = (pNbr[Right] ? xspline[Right]->t(iMline[Right]).real() : 0);
    Mline[Down]  = (pNbr[Down]  ? yspline[Down] ->t(iMline[Down]) .real() : 0);
    Mline[Up]    = (pNbr[Up]    ? yspline[Up]   ->t(iMline[Up])   .real() : 0);
    
    // position of the surrogate source
    std::array<std::size_t,nNbrs> iSrc;
    iSrc[Left]  = (pNbr[Left]  ? iR1x[Left]              : 0);
    iSrc[Right] = (pNbr[Right] ? iR2x[Right] - order - 1 : 0);
    iSrc[Down]  = (pNbr[Down]  ? iR1y[Down]              : 0);
    iSrc[Up]    = (pNbr[Up]    ? iR2y[Up]    - order - 1 : 0);
    
    // each of the edges needs also the overlap matrix with the support
    // between the complex grid and the mathing line
    std::array<SymBandMatrix<Complex>,nNbrs> Swgt;
    if (pNbr[Left] != nullptr)
    {
        Swgt[Left] = SymBandMatrix<Complex>(Nxspline[Left], order + 1).populate
        (
            [&](int m, int n)
            {
                return rint[Left]->computeM(rint[Left]->bspline_x(), rint[Left]->gaussleg_x(), 0, m, n, iR1x[Left], iMline[Left], false);
            }
        );
    }
    if (pNbr[Right] != nullptr)
    {
        Swgt[Right] = SymBandMatrix<Complex>(Nxspline[Right], order + 1).populate
        (
            [&](int m, int n)
            {
                return rint[Right]->computeM(rint[Right]->bspline_x(), rint[Right]->gaussleg_x(), 0, m, n, iMline[Right], iR2x[Right], false);
            }
        );
    }
    if (pNbr[Down] != nullptr)
    {
        Swgt[Down] = SymBandMatrix<Complex>(Nyspline[Down], order + 1).populate
        (
            [&](int m, int n)
            {
                return rint[Down]->computeM(rint[Down]->bspline_y(), rint[Down]->gaussleg_y(), 0, m, n, iR1y[Down], iMline[Down], false);
            }
        );
    }
    if (pNbr[Up] != nullptr)
    {
        Swgt[Up] = SymBandMatrix<Complex>(Nyspline[Up], order + 1).populate
        (
            [&](int m, int n)
            {
                return rint[Up]->computeM(rint[Up]->bspline_y(), rint[Up]->gaussleg_y(), 0, m, n, iMline[Up], iR2y[Up], false);
            }
        );
    }
    
    // for all angular states
    for (unsigned ill = 0; ill < ang_->states().size(); ill++)
    {
        // angular momenta of the two electrons
        int l1 = ang_->states()[ill].first;
        int l2 = ang_->states()[ill].second;
        
        //
        // Construct matrices of the free equations for all four grids
        //
        
            std::array<BlockSymBandMatrix<Complex>,nNbrs> A;
            
            for (int nbr = 0; nbr < nNbrs; nbr++) if (pNbr[nbr] != nullptr)
            {
                A[nbr] = kron(Complex(E_) * rint[nbr]->S_x(), rint[nbr]->S_y())
                    - kron(0.5_z * rint[nbr]->D_x(), rint[nbr]->S_y())
                    - kron(0.5_z * rint[nbr]->S_x(), rint[nbr]->D_y())
                    - kron(Complex(0.5 * l1 * (l1 + 1)) * rint[nbr]->Mm2_x(), rint[nbr]->S_y())
                    - kron(Complex(0.5 * l2 * (l2 + 1)) * rint[nbr]->S_x(), rint[nbr]->Mm2_y())
                    + kron(rint[nbr]->Mm1_x(), rint[nbr]->S_y())
                    + kron(rint[nbr]->S_x(), rint[nbr]->Mm1_y());
            
                for (int lambda = 0; lambda <= ang_->maxlambda(); lambda++)
                    A[nbr] -= Complex(ang_->f(ill, ill, lambda)) * rint[nbr]->R_tr_dia(lambda);
            }
        
        //
        // Own surrogate source coupling matrices (-C)
        //
        
            // surrogate coupling blocks, -C are diagonal blocks, -X offdiagonal blocks
            std::array<std::array<CooMatrix<LU_int_t,Complex>,nNbrs>,nNbrs> X;
            
            // for all edges
            for (int nbr = 0; nbr < nNbrs; nbr++) if (pNbr[nbr] != nullptr)
            {
                // resize the matrix according to the grid
                X[nbr][nbr].resize(Nxspline[nbr] * Nyspline[nbr], ssrcsize[nbr]);
                
                // vertical matching line
                if (nbr == Left or nbr == Right)
                {
                    // for all B-spline elements of the surrogate source
                    for (std::size_t i = iSrc[nbr]; i <= iSrc[nbr]; i++) // 1 B-spline wide
                    for (std::size_t j = iBeg[nbr]; j < iEnd[nbr]; j++)
                    {
                        // for all solution B-spline that overlap this source B-spline
                        for (std::size_t a = (i < order ? 0 : i - order); a <= i + order and a < Nxspline[nbr]; a++)
                        for (std::size_t b = (j < order ? 0 : j - order); b <= j + order and b < Nyspline[nbr]; b++)
                            X[nbr][nbr].add(a * Nyspline[nbr] + b, j - iBeg[nbr], -rint[nbr]->S_x(a,i) * rint[nbr]->S_y(b,j));
                    }
                }
                
                // horizontal matching line
                if (nbr == Down or nbr == Up)
                {
                    // for all B-spline elements of the surrogate source
                    for (std::size_t i = iBeg[nbr]; i < iEnd[nbr]; i++)
                    for (std::size_t j = iSrc[nbr]; j <= iSrc[nbr]; j++) // 1 B-spline wide
                    {
                        // for all solution B-spline that overlap this source B-spline
                        for (std::size_t a = (i < order ? 0 : i - order); a <= i + order and a < Nxspline[nbr]; a++)
                        for (std::size_t b = (j < order ? 0 : j - order); b <= j + order and b < Nyspline[nbr]; b++)
                            X[nbr][nbr].add(a * Nyspline[nbr] + b, i - iBeg[nbr], -rint[nbr]->S_x(a,i) * rint[nbr]->S_y(b,j));
                    }
                }
            }
        
        //
        // Neighbour surrogate source coupling matrices (-X)
        //
        
            // for all edges and their neighbours
            for (int nbr = 0; nbr < nNbrs; nbr++) if (pNbr[nbr] != nullptr)
            for (int mbr = 0; mbr < nNbrs; mbr++) if (pNbr[mbr] != nullptr)
            {
                // coupling to self already done in matrix -C
                if (nbr == mbr)
                    continue;
                
                // coupling to opposite edge is null
                if (nbr == rev(mbr))
                    continue;
                
                // resize the matrix according to the grid
                X[nbr][mbr].resize(Nxspline[nbr] * Nyspline[nbr], ssrcsize[mbr]);
                
                // vertical matching line
                if (nbr == Left or nbr == Right)
                {
                    // for all relevant B-splines of the neighbour surrogate source
                    for (std::size_t k = iMline[nbr] - order; k < iMline[nbr]; k++)
                    for (std::size_t l = iSrc[mbr]; l <= iSrc[mbr]; l++) // 1 B-spline wide
                    {
                        // for all solution B-splines that overlap this source B-spline
                        for (std::size_t i = (k < order ? 0 : k - order); i <= k + order and i < Nxspline[nbr]; i++)
                        for (std::size_t j = (l < order ? 0 : l - order); j <= l + order and j < Nyspline[nbr]; j++)
                            X[nbr][mbr].add(i * Nyspline[nbr] + j, k - iBeg[mbr], -Swgt[nbr](i,k) * rint[nbr]->S_y(j,l));
                    }
                }
                
                // horizontal matching line
                if (nbr == Down or nbr == Up)
                {
                    // for all relevant B-splines of the neighbour surrogate source
                    for (std::size_t k = iSrc[mbr]; k <= iSrc[mbr]; k++) // 1 B-spline wide
                    for (std::size_t l = iMline[nbr] - order; l < iMline[nbr]; l++)
                    {
                        // for all solution B-splines that overlap this source B-spline
                        for (std::size_t i = (k < order ? 0 : k - order); i <= k + order and i < Nxspline[nbr]; i++)
                        for (std::size_t j = (l < order ? 0 : l - order); j <= l + order and j < Nyspline[nbr]; j++)
                            X[nbr][mbr].add(i * Nyspline[nbr] + j, l - iBeg[mbr], -rint[nbr]->S_x(i,k) * Swgt[nbr](j,l));
                    }
                }
            }
        
        //
        // Boundary field coupling matrices (D)
        //
        
            std::array<CooMatrix<LU_int_t,Complex>,nNbrs> D;
            
            for (int nbr = 0; nbr < nNbrs; nbr++) if (pNbr[nbr] != nullptr)
            {
                D[nbr].resize(ssrcsize[nbr], Nxspline[nbr] * Nyspline[nbr]);
                
                // vertical matching line
                if (nbr == Left or nbr == Right)
                {
                    for (std::size_t j = iBeg[nbr]; j < iEnd[nbr]; j++)
                    for (std::size_t k = iMline[nbr] - order; k < iMline[nbr]; k++)
                        D[nbr].add(j - iBeg[nbr], k * Nyspline[nbr] + j, xspline[nbr]->bspline(k, iMline[nbr], order, Mline[nbr]));
                }
                
                // horizontal matching line
                if (nbr == Down or nbr == Up)
                {
                    for (std::size_t i = iBeg[nbr]; i < iEnd[nbr]; i++)
                    for (std::size_t l = iMline[nbr] - order; l < iMline[nbr]; l++)
                        D[nbr].add(i - iBeg[nbr], i * Nyspline[nbr] + l, yspline[nbr]->bspline(l, iMline[nbr], order, Mline[nbr]));
                }
            }
        
        //
        // Merge the blocks
        //
        
            NumberArray<LU_int_t> I, J; cArray V;
            std::size_t A_rows = 0, A_cols = 0, D_rows = 0, D_cols = 0, C_rows = 0, C_cols = 0;
            
            for (int nbr = 0; nbr < nNbrs; nbr++) if (pNbr[nbr] != nullptr)
            {
                // append elements of A
                CooMatrix<LU_int_t,Complex> A_coo = A[nbr].tocoo<LU_int_t>();
                I.reserve(I.size() + A_coo.i().size());
                J.reserve(J.size() + A_coo.j().size());
                V.reserve(V.size() + A_coo.v().size());
                for (std::size_t i = 0; i < A_coo.i().size(); i++)
                {
                    I.push_back(A_coo.i()[i] + A_rows);
                    J.push_back(A_coo.j()[i] + A_cols);
                    V.push_back(A_coo.v()[i]);
                }
                
                // update matrix size
                A_rows += Nxspline[nbr] * Nyspline[nbr]; D_rows = A_rows;
                A_cols += Nxspline[nbr] * Nyspline[nbr]; C_cols = A_cols;
            }
            
            for (int nbr = 0; nbr < nNbrs; nbr++) if (pNbr[nbr] != nullptr)
            {
                // reset block column
                C_cols = A_cols;
                
                // for all C and X blocks in this blocks row
                for (int mbr = 0; mbr < nNbrs; mbr++) if (pNbr[mbr] != nullptr)
                {
                    // append elements of -C and -X, if available
                    if (nbr == mbr) /// DEBUG
                    //if (nbr != rev(mbr))
                    {
                        I.reserve(I.size() + X[nbr][mbr].i().size());
                        J.reserve(J.size() + X[nbr][mbr].j().size());
                        V.reserve(V.size() + X[nbr][mbr].v().size());
                        for (std::size_t i = 0; i < X[nbr][mbr].i().size(); i++)
                        {
                            I.push_back(X[nbr][mbr].i()[i] + C_rows);
                            J.push_back(X[nbr][mbr].j()[i] + C_cols);
                            V.push_back(X[nbr][mbr].v()[i]);
                        }
                    }
                    
                    // move to next block column
                    C_cols += ssrcsize[mbr];
                }
                
                // move to next block row
                C_rows += Nxspline[nbr] * Nyspline[nbr];
            }
            
            for (int nbr = 0; nbr < nNbrs; nbr++) if (pNbr[nbr] != nullptr)
            {
                // append elements of D
                I.reserve(I.size() + D[nbr].i().size());
                J.reserve(J.size() + D[nbr].j().size());
                V.reserve(V.size() + D[nbr].v().size());
                for (std::size_t i = 0; i < D[nbr].i().size(); i++)
                {
                    I.push_back(D[nbr].i()[i] + D_rows);
                    J.push_back(D[nbr].j()[i] + D_cols);
                    V.push_back(D[nbr].v()[i]);
                }
                
                // update matrix size
                D_rows += ssrcsize[nbr];
                D_cols += Nxspline[nbr] * Nyspline[nbr];
            }
            
            // construct the final matrix
            CooMatrix<LU_int_t,Complex> A_coo (D_rows, C_cols, std::move(I), std::move(J), std::move(V));
        
        //
        // Righ-hand side
        //
        
            cArray rhs (A_coo.rows());
            
            // offset rhs pointer to the start of surrogate sources segment
            Complex * rhs_ptr = rhs.data();
            for (int nbr = 0; nbr < nNbrs; nbr++) if (pNbr[nbr] != nullptr)
                rhs_ptr += Nxspline[nbr] * Nyspline[nbr];
            
            // evaluate the neighbour fields
            for (int nbr = 0; nbr < nNbrs; nbr++) if (pNbr[nbr] != nullptr)
            {
                // panel's and neighbour's boundary field
                cArray const & ownf = panel->outf[nbr][ill];
                cArray const & nbrf = pNbr[nbr]->outf[rev(nbr)][ill];
                
                // copy neighbour field difference to right-hand side
                for (std::size_t i = iBeg[nbr]; i < iEnd[nbr]; i++)
                    *rhs_ptr++ = nbrf[i] - ownf[i];
            }
        
        //
        // Solve the equations
        //
            
            // convert COO matrix to CSR
            CsrMatrix<LU_int_t,Complex> A_csr = A_coo.tocsr();
            std::unique_ptr<LUft> A_lu;
            
            // calculate the LU decomposition
            A_lu.reset(LUft::Choose("umfpack")); // TODO : Make runtime selectable.
            A_lu->factorize(A_csr);
            
            // solve the system by application of the LU
            cArray solution (A_coo.cols());
            A_lu->solve(rhs, solution, 1);
            
            // DEBUG
            static int ble = 0; ble++;
            Complex * ptr = solution.data();
            for (int nbr = 0; nbr < nNbrs; nbr++) if (pNbr[nbr] != nullptr)
            {
                std::ofstream ofs (format("dom-%d-[%d]-ispsi-%d.vtk", ble, nbr, ill));
                writeVTK_points
                (
                    ofs,
                    Bspline::zip
                    (
                        *xspline[nbr],
                        *yspline[nbr],
                        cArrayView(Nxspline[nbr] * Nyspline[nbr], ptr),
                        linspace(xspline[nbr]->Rmin(), xspline[nbr]->Rmax(), 1001),
                        linspace(yspline[nbr]->Rmin(), yspline[nbr]->Rmax(), 1001)
                    ),
                    linspace(xspline[nbr]->Rmin(), xspline[nbr]->Rmax(), 1001),
                    linspace(yspline[nbr]->Rmin(), yspline[nbr]->Rmax(), 1001),
                    rArray{0}
                );
                ptr += Nxspline[nbr] * Nyspline[nbr];
            }
        
        //
        // Update the surrogate sources
        //
            
            // offset solution pointer to the start of surrogate sources segment
            Complex const * sol = solution.data();
            for (int nbr = 0; nbr < nNbrs; nbr++) if (pNbr[nbr] != nullptr)
                sol += Nxspline[nbr] * Nyspline[nbr];
            
            // extract individual surrogate sources
            for (int nbr = 0; nbr < nNbrs; nbr++) if (pNbr[nbr] != nullptr)
            for (std::size_t i = iBeg[nbr]; i < iEnd[nbr]; i++)
                panel->ssrc[nbr][ill][i] += *sol++;
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
