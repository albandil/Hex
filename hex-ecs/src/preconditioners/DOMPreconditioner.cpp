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

#include "hex-csrmatrix.h"
#include "hex-itersolve.h"

// --------------------------------------------------------------------------------- //

#include "radial.h"

#include "DOMPreconditioner.h"

// --------------------------------------------------------------------------------- //

DOMPreconditioner::DOMPreconditioner
(
    CommandLine  const & cmd,
    InputFile    const & inp,
    Parallel     const & par,
    AngularBasis const & ang,
    Bspline const & bspline_x_inner,
    Bspline const & bspline_x_full,
    Bspline const & bspline_y_inner,
    Bspline const & bspline_y_full
) : NoPreconditioner
    (
        cmd, inp, par, ang,
        bspline_x_inner, bspline_x_full,
        bspline_y_inner, bspline_y_full
    )
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
    /// DEBUG : fixed division into sub-panels
    int xpanels = 2, ypanels = 2;
    
    // B-spline parameters
    int order = rad().bspline_inner_x().order();
    Real theta = rad().bspline_inner_x().ECStheta();
    
    // construct sub-domain bases
    std::vector<PanelSolution> p;
    for (int ixpanel = 0; ixpanel < xpanels; ixpanel++)
    for (int iypanel = 0; iypanel < ypanels; iypanel++)
    {
        // knot sub-sequences
        rArray rxknots, cxknots1, cxknots2;
        rArray ryknots, cyknots1, cyknots2;
        
        // calculate the knot sub-sequences
        knotSubsequence(ixpanel, xpanels, rad().bspline_inner_x(), rxknots, cxknots1, cxknots2);
        knotSubsequence(iypanel, ypanels, rad().bspline_inner_x(), ryknots, cyknots1, cyknots2);
        
        // create the B-spline bases for the panel
        p.emplace_back
        (
            order,
            theta,
            cxknots1, rxknots, cxknots2,
            cyknots1, ryknots, cyknots2,
            cxknots1, rxknots, cxknots2,
            cyknots1, ryknots, cyknots2,
            ang_->states().size()
        );
    }
    
    // interpolate the residual into sub-domains
    splitResidual(r, p);
    
    // find the solution on sub-domains
    int cycles = 1;
    for (int cycle = 0; cycle < cycles; cycle++)
    for (int ixpanel = 0; ixpanel < xpanels; ixpanel++)
    for (int iypanel = 0; iypanel < ypanels; iypanel++)
        solvePanel(xpanels, p, ixpanel, iypanel);
    
    // interpolate the solution from sub-domains
    collectSolution(z, p);
}

void DOMPreconditioner::finish ()
{
    NoPreconditioner::finish();
}

void DOMPreconditioner::solvePanel (int n, std::vector<PanelSolution> & p, int i, int j) const
{
    // get reference to the current panel and its neighbours
    PanelSolution * pCentre = &p[i * n + j];
    PanelSolution * pNbr [nNbrs] =
    {
        i     == 0 ? nullptr : &p[(i - 1) * n + j],
        i + 1 == n ? nullptr : &p[(i + 1) * n + j],
        j     == 0 ? nullptr : &p[i * n + (j - 1)],
        j + 1 == n ? nullptr : &p[i * n + (j + 1)]
    };
    
    // B-spline order
    int order = pCentre->xspline_inner.order();
    
    // construct the surrogate sources for all boundaries shared with other panels
    for (int nbr = 0; nbr < nNbrs; nbr++)
        surrogateSource(pCentre, nbr, pNbr[nbr]);
    
    // create the preconditioner object
    PreconditionerBase * prec = PreconditionerBase::Choose
    (
        "ILU",
        *cmd_, *inp_, *par_, *ang_,
        pCentre->xspline_inner, pCentre->xspline_full,
        pCentre->yspline_inner, pCentre->yspline_full
    );
    
    // construct the matrix of the equations etc.
    prec->setup();
    prec->update(E_);
    
    // get radial integrals
    RadialIntegrals const & rad = dynamic_cast<NoPreconditioner*>(prec)->rad();
    
    // B-spline bases
    Bspline const & xspline_inner = pCentre->xspline_inner;
    Bspline const & yspline_inner = pCentre->yspline_inner;
    
    // number of B-splines in both directions of this panel
    int Nxspline = xspline_inner.Nspline();
    int Nyspline = yspline_inner.Nspline();
    
    // contruct the right-hand side
    cBlockArray & psi = pCentre->z;
    cBlockArray   chi = pCentre->r;
    
    // add surrogate sources
    if (pNbr[Left] != nullptr)
    {
        int iR1 = pCentre->xspline_inner.iR1();
        
        for (unsigned ill = 0; ill < chi.size(); ill++)
        for (int i = iR1 - 2 * order - 1; i < iR1; i++)
        for (int j = 0; j < Nyspline; j++)
            chi[ill][i * Nyspline + j] += rad.S_inner_x(iR1 - order - 1, i) * pNbr[Left]->ssrc[Right][ill][j];
    }
    if (pNbr[Right] != nullptr)
    {
        int iR2 = pCentre->xspline_inner.iR2();
        
        for (unsigned ill = 0; ill < chi.size(); ill++)
        for (int i = iR2 - order; i < iR2 + order + 1; i++)
        for (int j = 0; j < Nyspline; j++)
            chi[ill][i * Nyspline + j] += rad.S_inner_x(iR2, i) * pNbr[Right]->ssrc[Left][ill][j];
    }
    if (pNbr[Down] != nullptr)
    {
        int iR1 = pCentre->yspline_inner.iR1();
        
        for (unsigned ill = 0; ill < chi.size(); ill++)
        for (int i = 0; i < Nxspline; i++)
        for (int j = iR1 - 2 * order - 1; j < iR1; j++)
            chi[ill][i * Nyspline + j] += rad.S_inner_y(iR1 - order - 1, j) * pNbr[Down]->ssrc[Up][ill][i];
    }
    if (pNbr[Up] != nullptr)
    {
        int iR2 = pCentre->yspline_inner.iR2();
        
        for (unsigned ill = 0; ill < chi.size(); ill++)
        for (int i = 0; i < Nxspline; i++)
        for (int j = iR2 - order; j < iR2 + order + 1; j++)
            chi[ill][i * Nyspline + j] += rad.S_inner_y(iR2, j) * pNbr[Up]->ssrc[Down][ill][i];
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
        for (std::size_t j = 0; j < A[i].size(); j++)
            A[i] = a * A[i] + b * B[i];
    };
    auto new_array = [pCentre](std::size_t N, std::string name) -> cBlockArray
    {
        cBlockArray A (N);
        for (cArray & a : A)
            a.resize(pCentre->xspline_inner.Nspline() * pCentre->yspline_inner.Nspline());
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
    CG.solve(chi, psi, cmd_->itertol, 0, 1000);
    
    // evaluate the outgoing field by subtracting the incoming field
    if (pNbr[Left] != nullptr)
    {
        int iR1 = pCentre->xspline_inner.iR1();
        Real R1 = pCentre->xspline_inner.R1();
        
        for (unsigned ill = 0; ill < chi.size(); ill++)
        for (int i = iR1 - 2 * order - 1; i < iR1; i++)
        for (int j = 0; j < Nyspline; j++)
            pCentre->outf[Right][ill][j] = psi[ill][i * Nyspline + j] * xspline_inner.bspline(i, iR1, order, R1) - pNbr[Left]->outf[Right][ill][j];
    }
    if (pNbr[Right] != nullptr)
    {
        int iR2 = pCentre->xspline_inner.iR2();
        Real R2 = pCentre->xspline_inner.R2();
        
        for (unsigned ill = 0; ill < chi.size(); ill++)
        for (int i = iR2 - order; i < iR2 + order + 1; i++)
        for (int j = 0; j < Nyspline; j++)
            pCentre->outf[Left][ill][j] = psi[ill][i * Nyspline + j] * xspline_inner.bspline(i, iR2, order, R2) - pNbr[Right]->outf[Left][ill][j];
    }
    if (pNbr[Down] != nullptr)
    {
        int iR1 = pCentre->yspline_inner.iR1();
        Real R1 = pCentre->yspline_inner.R1();
        
        for (unsigned ill = 0; ill < chi.size(); ill++)
        for (int i = 0; i < Nxspline; i++)
        for (int j = iR1 - 2 * order - 1; j < iR1; j++)
            pCentre->outf[Down][ill][i] = psi[ill][i * Nyspline + j] * yspline_inner.bspline(j, iR1, order, R1) - pNbr[Down]->outf[Up][ill][i];
    }
    if (pNbr[Up] != nullptr)
    {
        int iR2 = pCentre->yspline_inner.iR2();
        Real R2 = pCentre->yspline_inner.R2();
        
        for (unsigned ill = 0; ill < chi.size(); ill++)
        for (int i = 0; i < Nxspline; i++)
        for (int j = iR2 - order; j < iR2 + order + 1; j++)
            pCentre->outf[Up][ill][i] = psi[ill][i * Nyspline + j] * yspline_inner.bspline(j, iR2, order, R2) - pNbr[Up]->outf[Down][ill][i];
    }
}

cArray DOMPreconditioner::surrogateSource (PanelSolution * panel, int direction, PanelSolution * neighbour) const
{
    if (direction == Left)
    {
        // original B-spline basis
        Bspline const & org_xspline = panel->xspline_inner;
        
        // short B-spline basis
        Bspline xspline
        (
            org_xspline.order(),
            org_xspline.ECStheta(),
            org_xspline.cknots1(),
            org_xspline.rknots().slice(0, org_xspline.order() + 1),
            org_xspline.cknots2() - org_xspline.cknots2().front() + org_xspline.rknots().front(org_xspline.order())
        );
        Bspline const & yspline = panel->yspline_inner;
        
        // calculate radial integrals
        RadialIntegrals rint (xspline, yspline, xspline, yspline, 0);
        
        // for all angular states
        for (unsigned ill = 0; ill < ang_->states().size(); ill++)
        {
            // angular momenta
            int l1 = ang_->states()[ill].first;
            int l2 = ang_->states()[ill].second;
        
            // construct the matrix of the free equations
            BlockSymBandMatrix<Complex> A =
                kron(Complex(E_) * rint.S_inner_x(), rint.S_inner_y())
                - kron(0.5_z * rint.D_inner_x(), rint.S_inner_x())
                - kron(0.5_z * rint.S_inner_x(), rint.D_inner_x())
                - kron(Complex(0.5 * l1 * (l1 + 1)) * rint.Mm2_inner_x(), rint.S_inner_x())
                - kron(Complex(0.5 * l2 * (l2 + 1)) * rint.S_inner_x(), rint.Mm2_inner_x());
            
            // add additional matrix elements
            CooMatrix<LU_int_t,Complex> A_coo = A.tocoo<LU_int_t>();
            A_coo.resize
            (
                A_coo.rows() + yspline.Nspline(),
                A_coo.cols() + yspline.Nspline()
            );
            for (int i = xspline.iR1() - 2*xspline.order() - 1; i < xspline.iR1(); i++)
            for (int j = 0; j < yspline.Nspline(); j++)
            {
                A_coo.add
                (
                    i * yspline.Nspline() + j,
                    xspline.Nspline() * yspline.Nspline() + j,
                    -rint.S_inner_x(i, xspline.iR1() - 1)
                );
            }
            for (int j = 0; j < yspline.Nspline(); j++)
            for (int k = xspline.iR1() - xspline.order() - 1; k < xspline.iR1(); k++)
            {
                A_coo.add
                (
                    xspline.Nspline() * yspline.Nspline() + j,
                    k * yspline.Nspline() + j,
                    xspline.bspline(k, xspline.iR1(), xspline.order(), xspline.R1())
                );
            }
            
            // calculate the right-hand side
            cArray rhs (xspline.Nspline() * yspline.Nspline() + yspline.Nspline());
            for (int i = neighbour->xspline_inner.iR2() - neighbour->xspline_inner.order() - 1; i < neighbour->xspline_inner.iR2(); j++)
            for (int j = 0; j < neighbour->yspline_inner.Nspline(); j++)
            {
                rhs[xspline.Nspline() * yspline.Nspline() + j] += neighbour->outf[Right][ill][j]
                    * neighbour->xspline_inner.bspline(i, neighbour->xspline_inner.iR2(), neighbour->xspline_inner.order(), neighbour->xspline_inner.R2());
            }
            
            // factorize and solve
            cArray solution (rhs.size());
            CsrMatrix<LU_int_t,Complex> A_csr = A_coo.tocsr();
            std::unique_ptr<LUft> A_lu;
            A_lu.reset(LUft::Choose("umfpack"));
            A_lu->factorize(A_csr);
            A_lu->solve(rhs, solution, 1);
            
            // store the surrogate source
            for (int j = 0; j < neighbour->yspline_inner.Nspline(); j++)
                panel->ssrc[Left][ill][j] = solution[xspline.Nspline() * yspline.Nspline() + j];
        }
    }
    if (direction == Right)
    {
        // TODO
    }
    if (direction == Down)
    {
        // TODO
    }
    if (direction == Up)
    {
        // TODO
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
    
    // new knot sub-sequence
    rknots = inp_->rknots.slice(iRa, iRb + 1);
    cknots1 = inp_->cknots - inp_->cknots.back() + rknots.front();
    cknots2 = inp_->cknots + rknots.back();
}

void DOMPreconditioner::splitResidual (cBlockArray const & r, std::vector<PanelSolution> & panels) const
{
    for (int ixpanel = 0; ixpanel < 2; ixpanel++)
    for (int iypanel = 0; iypanel < 2; iypanel++)
    for (unsigned ill = 0; ill < r.size(); ill++)
    {
        PanelSolution & p = panels[ixpanel * 2 + iypanel];
        
        std::size_t N = p.xspline_inner.Nspline() * p.yspline_inner.Nspline();
        
        cArray rab0 = kron_dot(p.SaF, p.SbF, r[ill]), rab1(N), rab2(N);
        p.lu_Sbb->solve(rab0, rab1, p.xspline_inner.Nspline());
        transpose(rab1);
        p.lu_Saa->solve(rab1, rab2, p.yspline_inner.Nspline());
        transpose(rab2);
        
        p.r[ill] = rab2;
    }
}

void DOMPreconditioner::collectSolution (cBlockArray & z, std::vector<PanelSolution> & panels) const
{
    // TODO
}

DOMPreconditioner::PanelSolution::PanelSolution
(
    int order,
    Real theta,
    rArray cxspline1_inner, rArray rxspline_inner, rArray cxspline2_inner,
    rArray cyspline1_inner, rArray ryspline_inner, rArray cyspline2_inner,
    rArray cxspline1_full,  rArray rxspline_full,  rArray cxspline2_full,
    rArray cyspline1_full,  rArray ryspline_full,  rArray cyspline2_full,
    int Nang
) : xspline_inner (order, theta, cxspline1_inner, rxspline_inner, cxspline2_inner),
    yspline_inner (order, theta, cyspline1_inner, ryspline_inner, cyspline2_inner),
    xspline_full (order, theta, cxspline1_full, rxspline_full, cxspline2_full),
    yspline_full (order, theta, cyspline1_full, ryspline_full, cyspline2_full),
    r (Nang), z (Nang)
{
    for (int nbr = 0; nbr < nNbrs; nbr++)
    {
        ssrc[nbr] = cBlockArray(Nang);
        outf[nbr] = cBlockArray(Nang);
    }
}

// --------------------------------------------------------------------------------- //

addClassToParentRunTimeSelectionTable(PreconditionerBase, DOMPreconditioner)

// --------------------------------------------------------------------------------- //
