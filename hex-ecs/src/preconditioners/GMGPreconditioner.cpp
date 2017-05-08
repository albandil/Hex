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

#include "hex-matrix.h"

// --------------------------------------------------------------------------------- //

#include "GMGPreconditioner.h"

// --------------------------------------------------------------------------------- //

std::string GMGPreconditioner::description () const
{
    return "Geometric multigrid.";
}

template <class T> NumberArray<T> dither (NumberArray<T> const & arr)
{
    NumberArray<T> brr;
    
    // copy al zeros, the last knot and also every second knot from arr to brr
    for (std::size_t i = 0; i < arr.size(); i++)
    {
        if (arr[i] == 0 or i + 1 == arr.size() or i % 2 == 0)
            brr.push_back(arr[i]);
    }
    
    return brr;
}

GMGPreconditioner::GMGPreconditioner
(
    CommandLine  const & cmd,
    InputFile    const & inp,
    Parallel     const & par,
    AngularBasis const & ang,
    Bspline const & bspline_x_inner,
    Bspline const & bspline_x_full,
    Bspline const & bspline_y_inner,
    Bspline const & bspline_y_full
) : GMGPreconditioner
    (
        cmd, inp, par, ang,
        bspline_x_inner, bspline_x_full,
        bspline_y_inner, bspline_y_full,
        cmd.multigrid_depth
    )
{
    // nothing
}

GMGPreconditioner::GMGPreconditioner
(
    CommandLine  const & cmd,
    InputFile    const & inp,
    Parallel     const & par,
    AngularBasis const & ang,
    Bspline const & bspline_x_inner,
    Bspline const & bspline_x_full,
    Bspline const & bspline_y_inner,
    Bspline const & bspline_y_full
    int level
) : NoPreconditioner
    (
        cmd, inp, par, ang,
        bspline_x_inner, bspline_x_full,
        bspline_y_inner, bspline_y_full
    ),
    level_(level),
    bspline_inner_fine_(bspline_x_inner),
    bspline_full_fine_(bspline_x_full),
    bspline_inner_coarse_
    (
        bspline_x_inner.order(),
        dither(bspline_x_inner.rknots()),
        bspline_x_inner.ECStheta(),
        dither(bspline_x_inner.cknots1())
    ),
    bspline_full_coarse_
    (
        bspline_x_full.order(),
        dither(bspline_x_full.rknots()),
        bspline_x_full.ECStheta(),
        dither(bspline_x_full.cknots1())
    )
{
    // setup subgrid
    if (level_ > 0)
    {
        subgrid_ = new GMGPreconditioner
        (
            par, inp, ll,
            bspline_inner_coarse_,
            bspline_full_coarse_,
            cmd, level_ - 1
        );
    }
    else
    {
        subgrid_ = Preconditioners::choose
        (
            cmd.multigrid_coarse_prec,
            par, inp, ll,
            bspline_inner_fine_,
            bspline_full_fine_,
            cmd
        );
    }
}

GMGPreconditioner::~GMGPreconditioner ()
{
    delete subgrid_;
}

void GMGPreconditioner::setup ()
{
    // setup subgrid
    subgrid_->setup();
    
    // setup parent
    NoPreconditioner::setup();
    
    if (level_ > 0)
    {
        std::cout << "Setting up GMG preconditioner level " << level_ << " (" << bspline_full_coarse_.Nspline() << " B-splines)" << std::endl << std::endl;
        
        // create integrator
        GaussLegendre g;
        
        // calculate overlap matrix between the fine and coarse bases
        ColMatrix<Complex> Sfc (bspline_inner_fine_.Nspline(), bspline_inner_coarse_.Nspline());
        for (int i = 0; i < bspline_inner_fine_.Nspline(); i++)
        for (int j = 0; j < bspline_inner_coarse_.Nspline(); j++)
            Sfc(i,j) = rad_.computeS12(g, bspline_inner_fine_, bspline_inner_coarse_, i, j);
        ColMatrix<Complex> Scf (Sfc.T());
        
        // factorize inner overlap matrix for both the fine and coarse basis
        SymBandMatrix<Complex> S_fine = this->rad_.S_inner();
        SymBandMatrix<Complex> S_coarse = dynamic_cast<GMGPreconditioner*>(subgrid_)->rad().S_inner();
        
        CooMatrix<LU_int_t,Complex> S_fine_coo = S_fine.tocoo<LU_int_t>();
        CooMatrix<LU_int_t,Complex> S_coarse_coo = S_coarse.tocoo<LU_int_t>();
        
        CsrMatrix<LU_int_t,Complex> S_fine_csr = S_fine_coo.tocsr();
        CsrMatrix<LU_int_t,Complex> S_coarse_csr = S_coarse_coo.tocsr();
        
        std::shared_ptr<LUft<LU_int_t,Complex>> lu_S_inner_fine = S_fine_csr.factorize();
        std::shared_ptr<LUft<LU_int_t,Complex>> lu_S_inner_coarse = S_coarse_csr.factorize();
        
        // create restrictors and prolongators
        ColMatrix<Complex> restrictor_inner (bspline_inner_coarse_.Nspline(), bspline_inner_fine_.Nspline());
        ColMatrix<Complex> prolongator_inner (bspline_inner_fine_.Nspline(), bspline_inner_coarse_.Nspline());
        lu_S_inner_coarse->solve(Scf.data(), restrictor_inner.data(), bspline_inner_fine_.Nspline());
        lu_S_inner_fine->solve(Sfc.data(), prolongator_inner.data(), bspline_inner_coarse_.Nspline());
        restrictor_inner_ = RowMatrix<Complex>(restrictor_inner);
        prolongator_inner_ = RowMatrix<Complex>(prolongator_inner);
    }
}

void GMGPreconditioner::update (Real E)
{
    // update subgrid
    subgrid_->update(E);
    
    // update parent
    NoPreconditioner::update(E);
    
    if (level_ > 0)
    {
        // calculate full matrix diagonal for use in Gauss-Seidel iterations
        D.resize(ang_.states().size());
        for (unsigned ill = 0; ill < ang_.states().size(); ill++)
        {
            D[ill].resize(bspline_inner_fine_.Nspline() * bspline_inner_fine_.Nspline());
            for (int i = 0; i < bspline_inner_fine_.Nspline(); i++)
            for (int j = 0; j < bspline_inner_fine_.Nspline(); j++)
            {
                D[ill][i * bspline_inner_fine_.Nspline() + j] =
                    E_ * rad_.S_inner()(i,i) * rad_.S_inner()(j,j)
                    - 0.5 * rad_.D_inner()(i,i) * rad_.S_inner()(j,j)
                    - 0.5 * rad_.S_inner()(i,i) * rad_.D_inner()(j,j)
                    - rad_.Mm2_inner()(i,i) * rad_.S_inner()(j,j)
                    - rad_.S_inner()(i,i) * rad_.Mm2_inner()(j,j)
                    + rad_.Mm1_tr_inner()(i,i) * rad_.S_inner()(j,j)
                    + rad_.S_inner()(i,i) * rad_.Mm1_tr_inner()(j,j);
                
                for (int lambda = 0; lambda <= rad_.maxlambda(); lambda++)
                {
                    D[ill][i * bspline_inner_fine_.Nspline() + j]
                        += ang_.f(lambda,ill,ill) * rad_.computeR(lambda, i, j, i, j);
                }
            }
        }
    }
}

void GMGPreconditioner::precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const
{
    // coarse arrays
    BlockArray<Complex> rn (r.size()), zn (z.size());
    
    // shorthands
    std::size_t Nspline_inner_fine   = bspline_inner_fine_.Nspline();
    std::size_t Nspline_inner_coarse = bspline_inner_coarse_.Nspline();
    std::size_t Nspline_full_fine    = bspline_full_fine_.Nspline();
    std::size_t Nspline_full_coarse  = bspline_full_coarse_.Nspline();
    std::size_t Nspline_outer_fine   = Nspline_full_fine   - Nspline_inner_fine;
    std::size_t Nspline_outer_coarse = Nspline_full_coarse - Nspline_inner_coarse;
    
    // restriction: down-sample the residual
    if (level_ > 0)
    {
//         std::cout << "before restriction r[0].norm() = " << r[0].norm() << std::endl;
//         std::cout << "r[0].size() = " << r[0].size() << ", bspline_inner_fine_.Nspline()^2 = " << bspline_inner_fine_.Nspline()*bspline_inner_fine_.Nspline() << std::endl;
//         std::ofstream ofs (format("r-%.4x.vtk", bspline_inner_fine_.hash()));
//         rArray grid = linspace(0., 100., 1001);
//         writeVTK_points
//         (
//             ofs,
//             bspline_inner_fine_.zip(r[0], grid, grid),
//             grid,
//             grid,
//             rArray { 0. }
//         );
//         ofs.close();
        
        for (unsigned ill = 0; ill < r.size(); ill++)
        {
            // get number of channels
            int nc = (Nspline_inner_fine == Nspline_full_fine ? 0 : (r[ill].size() - Nspline_inner_fine * Nspline_inner_fine) / Nspline_outer_fine);
            
            // allocate the memory
            rn[ill].resize(Nspline_inner_coarse * Nspline_inner_coarse + nc * Nspline_outer_coarse);
            zn[ill].resize(Nspline_inner_coarse * Nspline_inner_coarse + nc * Nspline_outer_coarse);
            
            // down-sample the inner region
            cArrayView
            (
                rn[ill],
                0,
                Nspline_inner_coarse * Nspline_inner_coarse
            ) = kron_dot
            (
                restrictor_inner_,
                restrictor_inner_,
                cArrayView(r[ill], 0, Nspline_inner_fine * Nspline_inner_fine)
            );
            
            // down-sample the outer region
            for (int ic = 0; ic < nc; ic++)
            {
                cArrayView
                (
                    rn[ill],
                    Nspline_inner_coarse * Nspline_inner_coarse + ic * Nspline_outer_coarse,
                    Nspline_outer_coarse
                ) = restrictor_outer_ * cArrayView
                (
                    r[ill],
                    Nspline_inner_fine * Nspline_inner_fine + ic * Nspline_outer_fine,
                    Nspline_outer_fine
                );
            }
        }
        
//         std::cout << "rn[0].norm() = " << rn[0].norm() << std::endl;
//         std::cout << "rn[0].size() = " << rn[0].size() << ", bspline_inner_coarse_.Nspline()^2 = " << bspline_inner_coarse_.Nspline()*bspline_inner_coarse_.Nspline() << std::endl;
//         ofs.open(format("rn-%.4x.vtk", bspline_inner_coarse_.hash()));
//         writeVTK_points
//         (
//             ofs,
//             bspline_inner_coarse_.zip(rn[0], grid, grid),
//             grid,
//             grid,
//             rArray { 0. }
//         );
        
        /// VERIFY
//         kron_dot(S_fine, S_fine, r[0]);
//         kron_dot(S_coar, S_coar, rn[0]);
        
//         ofs.close();
//         std::exit(0);
    }
    
    // solution: precondition by sub-grid
    if (level_ > 0)
    {
        std::cout << "Subgrid" << std::endl;
        subgrid_->precondition(rn, zn);
    }
    else
    {
        std::cout << "r[0].norm() = " << r[0].norm() << std::endl;
        std::ofstream ofs (format("r-%.4x.vtk", bspline_inner_fine_.hash()));
        rArray grid = linspace(0., 100., 1001);
        writeVTK_points
        (
            ofs,
            bspline_inner_fine_.zip(r[0], grid, grid),
            grid,
            grid,
            rArray { 0. }
        );
        ofs.close();
        
        std::cout << "Inner preconditioner" << std::endl;
        subgrid_->precondition(r, z);
        
        std::cout << "z[0].norm() = " << z[0].norm() << std::endl;
        ofs.open(format("z-%.4x.vtk", bspline_inner_fine_.hash()));
        writeVTK_points
        (
            ofs,
            bspline_inner_fine_.zip(z[0], grid, grid),
            grid,
            grid,
            rArray { 0. }
        );
        ofs.close();
        std::exit(1);
    }
    
    // prolongation: interpolate the solution
    if (level_ > 0)
    {
        std::cout << "zn[0].norm() = " << zn[0].norm() << std::endl;
        std::ofstream ofs (format("zn-%.4x.vtk", bspline_inner_coarse_.hash()));
        rArray grid = linspace(0., 100., 1001);
        writeVTK_points
        (
            ofs,
            bspline_inner_coarse_.zip(zn[0], grid, grid),
            grid,
            grid,
            rArray { 0. }
        );
        ofs.close();
        
        for (unsigned ill = 0; ill < r.size(); ill++)
        {
            // get number of channels
            int nc = (Nspline_inner_fine == Nspline_full_fine ? 0 : (r[ill].size() - Nspline_inner_fine * Nspline_inner_fine) / Nspline_outer_fine);
            
            // interpolate the inner region
            cArrayView
            (
                z[ill],
                0,
                Nspline_inner_fine * Nspline_inner_fine
            ) = kron_dot
            (
                prolongator_inner_,
                prolongator_inner_,
                cArrayView(zn[ill], 0, Nspline_inner_coarse * Nspline_inner_coarse)
            );
            
            // interpolate the outer region
            for (int ic = 0; ic < nc; ic++)
            {
                cArrayView
                (
                    z[ill],
                    Nspline_inner_fine * Nspline_inner_fine + ic * Nspline_outer_fine,
                    Nspline_outer_fine
                ) = prolongator_outer_ * cArrayView
                (
                    zn[ill],
                    Nspline_inner_coarse * Nspline_inner_coarse + ic * Nspline_outer_coarse,
                    Nspline_outer_coarse
                );
            }
        }
        
        std::cout << "z[0].norm() = " << z[0].norm() << std::endl;
        ofs.open(format("z-%.4x.vtk", bspline_inner_fine_.hash()));
        writeVTK_points
        (
            ofs,
            bspline_inner_fine_.zip(z[0], grid, grid),
            grid,
            grid,
            rArray { 0. }
        );
        ofs.close();
        std::exit(1);
    }
    
    std::cout << "z[0].norm() = " << z[0].norm() << std::endl;
    std::ofstream ofs (format("z-%.4x.vtk", bspline_inner_fine_.hash()));
    rArray grid = linspace(0., 100., 1001);
    writeVTK_points
    (
        ofs,
        bspline_inner_fine_.zip(z[0], grid, grid),
        grid,
        grid,
        rArray { 0. }
    );
    ofs.close();
    
    // correct high-frequency error using Gauss-Seidel iterations
    if (level_ > 0)
    {
        int nSmoothCycles = 0;
//     std::size_t I = 0, J = 0;
        Complex elem;
        for (int cycle = 0; cycle < nSmoothCycles; cycle++)
        {
            BlockArray<Complex> w = z;
            for (unsigned i = 0; i < zn.size(); i++)
            {
                std::cout << "+@" << cycle << ": z[" << i << "].norm() = " << w[i].norm() << ", D[" << i << "].norm() = " << D[i].norm() << std::endl;
            }
            this->multiply(w, z, MatrixSelection::StrictUpper);
            this->multiply(z, w, MatrixSelection::StrictLower);
            for (unsigned ill = 0; ill < ang_.states().size(); ill++)
            {
                z[ill] -= w[ill];
                z[ill] /= D[ill];
                std::cout << "-@" << cycle << ": z[" << ill << "].norm() = " << z[ill].norm() << std::endl;
            }
            
            /*
            // forward : multiply by the strict upper triangle
            for (int ill = 0; ill < (int)ang_.states().size(); ill++)
            for (std::size_t i = 0; i < Nspline_inner_fine; i++)
            for (std::size_t j = 0; j < Nspline_inner_fine; j++)
            {
                elem = 0;
                J = 0;
                
                for (int illp = 0; illp < (int)ang_.states().size(); illp++)
                for (std::size_t k = 0; k < Nspline_inner_fine; k++)
                for (std::size_t l = 0; l < Nspline_inner_fine; l++)
                {
                    if (J > I)
                    {
                        elem += calc_matrix_elem(ill, illp, i, j, k, l) * z[illp][k * Nspline_inner_fine + l];
                    }
                    J++;
                }
                
                z[ill][i * Nspline_inner_fine + j] = elem;
                I++;
            }
            
            // backward : solve the lower triangle
            for (int ill = 0; ill < (int)ang_.states().size(); ill++)
            for (std::size_t i = 0; i < Nspline_inner_fine; i++)
            for (std::size_t j = 0; j < Nspline_inner_fine; j++)
            {
                elem = 0;
                J = 0;
                
                for (int illp = 0; illp < (int)ang_.states().size(); illp++)
                for (std::size_t k = 0; k < Nspline_inner_fine; k++)
                for (std::size_t l = 0; l < Nspline_inner_fine; l++)
                {
                    if (J < I)
                    {
                        elem += calc_matrix_elem(ill, illp, i, j, k, l) * z[illp][k * Nspline_inner_fine + l];
                    }
                    J++;
                }
                
                z[ill][i * Nspline_inner_fine + j] -= elem;
                z[ill][i * Nspline_inner_fine + j] /= elem += calc_matrix_elem(ill, ill, i, j, i, j);
                
                I++;
            }*/
        }
    }
}

// --------------------------------------------------------------------------------- //

addClassToParentRunTimeSelectionTable(PreconditionerBase, GMGPreconditioner)

// --------------------------------------------------------------------------------- //
