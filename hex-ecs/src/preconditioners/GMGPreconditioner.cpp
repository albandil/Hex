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

#include "GMGPreconditioner.h"

const std::string GMGPreconditioner::prec_name = "GMG";
const std::string GMGPreconditioner::prec_description = "Geometric multigrid.";

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
    Parallel const & par,
    InputFile const & inp,
    AngularBasis const & ll,
    Bspline const & bspline_inner,
    Bspline const & bspline_full,
    CommandLine const & cmd
) : GMGPreconditioner(par, inp, ll, bspline_inner, bspline_full, cmd, cmd.multigrid_depth)
{
    // nothing
}

GMGPreconditioner::GMGPreconditioner
(
    Parallel const & par,
    InputFile const & inp,
    AngularBasis const & ll,
    Bspline const & bspline_inner,
    Bspline const & bspline_full,
    CommandLine const & cmd,
    int level
) : NoPreconditioner(par, inp, ll, bspline_inner, bspline_full, cmd),
    level_(level),
    bspline_inner_fine_(bspline_inner),
    bspline_inner_coarse_
    (
        bspline_inner.order(),
        dither(bspline_inner.rknots()),
        bspline_inner.ECStheta(),
        dither(bspline_inner.cknots())
    ),
    bspline_full_fine_(bspline_full),
    bspline_full_coarse_
    (
        bspline_full.order(),
        dither(bspline_full.rknots()),
        bspline_full.ECStheta(),
        dither(bspline_full.cknots())
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
            bspline_inner_coarse_,
            bspline_full_coarse_,
            cmd
        );
    }
}

void GMGPreconditioner::setup ()
{
    // setup subgrid
    subgrid_->setup();
    
    if (level_ > 0)
    {
        std::cout << "Setting up GMG preconditioner level " << level_ << " (" << bspline_full_coarse_.Nspline() << " B-splines)" << std::endl << std::endl;
        
        // update parent
        NoPreconditioner::setup();
        
        // create integrator
        GaussLegendre g;
        
        // calculate overlap matrix between the fine and coarse bases
        ColMatrix<Complex> Sfc (bspline_inner_fine_.Nspline(), bspline_inner_coarse_.Nspline());
        for (int i = 0; i < bspline_inner_fine_.Nspline(); i++)
        for (int j = 0; j < bspline_inner_coarse_.Nspline(); j++)
            Sfc(i,j) = rad_.computeS12(g, bspline_inner_fine_, bspline_inner_coarse_, i, j);
        ColMatrix<Complex> Scf (Sfc);
        
        // factorize inner overlap matrix for both the fine and coarse basis
        SymBandMatrix<Complex> S_fine = this->rad_.S_inner();
        SymBandMatrix<Complex> S_coarse = dynamic_cast<GMGPreconditioner*>(subgrid_)->rad().S_inner();
        std::cout << "S_coarse.data().size() = " << S_coarse.data().size() << std::endl;
        std::cout << "S_coarse.data().norm() = " << S_coarse.data().norm() << std::endl;
        
        CooMatrix<LU_int_t,Complex> S_fine_coo = S_fine.tocoo<LU_int_t>();
        CooMatrix<LU_int_t,Complex> S_coarse_coo = S_coarse.tocoo<LU_int_t>();
        std::cout << "S_coarse_coo.v().size() = " << S_coarse_coo.v().size() << std::endl;
        
        CsrMatrix<LU_int_t,Complex> S_fine_csr = S_fine_coo.tocsr();
        CsrMatrix<LU_int_t,Complex> S_coarse_csr = S_coarse_coo.tocsr();
        std::cout << "S_coarse_csr.x().size() = " << S_coarse_csr.x().size() << std::endl;
        
        std::cout << "a\n";
        std::shared_ptr<LUft<LU_int_t,Complex>> lu_S_inner_fine = S_fine_csr.factorize();
        std::cout << "b\n";
        std::shared_ptr<LUft<LU_int_t,Complex>> lu_S_inner_coarse = S_coarse_csr.factorize();
        std::cout << "c\n";
        
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
    
    if (level_ > 0)
    {
        // update parent
        NoPreconditioner::update(E);
        
        // calculate full matrix diagonal for use in Gauss-Seidel iterations
        D.resize(ang_.states().size());
        for (unsigned ill = 0; ill < ang_.states().size(); ill++)
        {
            D[ill].resize(bspline_inner_fine_.Nspline() * bspline_inner_fine_.Nspline());
            for (int i = 0; i < bspline_inner_fine_.Nspline(); i++)
            for (int j = 0; j < bspline_inner_fine_.Nspline(); j++)
            {
                D[ill][i * bspline_inner_fine_.Nspline() + j] =
                    E_ * rad_.S_inner()(i,j) * rad_.S_inner()(i,j)
                    - 0.5 * rad_.D_inner()(i,j) * rad_.S_inner()(i,j)
                    - 0.5 * rad_.S_inner()(i,j) * rad_.D_inner()(i,j)
                    - rad_.Mm2_inner()(i,j) * rad_.S_inner()(i,j)
                    - rad_.S_inner()(i,j) * rad_.Mm2_inner()(i,j)
                    + rad_.Mm1_tr_inner()(i,j) * rad_.S_inner()(i,j)
                    + rad_.S_inner()(i,j) * rad_.Mm1_tr_inner()(i,j);
                
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
        for (unsigned ill = 0; ill < r.size(); ill++)
        {
            // get number of channels
            int nc = (r[ill].size() - Nspline_inner_fine * Nspline_inner_fine) / Nspline_outer_fine;
            
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
    }
    
    // solution: precondition by sub-grid
    subgrid_->precondition(rn, zn);
    
    // prolongation: interpolate the solution
    if (level_ > 0)
    {
        for (unsigned ill = 0; ill < r.size(); ill++)
        {
            // get number of channels
            int nc = (r[ill].size() - Nspline_inner_fine * Nspline_inner_fine) / Nspline_outer_fine;
            
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
    }
    
    // correct high-frequency error using Gauss-Seidel iterations
    int nSmoothCycles = 5;
//     std::size_t I = 0, J = 0;
    Complex elem;
    for (int cycle = 0; cycle < nSmoothCycles; cycle++)
    {
        BlockArray<Complex> w = z;
        this->multiply(w, z, MatrixSelection::StrictUpper);
        this->multiply(z, w, MatrixSelection::StrictLower);
        for (unsigned ill = 0; ill < ang_.states().size(); ill++)
        {
            z[ill] -= w[ill];
            z[ill] /= D[ill];
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
