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

#include <iostream>

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"
#include "hex-csrmatrix.h"
#include "hex-densematrix.h"
#include "hex-misc.h"
#include "hex-openmp.h"

// --------------------------------------------------------------------------------- //

#include "gauss.h"
#include "parallel.h"
#include "radial.h"

// --------------------------------------------------------------------------------- //

#include "NoPreconditioner.h"
#include "CGPreconditioner.h"

// --------------------------------------------------------------------------------- //

NoPreconditioner::NoPreconditioner ()
  : PreconditionerBase(),
    E_(0), cmd_(nullptr), par_(nullptr), inp_(nullptr), ang_(nullptr),
    rad_inner_(nullptr), rad_full_(nullptr), rad_panel_(nullptr)
{
    // nothing to do
}

NoPreconditioner::NoPreconditioner
(
    CommandLine  const & cmd,
    InputFile    const & inp,
    Parallel     const & par,
    AngularBasis const & ang,
    Bspline const & bspline_inner,
    Bspline const & bspline_full,
    Bspline const & bspline_panel_x,
    Bspline const & bspline_panel_y
) : PreconditionerBase(),
    E_(0), cmd_(&cmd), par_(&par), inp_(&inp), ang_(&ang),
    A_blocks_ (ang.states().size() * ang.states().size()),
    B1_blocks_(ang.states().size() * ang.states().size()),
    B2_blocks_(ang.states().size() * ang.states().size()),
    Cu_blocks_(ang.states().size() * ang.states().size()),
    Cl_blocks_(ang.states().size() * ang.states().size()),
    rad_inner_(new RadialIntegrals(bspline_inner,   bspline_inner,   ang.maxlambda() + 1)),
    rad_full_ (new RadialIntegrals(bspline_full,    bspline_full,    ang.maxlambda() + 1)),
    rad_panel_(new RadialIntegrals(bspline_panel_x, bspline_panel_y, ang.maxlambda() + 1))
{
    // nothing to do
}

NoPreconditioner::~NoPreconditioner ()
{
    if (rad_inner_) delete rad_inner_;
    if (rad_full_)  delete rad_full_;
    if (rad_panel_) delete rad_panel_;
}

std::string NoPreconditioner::description () const
{
    return "\"Preconditioning\" by the identity matrix.";
}

void NoPreconditioner::setup ()
{
    // shorthands
    Bspline const & bspline_inner = rad_inner_->bspline();
    Bspline const & bspline_full  = rad_full_ ->bspline();
    Bspline const & bspline_x     = rad_panel_->bspline_x();
    Bspline const & bspline_y     = rad_panel_->bspline_y();
    
    // calculate all radial integrals in inner region
    rad_inner_->setupOneElectronIntegrals(*par_, *cmd_);
    rad_inner_->setupTwoElectronIntegrals(*par_, *cmd_);
    
    // calculate one-electron integrals in full domain
    rad_full_->setupOneElectronIntegrals(*par_, *cmd_);
    
    // calculate/copy all radial integrals on panel
    rad_panel_->setupOneElectronIntegrals(*par_, *cmd_);
    rad_panel_->subsetTwoElectronIntegrals(*par_, *cmd_, *rad_inner_);
    
    // number of basis B-splines
    std::size_t Nspline_inner = bspline_inner.Nspline();
    std::size_t Nspline_full  = bspline_full .Nspline();
    
    // number of panel B-splines
    std::size_t Nspline_x = bspline_x.Nspline();
    std::size_t Nspline_y = bspline_y.Nspline();
    
    // is this panel the full solution domain?
    bool full_domain = (bspline_x.hash() == bspline_full.hash() and bspline_y.hash() == bspline_y.hash());
    
    // angular momenta of electrons needed by angular blocks
    std::array<std::vector<int>,2> ells;
    
    // assemble all electrons' angular momenta
    for (std::pair<int,int> ll : ang_->states()) ells[0].push_back(ll.first);
    for (std::pair<int,int> ll : ang_->states()) ells[1].push_back(ll.second);
    
    // sort electrons' angular momenta in ascending order
    std::sort(ells[0].begin(), ells[0].end());
    std::sort(ells[1].begin(), ells[1].end());
    
    // remove duplicates
    ells[0].resize(std::unique(ells[0].begin(), ells[0].end()) - ells[0].begin());
    ells[1].resize(std::unique(ells[1].begin(), ells[1].end()) - ells[1].begin());
    
    // one-electron hamiltonian diagonalization and other data
    Hl_[0].resize(ells[0].back() + 1);
    Hl_[1].resize(ells[1].back() + 1);
    
    // channel functions (pseudostates), their overlaps and energies
    Xp_[0].resize(ells[0].back() + 1);  Xp_[0].fill(cArrays());
    Sp_[0].resize(ells[0].back() + 1);  Sp_[0].fill(cArrays());
    Eb_[0].resize(ells[0].back() + 1);  Eb_[0].fill(cArray());
    Xp_[1].resize(ells[1].back() + 1);  Xp_[1].fill(cArrays());
    Sp_[1].resize(ells[1].back() + 1);  Sp_[1].fill(cArrays());
    Eb_[1].resize(ells[1].back() + 1);  Eb_[1].fill(cArray());
    
    // If the preconditioner is called for full domain, diagonalize the overlap matrix and
    // find the eigenvectors of the one-electron hamiltonian. It will be needed to construct
    // matrix of the system.
    
    if (full_domain)
    {
        std::cout << "Setting up the hydrogen eigenstates..." << std::endl << std::endl;
        
        cArray D (Nspline_inner);
        ColMatrix<Complex> CR (Nspline_inner, Nspline_inner);
        ColMatrix<Complex> invCR (Nspline_inner, Nspline_inner);
        ColMatrix<Complex> invsqrtS (Nspline_inner, Nspline_inner);
        
        std::cout << "\t- inner region basis overlap matrix diagonalization" << std::endl;
        Timer timer;
        
        ColMatrix<Complex> S = rad_inner_->S_x().torow().T();
        S.diagonalize(D, nullptr, &CR);
        CR.invert(invCR);
        
        // Now S = CR * (D * CR⁻¹)
        std::cout << "\t\t- time: " << timer.nice_time() << std::endl;
        for (std::size_t i = 0; i < (std::size_t)Nspline_inner * (std::size_t)Nspline_inner; i++)
            invCR.data()[i] *= D[i % Nspline_inner];
        
        // S = S - CR * invCR
        blas::gemm(-1., CR, invCR, 1., S);
        std::cout << "\t\t- residual: " << S.data().norm() << std::endl;
        
        // compute √S⁻¹
        for (std::size_t i = 0; i < (std::size_t)Nspline_inner * (std::size_t)Nspline_inner; i++)
            invCR.data()[i] /= std::pow(D.data()[i % Nspline_inner], 1.5);
        blas::gemm(1., CR, invCR, 0., invsqrtS);
        
        std::cout << std::endl;
        
        //
        // now calculate the eigenstates of the inner one-electron hamiltonian
        // (these will be used in the asymptotic expansion in the outer region)
        //
        
        // diagonalize the one-electron hamiltonians for both the atomic and projectile particle
        for (int i = 0; i < 2; i++)
        {
            // particle charge (first is electron, second is either electron or positron)
            Real Z = (i == 0 ? -1.0 : inp_->Zp);
            
            // for all one-electron angular momenta
            bool written = false;
            for (int l : ells[i])
            {
                // compose name of the file containing the hamiltonian data
                Hl_[i][l].hdflink(format("Hl%+g-%d-%.4x.hdf", Z, l, bspline_inner.hash()).c_str());
                
                // do not calculate if this work is supposed to be done by someone else
                if (cmd_->shared_scratch and not (par_->isMyGroupWork(l) and par_->IamGroupMaster()))
                    continue;
                
                // check if the file already exists; skip calculation in that case
                if (Hl_[i][l].hdfcheck())
                    continue;
                
                written = true;
                std::cout << "\t- inner region one-electron Hamiltonian matrix diagonalization (Z = " << Z << ", l = " << l << ")" << std::endl;
                timer.reset();
                
                // compose the symmetrical one-electron hamiltonian
                ColMatrix<Complex> tHl;
                if (i == 0)
                    tHl = (Complex(0.5) * rad_inner_->D_x() + Complex(Z) * rad_inner_->Mm1_tr_x() + Complex(0.5*l*(l+1)) * rad_inner_->Mm2_x()).torow().T();
                else
                    tHl = (Complex(0.5) * rad_inner_->D_y() + Complex(Z) * rad_inner_->Mm1_tr_y() + Complex(0.5*l*(l+1)) * rad_inner_->Mm2_y()).torow().T();
                
                // symmetrically transform by inverse square root of the overlap matrix, tHl <- invsqrtS * tHl * invsqrtS
                blas::gemm(1., invsqrtS, tHl, 0., S);
                blas::gemm(1., S, invsqrtS, 0., tHl);
                
                // diagonalize the transformed matrix
                tHl.diagonalize(D, nullptr, &CR);
                CR.invert(invCR);
                
                // store the KPA preconditioner data
                Hl_[i][l].Dl = D;
                Hl_[i][l].invsqrtS_Cl = RowMatrix<Complex>(Nspline_inner, Nspline_inner);
                Hl_[i][l].invCl_invsqrtS = RowMatrix<Complex>(Nspline_inner, Nspline_inner);
                blas::gemm(1., invsqrtS, CR, 0., Hl_[i][l].invsqrtS_Cl);
                blas::gemm(1., invCR, invsqrtS, 0., Hl_[i][l].invCl_invsqrtS);
                
                // Now Hl = ClR * D * ClR⁻¹
                std::cout << "\t\t- time: " << timer.nice_time() << std::endl;
                for (std::size_t i = 0; i < (std::size_t)Nspline_inner * (std::size_t)Nspline_inner; i++)
                    invCR.data()[i] *= D[i % Nspline_inner];
                
                // Hl <- Hl - CR * invCR
                blas::gemm(-1., CR, invCR, 1., tHl);
                std::cout << "\t\t- residual: " << tHl.data().norm() << std::endl;
                
                // copy the eigenvectors as columns
                // - already normalized by xGEEV to "Euclidean norm equal to 1 and largest component real"
                Hl_[i][l].Cl = ColMatrix<Complex>(Hl_[i][l].invsqrtS_Cl);
                
                // write to disk and abandon for now
                Hl_[i][l].hdfsave();
                Hl_[i][l].drop();
            }
            if (written)
            {
                std::cout << std::endl;
            }
        }
        
        // wait for all processes so that all Hlxxx files get written
        par_->wait();
    }
    
    // load the requested one-electron eigenstates from disk files
    for (int i = 0; i < 2; i++)
    {
        // particle charge (first is electron, second is either electron or positron)
        Real Z = (i == 0 ? -1.0 : inp_->Zp);
        
        // for all one-electron angular momenta needed by this particle
        bool written = false;
        for (int l : ells[i])
        {
            // compose name of the file containing the hamiltonian data
            Hl_[i][l].hdflink(format("Hl%+g-%d-%.4x.hdf", Z, l, bspline_inner.hash()).c_str());
            
            // load the factorization file
            if (not Hl_[i][l].hdfload())
                HexException("Failed to load one-electron diagonalization file %s for Z = %g, l = %d.", Hl_[i][l].filename.c_str(), Z, l);
            
            written = true;
            std::cout << "\t- one-electron Hamiltonian data loaded (Z = " << Z << ", l = " << l << ")" << std::endl;
            
            // get sorted energies (ascending real parts)
            std::vector<int> indices (Nspline_inner);
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(), [=](int a, int b){ return Hl_[i][l].Dl[a].real() < Hl_[i][l].Dl[b].real(); });
            
            // get maximal element that has accurate bound state energy
            int max_nr = -1;
            for (unsigned nr = 0; nr < Nspline_inner; nr++)
            {
                Real E = Hl_[i][l].Dl[indices[nr]].real();
                Real E0 = -0.5 / ((nr + l + 1) * (nr + l + 1));
                
                if (E < 0 and std::abs(E0 - E) < 1e-3 * std::abs(E0))
                    max_nr = nr;
                else
                    break;
            }
            if (max_nr >= 0)
            {
                std::cout << "\t\t- bound states with energy within 0.1 % from exact value: " << l + 1 << " <= n <= " << max_nr + l + 1 << std::endl;
            }
            
            // get all valid asymptotic states
            for (unsigned nr = 0; nr < Nspline_inner; nr++)
            {
                // bound energy of the state (a.u. -> Ry)
                Complex Eb = 2.0_r * Hl_[i][l].Dl[indices[nr]];
                
                // maximal bound state energy allowed in the asymptotic expantion (or at least the bound energy allowed by impact channels)
                Real Em = std::max(inp_->channel_max_E, inp_->max_Etot);
                
                // add all requested channels
                if (Eb.real() <= Em)
                {
                    // bound energy of the state (a.u. -> Ry)
                    Complex Eb = 2.0_r * Hl_[i][l].Dl[indices[nr]];
                    
                    // maximal bound state energy allowed in the asymptotic expantion (or at least the bound energy allowed by impact channels)
                    Real Em = std::max(inp_->channel_max_E, inp_->max_Etot);
                    
                    // add all requested channels
                    if (Eb.real() <= Em)
                    {
                        Xp_[i][l].push_back(Hl_[i][l].Cl.col(indices[nr]));
                        Sp_[i][l].push_back(rad_inner_->S().dot(Xp_[i][l][nr]));
                        Eb_[i][l].push_back(Eb);
                        
                        // Adjust the overall sign of the eigenvector so that the result is compatible with the
                        // sign convention of GSL's function gsl_sf_hydrogenicR (used in previous versions of hex-ecs).
                        // That is, the radial function should increase from origin to positive values, then turn back
                        // and (potentially) dive through zero.
                        
                        if (Xp_[i][l].back().front().real() < 0.0_r)
                        {
                            Xp_[i][l].back() = -Xp_[i][l].back();
                            Sp_[i][l].back() = -Sp_[i][l].back();
                        }
                    }
                    else
                    {
                        // there will be no more states (the energy rises with 'nr')
                        break;
                    }
                }
                else
                {
                    // there will be no more states (the energy rises with 'nr')
                    break;
                }
            }
            if (Xp_[i][l].size() >= 1)
            {
                std::cout << "\t\t- eigenstates used in asymptotic (outer) domain: " << l + 1 << " <= n <= " << l + Xp_[i][l].size() << std::endl;
            }
            
            // unload the factorization file
            Hl_[i][l].drop();
        }
        if (written)
        {
            std::cout << std::endl;
        }
    }
}

std::pair<int,int> NoPreconditioner::bstates (Real E, int l1, int l2) const
{
    std::pair<int,int> nstates = { 0, 0 };
    
    if (l1 < (int)Eb_[0].size())
    {
        for (Complex Eb : Eb_[0][l1])
        {
            if (Eb.real() <= E)
            {
                nstates.first++;
            }
        }
    }
    
    if (l2 < (int)Eb_[1].size())
    {
        for (Complex Eb : Eb_[1][l2])
        {
            if (Eb.real() <= E)
            {
                nstates.second++;
            }
        }
    }
    
    return nstates;
}

BlockSymBandMatrix<Complex> NoPreconditioner::calc_A_block (int ill, int illp, bool twoel) const
{
    // number of panel B-splines common with the inner basis
    int Nspline_x_inner = rad_panel_->bspline_x().Nspline();
    int Nspline_y_inner = rad_panel_->bspline_y().Nspline();
    
    // angular momenta
    int l1 = ang_->states()[ill].first;
    int l2 = ang_->states()[ill].second;
    
    BlockSymBandMatrix<Complex> A
    (
        Nspline_x_inner, inp_->order + 1, // block structure
        Nspline_y_inner, inp_->order + 1, // nested structure
        not cmd_->outofcore,              // keep in memory?
        format                            // scratch disk file name for this block and panels
        (
            "blk-A-%d-%d-%04x-%04x.ooc",
            ill, illp,
            rad_panel_->bspline_x().hash(),
            rad_panel_->bspline_y().hash()
        )
    );
    
    // for all sub-blocks
    # pragma omp parallel for if (!cmd_->outofcore)
    for (int i = 0; i < Nspline_x_inner; i++)
    for (int d = 0; d <= inp_->order; d++)
    if (i + d < Nspline_x_inner)
    {
        int k = i + d;
        
        SymBandMatrix<Complex> subblock (Nspline_y_inner, inp_->order + 1);
        
        // one-electron part
        if (ill == illp)
        {
            subblock.populate
            (
                [&](int j, int l)
                {
                    return E_ * rad_panel_->S_x()(i,k) * rad_panel_->S_y()(j,l)
                            - (0.5_z * rad_panel_->D_x()(i,k) - rad_panel_->Mm1_tr_x()(i,k)) * rad_panel_->S_y()(j,l)
                            - 0.5_r * l1 * (l1 + 1) * rad_panel_->Mm2_x()(i,k) * rad_panel_->S_y()(j,l)
                            - rad_panel_->S_x()(i,k) * (0.5_z * rad_panel_->D_y()(j,l) + inp_->Zp * rad_panel_->Mm1_tr_y()(j,l))
                            - 0.5_r * l2 * (l2 + 1) * rad_panel_->S_x()(i,k) * rad_panel_->Mm2_y()(j,l);
                }
            );
        }
        
        // two-electron part
        if(twoel)
        for (int lambda = 0; lambda <= rad_panel_->maxlambda(); lambda++) if (ang_->f(ill,illp,lambda) != 0)
        {
            // calculate two-electron term
            if (not cmd_->lightweight_radial_cache)
            {
                // use precomputed block from scratch file or from memory
                subblock.data() += inp_->Zp * ang_->f(ill,illp,lambda) * rad_panel_->R_tr_dia(lambda).getBlock(i * (inp_->order + 1) + d).slice(0, Nspline_y_inner * (inp_->order + 1));
            }
            else
            {
                // compute the data anew
                subblock.data() += inp_->Zp * ang_->f(ill,illp,lambda) * rad_full_->calc_R_tr_dia_block(lambda, i, k).data().slice(0, Nspline_y_inner * (inp_->order + 1));
            }
        }
        
        // save block
        A.setBlock(i * (inp_->order + 1) + d, subblock.data());
    }
    
    if (cmd_->outofcore)
        A.drop();
    
    return A;
}

void NoPreconditioner::update (Real E)
{
    // shorthands
    int order = inp_->order;
    int Nang = ang_->states().size();
    
    // global basis
    int Nspline_full = rad_full_->bspline().Nspline();
    int Nspline_inner = rad_inner_->bspline().Nspline();
    int Nspline_outer = Nspline_full - Nspline_inner;
    
    // panel x basis
    int Nspline_x = rad_panel_->bspline_x().Nspline();
    int Nspline_x_inner = rad_panel_->bspline_x().knot(rad_inner_->bspline_x().R2());
    int Nspline_x_outer = Nspline_x - Nspline_x_inner;
    
    // panel y basis
    int Nspline_y = rad_panel_->bspline_y().Nspline();
    int Nspline_y_inner = rad_panel_->bspline_y().knot(rad_inner_->bspline_y().R2());
    int Nspline_y_outer = Nspline_y - Nspline_y_inner;
    
    // size of the A-matrix
    std::size_t A_size = std::size_t(Nspline_x_inner) * std::size_t(Nspline_y_inner);
    
    // update energy
    E_ = E;
    
    std::cout << "\tUpdate the common preconditioner base" << std::endl;
    
    std::cout << "\tPrecompute matrix blocks ... " << std::flush;
    Timer t;
    
    // outer one-electron overlap matrix
    SymBandMatrix<Complex> S_outer (Nspline_outer, order + 1);
    S_outer.populate([&](int m, int n) { return rad_full_->S()(Nspline_inner + m, Nspline_inner + n); });
    
    // outer one-electron derivative matrix
    SymBandMatrix<Complex> D_outer (Nspline_outer, order + 1);
    D_outer.populate([&](int m, int n) { return rad_full_->D()(Nspline_inner + m, Nspline_inner + n); });
    
    // outer one-electron centrifugal moment matrix
    SymBandMatrix<Complex> Mm2_outer (Nspline_outer, order + 1);
    Mm2_outer.populate([&](int m, int n) { return rad_full_->Mm2()(Nspline_inner + m, Nspline_inner + n); });
    
    // inner one-electron multipole moment matrices
    std::vector<SymBandMatrix<Complex>> Mtr_L_inner;
    for (int lambda = 0; lambda <= rad_full_->maxlambda(); lambda++)
    {
        Mtr_L_inner.push_back(SymBandMatrix<Complex>(Nspline_inner, order + 1));
        Mtr_L_inner.back().populate
        (
            [ & ] (int m, int n)
            {
                return rad_inner_->Mtr_L(lambda)(m, n)
                     * special::pow_int(rad_full_->bspline().t(std::min(m,n) + order + 1).real(), lambda);
            }
        );
    }
    
    // outer one-electron multipole moment matrices
    std::vector<SymBandMatrix<Complex>> Mtr_mLm1_outer;
    for (int lambda = 0; lambda <= rad_full_->maxlambda(); lambda++)
    {
        Mtr_mLm1_outer.push_back(SymBandMatrix<Complex>(Nspline_outer, order + 1));
        Mtr_mLm1_outer.back().populate
        (
            [ & ] (int m, int n)
            {
                return rad_full_->Mtr_mLm1(lambda)(Nspline_inner + m, Nspline_inner + n)
                     * special::pow_int(rad_full_->bspline().t(Nspline_inner + std::min(m,n) + order + 1).real(), -lambda-1);
            }
        );
    }
    
    // calculate number of asymptotic channels (i.e. bound states of the other particle)
    Nchan_.resize(Nang);
    for (int ill = 0; ill < Nang; ill++)
    {
        int l1 = ang_->states()[ill].first;
        int l2 = ang_->states()[ill].second;
        
        std::pair<int,int> Nbound = bstates(std::max(2*E, inp_->channel_max_E), l1, l2);
        
        Nchan_[ill].first  = Nbound.second;
        Nchan_[ill].second = Nbound.first;
    }
    
    // setup blocks
    for (int ill = 0; ill < Nang; ill++) if (par_->isMyGroupWork(ill))
    for (int illp = 0; illp < Nang; illp++)
    {
        // angular momenta
        int l1 = ang_->states()[ill].first;
        int l2 = ang_->states()[ill].second;
        int l1p = ang_->states()[illp].first;
        int l2p = ang_->states()[illp].second;
        
        // get number of asymptotic bound channels
        int Nchan1 = Nchan_[ill].first;     // # r1 -> inf, l2 bound
        int Nchan2 = Nchan_[ill].second;    // # r2 -> inf, l1 bound
        int Nchan1p = Nchan_[illp].first;   // # r1 -> inf, l2p bound
        int Nchan2p = Nchan_[illp].second;  // # r2 -> inf, l1p bound
        
        // initialize diagonal block of the inner problem
        // - do not precompute off-diagonal blocks in lightweight mode
        if (not cmd_->lightweight_full or ill == illp)
            A_blocks_[ill * Nang + illp] = calc_A_block(ill, illp);
        
        // create inner-outer coupling blocks
        Cu_blocks_[ill * Nang + illp] = CooMatrix<LU_int_t,Complex>
        (
            A_size + Nchan1  * Nspline_x_outer + Nchan2  * Nspline_y_outer,
            A_size + Nchan1p * Nspline_x_outer + Nchan2p * Nspline_y_outer
        );
        Cl_blocks_[ill * Nang + illp] = CooMatrix<LU_int_t,Complex>
        (
            A_size + Nchan1  * Nspline_x_outer + Nchan2  * Nspline_y_outer,
            A_size + Nchan1p * Nspline_x_outer + Nchan2p * Nspline_y_outer
        );
        
        // setup stretched inner-outer problem
        if (not inp_->inner_only)
        {
            // outer problem matrix : r2 -> inf, r1 bound
            B2_blocks_[ill * Nang + illp].resize(Nchan2 * Nchan2p);
            for (int m = 0; m < Nchan2; m++)
            for (int n = 0; n < Nchan2p; n++)
            {
                SymBandMatrix<Complex> subblock (Nspline_y_outer, order + 1);
                
                // channel-diagonal contribution
                if (ill == illp and m == n)
                {
                    subblock += (E_ + 1.0_z / (2.0_z * (l1 + m + 1.0_r) * (l1 + m + 1.0_r))) * S_outer
                             - 0.5_z * D_outer
                             - 0.5_z * (l2 * (l2 + 1.0_r)) * Mm2_outer;
                }
                
                // channel-offdiagonal contribution
                for (int lambda = 1; lambda <= rad_full_->maxlambda(); lambda++) if (ang_->f(ill,illp,lambda) != 0.0_r)
                    subblock += inp_->Zp * ang_->f(ill,illp,lambda) * (Xp_[0][l1p][n] | Mtr_L_inner[lambda].dot(Xp_[0][l1][m])) * Mtr_mLm1_outer[lambda];
                
                // use the block
                B2_blocks_[ill * Nang + illp][m * Nchan2p + n].hdflink(format("blk-B2-%d-%d-%d-%d.ooc", ill, illp, m, n));
                B2_blocks_[ill * Nang + illp][m * Nchan2p + n] = std::move(subblock);
                if (cmd_->outofcore)
                {
                    B2_blocks_[ill * Nang + illp][m * Nchan2p + n].hdfsave();
                    B2_blocks_[ill * Nang + illp][m * Nchan2p + n].drop();
                }
            }
            
            // outer problem matrix : r1 -> inf, r2 bound
            B1_blocks_[ill * Nang + illp].resize(Nchan1 * Nchan1p);
            for (int m = 0; m < Nchan1; m++)
            for (int n = 0; n < Nchan1p; n++)
            {
                SymBandMatrix<Complex> subblock (Nspline_y_outer, order + 1);
                
                // channel-diagonal contribution
                if (ill == illp and m == n)
                {
                    subblock += (E_ + 1.0_z / (2.0_z * (l2 + m + 1.0_r) * (l2 + m + 1.0_r))) * S_outer
                             - 0.5_z * D_outer
                             - 0.5_z * (l1 * (l1 + 1.0_r)) * Mm2_outer;
                }
                
                // channel-offdiagonal contribution
                for (int lambda = 1; lambda <= rad_full_->maxlambda(); lambda++) if (ang_->f(ill,illp,lambda) != 0.0_r)
                    subblock += inp_->Zp * ang_->f(ill,illp,lambda) * (Xp_[1][l2p][n] | Mtr_L_inner[lambda].dot(Xp_[1][l2][m])) * Mtr_mLm1_outer[lambda];
                
                // use the block
                B1_blocks_[ill * Nang + illp][m * Nchan1p + n].hdflink(format("blk-B1-%d-%d-%d-%d.ooc", ill, illp, m, n));
                B1_blocks_[ill * Nang + illp][m * Nchan1p + n] = std::move(subblock);
                if (cmd_->outofcore)
                {
                    B1_blocks_[ill * Nang + illp][m * Nchan1p + n].hdfsave();
                    B1_blocks_[ill * Nang + illp][m * Nchan1p + n].drop();
                }
            }
            
            // transition area r2 > r1, upper : psi_kl expressed in terms of F_nl for 'l' out of inner area
            for (int i = 0; i < Nspline_inner; i++)
            for (int j = 0; j < Nspline_inner; j++) // *
            for (int k = std::max(0, i - order); k <= std::min(i + order, Nspline_inner - 1); k++)
            for (int l = Nspline_inner; l <= j + order; l++)
            for (int n = 0; n < Nchan2p; n++)
            {
                std::size_t row = i * Nspline_inner + j;
                std::size_t col = A_size + (Nchan1p + n) * Nspline_outer + (l - Nspline_inner);
                
                Complex elem = 0; // A_ij,kl Xp_nk
                
                if (ill == illp)
                {
                    elem += E_ * rad_full_->S_x()(i,k) * rad_full_->S_y()(j,l)
                         - 0.5_r * rad_full_->D_x()(i,k) * rad_full_->S_y()(j,l)
                         - 0.5_r * rad_full_->S_x()(i,k) * rad_full_->D_y()(j,l)
                         - 0.5_r * (l1 * (l1 + 1.0_r)) * rad_full_->Mm2_x()(i,k) * rad_full_->S_y()(j,l)
                         - 0.5_r * (l2 * (l2 + 1.0_r)) * rad_full_->S_x()(i,k) * rad_full_->Mm2_y()(j,l)
                         + rad_full_->Mm1_tr_x()(i,k) * rad_full_->S_y()(j,l)
                         - inp_->Zp * rad_full_->S_x()(i,k) * rad_full_->Mm1_tr_y()(j,l);
                }
                
                Real r1 = rad_full_->bspline_x().t(std::min(i,k) + order + 1).real();
                Real r2 = rad_full_->bspline_y().t(std::min(j,l) + order + 1).real();
                
                for (int lambda = 0; lambda <= rad_full_->maxlambda(); lambda++) if (ang_->f(ill,illp,lambda) != 0)
                {
                    Real multipole = special::pow_int(r1/r2, lambda) / r2;
                    elem += inp_->Zp * ang_->f(ill,illp,lambda) * multipole * rad_full_->Mtr_L_x(lambda)(i,k) * rad_full_->Mtr_mLm1_y(lambda)(j,l);
                }
                
                Cu_blocks_[ill * Nang + illp].add(row, col, Xp_[0][l1p][n][k] * elem);
            }
            
            // transition area r1 > r2, upper : psi_kl expressed in terms of F_nk for 'k' out of inner area
            for (int i = 0; i < Nspline_inner; i++)
            for (int j = 0; j < Nspline_inner; j++) // *
            for (int k = Nspline_inner; k <= i + order; k++)
            for (int l = std::max(0, j - order); l <= std::min(j + order, Nspline_inner - 1); l++)
            for (int n = 0; n < Nchan1p; n++)
            {
                std::size_t row = i * Nspline_inner + j;
                std::size_t col = A_size + n * Nspline_outer + (k - Nspline_inner);
                
                Complex elem = 0; // A_ij,kl Xp_nl
                
                if (ill == illp)
                {
                    elem += E_ * rad_full_->S_x()(i,k) * rad_full_->S_y()(j,l)
                         - 0.5_r * rad_full_->D_x()(i,k) * rad_full_->S_y()(j,l)
                         - 0.5_r * rad_full_->S_x()(i,k) * rad_full_->D_y()(j,l)
                         - 0.5_r * (l1 * (l1 + 1.0_r)) * rad_full_->Mm2_x()(i,k) * rad_full_->S_y()(j,l)
                         - 0.5_r * (l2 * (l2 + 1.0_r)) * rad_full_->S_x()(i,k) * rad_full_->Mm2_y()(j,l)
                         + rad_full_->Mm1_tr_x()(i,k) * rad_full_->S_y()(j,l)
                         - inp_->Zp * rad_full_->S_x()(i,k) * rad_full_->Mm1_tr_y()(j,l);
                }
                
                Real r1 = rad_full_->bspline_x().t(std::min(i,k) + order + 1).real();
                Real r2 = rad_full_->bspline_y().t(std::min(j,l) + order + 1).real();
                
                for (int lambda = 0; lambda <= rad_full_->maxlambda(); lambda++) if (ang_->f(ill,illp,lambda) != 0)
                {
                    Real multipole = special::pow_int(r2/r1, lambda) / r1;
                    elem += inp_->Zp * ang_->f(ill,illp,lambda) * multipole * rad_full_->Mtr_mLm1_x(lambda)(i,k) * rad_full_->Mtr_L_y(lambda)(j,l);
                }
                
                Cu_blocks_[ill * Nang + illp].add(row, col, Xp_[1][l2p][n][l] * elem);
            }
            
            // transition area r2 > r1, lower : F_nl expressed in terms of psi_kl for 'l' out of outer area
            for (int m = 0; m < Nchan2; m++)
            for (int n = 0; n < Nchan2p; n++)
            for (int j = Nspline_inner; j < Nspline_full; j++) // *
            for (int k = 0; k < Nspline_inner; k++)
            for (int l = j - order; l < Nspline_inner; l++)
            {
                std::size_t row = A_size + (Nchan1 + m) * Nspline_outer + (j - Nspline_inner);
                std::size_t col = k * Nspline_inner + l;
                
                Complex elem = 0; // B_mj,nl Sp_nk
                
                if (ill == illp and m == n)
                {
                    elem += (E_ + 0.5_r / ((n + l1 + 1) * (n + l1 + 1))) * rad_full_->S_y()(j,l)
                         - 0.5_r * rad_full_->D_y()(j,l)
                         - 0.5_r * (l2 * (l2 + 1.0_r)) * rad_full_->Mm2_y()(j,l);
                }
                
                Real r2 = rad_full_->bspline_y().t(std::min(j,l) + order + 1).real();
                
                for (int lambda = 1; lambda <= rad_full_->maxlambda(); lambda++) if (ang_->f(ill,illp,lambda) != 0)
                {
                    Real multipole = special::pow_int(1/r2, lambda + 1);
                    elem += inp_->Zp * ang_->f(ill,illp,lambda) * multipole * rad_full_->Mtr_mLm1_x(lambda)(j,l) * (Xp_[0][l1p][n] | Mtr_L_inner[lambda].dot(Xp_[0][l1][m]));
                }
                
                Cl_blocks_[ill * Nang + illp].add(row, col, Sp_[0][l1p][n][k] * elem);
            }
            
            // transition area r1 > r2, lower : F_nk expressed in terms of psi_kl for 'k' out of outer area
            for (int m = 0; m < Nchan1; m++)
            for (int n = 0; n < Nchan1p; n++)
            for (int i = Nspline_inner; i < Nspline_full; i++) // *
            for (int k = i - order; k < Nspline_inner; k++)
            for (int l = 0; l < Nspline_inner; l++)
            {
                std::size_t row = A_size + m * Nspline_outer + (i - Nspline_inner);
                std::size_t col = k * Nspline_inner + l;
                
                Complex elem = 0; // B_mi,nk Sp_nl
                
                if (ill == illp and m == n)
                {
                    elem += (E_ + 0.5_r / ((n + l2 + 1) * (n + l2 + 1))) * rad_full_->S_x()(i,k)
                         - 0.5_r * rad_full_->D_x()(i,k)
                         - 0.5_r * (l1 * (l1 + 1.0_r)) * rad_full_->Mm2_x()(i,k);
                }
                
                Real r1 = rad_full_->bspline_x().t(std::min(i,k) + order + 1).real();
                
                for (int lambda = 1; lambda <= rad_full_->maxlambda(); lambda++) if (ang_->f(ill,illp,lambda) != 0)
                {
                    Real multipole = special::pow_int(1/r1, lambda + 1);
                    elem += inp_->Zp * ang_->f(ill,illp,lambda) * multipole * rad_full_->Mtr_mLm1_x(lambda)(i,k) * (Xp_[1][l2p][n] | Mtr_L_inner[lambda].dot(Xp_[1][l2][m]));
                }
                
                Cl_blocks_[ill * Nang + illp].add(row, col, Sp_[1][l2p][n][l] * elem);
            }
        }
    }
    
    std::cout << "done after " << t.nice_time() << std::endl;
    par_->wait();
}

void NoPreconditioner::rhs (BlockArray<Complex> & chi, int ie, int instate) const
{
    // shorthands
    int ni = std::get<0>(inp_->instates[instate]);
    int li = std::get<1>(inp_->instates[instate]);
    int mi = std::get<2>(inp_->instates[instate]);
    
    // shorthands
    int order = rad_inner_->bspline().order();
    std::size_t Nspline_inner = rad_inner_->bspline().Nspline();
    std::size_t Nspline_full  = rad_full_->bspline().Nspline();
    std::size_t Nspline_outer = Nspline_full - Nspline_inner;
    
    // impact momentum
    rArray ki = { std::sqrt(inp_->Etot[ie] + 1.0_r/(ni*ni)) };
    
    // get the initial bound pseudo-state B-spline expansion
    cArray Xp = Hl_[0][li].readPseudoState(li, ni - li - 1);
    
    // (anti)symmetrization
    Real Sign = ((ang_->S() + ang_->Pi()) % 2 == 0) ? 1. : -1.;
    
    // inner region multipole matrix
    std::vector<SymBandMatrix<Complex>> Mtr_L_inner;
    for (int lambda = 0; lambda <= rad_full_->maxlambda(); lambda++)
    {
        Mtr_L_inner.push_back(SymBandMatrix<Complex>(Nspline_inner, order + 1));
        Mtr_L_inner.back().populate
        (
            [ & ] (int m, int n)
            {
                return rad_inner_->Mtr_L_x(lambda)(m, n) * special::pow_int(rad_full_->bspline().t(std::min(m,n) + order + 1).real(), lambda);
            }
        );
    }
    
    // for all segments constituting the RHS
    for (unsigned ill = 0; ill < ang_->states().size(); ill++) if (par_->isMyGroupWork(ill))
    {
        int l1 = ang_->states()[ill].first;
        int l2 = ang_->states()[ill].second;
        
        // get number of open channels in the outer region
        int Nchan1 = Nchan_[ill].first;     // r1 -> inf, l2 bound
        int Nchan2 = Nchan_[ill].second;    // r2 -> inf, l1 bound
        
        // setup storage
        cArray chi_block (Nspline_inner * Nspline_inner + (Nchan1 + Nchan2) * Nspline_outer);
        
        // for all allowed angular momenta (by momentum composition) of the projectile
        for (int l = std::abs(li - ang_->L()); l <= li + ang_->L(); l++)
        {
            // skip wrong parity
            if ((ang_->L() + li + l) % 2 != ang_->Pi())
                continue;
            
            // compute energy- and angular momentum-dependent prefactor
            Complex prefactor = std::pow(1.0_i,l)
                              * std::sqrt(4.0_r * special::constant::pi * (2 * l + 1))
                              * (Real)special::ClebschGordan(li,mi, l,0, inp_->L,mi) / ki[0];
            
            // skip non-contributing terms
            if (prefactor == 0.0_r)
                continue;
            
            // calculate angular integrals
            rArray f1 (rad_full_->maxlambda() + 1), f2 (rad_full_->maxlambda() + 1);
            for (int lambda = 0; lambda <= rad_full_->maxlambda(); lambda++)
            {
                f1[lambda] = special::computef(lambda, l1, l2, li, l, inp_->L);
                f2[lambda] = special::computef(lambda, l1, l2, l, li, inp_->L);
                
                // abort if any of the coefficients is non-number (factorial overflow etc.)
                if (not std::isfinite(f1[lambda]))
                    HexException("Invalid result of computef(%d,%d,%d,%d,%d,%d)\n", lambda,l1,l2,li,l,inp_->L);
                if (not std::isfinite(f2[lambda]))
                    HexException("Invalid result of computef(%d,%d,%d,%d,%d,%d)\n", lambda,l1,l2,l,li,inp_->L);
            }
            
            // quadrature degree
            int points = order + li + l + 1;
            
            // prepare quadrature nodes and weights
            rad_full_->gaussleg().precompute_nodes_and_weights(points);
            
            // precompute quadrature nodes and weights
            cArray xs ((rad_full_->bspline().Nreknot() - 1) * points), xws ((rad_full_->bspline().Nreknot() - 1) * points);
            # pragma omp parallel for
            for (int ixknot = 0; ixknot < rad_full_->bspline().Nreknot() - 1; ixknot++)
                rad_full_->gaussleg().scaled_nodes_and_weights(points, rad_full_->bspline().t(ixknot), rad_full_->bspline().t(ixknot + 1), &xs[ixknot * points], &xws[ixknot * points]);
            
            // precompute B-splines
            cArray B_x (rad_full_->bspline().Nspline() * (order + 1) * points);
            # pragma omp parallel for
            for (int ixspline = 0; ixspline < rad_full_->bspline().Nspline(); ixspline++)
            for (int ixknot = ixspline; ixknot <= ixspline + order and ixknot < rad_full_->bspline().Nreknot() - 1; ixknot++)
                rad_full_->bspline().B(ixspline, ixknot, points, &xs[ixknot * points], &B_x[(ixspline * (order + 1) + ixknot - ixspline) * points]);
            
            // precompute radial functions and Riccati-Bessel functions
            rArray Pi_x = realpart(rad_inner_->bspline().zip(Xp, realpart(xs)));
            rArray ji_x (xs.size());
            # pragma omp parallel for
            for (unsigned ix = 0; ix < xs.size(); ix++)
                ji_x[ix] = special::ric_j(l, ki[0] * xs[ix].real());
            
            // precompute integral moments
            cArrays M_L_P (rad_full_->maxlambda() + 1), M_mLm1_P (rad_full_->maxlambda() + 1);
            cArrays M_L_j (rad_full_->maxlambda() + 1), M_mLm1_j (rad_full_->maxlambda() + 1);
            for (int lambda = 0; lambda <= rad_full_->maxlambda(); lambda++)
            {
                M_L_P[lambda].resize(Nspline_full); M_mLm1_P[lambda].resize(Nspline_full);
                M_L_j[lambda].resize(Nspline_full); M_mLm1_j[lambda].resize(Nspline_full);
                
                # pragma omp parallel for
                for (int ispline = 0; ispline < (int)Nspline_full; ispline++)
                {
                    for (int iknot = ispline; iknot < std::min(ispline + order + 1, rad_full_->bspline().Nreknot() - 1); iknot++)
                    if (rad_full_->bspline().t(iknot).real() != rad_full_->bspline().t(iknot + 1).real())
                    {
                        for (int ipoint = 0; ipoint < points; ipoint++)
                        {
                            Real x = xs[iknot * points + ipoint].real();
                            Real t = rad_full_->bspline().t(ispline + order + 1).real();
                            
                            Real L = special::pow_int(x / t, lambda);
                            Real mLm1 = special::pow_int(x / t, -lambda-1);
                            
                            Complex w = xws[iknot * points + ipoint];
                            Complex B = B_x[(ispline * (order + 1) + iknot - ispline) * points + ipoint];
                            
                            Real P = Pi_x[iknot * points + ipoint];
                            Real j = ji_x[iknot * points + ipoint];
                            
                            Real d = 1;//damp(x, 0, distance);
                            
                            M_L_P[lambda][ispline] += w * B * L * P * d;
                            M_L_j[lambda][ispline] += w * B * L * j * d;
                            M_mLm1_P[lambda][ispline] += w * B * mLm1 * P * d;
                            M_mLm1_j[lambda][ispline] += w * B * mLm1 * j * d;
                        }
                    }
                }
            }
            
            // precompute hydrogen multipoles
            cArrays rho_l1 (Nchan2), rho_l2 (Nchan1);
            for (int ichan1 = 0; ichan1 < Nchan1; ichan1++)
            {
                // this is needed only for exchange contribution
                rho_l2[ichan1].resize(rad_full_->maxlambda() + 1);
                for (int lambda = 1; lambda <= rad_full_->maxlambda(); lambda++)
                    rho_l2[ichan1][lambda] = (Xp_[1][l2][ichan1] | Mtr_L_inner[lambda].dot(Xp));
            }
            for (int ichan2 = 0; ichan2 < Nchan2; ichan2++)
            {
                // this is needed only for direct contribution
                rho_l1[ichan2].resize(rad_full_->maxlambda() + 1);
                for (int lambda = 1; lambda <= rad_full_->maxlambda(); lambda++)
                    rho_l1[ichan2][lambda] = (Xp_[0][l1][ichan2] | Mtr_L_inner[lambda].dot(Xp));
            }
            
            // for all B-spline pairs (elements of the right-hand side)
            # pragma omp parallel for schedule(dynamic,Nspline_inner)
            for (std::size_t ispline = 0; ispline < Nspline_inner * Nspline_inner + (Nchan1 + Nchan2) * Nspline_outer; ispline++)
            {
                // contributions to the element of the right-hand side
                Complex contrib_direct = 0, contrib_exchange = 0;
                
                // determine inner/outer region
                if (ispline >= Nspline_inner * Nspline_inner and ispline < Nspline_inner * Nspline_inner + Nchan1 * Nspline_outer)
                {
                    // r1 > Ra, r2 < Ra
                    int ixspline = (ispline - Nspline_inner * Nspline_inner) % Nspline_outer + Nspline_inner;
                    int ichan1 = (ispline - Nspline_inner * Nspline_inner) / Nspline_outer;
                    
                    // calculate the exchange contribution
                    Real x = rad_full_->bspline().t(ixspline + order + 1).real(), multipole = x;
                    for (int lambda = 1; lambda <= rad_full_->maxlambda(); lambda++)
                    {
                        multipole *= x;
                        
                        if (f2[lambda] != 0 and inp_->exchange)
                            contrib_exchange += f2[lambda] * rho_l2[ichan1][lambda] * M_mLm1_j[lambda][ixspline] / multipole;
                    }
                }
                else if (ispline >= Nspline_inner * Nspline_inner + Nchan1 * Nspline_outer)
                {
                    // r1 < Ra, r2 > Ra
                    int iyspline = (ispline - Nspline_inner * Nspline_inner - Nchan1 * Nspline_outer) % Nspline_outer + Nspline_inner;
                    int ichan2 = (ispline - Nspline_inner * Nspline_inner - Nchan1 * Nspline_outer) / Nspline_outer;
                    
                    // calculate the direct contribution
                    Real y = rad_full_->bspline().t(iyspline + order + 1).real(), multipole = y;
                    for (int lambda = 1; lambda <= rad_full_->maxlambda(); lambda++)
                    {
                        multipole *= y;
                        
                        if (f1[lambda] != 0)
                            contrib_direct += f1[lambda] * rho_l1[ichan2][lambda] * M_mLm1_j[lambda][iyspline] / multipole;
                    }
                }
                else /* i.e. ispline < Nspline_inner * Nspline_inner */
                {
                    // r1 < Ra, r2 < Ra
                    int ixspline = ispline / Nspline_inner;
                    int iyspline = ispline % Nspline_inner;
                    
                    // non-overlapping B-splines ?
                    if (std::abs(ixspline - iyspline) > order)
                    {
                        // monopole contribution
                        if (ixspline > iyspline and f1[0] != 0)
                        {
                            contrib_direct   += f1[0] * (M_mLm1_P[0][ixspline] * M_L_j[0][iyspline] / rad_full_->bspline().t(ixspline + order + 1) - M_L_P[0][ixspline] * M_mLm1_j[0][iyspline] / rad_full_->bspline_y().t(iyspline + order + 1));
                            contrib_exchange += 0;
                        }
                        if (ixspline < iyspline and f2[0] != 0 and inp_->exchange)
                        {
                            contrib_direct   += 0;
                            contrib_exchange += f2[0] * (M_L_j[0][ixspline] * M_mLm1_P[0][iyspline] / rad_full_->bspline().t(iyspline + order + 1) - M_mLm1_j[0][ixspline] * M_L_P[0][iyspline] / rad_full_->bspline_x().t(ixspline + order + 1));
                        }
                        
                        // multipole contributions
                        Real x = rad_full_->bspline().t(ixspline + order + 1).real(), y = rad_full_->bspline().t(iyspline + order + 1).real();
                        Real multipole1 = 1 / x, multipole2 = 1 / y, y_over_x = y / x, x_over_y = x / y;
                        for (int lambda = 1; lambda <= rad_full_->maxlambda(); lambda++)
                        {
                            multipole1 *= y_over_x;
                            multipole2 *= x_over_y;
                            
                            if (ixspline > iyspline)
                            {
                                /* always */        contrib_direct   += f1[lambda] * M_mLm1_P[lambda][ixspline] * M_L_j[lambda][iyspline] * multipole1;
                                if (inp_->exchange) contrib_exchange += f2[lambda] * M_mLm1_j[lambda][ixspline] * M_L_P[lambda][iyspline] * multipole1;
                            }
                            if (ixspline < iyspline)
                            {
                                /* always */        contrib_direct   += f1[lambda] * M_L_P[lambda][ixspline] * M_mLm1_j[lambda][iyspline] * multipole2;
                                if (inp_->exchange) contrib_exchange += f2[lambda] * M_L_j[lambda][ixspline] * M_mLm1_P[lambda][iyspline] * multipole2;
                            }
                        }
                    }
                    else
                    {
                        // for all knots
                        for (int ixknot = ixspline; ixknot <= ixspline + order and ixknot < rad_full_->bspline().Nreknot() - 1; ixknot++) if (rad_full_->bspline().t(ixknot).real() != rad_full_->bspline_x().t(ixknot + 1).real())
                        for (int iyknot = iyspline; iyknot <= iyspline + order and iyknot < rad_full_->bspline().Nreknot() - 1; iyknot++) if (rad_full_->bspline().t(iyknot).real() != rad_full_->bspline_y().t(iyknot + 1).real())
                        {
                            // off-diagonal contribution
                            if (ixknot != iyknot)
                            {
                                // for all quadrature points
                                for (int ix = 0; ix < points; ix++)
                                for (int iy = 0; iy < points; iy++)
                                {
                                    // radii
                                    Real rx = xs[ixknot * points + ix].real(), ry = xs[iyknot * points + iy].real(), rmin = std::min(rx,ry), rmax = std::max(rx,ry);
                                    
                                    // evaluated functions
                                    Complex Bx = B_x[(ixspline * (order + 1) + ixknot - ixspline) * points + ix];
                                    Complex By = B_x[(iyspline * (order + 1) + iyknot - iyspline) * points + iy];
                                    Real Pix = Pi_x[ixknot * points + ix];
                                    Real Piy = Pi_x[iyknot * points + iy];
                                    Real jix = ji_x[ixknot * points + ix];
                                    Real jiy = ji_x[iyknot * points + iy];
                                    Complex wx = xws[ixknot * points + ix];
                                    Complex wy = xws[iyknot * points + iy];
                                    
                                    // damp factor
                                    Real dampfactor = 1;//damp(rx, ry, distance);
                                    
                                    // monopole contribution
                                    if (rx > ry and li == l1 and l == l2)                    contrib_direct   += Bx * By * (1.0_r/rx - 1.0_r/ry) * Pix * jiy * dampfactor * wx * wy;
                                    if (ry > rx and li == l2 and l == l1 and inp_->exchange) contrib_exchange += Bx * By * (1.0_r/ry - 1.0_r/rx) * jix * Piy * dampfactor * wx * wy;
                                    
                                    // higher multipoles contribution
                                    Real multipole = 1 / rmax, rmin_over_rmax = rmin / rmax;
                                    for (int lambda = 1; lambda <= rad_full_->maxlambda(); lambda++)
                                    {
                                        multipole *= rmin_over_rmax;
                                        if (f1[lambda] != 0)                    contrib_direct   += f1[lambda] * Bx * By * multipole * Pix * jiy * dampfactor * wx * wy;
                                        if (f2[lambda] != 0 and inp_->exchange) contrib_exchange += f2[lambda] * Bx * By * multipole * jix * Piy * dampfactor * wx * wy;
                                    }
                                }
                            }
                            // diagonal contribution: needs to be integrated more carefully
                            else if (ixknot < rad_full_->bspline_x().Nreknot() - 1)
                            {
                                // for all quadrature points from the triangle x < y
                                for (int ix = 0; ix < points; ix++)
                                {
                                    cArray ys (points), yws (points), B_y (points);
                                    rad_full_->gaussleg().scaled_nodes_and_weights(points, xs[ixknot * points + ix], rad_full_->bspline().t(iyknot + 1), &ys[0], &yws[0]);
                                    rad_full_->bspline().B(iyspline, iyknot, points, &ys[0], &B_y[0]);
                                    
                                    for (int iy = 0; iy < points; iy++)
                                    {
                                        // radii
                                        Real rx = xs[ixknot * points + ix].real(), ry = ys[iy].real(), rmin = std::min(rx,ry), rmax = std::max(rx,ry);
                                        
                                        // evaluated functions
                                        Complex Bx = B_x[(ixspline * (order + 1) + ixknot - ixspline) * points + ix];
                                        Complex By = B_y[iy];
                                        Real Pix = Pi_x[ixknot * points + ix];
                                        Real Piy = rad_inner_->bspline().eval(Xp, ry).real();
                                        Real jix = ji_x[ixknot * points + ix];
                                        Real jiy = special::ric_j(l, ki[0] * ry);
                                        Complex wx = xws[ixknot * points + ix];
                                        Complex wy = yws[iy];
                                        
                                        // damp factor
                                        Real dampfactor = 1;//damp(rx, ry, distance);
                                        
                                        // monopole contribution
                                        if (rx > ry and li == l1 and l == l2)                    contrib_direct   += Bx * By * (1.0_r/rx - 1.0_r/ry) * Pix * jiy * dampfactor * wx * wy;
                                        if (ry > rx and li == l2 and l == l1 and inp_->exchange) contrib_exchange += Bx * By * (1.0_r/ry - 1.0_r/rx) * jix * Piy * dampfactor * wx * wy;
                                        
                                        // higher multipoles contribution
                                        Real multipole = 1 / rmax, rmin_over_rmax = rmin / rmax;
                                        for (int lambda = 1; lambda <= rad_full_->maxlambda(); lambda++)
                                        {
                                            multipole *= rmin_over_rmax;
                                            if (f1[lambda] != 0)                    contrib_direct   += f1[lambda] * Bx * By * multipole * Pix * jiy * dampfactor * wx * wy;
                                            if (f2[lambda] != 0 and inp_->exchange) contrib_exchange += f2[lambda] * Bx * By * multipole * jix * Piy * dampfactor * wx * wy;
                                        }
                                    }
                                }
                                
                                // for all quadrature points from the triangle x > y
                                for (int ix = 0; ix < points; ix++)
                                {
                                    cArray ys (points), yws (points), B_y (points);
                                    rad_full_->gaussleg_x().scaled_nodes_and_weights(points, rad_full_->bspline().t(iyknot), xs[ixknot * points + ix], &ys[0], &yws[0]);
                                    rad_full_->bspline().B(iyspline, iyknot, points, &ys[0], &B_y[0]);
                                    
                                    for (int iy = 0; iy < points; iy++)
                                    {
                                        // radii
                                        Real rx = xs[ixknot * points + ix].real(), ry = ys[iy].real(), rmin = std::min(rx,ry), rmax = std::max(rx,ry);
                                        
                                        // evaluated functions
                                        Complex Bx = B_x[(ixspline * (order + 1) + ixknot - ixspline) * points + ix];
                                        Complex By = B_y[iy];
                                        Real Pix = Pi_x[ixknot * points + ix];
                                        Real Piy = rad_inner_->bspline().eval(Xp, ry).real();
                                        Real jix = ji_x[ixknot * points + ix];
                                        Real jiy = special::ric_j(l, ki[0] * ry);
                                        Complex wx = xws[ixknot * points + ix];
                                        Complex wy = yws[iy];
                                        
                                        // damp factor
                                        Real dampfactor = 1;//damp(rx, ry, distance);
                                        
                                        // monopole contribution
                                        if (rx > ry and li == l1 and l == l2)                    contrib_direct   += Bx * By * (1.0_r/rx - 1.0_r/ry) * Pix * jiy * dampfactor * wx * wy;
                                        if (ry > rx and li == l2 and l == l1 and inp_->exchange) contrib_exchange += Bx * By * (1.0_r/ry - 1.0_r/rx) * jix * Piy * dampfactor * wx * wy;
                                        
                                        // higher multipoles contribution
                                        Real multipole = 1 / rmax, rmin_over_rmax = rmin / rmax;
                                        for (int lambda = 1; lambda <= rad_full_->maxlambda(); lambda++)
                                        {
                                            multipole *= rmin_over_rmax;
                                            if (f1[lambda] != 0)                    contrib_direct   += f1[lambda] * Bx * By * multipole * Pix * jiy * dampfactor * wx * wy;
                                            if (f2[lambda] != 0 and inp_->exchange) contrib_exchange += f2[lambda] * Bx * By * multipole * jix * Piy * dampfactor * wx * wy;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                
                // update element of the right-hand side
                if (inp_->Zp > 0)
                {
                    chi_block[ispline] += -prefactor * contrib_direct;
                }
                else
                {
                    chi_block[ispline] += prefactor * (contrib_direct + Sign * contrib_exchange) / special::constant::sqrt_two;
                }
            }
        }
        
        // use the calculated block
        chi[ill] = chi_block;
        
        /// v --- v DEBUG
        std::ofstream out (format("chi-%d.vtk", ill));
        writeVTK_points
        (
            out,
            Bspline::zip
            (
                rad_full().bspline_x(),
                rad_full().bspline_y(),
                chi_block,
                linspace(rad_full().bspline_x().Rmin(), rad_full().bspline_x().Rmax(), 1001),
                linspace(rad_full().bspline_y().Rmin(), rad_full().bspline_y().Rmax(), 1001)
            ),
            linspace(rad_full().bspline_x().Rmin(), rad_full().bspline_x().Rmax(), 1001),
            linspace(rad_full().bspline_y().Rmin(), rad_full().bspline_y().Rmax(), 1001),
            rArray{ 0. }
        );
        /// ^ --- ^ DEBUG
        
        // optionally transfer to disk
        if (not chi.inmemory())
        {
            chi.hdfsave(ill);
            chi[ill].drop();
        }
    }
}

void NoPreconditioner::multiply (BlockArray<Complex> const & p, BlockArray<Complex> & q, MatrixSelection::Selection tri) const
{
    // full basis
    unsigned Nspline_inner = rad_inner_->bspline_x().Nspline();
    unsigned Nspline_full  = rad_full_->bspline_x().Nspline();
    unsigned Nspline_outer = Nspline_full - Nspline_inner;
    
    // panel x basis
    unsigned Nspline_x_full  = rad_panel_->bspline_x().Nspline();
    unsigned Nspline_x_inner = Nspline_x_full; // FIXME
    unsigned Nspline_x_outer = Nspline_x_full - Nspline_x_inner;
    
    // panel y basis
    unsigned Nspline_y_full  = rad_panel_->bspline_y().Nspline();
    unsigned Nspline_y_inner = Nspline_y_full; // FIXME
    unsigned Nspline_y_outer = Nspline_y_full - Nspline_y_inner;
    
    // number of diagonal blocks
    unsigned Nang = ang_->states().size();
    
    // make sure no process is playing with the data
    par_->wait();
    
    // TODO : It is slightly more subtle to do this efficiently in out-of-core mode, so we are just
    //        loading the destination vectors here. But it would be possible to rewrite the code to load
    //        only the necessary pieces of those vectors on the fly.
    
    // de-const-ed source vector reference
    BlockArray<Complex> & v = const_cast<BlockArray<Complex> &>(p);
    
    // in distributed case we first need to collect the whole source vector
    for (unsigned ill = 0; ill < Nang; ill++)
    {
        // calculate rank of the process that owns this source vector segment (use group master process)
        int owner = (ill % par_->Ngroup()) * par_->groupsize();
        
        // load the source segment from disk, if necessary
        if (par_->iproc() == owner and cmd_->outofcore)
            v.hdfload(ill);
        
        // broadcast the source segment from the owner process to all others
        par_->bcast(owner, v[ill]);
        
        // also load the destination segment
        if (par_->isMyGroupWork(ill) and cmd_->outofcore)
            q.hdfload(ill);
    }
    
    // for all angular block rows
    for (unsigned ill = 0; ill < Nang; ill++) if (par_->isMyGroupWork(ill))
    {
        std::memset(q[ill].data(), 0, q[ill].size() * sizeof(Complex));
        
        // for all angular blocks in a block row; only executed by one of the processes in a process group
        for (unsigned illp = 0; illp < Nang; illp++) if (par_->igroupproc() == (int)illp % par_->groupsize())
        {
            // determine which part of the block should be considered non-zero for a particlar selection
            MatrixSelection::Selection selection = tri;
            if (ill < illp) selection = (tri & MatrixSelection::StrictUpper ? MatrixSelection::Both : MatrixSelection::None);
            if (ill > illp) selection = (tri & MatrixSelection::StrictLower ? MatrixSelection::Both : MatrixSelection::None);
            
            // near-origin part multiplication
            if (cmd_->lightweight_full)
            {
                // only one-electron contribution; the rest is below
                calc_A_block(ill, illp, false).dot
                (
                    1.0_z, cArrayView(p[illp], 0, Nspline_x_inner * Nspline_y_inner),
                    1.0_z, cArrayView(q[ill],  0, Nspline_x_inner * Nspline_y_inner),
                    true,
                    selection
                );
            }
            else
            {
                // read matrix from disk
                if (cmd_->outofcore and cmd_->wholematrix)
                    const_cast<BlockSymBandMatrix<Complex> &>(A_blocks_[ill * Nang + illp]).hdfload();
                
                // full diagonal block multiplication
                A_blocks_[ill * Nang + illp].dot
                (
                    1.0_z, cArrayView(p[illp], 0, Nspline_x_inner * Nspline_y_inner),
                    1.0_z, cArrayView(q[ill],  0, Nspline_x_inner * Nspline_y_inner),
                    true,
                    selection
                );
                
                // release memory
                if (cmd_->outofcore and cmd_->wholematrix)
                    const_cast<BlockSymBandMatrix<Complex> &>(A_blocks_[ill * Nang + illp]).drop();
            }
            
            // channel expansion part multiplication
            if (not inp_->inner_only)
            {
                int Nchan1 = Nchan_[ill].first;     // # r1 -> inf; l2 bound
                int Nchan2 = Nchan_[ill].second;    // # r2 -> inf; l1 bound
                int Nchan1p = Nchan_[illp].first;   // # r1 -> inf; l2p bound
                int Nchan2p = Nchan_[illp].second;  // # r2 -> inf; l1p bound
                
                // r1 -> inf
                # pragma omp parallel for
                for (int m = 0; m < Nchan1; m++)
                for (int n = 0; n < Nchan1p; n++)
                {
                    // read matrix from disk
                    if (cmd_->outofcore)
                        const_cast<SymBandMatrix<Complex>&>(B1_blocks_[ill * Nang + illp][m * Nchan1p + n]).hdfload();
                    
                    // multiply
                    B1_blocks_[ill * Nang + illp][m * Nchan1p + n].dot
                    (
                        1.0_z, cArrayView(p[illp], Nspline_x_inner * Nspline_y_inner + n * Nspline_x_outer, Nspline_x_outer),
                        1.0_z, cArrayView(q[ill],  Nspline_x_inner * Nspline_y_inner + m * Nspline_x_outer, Nspline_x_outer)
                    );
                    
                    // release memory
                    if (cmd_->outofcore)
                        const_cast<SymBandMatrix<Complex>&>(B1_blocks_[ill * Nang + illp][m * Nchan1p + n]).drop();
                }
                
                // r2 -> inf
                # pragma omp parallel for
                for (int m = 0; m < Nchan2; m++)
                for (int n = 0; n < Nchan2p; n++)
                {
                    // read matrix from disk
                    if (cmd_->outofcore)
                        const_cast<SymBandMatrix<Complex>&>(B2_blocks_[ill * Nang + illp][m * Nchan2p + n]).hdfload();
                    
                    // multiply
                    B2_blocks_[ill * Nang + illp][m * Nchan2p + n].dot
                    (
                        1.0_z, cArrayView(p[illp], Nspline_x_inner * Nspline_y_inner + Nchan1p * Nspline_x_outer + n * Nspline_y_outer, Nspline_y_outer),
                        1.0_z, cArrayView(q[ill],  Nspline_x_inner * Nspline_y_inner + Nchan1  * Nspline_x_outer + m * Nspline_y_outer, Nspline_y_outer)
                    );
                    
                    // release memory
                    if (cmd_->outofcore) const_cast<SymBandMatrix<Complex>&>(B2_blocks_[ill * Nang + illp][m * Nchan2p + n]).drop();
                }
                
                // multiply by coupling matrices
                Cu_blocks_[ill * Nang + illp].dot(1.0_z, p[illp], 1.0_z, q[ill]);
                Cl_blocks_[ill * Nang + illp].dot(1.0_z, p[illp], 1.0_z, q[ill]);
            }
        }
    }
    
    // lightweight-full off-diagonal contribution
    if (cmd_->lightweight_full)
    {
        OMP_CREATE_LOCKS(Nang * Nspline_x_inner);
        
        int maxlambda = rad_full_->maxlambda();
        
        # pragma omp parallel for collapse (3) schedule (dynamic,1)
        for (int lambda = 0; lambda <= maxlambda; lambda++)
        for (unsigned i = 0; i < Nspline_x_inner; i++)
        for (unsigned d = 0; d <= (unsigned)inp_->order; d++)
        if (i + d < Nspline_x_inner)
        {
            unsigned k = i + d;
            
            SymBandMatrix<Complex> R = rad_full_->calc_R_tr_dia_block(lambda, i, k);
            
            for (unsigned ill = 0; ill < Nang; ill++) if (par_->isMyGroupWork(ill))
            for (unsigned illp = 0; illp < Nang; illp++) if (par_->igroupproc() == (int)illp % par_->groupsize())
            if (Real f = ang_->f(ill, illp, lambda))
            {
                // diagonal blocks
                if (d == 0)
                {
                    OMP_LOCK_LOCK(ill * Nspline_x_inner + i);
                    
                    R.dot
                    (
                        inp_->Zp * f, cArrayView(p[illp], i * Nspline_y_inner, Nspline_y_inner),
                        1.0_z,        cArrayView(q[ill],  i * Nspline_y_inner, Nspline_y_inner),
                        ill == illp ? tri : (ill < illp ? (tri & MatrixSelection::StrictUpper ? MatrixSelection::Both : MatrixSelection::None) :
                                                          (tri & MatrixSelection::StrictLower ? MatrixSelection::Both : MatrixSelection::None))
                    );
                    
                    OMP_UNLOCK_LOCK(ill * Nspline_x_inner + i);
                }
                
                // off-diagonal blocks
                else
                {
                    if (tri & MatrixSelection::StrictUpper)
                    {
                        OMP_LOCK_LOCK(ill * Nspline_x_inner + i);
                        
                        R.dot
                        (
                            inp_->Zp * f, cArrayView(p[illp], k * Nspline_y_inner, Nspline_y_inner),
                            1.0_z,        cArrayView(q[ill],  i * Nspline_y_inner, Nspline_y_inner)
                        );
                        
                        OMP_UNLOCK_LOCK(ill * Nspline_inner + i);
                    }
                    
                    if (tri & MatrixSelection::StrictLower)
                    {
                        OMP_LOCK_LOCK(ill * Nspline_inner + k);
                        
                        R.dot
                        (
                            inp_->Zp * f, cArrayView(p[illp], i * Nspline_y_inner, Nspline_y_inner),
                            1.0_z,        cArrayView(q[ill],  k * Nspline_y_inner, Nspline_y_inner)
                        );
                        
                        OMP_UNLOCK_LOCK(ill * Nspline_x_inner + k);
                    }
                }
            }
        }
        
        OMP_DELETE_LOCKS();
    }
    
    // release source vectors
    for (unsigned ill = 0; ill < Nang; ill++) if (cmd_->outofcore)
        v[ill].drop();
    
    // synchronize and release the result vectors
    for (unsigned ill = 0; ill < Nang; ill++) if (par_->isMyGroupWork(ill))
    {
        // synchronize (sum) across the group
        par_->syncsum_g(q[ill].data(), q[ill].size());
        
        // constrain
        if (const CGPreconditioner * cgprec = dynamic_cast<const CGPreconditioner*>(this))
            cgprec->CG_constrain(q[ill]);
        
        // release memory
        if (cmd_->outofcore)
        {
            if (par_->IamGroupMaster())
                q.hdfsave(ill);
            
            q[ill].drop();
        }
    }
}

void NoPreconditioner::precondition (const BlockArray< Complex >& r, BlockArray< Complex >& z) const
{
    z = r;
}

void NoPreconditioner::finish ()
{
    A_blocks_ .resize(0);
    B1_blocks_.resize(0);
    B2_blocks_.resize(0);
    Cu_blocks_.resize(0);
    Cl_blocks_.resize(0);
    
    Xp_[0].resize(0);
    Sp_[0].resize(0);
    Eb_[0].resize(0);
    
    Xp_[1].resize(0);
    Sp_[1].resize(0);
    Eb_[1].resize(0);
}

// --------------------------------------------------------------------------------- //

addClassToParentRunTimeSelectionTable(PreconditionerBase, NoPreconditioner)

// --------------------------------------------------------------------------------- //
