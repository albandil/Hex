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
#include "KPAPreconditioner.h"

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
    block_rank_(ang.states().size()),
    rad_inner_(new RadialIntegrals(bspline_inner,   bspline_inner,   ang.maxlambda() + 1)),
    rad_full_ (new RadialIntegrals(bspline_full,    bspline_full,    ang.maxlambda() + 1)),
    rad_panel_(new RadialIntegrals(bspline_panel_x, bspline_panel_y, ang.maxlambda() + 1)),
    luS_(LUft::Choose("lapack"))
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
    if (verbose_) std::cout << "[ Inner region radial integrals ]" << std::endl << std::endl;
    rad_inner_->verbose(verbose_);
    rad_inner_->setupOneElectronIntegrals(*par_, *cmd_);
    rad_inner_->setupTwoElectronIntegrals(*par_, *cmd_);
    
    // calculate one-electron integrals in full domain
    if (verbose_) std::cout << "[ Full-domain radial integrals ]" << std::endl << std::endl;
    rad_full_->verbose(verbose_);
    rad_full_->setupOneElectronIntegrals(*par_, *cmd_);
    
    // calculate all radial integrals on panel
    if (verbose_) std::cout << "[ Current panel full radial integrals ]" << std::endl << std::endl;
    rad_panel_->verbose(verbose_);
    rad_panel_->setupOneElectronIntegrals(*par_, *cmd_);
    rad_panel_->setupTwoElectronIntegrals(*par_, *cmd_);
    
    // number of global inner basis B-splines
    std::size_t Nspline_inner = bspline_inner.Nspline();
    
    // is this panel the full solution domain?
    bool full_domain = bspline_x.hash() == bspline_full.hash()
                   and bspline_y.hash() == bspline_full.hash();
    
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
    
    // Diagonalize the global overlap matrix and find the eigenvectors of the one-electron hamiltonian.
    // The diagonalized Hamiltonians are used by the KPA preconditioner and the eigenstates
    // for the initial and final state and also to match the inner and outer problems.
    
    if (dynamic_cast<KPACGPreconditioner*>(this) or not cmd_->analytic_eigenstates)
    for (int i = 0; i < 2; i++)
    {
        // shorthands
        RadialIntegrals const * rint = full_domain ? rad_inner_ : rad_panel_;
        Bspline const & bspline = (i == 0 ? rint->bspline_x() : rint->bspline_y());
        std::size_t Nspline = bspline.Nspline();
        
        // particle charge (first is electron, second is either electron or positron)
        Real Z = (i == 0 ? -1.0 : inp_->Zp);
        
        // link all factorizations to the corresponding disk files
        for (int l : ells[i])
            Hl_[i][l].hdflink(format("Hl%+g-%d-%.4x.hdf", Z, l, bspline.hash()).c_str());
        
        // if all diagonalizations exist, do not do anything
        if (std::all_of(ells[i].begin(), ells[i].end(), [&](int l){ return Hl_[i][l].hdfcheck(); }))
            continue;
        
        if (verbose_) std::cout << "Setting up the hydrogen eigenstates for electron #" << i + 1 << " ..." << std::endl << std::endl;
        
        cArray D (Nspline);
        ColMatrix<Complex> CR (Nspline, Nspline);
        ColMatrix<Complex> invCR (Nspline, Nspline);
        ColMatrix<Complex> invsqrtS (Nspline, Nspline);
        
        if (verbose_) std::cout << "\t- basis overlap matrix diagonalization" << std::endl;
        Timer timer;
        
        SymBandMatrix<Complex> const & S_sym = (i == 0 ? rint->S_x() : rint->S_y());
        ColMatrix<Complex> S = std::move(S_sym.torow().T());
        S.diagonalize(D, nullptr, &CR);
        CR.invert(invCR);
        
        // Now S = CR * (D * CR⁻¹)
        if (verbose_) std::cout << "\t\t- time: " << timer.nice_time() << std::endl;
        for (std::size_t j = 0; j < Nspline * Nspline; j++)
            invCR.data()[j] *= D[j % Nspline];
        
        // S = S - CR * invCR
        blas::gemm(-1., CR, invCR, 1., S);
        if (verbose_) std::cout << "\t\t- residual: " << S.data().norm() << std::endl;
        
        // compute √S⁻¹
        for (std::size_t j = 0; j < Nspline * Nspline; j++)
            invCR.data()[j] /= std::pow(D.data()[j % Nspline], 1.5);
        blas::gemm(1., CR, invCR, 0., invsqrtS);
        
        if (verbose_) std::cout << std::endl;
        
        // for all one-electron angular momenta
        bool written = false;
        for (int l : ells[i])
        {
            // do not calculate if this work is supposed to be done by someone else
            if (cmd_->shared_scratch and not (par_->isMyGroupWork(l) and par_->IamGroupMaster()))
                continue;
            
            // check if the file already exists; skip calculation in that case
            if (Hl_[i][l].hdfcheck())
                continue;
            
            written = true;
            if (verbose_) std::cout << "\t- one-electron Hamiltonian matrix diagonalization (Z = " << Z << ", l = " << l << ")" << std::endl;
            timer.reset();
            
            // compose the symmetrical one-electron hamiltonian
            ColMatrix<Complex> tHl;
            if (i == 0)
                tHl = (Complex(0.5) * rint->D_x() + Complex(inp_->Za * Z) * rint->Mm1_x() + Complex(0.5*l*(l+1)) * rint->Mm2_x()).torow().T();
            else
                tHl = (Complex(0.5) * rint->D_y() + Complex(inp_->Za * Z) * rint->Mm1_y() + Complex(0.5*l*(l+1)) * rint->Mm2_y()).torow().T();
            
            // symmetrically transform by inverse square root of the overlap matrix, tHl <- invsqrtS * tHl * invsqrtS
            blas::gemm(1., invsqrtS, tHl, 0., S);
            blas::gemm(1., S, invsqrtS, 0., tHl);
            
            // diagonalize the transformed matrix
            tHl.diagonalize(D, nullptr, &CR);
            CR.invert(invCR);
            
            // store the KPA preconditioner data
            Hl_[i][l].Dl = D;
            Hl_[i][l].invsqrtS_Cl = std::move(RowMatrix<Complex>(Nspline, Nspline));
            Hl_[i][l].invCl_invsqrtS = std::move(RowMatrix<Complex>(Nspline, Nspline));
            blas::gemm(1., invsqrtS, CR, 0., Hl_[i][l].invsqrtS_Cl);
            blas::gemm(1., invCR, invsqrtS, 0., Hl_[i][l].invCl_invsqrtS);
            
            // Now Hl = ClR * D * ClR⁻¹
            if (verbose_) std::cout << "\t\t- time: " << timer.nice_time() << std::endl;
            for (std::size_t j = 0; j < Nspline * Nspline; j++)
                invCR.data()[j] *= D[j % Nspline];
            
            // Hl <- Hl - CR * invCR
            blas::gemm(-1., CR, invCR, 1., tHl);
            if (verbose_) std::cout << "\t\t- residual: " << tHl.data().norm() << std::endl;
            
            // copy the eigenvectors as columns
            // - already normalized by xGEEV to "Euclidean norm equal to 1 and largest component real"
            Hl_[i][l].Cl = std::move(ColMatrix<Complex>(Hl_[i][l].invsqrtS_Cl));
            
            // write to disk and abandon for now
            Hl_[i][l].hdfsave();
            Hl_[i][l].drop();
        }
        if (written)
        {
            if (verbose_) std::cout << std::endl;
        }
        
        // wait for all processes so that all Hlxxx files get written
        par_->wait();
    }
    
    // LU decomposition of the overlap matrix
    CsrMatrix<LU_int_t,Complex> csr_S = rad_inner_->S().tocoo<LU_int_t>().tocsr();
    luS_->factorize(csr_S);
    
    // load the requested one-electron eigenstates from disk files (only needed for full domain, where the RHS is constructed)
    if (full_domain) for (int i = 0; i < 2; i++)
    {
        if (verbose_) std::cout << "Loading the hydrogen eigenstates for electron #" << i + 1 << " ..." << std::endl;
        
        // particle charge (first is electron, second is either electron or positron)
        Real Z = (i == 0 ? -1.0 : inp_->Zp);
        
        // for all one-electron angular momenta needed by this particle
        for (int l : ells[i])
        {
            std::vector<int> indices (Nspline_inner);
            
            // load the factorization file
            bool loaded = Hl_[i][l].hdfload();
            
            if (not loaded)
            {
                if (not cmd_->analytic_eigenstates)
                {
                    HexException("Failed to load one-electron diagonalization file %s for Z = %g, l = %d.", Hl_[i][l].filename.c_str(), Z, l);
                }
                else
                {
                    if (verbose_)
                        std::cout << "\t- one-electron Hamiltonian data skipped (Z = " << Z << ", l = " << l << ")" << std::endl;
                }
            }
            else
            {
                if (verbose_) std::cout << "\t- one-electron Hamiltonian data loaded (Z = " << Z << ", l = " << l << ")" << std::endl;
                
                // get sorted energies (ascending real parts)
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
                    if (verbose_) std::cout << "\t\t- bound states with energy within 0.1 % from exact value: " << l + 1 << " <= n <= " << max_nr + l + 1 << std::endl;
                }
            }
            
            // get all valid asymptotic states
            for (unsigned nr = 0; nr < Nspline_inner; nr++)
            {
                // bound energy of the state (Ry)
                Complex Eb = loaded ? 2.0_r * Hl_[i][l].Dl[indices[nr]] : -1.0_r / ((nr + l + 1) * (nr + l + 1));
                
                // add all requested channels
                if
                (
                    // in inner-region-only calculation we need just all required bound channels
                        (inp_->inner_only and Eb.real() < inp_->max_Ebound) or
                    // in channel-reduced calculation we need also all states up to the total energy of the system (or higher, if requested)
                        (not inp_->inner_only and Eb.real() <= mmax(inp_->channel_max_E, inp_->max_Etot, inp_->max_Ebound))
                )
                {
                    if (loaded)
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
                    else if (Eb.real() >= 0)
                    {
                        // This is a positive-energy pseudo-bound state, but as we are supposed to use only
                        // analytic eigenstates, there is no way how to get it. So we will terminate here and
                        // consider no more scattering channels.
                        
                        break;
                    }
                    else
                    {
                        // Evaluate the analytic formula for the bound state. This is useful for inner-region-only calculation,
                        // where only the initial and final states need to be calculated, and for channel-reduced calculation
                        // sufficiently below zero, where not many states exist. Otherwise this will blow the memory of the
                        // computer.
                        
                        Sp_[i][l].push_back(RadialIntegrals::overlapP(rad_inner_->bspline(), rad_inner_->gaussleg(), inp_->Za, nr + l + 1, l));
                        Xp_[i][l].push_back(luS_->solve(Sp_[i][l].back(), 1));
                        Eb_[i][l].push_back(Eb);
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
                if (verbose_) std::cout << "\t\t- asymptotical scattering channels: " << l + 1 << " <= n <= " << l + Xp_[i][l].size() << std::endl;
            }
            
            // unload the factorization file
            Hl_[i][l].drop();
        }
        if (verbose_) std::cout << std::endl;
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
    // inner B-spline count
    int Nspline_x_inner = rad_panel_->bspline_x().hash() == rad_full_->bspline().hash() ? rad_inner_->bspline().Nspline() : rad_panel_->bspline_x().Nspline();
    int Nspline_y_inner = rad_panel_->bspline_y().hash() == rad_full_->bspline().hash() ? rad_inner_->bspline().Nspline() : rad_panel_->bspline_y().Nspline();
    
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
                    Complex Sx = rad_panel_->S_x()(i,k), Dx = rad_panel_->D_x()(i,k), Mm1x = rad_panel_->Mm1_x()(i,k), Mm2x = rad_panel_->Mm2_x()(i,k);
                    Complex Sy = rad_panel_->S_y()(j,l), Dy = rad_panel_->D_y()(j,l), Mm1y = rad_panel_->Mm1_y()(j,l), Mm2y = rad_panel_->Mm2_y()(j,l);
                    
                    return E_ * Sx * Sy - 0.5_r * Dx * Sy + inp_->Za *            Mm1x * Sy - 0.5_r * l1 * (l1 + 1) * Mm2x * Sy
                                        - 0.5_r * Dy * Sx - inp_->Za * inp_->Zp * Mm1y * Sx - 0.5_r * l2 * (l2 + 1) * Mm2y * Sx;
                }
            );
        }
        
        // two-electron part
        if(twoel)
        for (int lambda = 0; lambda <= ang_->maxlambda(); lambda++) if (ang_->f(ill,illp,lambda) != 0)
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
                subblock.data() += inp_->Zp * ang_->f(ill,illp,lambda) * rad_panel_->calc_R_tr_dia_block(lambda, i, k).data().slice(0, Nspline_y_inner * (inp_->order + 1));
            }
        }
        
        // save block
        A.setBlock(i * (inp_->order + 1) + d, subblock.data());
    }
    
    if (cmd_->outofcore)
        A.drop();
    
    return A;
}

CooMatrix<LU_int_t, Complex> NoPreconditioner::calc_full_block (int ill, int illp) const
{
    // number of asymptotic channels
    int Nchan1 = Nchan_[ill].first;
    int Nchan2 = Nchan_[ill].second;
    int Nchan1p = Nchan_[illp].first;
    int Nchan2p = Nchan_[illp].second;
    
    // B-spline bases
    Bspline const & bspline_full  = rad_full_ ->bspline();
    Bspline const & bspline_inner = rad_inner_->bspline();
    Bspline const & bspline_x     = rad_panel_->bspline_x();
    Bspline const & bspline_y     = rad_panel_->bspline_y();
    
    // panel x basis
    LU_int_t Nspline_full_x  = bspline_x.Nspline();
    LU_int_t Nspline_inner_x = bspline_x.hash() == bspline_full.hash() ? bspline_inner.Nspline() : bspline_x.Nspline();
    LU_int_t Nspline_outer_x = Nspline_full_x - Nspline_inner_x;
    
    // panel y basis
    LU_int_t Nspline_full_y  = bspline_y.Nspline();
    LU_int_t Nspline_inner_y = bspline_y.hash() == bspline_full.hash() ? bspline_inner.Nspline() : bspline_y.Nspline();
    LU_int_t Nspline_outer_y = Nspline_full_y - Nspline_inner_y;
    
    // angular block
    int iang = ill * ang_->states().size() + illp;
    
    // convert inner region matrix block to COO matrix
    CooMatrix<LU_int_t,Complex> coo_block;
    if (cmd_->lightweight_full)
        coo_block = std::move(calc_A_block(ill, illp).tocoo<LU_int_t>());
    else
        coo_block = std::move(A_blocks_[iang].tocoo<LU_int_t>());
    
    // add the A-block
    coo_block.resize
    (
        Nspline_inner_x * Nspline_inner_y + Nchan1  * Nspline_outer_x + Nchan2  * Nspline_outer_y,
        Nspline_inner_x * Nspline_inner_y + Nchan1p * Nspline_outer_x + Nchan2p * Nspline_outer_y
    );
    
    if (not inp_->inner_only)
    {
        // add the outer region C-blocks
        coo_block += Cu_blocks_[iang];
        coo_block += Cl_blocks_[iang];
        
        // reserve memory for B-blocks; otherwise the += operation below would perpetually
        // reallocate the arrays, making the addition of new elements rather slow, especially
        // for a large number of channels
        std::size_t addsize = Nchan1 * Nchan1p * Nspline_outer_x * (2*inp_->order + 1)
                            + Nchan2 * Nchan2p * Nspline_outer_y * (2*inp_->order + 1);
        const_cast<NumberArray<LU_int_t>&>(coo_block.i()).reserve(coo_block.i().size() + addsize);
        const_cast<NumberArray<LU_int_t>&>(coo_block.j()).reserve(coo_block.j().size() + addsize);
        const_cast<NumberArray<Complex>& >(coo_block.v()).reserve(coo_block.v().size() + addsize);
        
        // add the B-blocks
        for (int m = 0; m < Nchan1; m++)
        for (int n = 0; n < Nchan1p; n++)
        {
            if (cmd_->outofcore) const_cast<SymBandMatrix<Complex>&>(B1_blocks_[iang][m * Nchan1p + n]).hdfload();
            CooMatrix<LU_int_t,Complex> B_coo_small = B1_blocks_[iang][m * Nchan1p + n].tocoo<LU_int_t>();
            if (cmd_->outofcore) const_cast<SymBandMatrix<Complex>&>(B1_blocks_[iang][m * Nchan1p + n]).drop();
            CooMatrix<LU_int_t,Complex> B_coo_large
            (
                coo_block.rows(), coo_block.cols(),
                B_coo_small.i() + Nspline_inner_x * Nspline_inner_y + m * Nspline_outer_x,
                B_coo_small.j() + Nspline_inner_x * Nspline_inner_y + n * Nspline_outer_x,
                B_coo_small.v()
            );
            coo_block += B_coo_large;
        }
        for (int m = 0; m < Nchan2; m++)
        for (int n = 0; n < Nchan2p; n++)
        {
            if (cmd_->outofcore) const_cast<SymBandMatrix<Complex>&>(B2_blocks_[iang][m * Nchan2p + n]).hdfload();
            CooMatrix<LU_int_t,Complex> B_coo_small = B2_blocks_[iang][m * Nchan2p + n].tocoo<LU_int_t>();
            if (cmd_->outofcore) const_cast<SymBandMatrix<Complex>&>(B2_blocks_[iang][m * Nchan2p + n]).drop();
            CooMatrix<LU_int_t,Complex> B_coo_large
            (
                coo_block.rows(), coo_block.cols(),
                B_coo_small.i() + Nspline_inner_x * Nspline_inner_y + Nchan1 * Nspline_outer_x + m * Nspline_outer_y,
                B_coo_small.j() + Nspline_inner_x * Nspline_inner_y + Nchan1p * Nspline_outer_x + n * Nspline_outer_y,
                B_coo_small.v()
            );
            coo_block += B_coo_large;
        }
    }
    
    return coo_block;
}

void NoPreconditioner::update (Real E)
{
    // shorthands
    int order = inp_->order;
    int Nang = ang_->states().size();
    
    // B-spline bases
    Bspline const & bspline_full  = rad_full_ ->bspline();
    Bspline const & bspline_inner = rad_inner_->bspline();
    Bspline const & bspline_x     = rad_panel_->bspline_x();
    Bspline const & bspline_y     = rad_panel_->bspline_y();
    
    // global basis
    int Nspline_full  = bspline_full.Nspline();
    int Nspline_inner = bspline_inner.Nspline();
    int Nspline_outer = Nspline_full - Nspline_inner;
    
    // panel x basis
    int Nspline_x_full  = bspline_x.Nspline();
    int Nspline_x_inner = bspline_x.hash() == bspline_full.hash() ? bspline_inner.Nspline() : bspline_x.Nspline();
    int Nspline_x_outer = Nspline_x_full - Nspline_x_inner;
    
    // panel y basis
    int Nspline_y_full  = bspline_y.Nspline();
    int Nspline_y_inner = bspline_y.hash() == bspline_full.hash() ? bspline_inner.Nspline() : bspline_y.Nspline();
    int Nspline_y_outer = Nspline_y_full - Nspline_y_inner;
    
    // size of the A-matrix
    std::size_t A_size = std::size_t(Nspline_x_inner) * std::size_t(Nspline_y_inner);
    
    // update energy
    E_ = E;
    
    if (verbose_) std::cout << "\tUpdate the common preconditioner base" << std::endl;
    
    if (verbose_) std::cout << "\tPrecompute matrix blocks ... " << std::flush;
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
        
        block_rank_[ill] = A_size + Nchan_[ill].first * Nspline_x_outer + Nchan_[ill].second * Nspline_y_outer;
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
        if (not cmd_->lightweight_full)
            A_blocks_[ill * Nang + illp] = std::move(calc_A_block(ill, illp));
        
        // create inner-outer coupling blocks
        Cu_blocks_[ill * Nang + illp] = CooMatrix<LU_int_t,Complex>
        (
            block_rank_[ill],
            block_rank_[illp]
        );
        Cl_blocks_[ill * Nang + illp] = CooMatrix<LU_int_t,Complex>
        (
            block_rank_[ill],
            block_rank_[illp]
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
                         + inp_->Za * rad_full_->Mm1_x()(i,k) * rad_full_->S_y()(j,l)
                         - inp_->Za * inp_->Zp * rad_full_->S_x()(i,k) * rad_full_->Mm1_y()(j,l);
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
                         + inp_->Za * rad_full_->Mm1_x()(i,k) * rad_full_->S_y()(j,l)
                         - inp_->Za * inp_->Zp * rad_full_->S_x()(i,k) * rad_full_->Mm1_y()(j,l);
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
            for (int l = j - order; l < Nspline_inner; l++)
            {
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
                
                for (int k = 0; k < Nspline_inner; k++)
                {
                    std::size_t row = A_size + (Nchan1 + m) * Nspline_outer + (j - Nspline_inner);
                    std::size_t col = k * Nspline_inner + l;
                    
                    Cl_blocks_[ill * Nang + illp].add(row, col, Sp_[0][l1p][n][k] * elem);
                }
            }
            
            // transition area r1 > r2, lower : F_nk expressed in terms of psi_kl for 'k' out of outer area
            for (int m = 0; m < Nchan1; m++)
            for (int n = 0; n < Nchan1p; n++)
            for (int i = Nspline_inner; i < Nspline_full; i++) // *
            for (int k = i - order; k < Nspline_inner; k++)
            {
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
                
                for (int l = 0; l < Nspline_inner; l++)
                {
                    std::size_t row = A_size + m * Nspline_outer + (i - Nspline_inner);
                    std::size_t col = k * Nspline_inner + l;
                    
                    Cl_blocks_[ill * Nang + illp].add(row, col, Sp_[1][l2p][n][l] * elem);
                }
            }
        }
    }
    
    if (verbose_) std::cout << "done after " << t.nice_time() << std::endl;
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
    Real ki = std::sqrt(inp_->Etot[ie] + 1.0_r/(ni*ni));
    
    // get the initial bound pseudo-state B-spline expansion
    cArray Xp = cmd_->analytic_eigenstates ?
        luS_->solve(RadialIntegrals::overlapP(rad_inner_->bspline(), rad_inner_->gaussleg(), inp_->Za, ni, li), 1) :
        Hl_[0][li].readPseudoState(li, ni - li - 1);
    
    // calculate Ricatti-Bessel B-spline expansions
    cArray XJ = luS_->solve(RadialIntegrals::overlapj(rad_inner_->bspline(), rad_inner_->gaussleg(), ang_->maxell(), rArray{ ki }, cmd_->fast_bessel), ang_->maxell() + 1);
    
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
            // Ricatti-Bessel function B-spline overlaps for this angular momentum (only for the inner region)
            cArrayView Xj (XJ, l * Nspline_inner, Nspline_inner);
            
            // skip wrong parity
            if ((ang_->L() + li + l) % 2 != ang_->Pi())
                continue;
            
            // compute energy- and angular momentum-dependent prefactor
            Complex prefactor = std::pow(1.0_i,l)
                              * std::sqrt(4.0_r * special::constant::pi * (2 * l + 1))
                              * (Real)special::ClebschGordan(li,mi, l,0, inp_->L,mi) / ki;
            
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
            
            // precompute hydrogen radial functions
            rArray Pi_x = realpart(rad_inner_->bspline().zip(Xp, realpart(xs)));
            
            // evaluate projectile radial functions (Ricatti-Bessel or Coulomb)
            rArray ji_x (xs.size());
            # pragma omp parallel for
            for (unsigned ix = 0; ix < xs.size(); ix++)
                ji_x[ix] = special::ric_j(l, ki * xs[ix].real()); // FIXME : Coulomb for charge (inp_->Za - 1)
            
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
                    
                    // skip B-splines too close to the complex region
                    if (ixspline + order + 1 > rad_full().bspline_x().iR2())
                        continue;
                    
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
                    
                    // skip B-splines too close to the complex region
                    if (iyspline + order + 1 > rad_full().bspline_y().iR2())
                        continue;
                    
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
                    
                    // skip B-splines too close to the complex region
                    if (ixspline + order + 1 > rad_full().bspline_x().iR2() or iyspline + order + 1 > rad_full().bspline_y().iR2())
                        continue;
                    
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
                                    
                                    Real rx = xs[ixknot * points + ix].real();
                                    Complex Bx = B_x[(ixspline * (order + 1) + ixknot - ixspline) * points + ix];
                                    Real Pix = Pi_x[ixknot * points + ix];
                                    Real jix = ji_x[ixknot * points + ix];
                                    Complex wx = xws[ixknot * points + ix];
                                    
                                    for (int iy = 0; iy < points; iy++)
                                    {
                                        Real ry = ys[iy].real();
                                        Complex By = B_y[iy];
                                        Complex wy = yws[iy];
                                        
                                        Real Piy = rad_inner_->bspline().eval(Xp, ry).real();
                                        Real jiy = rad_inner_->bspline().R2() == rad_full_->bspline_y().R2()
                                                 ? rad_inner_->bspline().eval(Xj, ry).real() // if only inner region - just evaluate the precomputed B-spline expansion
                                                 : special::ric_j(l, ki * ry); // if also outer region - evaluate the wave now
                                        
                                        // damp factor
                                        Real dampfactor = 1;//damp(rx, ry, distance);
                                        
                                        // monopole contribution
                                        if (rx > ry and li == l1 and l == l2)                    contrib_direct   += Bx * By * (1.0_r/rx - 1.0_r/ry) * Pix * jiy * dampfactor * wx * wy;
                                        if (ry > rx and li == l2 and l == l1 and inp_->exchange) contrib_exchange += Bx * By * (1.0_r/ry - 1.0_r/rx) * jix * Piy * dampfactor * wx * wy;
                                        
                                        // higher multipoles contribution
                                        Real rmin = std::min(rx,ry), rmax = std::max(rx,ry);
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
                                    
                                    Real rx = xs[ixknot * points + ix].real();
                                    Complex Bx = B_x[(ixspline * (order + 1) + ixknot - ixspline) * points + ix];
                                    Real Pix = Pi_x[ixknot * points + ix];
                                    Real jix = ji_x[ixknot * points + ix];
                                    Complex wx = xws[ixknot * points + ix];
                                    
                                    for (int iy = 0; iy < points; iy++)
                                    {
                                        Real ry = ys[iy].real();
                                        Complex By = B_y[iy];
                                        Complex wy = yws[iy];
                                        
                                        Real Piy = rad_inner_->bspline().eval(Xp, ry).real();
                                        Real jiy = rad_inner_->bspline().R2() == rad_full_->bspline_y().R2()
                                                 ? rad_inner_->bspline().eval(Xj, ry).real() // if only inner region - just evaluate the precomputed B-spline expansion
                                                 : special::ric_j(l, ki * ry); // if also outer region - evaluate the wave now
                                        
                                        // damp factor
                                        Real dampfactor = 1;//damp(rx, ry, distance);
                                        
                                        // monopole contribution
                                        if (rx > ry and li == l1 and l == l2)                    contrib_direct   += Bx * By * (1.0_r/rx - 1.0_r/ry) * Pix * jiy * dampfactor * wx * wy;
                                        if (ry > rx and li == l2 and l == l1 and inp_->exchange) contrib_exchange += Bx * By * (1.0_r/ry - 1.0_r/rx) * jix * Piy * dampfactor * wx * wy;
                                        
                                        // higher multipoles contribution
                                        Real rmin = std::min(rx,ry), rmax = std::max(rx,ry);
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
        chi[ill] = std::move(chi_block);
        
        // optionally transfer to disk
        if (not chi.inmemory())
        {
            if (not cmd_->shared_scratch and par_->IamGroupMaster())
                chi.hdfsave(ill);
            
            chi.drop(ill);
        }
    }
}

void NoPreconditioner::multiply (BlockArray<Complex> const & p, BlockArray<Complex> & q, MatrixSelection::Selection tri) const
{
    // shorthands
    unsigned order = inp_->order;
    unsigned Nang = ang_->states().size();
    
    // B-spline bases
    Bspline const & bspline_full  = rad_full_ ->bspline();
    Bspline const & bspline_inner = rad_inner_->bspline();
    Bspline const & bspline_x     = rad_panel_->bspline_x();
    Bspline const & bspline_y     = rad_panel_->bspline_y();
    
    // panel x basis
    std::size_t Nspline_x_full  = bspline_x.Nspline();
    std::size_t Nspline_x_inner = bspline_x.hash() == bspline_full.hash() ? bspline_inner.Nspline() : bspline_x.Nspline();
    std::size_t Nspline_x_outer = Nspline_x_full - Nspline_x_inner;
    
    // panel y basis
    std::size_t Nspline_y_full  = bspline_y.Nspline();
    std::size_t Nspline_y_inner = bspline_y.hash() == bspline_full.hash() ? bspline_inner.Nspline() : bspline_y.Nspline();
    std::size_t Nspline_y_outer = Nspline_y_full - Nspline_y_inner;
    
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
    
    // get number of right-hand sides (initial states)
    unsigned Nini = v[0].size() / block_rank_[0];
    
    // in parallel, only the master process and its fellows from the MPI group have the complete segment p[0]
    // so to get correct number of right-hand sides we need to let the process broadcast the number;
    par_->bcast(0, &Nini, 1);
    
    // for all initial states (right-hand sides)
    for (unsigned ini = 0; ini < Nini; ini++)
    {
        // for all angular block rows
        for (unsigned ill = 0; ill < Nang; ill++) if (par_->isMyGroupWork(ill))
        {
            std::size_t offset = block_rank_[ill] * ini;
            cArrayView(q[ill], offset, block_rank_[ill]).fill(0.0_z);
            
            // for all angular blocks in a block row; only executed by one of the processes in a process group
            for (unsigned illp = 0; illp < Nang; illp++) if (par_->igroupproc() == (int)illp % par_->groupsize())
            {
                std::size_t offsetp = block_rank_[illp] * ini;
                
                // skip unwanted angular blocks
                if (ill == illp and not (tri & MatrixSelection::BlockDiagonal))    continue;
                if (ill <  illp and not (tri & MatrixSelection::BlockStrictUpper)) continue;
                if (ill >  illp and not (tri & MatrixSelection::BlockStrictLower)) continue;
                
                // determine which part of the block should be considered non-zero for a particlar selection
                MatrixSelection::Selection selection = tri;
                if (ill < illp) selection = (tri & MatrixSelection::StrictUpper ? MatrixSelection::Both : MatrixSelection::None);
                if (ill > illp) selection = (tri & MatrixSelection::StrictLower ? MatrixSelection::Both : MatrixSelection::None);
                
                // inner-region subset of the vectors
                cArrayView p_inner (p[illp], offsetp, Nspline_x_inner * Nspline_y_inner);
                cArrayView q_inner (q[ill],  offset,  Nspline_x_inner * Nspline_y_inner);
                
                // inner region part multiplication
                if (cmd_->lightweight_full)
                {
                    // multiply 'p' by the diagonal block (except for the two-electron term, which is done later)
                    if (ill == illp)
                    {
                        // get block angular momemnta
                        int l1 = ang_->states()[ill].first;
                        int l2 = ang_->states()[ill].second;
                        
                        // one-electron matrices
                        SymBandMatrix<Complex> const & Sx = rad_panel_->S_x();
                        SymBandMatrix<Complex> const & Sy = rad_panel_->S_y();
                        SymBandMatrix<Complex> Hx = 0.5_z * rad_panel_->D_x() + (0.5_z * (l1 * (l1 + 1.0_r))) * rad_panel_->Mm2_x() + Complex(inp_->Za *   -1.0_r) * rad_panel_->Mm1_x();
                        SymBandMatrix<Complex> Hy = 0.5_z * rad_panel_->D_y() + (0.5_z * (l2 * (l2 + 1.0_r))) * rad_panel_->Mm2_y() + Complex(inp_->Za * inp_->Zp) * rad_panel_->Mm1_y();
                        
                        // multiply 'p' by the diagonal block (except for the two-electron term)
                        kron_dot(1., q_inner, E_, p_inner, Sx, Sy, Nspline_x_inner, Nspline_y_inner);
                        kron_dot(1., q_inner, -1, p_inner, Hx, Sy, Nspline_x_inner, Nspline_y_inner);
                        kron_dot(1., q_inner, -1, p_inner, Sx, Hy, Nspline_x_inner, Nspline_y_inner);
                    }
                }
                else
                {
                    // read matrix from disk
                    if (cmd_->outofcore and cmd_->wholematrix)
                        const_cast<BlockSymBandMatrix<Complex> &>(A_blocks_[ill * Nang + illp]).hdfload();
                    
                    // full diagonal block multiplication
                    A_blocks_[ill * Nang + illp].dot(1.0_z, p_inner, 1.0_z, q_inner, true, selection);
                    
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
                            1.0_z, cArrayView(p[illp], offsetp + Nspline_x_inner * Nspline_y_inner + n * Nspline_x_outer, Nspline_x_outer),
                            1.0_z, cArrayView(q[ill],  offset  + Nspline_x_inner * Nspline_y_inner + m * Nspline_x_outer, Nspline_x_outer)
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
                            1.0_z, cArrayView(p[illp], offsetp + Nspline_x_inner * Nspline_y_inner + Nchan1p * Nspline_x_outer + n * Nspline_y_outer, Nspline_y_outer),
                            1.0_z, cArrayView(q[ill],  offset  + Nspline_x_inner * Nspline_y_inner + Nchan1  * Nspline_x_outer + m * Nspline_y_outer, Nspline_y_outer)
                        );
                        
                        // release memory
                        if (cmd_->outofcore) const_cast<SymBandMatrix<Complex>&>(B2_blocks_[ill * Nang + illp][m * Nchan2p + n]).drop();
                    }
                    
                    // multiply by coupling matrices
                    Cu_blocks_[ill * Nang + illp].dot(1.0_z, cArrayView(p[illp], offsetp, block_rank_[illp]), 1.0_z, cArrayView(q[ill], offset, block_rank_[ill]));
                    Cl_blocks_[ill * Nang + illp].dot(1.0_z, cArrayView(p[illp], offsetp, block_rank_[illp]), 1.0_z, cArrayView(q[ill], offset, block_rank_[ill]));
                }
            }
        }
        
        // lightweight-full off-diagonal contribution
	if (cmd_->lightweight_full)
        {
            int maxlambda = rad_panel_->maxlambda();
            
            cArray R;
            cArray Rp (Nspline_y_inner);
            
            # pragma omp parallel for schedule (dynamic, 1) firstprivate (R, Rp)
            for (std::size_t i = 0; i < Nspline_x_inner; i++)
            for (std::size_t k = (i < order ? 0 : i - order); k < Nspline_x_inner and k <= i + order; k++)
            for (int lambda = 0; lambda <= maxlambda; lambda++)
            {
                bool RReady = false;
                
                for (unsigned illp = 0; illp < Nang; illp++) if (par_->igroupproc() == (int)illp % par_->groupsize())
                {
                    bool RpReady = false;
                    
                    for (unsigned ill = 0; ill < Nang; ill++) if (par_->isMyGroupWork(ill))
                    if (Real f = ang_->f(ill, illp, lambda))
                    {
                        // skip unwanted angular blocks
                        if (ill == illp and not (tri & MatrixSelection::BlockDiagonal))    continue;
                        if (ill <  illp and not (tri & MatrixSelection::BlockStrictUpper)) continue;
                        if (ill >  illp and not (tri & MatrixSelection::BlockStrictLower)) continue;
                        
                        std::size_t chunk = p[ill].size() / Nini;
                        std::size_t chunkp = p[illp].size() / Nini;
                        std::size_t offset = chunk * ini;
                        std::size_t offsetp = chunkp * ini;
                        
                        // calculate sub-block R_ik^lambda
                        if (not RReady)
                        {
                            R = rad_panel_->calc_R_tr_dia_block(lambda, i, k).data().slice(0, Nspline_y_inner * (order + 1));
                            RReady = true;
                        }
                        
                        // calculate scalar products of R_ik^lambda with appropriate sub-segments of p[illp]
                        if (not RpReady)
                        {
                            // precompute Z * R_ijkl * p_jl
                            SymBandMatrix<Complex>::sym_band_dot
                            (
                                Nspline_y_inner, order + 1, R,
                                inp_->Zp, cArrayView(p[illp], offsetp + k * Nspline_y_inner, Nspline_y_inner),
                                0.0_z,    Rp
                            );
                            
                            RpReady = true;
                        }
                        
                        for (std::size_t r = 0; r < Nspline_y_inner; r++)
                            q[ill][offset + i * Nspline_y_inner + r] += f * Rp[r];
                    }
                }
            }
        }
    }
    
    // release source vectors
    for (unsigned ill = 0; ill < Nang; ill++) if (cmd_->outofcore)
        v.drop(ill);
    
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
            if (not cmd_->shared_scratch or par_->IamGroupMaster())
                q.hdfsave(ill);
            
            q.drop(ill);
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
