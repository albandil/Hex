/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2013                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>

#include "arrays.h"
#include "gauss.h"
#include "input.h"
#include "itersolve.h"
#include "misc.h"
#include "preconditioners.h"
#include "radial.h"

cArray IC (cArrayView const & A, lArrayView const & I, lArrayView const & P)
{
    // this will be returned
    cArray LD(A.size());
    
    // check lengths
    assert(A.size() == (size_t)P.back());
    assert(I.size() == (size_t)P.back());
    
    // current row
    int irow = 0;
    
    // for all elements of the output array
    for (int pos = 0; pos < (int)LD.size(); pos++)
    {
        // get column index of this element
        int icol = I[pos];
        
        // is this a diagonal?
        if (icol == irow)
        {
            // Compute an element of D
            // - start by copying corresponding coefficient from A
            LD[pos] = A[pos];
            
            // - continue by subtracting all existing contributions from the current row
            //   (loop over ELEMENTS)
            for (int ielem = P[irow]; ielem < pos; ielem++)
                LD[pos] -= LD[ielem] * LD[ielem] * LD[P[I[ielem]+1]-1];
            
            irow++;
        }
        else
        {
            // Compute an element of L
            // - start by copying corresponding coefficient from A
            LD[pos] = A[pos];
            
            // - continue by subtracting all existing contributions
            //   (loop over COLUMNS)
            int pos1 = P[irow], pos2 = P[icol];
            while (pos1 < P[irow+1] and I[pos1] < icol and pos2 < P[icol+1] and I[pos2] < icol)
            {
                if (I[pos1] < I[pos2])
                {
                    pos1++;
                }
                else if (I[pos1] > I[pos2])
                {
                    pos2++;
                }
                else
                {
                    LD[pos] -= LD[pos1] * LD[pos2] * LD[P[I[pos1]+1]-1];
                    pos1++; pos2++;
                }
            }
            
            // - finish by dividing by the diagonal element
            LD[pos] /= LD[P[icol+1]-1];
        }
    }
    
    return LD;
}

SymDiaMatrix DIC (SymDiaMatrix const & A)
{
    #define THRESHOLD 1e-3
    
    //
    // compute the preconditioned diagonal
    //
    
    // the preconditioned diagonal
    cArray D = A.main_diagonal();
    
    // pointer to the preconditioned diagonal data
    Complex * const restrict pD = D.data();
    
    // pointer to the A's diagonal labels
    int const * restrict pAd = A.diag().data();
    
    // array sizes
    register int Nrows = D.size();
    register int Ndiag = A.diag().size();
    
    // for all elements of the diagonal
    for (register int irow = 0; irow < Nrows; irow++)
    {
        // pointer to the A's raw concatenated diagonal data
        Complex const * restrict pA = A.data().data() + Nrows;
        
        // for all diagonals (except the main diagonal) constributing to DIC
        for (register int idiag = 1; idiag < Ndiag and pAd[idiag] <= irow; idiag++)
        {
            // update pivot
            pD[irow] -= pA[irow - pAd[idiag]] * pA[irow - pAd[idiag]] / pD[irow - pAd[idiag]];
            
            // if the pivod is too small, set it to one
            if (std::abs(pD[irow]) < THRESHOLD)
                pD[irow] = THRESHOLD;
            
            // move pA to the beginning of the next diagonal
            pA += Nrows - pAd[idiag];
        }
    }
    
    //
    // construct matrix of the preconditioner
    //
    
    // the preconditioner matrix, initialized to A
    SymDiaMatrix DIC(A);
    
    // pointer to DIC's data
    Complex * restrict pDIC = DIC.data().data();
    
    // pointer to DIC's diagonal labels (same as A.diag())
    int const * restrict pDICd = DIC.diag().data();
    
    // set the diagonal to the inverse of the DIC diagonal
    for (register int irow = 0; irow < Nrows; irow++)
        pDIC[irow] = 1. / pD[irow];
    
    // move on to the first non-main diagonal
    pDIC += Nrows;
    
    // for all non-main diagonals
    for (register int idiag = 1; idiag < Ndiag; idiag++)
    {
        // for all elements in the diagonal
        for (register int irow = 0; irow < Nrows - pDICd[idiag]; irow++)
        {
            // divide the off-diagonal element by the correct diagonal element
            pDIC[irow] *= pD[irow];
        }
        
        // move to next diagonal
        pDIC +=  Nrows - pDICd[idiag];
    }
    
    return DIC;
}

SymDiaMatrix SSOR (SymDiaMatrix const & A)
{
    // array sizes
    register int Nrows = A.size();
    register int Ndiag = A.diag().size();
    
    // the preconditioner matrix, initialized to A
    SymDiaMatrix SSOR(A);
    
    // pointer to SSOR's concatenated diagonal data
    Complex * restrict pSSOR = SSOR.data().data();
    
    // pointer to SSOR's main diagonal
    Complex const * const restrict pD = pSSOR;
    
    // pointer to SSOR's diagonal labels (same as A.diag())
    int const * restrict pSSORd = SSOR.diag().data();
    
    // set the diagonal to the inverse of the A's diagonal
    for (register int irow = 0; irow < Nrows; irow++)
        pSSOR[irow] = 1. / pSSOR[irow];
    
    // move on to the first non-main diagonal
    pSSOR += Nrows;
    
    // for all non-main diagonals
    for (register int idiag = 1; idiag < Ndiag; idiag++)
    {
        // for all elements in the diagonal
        for (register int irow = 0; irow < Nrows - pSSORd[idiag]; irow++)
        {
            // divide the off-diagonal element by the correct diagonal element
            pSSOR[irow] *= pD[irow];
        }
        
        // move to next diagonal
        pSSOR +=  Nrows - pSSORd[idiag];
    }
    
    return SSOR;
}


void NoPreconditioner::setup ()
{
    // compute large radial integrals
    s_rad_.setupOneElectronIntegrals();
    s_rad_.setupTwoElectronIntegrals(par_, inp_.L + 2 * inp_.levels);
    
    std::cout << "Creating Kronecker products... ";
    
    // Kronecker producs
    # pragma omp parallel sections
    {
        # pragma omp section
        S_kron_S_   = s_rad_.S().kron(s_rad_.S());
        # pragma omp section
        S_kron_Mm1_tr_ = s_rad_.S().kron(s_rad_.Mm1_tr());
        # pragma omp section
        S_kron_Mm2_ = s_rad_.S().kron(s_rad_.Mm2());
        # pragma omp section
        Mm1_tr_kron_S_ = s_rad_.Mm1_tr().kron(s_rad_.S());
        # pragma omp section
        Mm2_kron_S_ = s_rad_.Mm2().kron(s_rad_.S());
        # pragma omp section
        half_D_minus_Mm1_tr_ = 0.5 * s_rad_.D() - s_rad_.Mm1_tr();
    }
    # pragma omp parallel sections
    {
        # pragma omp section
        half_D_minus_Mm1_tr_kron_S_ = half_D_minus_Mm1_tr_.kron(s_rad_.S());
        # pragma omp section
        S_kron_half_D_minus_Mm1_tr_ = s_rad_.S().kron(half_D_minus_Mm1_tr_);
    }
    
    std::cout << "ok\n\n";
    
    // resize arrays
    dia_blocks_.resize(l1_l2_.size());
}

void NoPreconditioner::update (double E)
{
    std::cout << "\tPrecompute diagonal blocks... " << std::flush;
    
    // setup diagonal blocks
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++)
    {
        // skip computation of unwanted blocks for this process
        if (not par_.isMyWork(ill))
            continue;
        
        int l1 = l1_l2_[ill].first;
        int l2 = l1_l2_[ill].second;
        
        // one-electron parts
        SymDiaMatrix Hdiag =
            half_D_minus_Mm1_tr_kron_S_
            + (0.5*l1*(l1+1)) * Mm2_kron_S_
            + S_kron_half_D_minus_Mm1_tr_
            + (0.5*l2*(l2+1)) * S_kron_Mm2_;
        
        // two-electron part
        for (unsigned lambda = 0; lambda <= s_rad_.maxlambda(); lambda++)
        {
            Complex f = computef(lambda,l1,l2,l1,l2,inp_.L);
            if (f != 0.)
                Hdiag += f * s_rad_.R_tr_dia(lambda);
        }
        
        // finalize the matrix
        dia_blocks_[ill] = E * S_kron_S_ - Hdiag;
    }
    
    std::cout << "ok\n";
}

void NoPreconditioner::rhs (cArrayView chi, int ie, int instate) const
{
    // shorthands
    int li = std::get<1>(inp_.instates[instate]);
    int mi = std::get<2>(inp_.instates[instate]);
    
    // shorthands
    int Nspline = s_rad_.bspline().Nspline();
    
    // j-overlaps of shape [Nangmom × Nspline]
    cArray ji_overlaps = s_rad_.overlapj (
        inp_.maxell,
        inp_.ki.slice(ie, ie + 1), // use only one ki
        weightEdgeDamp(s_rad_.bspline())
    );
    
    // j-expansions
    cArray ji_expansion = s_rad_.S().tocoo().tocsr().solve(ji_overlaps);
    
    // compute P-overlaps and P-expansion
    cArray Pi_overlaps, Pi_expansion;
    Pi_overlaps = s_rad_.overlapP(inp_.ni, li, weightEndDamp(s_rad_.bspline()));
    Pi_expansion = s_rad_.S().tocoo().tocsr().solve(Pi_overlaps);
    
    // for all segments constituting the RHS
    # pragma omp parallel for
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++)
    {
        int l1 = l1_l2_[ill].first;
        int l2 = l1_l2_[ill].second;
        
        // setup storage
        cArrayView chi_block (chi, ill * Nspline * Nspline, Nspline * Nspline);
        chi_block.fill(0);
        
        // for all allowed angular momenta (by momentum composition) of the projectile
        for (int l = abs(li - inp_.L); l <= li + inp_.L; l++)
        {
            // skip wrong parity
            if ((inp_.L + li + l) % 2 != inp_.Pi)
                continue;
            
            // (anti)symmetrization
            int Sign = ((inp_.Spin + inp_.Pi) % 2 == 0) ? 1. : -1.;
            
            // compute energy- and angular momentum-dependent prefactor
            Complex prefactor = pow(Complex(0.,1.),l) * sqrt(2*M_PI*(2*l+1)) / Complex(inp_.ki[ie]); 
            prefactor *= ClebschGordan(li,mi,l,0,inp_.L,mi);
            if (prefactor == 0.)
                continue;
            
            // pick the correct Bessel function expansion
            cArrayView Ji_expansion (ji_expansion, l * Nspline, Nspline);
            
            // compute outer products of B-spline expansions
            cArray Pj1 = outer_product(Pi_expansion, Ji_expansion);
            cArray Pj2 = outer_product(Ji_expansion, Pi_expansion);
            
            // skip angular forbidden right hand sides
            for (unsigned lambda = 0; lambda <= s_rad_.maxlambda(); lambda++)
            {
                Complex f1 = computef(lambda, l1, l2, li, l, inp_.L);
                Complex f2 = computef(lambda, l1, l2, l, li, inp_.L);
                
                if (f1 != 0.)
                {
                    chi_block += (prefactor * f1) * s_rad_.R_tr_dia(lambda).dot(Pj1);
                }
                
                if (f2 != 0.)
                {
                    if (Sign > 0)
                        chi_block += (prefactor * f2) * s_rad_.R_tr_dia(lambda).dot(Pj2);
                    else
                        chi_block -= (prefactor * f2) * s_rad_.R_tr_dia(lambda).dot(Pj2);
                }
            }
            
            if (li == l1 and l == l2)
            {
                // direct contribution
                chi_block -= prefactor * S_kron_Mm1_tr_.dot(Pj1);
            }
            
            if (li == l2 and l == l1)
            {
                // exchange contribution with the correct sign
                if (Sign > 0)
                    chi_block -= prefactor * Mm1_tr_kron_S_.dot(Pj2);
                else
                    chi_block += prefactor * Mm1_tr_kron_S_.dot(Pj2);
            }
        }
    }
}

void NoPreconditioner::multiply (const cArrayView p, cArrayView q) const
{
    // shorthands
    int Nspline = s_rad_.bspline().Nspline();
    
    // clear all output segments that are going to be referenced by this process
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++)
        if (par_.isMyWork(ill))
            cArrayView(q, ill * Nspline * Nspline, Nspline * Nspline).fill(0);
    
    // multiply "q" by the matrix of the system
    # pragma omp parallel for schedule (dynamic,1) collapse(2)
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++)
    for (unsigned illp = 0; illp < l1_l2_.size(); illp++)
    if (par_.isMyWork(ill))
    {
        // row multi-index
        int l1 = l1_l2_[ill].first;
        int l2 = l1_l2_[ill].second;
        
        // column multi-index
        int l1p = l1_l2_[illp].first;
        int l2p = l1_l2_[illp].second;
        
        // product segment contribution
        cArray q_contrib (Nspline * Nspline);
        
        // copy-from segment of "p"
        cArrayView p_block(p, illp * Nspline * Nspline, Nspline * Nspline);
        
        // multiply by hamiltonian terms
        if (ill == illp)
        {
            // reuse the diagonal block
            q_contrib += dia_blocks_[ill].dot(p_block);
        }
        else
        {
            // compute the offdiagonal block
            for (unsigned lambda = 0; lambda <= s_rad_.maxlambda(); lambda++)
            {
                Complex f = computef(lambda, l1, l2, l1p, l2p, inp_.L);
                if (f != 0.)
                    q_contrib -= f * s_rad_.R_tr_dia(lambda).dot(p_block);
            }
        }
        
        // safely update shared output array "q"
        # pragma omp critical
        cArrayView (q, ill * Nspline * Nspline, Nspline * Nspline) += q_contrib;
    }
    
    // synchronize across processes
    par_.sync(q, Nspline * Nspline, l1_l2_.size());
}

void NoPreconditioner::precondition (const cArrayView r, cArrayView z) const
{
    z = r;
}

void JacobiPreconditioner::setup ()
{
    NoPreconditioner::setup();
    
    // resize attributes
    invd_.resize(l1_l2_.size());
}

void JacobiPreconditioner::update (double E)
{
    // update parent
    NoPreconditioner::update(E);
    
    // shorthands
    int Nspline = s_rad_.bspline().Nspline();
    
    # pragma omp parallel for
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++)
    {
        // resize data array
        invd_[ill] = cArray(Nspline * Nspline, 1.);
        
        // divide by the diagonal
        invd_[ill] /= dia_blocks_[ill].main_diagonal();
    }
}

void JacobiPreconditioner::precondition (const cArrayView r, cArrayView z) const
{
    // shorthands
    int Nspline = s_rad_.bspline().Nspline();
    
    # pragma omp parallel for
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyWork(ill))
    {
        // create segment views
        cArrayView rview (r, ill * Nspline * Nspline, Nspline * Nspline);
        cArrayView zview (z, ill * Nspline * Nspline, Nspline * Nspline);
        
        // apply preconditioner
        auto apply_preconditioner = [ & ](const cArrayView r, cArrayView z) -> void
        {
            z = invd_[ill] * r;
        };
        
        // multiply by matrix block
        auto matrix_multiply = [ & ](const cArrayView p, cArrayView q) -> void
        {
            q = dia_blocks_[ill].dot(p);
        };
            
        // solve using the CG solver
        cg_callbacks
        (
            rview,                  // rhs
            zview,                  // solution
            1e-11,                  // tolerance
            0,                      // min. iterations
            Nspline * Nspline,      // max. iteration
            apply_preconditioner,
            matrix_multiply,
            false
        );
    }
    
    // synchronize across processes
    par_.sync(z, Nspline * Nspline, l1_l2_.size());
}

void SSORPreconditioner::setup ()
{
    // setup parent
    NoPreconditioner::setup();
    
    // resize array
    SSOR_.resize(l1_l2_.size());
}

void SSORPreconditioner::update (double E)
{
    // update parent
    NoPreconditioner::update(E);
    
    // compute preconditioner matrix
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyWork(ill))
        SSOR_[ill] = SSOR(dia_blocks_[ill]);
}

void SSORPreconditioner::precondition (const cArrayView r, cArrayView z) const
{
    // shorthands
    int Nspline = s_rad_.bspline().Nspline();
    
    # pragma omp parallel for
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyWork(ill))
    {
        // create segment views
        cArrayView rview (r, ill * Nspline * Nspline, Nspline * Nspline);
        cArrayView zview (z, ill * Nspline * Nspline, Nspline * Nspline);
        
        // apply preconditioner
        auto apply_preconditioner = [ & ](const cArrayView r, cArrayView z) -> void
        {
            z = SSOR_[ill].upperSolve( SSOR_[ill].dot( SSOR_[ill].lowerSolve(r), diagonal ) );
        };
            
        // multiply by matrix block
        auto matrix_multiply = [ & ](const cArrayView p, cArrayView q) -> void
        {
            q = dia_blocks_[ill].dot(p);
        };
            
        // solve using the CG solver
        cg_callbacks
        (
            rview,                  // rhs
            zview,                  // solution
            1e-11,                  // tolerance
            0,                      // min. iterations
            Nspline * Nspline,      // max. iteration
            apply_preconditioner,
            matrix_multiply,
            false     
        );
    }
    
    // synchronize across processes
    par_.sync(z, Nspline * Nspline, l1_l2_.size());
}

void ILUPreconditioner::setup ()
{
    NoPreconditioner::setup();
    
    // resize arrays
    csr_blocks_.resize(l1_l2_.size());
    lu_.resize(l1_l2_.size());
}

void ILUPreconditioner::update (double E)
{
    // update parent
    NoPreconditioner::update(E);
    
    std::cout << "\t[" << par_.iproc() << "] Update preconditioner..." << std::flush;
    
    // for all diagonal blocks
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyWork(ill))
    {
        // start timer
        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
        
        // create CSR block
        csr_blocks_[ill] = dia_blocks_[ill].tocoo().tocsr();
        
        // factorize the block
        lu_[ill] = csr_blocks_[ill].factorize();
        
        // stop timer
        std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    
        // compute time usage
        std::chrono::seconds secs = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    
        // print info
        std::cout << "\b\b\b in " << secs.count() / 60 << ":" << std::setw(2) << std::setfill('0') << secs.count() % 60
                  << " (" << lu_[ill].size() / 1048576 << " MiB)\n";
    }
}

void ILUPreconditioner::multiply (const cArrayView p, cArrayView q) const
{
    NoPreconditioner::multiply(p, q);
}

void ILUPreconditioner::rhs (cArrayView chi, int ienergy, int instate) const
{
    NoPreconditioner::rhs(chi, ienergy, instate);
}

void ILUPreconditioner::precondition (const cArrayView r, cArrayView z) const
{
    // shorthands
    int Nspline = s_rad_.bspline().Nspline();
    
    // for all diagonal blocks
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++)
    {
        if (par_.isMyWork(ill))
        {
            // create copy-to view of "z"
            cArrayView zview (z, ill * Nspline * Nspline, Nspline * Nspline);
                    
            // create copy-from view of "r"
            cArrayView rview (r, ill * Nspline * Nspline, Nspline * Nspline);
            
            // apply preconditioner
            auto apply_preconditioner = [ & ](const cArrayView r, cArrayView z) -> void
            {
                z = lu_[ill].solve(r);
            };
            
            // multiply by matrix block
            auto matrix_multiply = [ & ](const cArrayView p, cArrayView q) -> void
            {
                q = dia_blocks_[ill].dot(p);
            };
            
            // solve using the CG solver
            cg_callbacks
            (
                rview,                  // rhs
                zview,                  // solution
                1e-11,                  // tolerance
                0,                      // min. iterations
                Nspline * Nspline,      // max. iteration
                apply_preconditioner,
                matrix_multiply,
                false     
            );
        }
    }
    
    // synchronize across processes
    par_.sync(z, Nspline * Nspline, l1_l2_.size());
}

void MultiLevelPreconditioner::setup ()
{
    // setup parent
    NoPreconditioner::setup();
    
    // compute radial integrals
    p_rad_.setupOneElectronIntegrals();
    p_rad_.setupTwoElectronIntegrals(par_, s_rad_.maxlambda());
    
    // precompute kronecker products
    p_half_D_minus_Mm1_tr_kron_S_ = (0.5 * p_rad_.D() - p_rad_.Mm1_tr()).kron(p_rad_.S());
    p_S_kron_half_D_minus_Mm1_tr_ = p_rad_.S().kron(0.5 * p_rad_.D() - p_rad_.Mm1_tr());
    p_Mm2_kron_S_ = p_rad_.Mm2().kron(p_rad_.S());
    p_S_kron_Mm2_ = p_rad_.S().kron(p_rad_.Mm2());
    
    // convert overlap matrices to CSR format
    s_csrS_ = s_rad_.S().tocoo().tocsr();
    p_csrS_ = p_rad_.S().tocoo().tocsr();
    
    // LU-factorize the overlap matrices
    s_luS_ = s_csrS_.factorize();
    p_luS_ = p_csrS_.factorize();
    
    // compute the transition overlap matrix
    computeSigma_();
}

void MultiLevelPreconditioner::update (double E)
{
    // update parent
    NoPreconditioner::update(E);
    
    // resize arrays
    p_csr_.resize(l1_l2_.size());
    p_lu_.resize(l1_l2_.size());
    
    // construct hamiltonian blocks
    for (size_t ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyWork(ill))
    {
        int l1 = l1_l2_[ill].first;
        int l2 = l1_l2_[ill].second;
        
        std::cout << "\t[" << par_.iproc() << "] Update preconditioner for block ("
                  << l1 << "," << l2 << ")..." << std::flush;
        
        // start timer
        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
        
        p_csr_[ill] = (
            E * p_rad_.S().kron(p_rad_.S())
            - p_half_D_minus_Mm1_tr_kron_S_
            - p_S_kron_half_D_minus_Mm1_tr_
            - (0.5 * l1 * (l1 + 1.)) * p_Mm2_kron_S_
            - (0.5 * l2 * (l2 + 1.)) * p_S_kron_Mm2_
        ).tocoo().tocsr();
        
        p_lu_[ill] = p_csr_[ill].factorize();
        
        // stop timer
        std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
        
        // compute time usage
        std::chrono::seconds secs = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    
        // print info
        std::cout << "\b\b\b in " << secs.count() / 60 << ":" << std::setw(2) << std::setfill('0') << secs.count() % 60
                  << " (" << p_lu_[ill].size() / 1048576 << " MiB)\n";
    }
}

void MultiLevelPreconditioner::computeSigma_()
{
    // sizes
    int Nss = s_bspline_.Nspline();
    int Nsp = p_bspline_.Nspline();
    int ords = s_bspline_.order();
    int ordp = p_bspline_.order();
    size_t N = Nss * Nsp;
    
    // allocate memory
    spSigma_ = RowMatrix (Nss, Nsp, cArray(N));
    psSigma_ = RowMatrix (Nsp, Nss, cArray(N));
    
    // integrator
    GaussLegendre gs(s_bspline_), gp(p_bspline_);
    
    // evaluate B-splines on their definition intervals
    int points = (ords + ordp)/2;
    cArray Bss(Nss * (ords + 1) * points), Bsp(Nsp * (ordp + 1) * points);
    for (int iss = 0; iss < Nss; iss++) // for all B-splines
    {
        for (int iknots = 0; iknots <= ords; iknots++) // for all definition intervals
        {
            s_bspline_.B (
                iss,
                iknots,
                points, 
                gs.p_points (
                    points,
                    s_bspline_.t(iknots),
                    s_bspline_.t(iknots+1)
                ).data(),
                Bss.begin() + (iss * (ords + 1) + iknots) * points
            );
        }
    }
    for (int isp = 0; isp < Nsp; isp++) // for all B-splines
    {
        for (int iknotp = 0; iknotp <= ordp; iknotp++) // for all definition intervals
        {
            p_bspline_.B (
                isp,
                iknotp,
                points,
                gp.p_points (
                    points,
                    p_bspline_.t(iknotp),
                    p_bspline_.t(iknotp+1)
                ).data(),
                Bsp.begin() + (isp * (ordp + 1) + iknotp) * points
            );
        }
    }
    
    // create temporary raw Sigma matrices
    ColMatrix spSigma(Nss, Nsp), SspSigma(Nss, Nsp);
    ColMatrix psSigma(Nsp, Nss), SpsSigma(Nsp, Nss);
    
    // for all B-spline pairs
    for (int iss = 0; iss < Nss; iss++)
    for (int isp = 0; isp < Nsp; isp++)
    {
        Complex elem = 0.;
        
        // get overlap interval boundaries (shifted by the respective order)
        // WARNING : Full multiplicity of zero knot (= order) is assumed!
        // FIXME : The intervals still could be optimized!
        int iknot_start = std::max(isp - ordp, iss - ords);
        int iknot_end   = std::min(isp + 1, iss + 1);
        
        // add contributions from all common intervals
        for (int iknot = iknot_start; iknot < iknot_end; iknot++)
        {
            // get quadrature weights
            // NOTE : Identical grids for (s) and (p) assumed (apart from the multiplicity).
            cArray ws = gs.p_weights(points, p_bspline_.t(iknot + ordp), p_bspline_.t(iknot + ordp + 1));
            
            // add contributions from all common quadrature points
            for (int ipt = 0; ipt < points; ipt++)
            {
                elem += ws[ipt]
                          * Bsp[(isp * (ordp + 1) + iknot + ordp) * points + ipt]
                          * Bss[(iss * (ords + 1) + iknot + ords) * points + ipt];
            }
        }
        
        // store spSigma and psSigma
        spSigma.data()[isp * Nsp + iss] = elem;
        psSigma.data()[iss * Nss + isp] = elem;
    }
    
    // compute inverse overlap matrices
    CsrMatrix csrSs = s_rad_.S().tocoo().tocsr();
    CsrMatrix csrSp = p_rad_.S().tocoo().tocsr();
    CsrMatrix::LUft luSs = csrSs.factorize();
    CsrMatrix::LUft luSp = csrSp.factorize();
    
    // multiply by S⁻¹(s)
    for (int icol = 0; icol < spSigma.cols(); icol++)
    {
        SspSigma.col(icol) = luSs.solve(spSigma.col(icol));
        SpsSigma.col(icol) = luSp.solve(psSigma.col(icol));
    }
    
    // save the matrix
    spSigma_ = RowMatrix(SspSigma);
    psSigma_ = RowMatrix(SpsSigma);
    
    /// DEBUG
    std::ofstream out;
    out.open("spSigma.txt"); spSigma_.write(out); out.close();
    out.open("psSigma.txt"); psSigma_.write(out); out.close();
}

void MultiLevelPreconditioner::rhs(cArrayView chi, int ienergy, int instate) const
{
    NoPreconditioner::rhs(chi, ienergy, instate);
}

void MultiLevelPreconditioner::multiply(const cArrayView p, cArrayView q) const
{
    NoPreconditioner::multiply(p, q);
}

void MultiLevelPreconditioner::precondition (const cArrayView rs, cArrayView zs) const
{
    // shorthands
    int Nss = s_bspline_.Nspline();
    int Nsp = p_bspline_.Nspline();
    size_t  N = Nss * Nsp;
    
    // for all work items
    # pragma omp parallel for schedule (dynamic, 1)
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++)
    {
        if (par_.isMyWork(ill))
        {
            // get segment view
            cArrayView rsview (rs, ill * N, (ill + 1) * N);
            cArrayView zsview (zs, ill * N, (ill + 1) * N);
            
            // reinterpret rs as a column matrix
            ColMatrix Rs (Nss, Nss, rsview);
            
            // convert R^s to preconditioner basis
            ColMatrix Rp = (psSigma_ * (psSigma_ * Rs).T()).T();
            
            // solve the low-order system of equations
            ColMatrix Zp (Nsp, Nsp);
            for (int icol = 0; icol < Rp.cols(); icol++)
                Zp.col(icol) = p_lu_[ill].solve(Rp.col(icol));
            
            // convert Z^p to solver basis
            ColMatrix Zs = (spSigma_ * (spSigma_ * Zp).T()).T();
            
            // use the segment
            zsview = Zs.data();
        }
    }
    
    // synchronize across processes
    par_.sync(zs, Nss * Nss, l1_l2_.size());
}

/*
        // block incomplete D-ILU
        std::vector<std::vector<SymDiaMatrix>> bcks(l1_l2_.size());
        std::vector<std::vector<CsrMatrix>> bcks_csr(l1_l2_.size());
        std::vector<std::vector<CsrMatrix::LUft>> bcks_lufts(l1_l2_.size());
        
        // sparse approximate inverse
        std::vector<SymDiaMatrix> spai(l1_l2_.size());
        
        // incomplete Cholesky factorization of the diagonal blocks
        //   L + Lt - I
        std::vector<SymDiaMatrix> icholL(l1_l2_.size());
        //   D⁻¹
        std::vector<cArray> icholD(l1_l2_.size());
        
        // diagonal incomplete Cholesky (DIC) preconditioner
        std::vector<SymDiaMatrix> DIC(l1_l2_.size());
 
            // diagonal incomplete Cholesky factorization
            if (preconditioner == dic_prec)
            {
                DIC[ill] = DIC_preconditioner(dia_blocks[ill]);
                cArray(DIC[ill].main_diagonal()).hdfsave(format("DIC-%d.hdf",ill));
            }
 
            if (preconditioner == bilu_prec)
            {
                // log output
                std::cout << "\n\t\t-> [" << par.iproc() << "] block D-ILU factorization "
                          << ill << " of (" << l1 << "," << l2 << ") block started\n";
                
                // allocate space
                std::vector<SymDiaMatrix> ibcks(Nspline);
                bcks[ill].resize(Nspline);
                bcks_csr[ill].resize(Nspline);
                bcks_lufts[ill].resize(Nspline);
                
                std::cout << "\t\t   [" << par.iproc() << "] ";
                std::cout << "computing blocks and SPAIs\n";
                
                // for all diagonal blocks
                for (int iblock = 0; iblock < Nspline; iblock++)
                {
                    // - compute the block
                    bcks[ill][iblock] = E * rad.S().main_diagonal()[iblock] * rad.S()
                            - 0.5 * rad.D().main_diagonal()[iblock] * rad.S()
                            - 0.5 * rad.S().main_diagonal()[iblock] * rad.D()
                            - 0.5 * l1 * (l1 + 1.) * rad.Mm2().main_diagonal()[iblock] * rad.S()
                            - 0.5 * l2 * (l2 + 1.) * rad.S().main_diagonal()[iblock] * rad.Mm2()
                            + rad.Mm1_tr().main_diagonal()[iblock] * rad.S()
                            + rad.S().main_diagonal()[iblock] * rad.Mm1_tr();
                    
                    // - invert using the SPAI
                    ibcks[iblock] = SymDiaMatrix (
                        Nspline,
                        iArray(1), // = (int[]){ 0 }
                        cArray(Nspline,1.)/bcks[ill][iblock].main_diagonal()
                    );
                }
                
                std::cout << "\t\t   [" << par.iproc() << "] ";
                std::cout << "preconditioning the diagonal\n";
                
                // for all diagonals of overlap matrix
                for (int idiag = 1; idiag <= order; idiag++)
                {
                    // for all elements of this diagonal
                    for (int irow = 0; irow < Nspline - idiag; irow++)
                    {
                        // - construct corresponding block of Kronecker product
                        SymDiaMatrix bck = E * rad.S().dptr(idiag)[irow] * rad.S()
                                - 0.5 * rad.D().dptr(idiag)[irow] * rad.S()
                                - 0.5 * rad.S().dptr(idiag)[irow] * rad.D()
                                - 0.5 * l1 * (l1 + 1.) * rad.Mm2().dptr(idiag)[irow] * rad.S()
                                - 0.5 * l2 * (l2 + 1.) * rad.S().dptr(idiag)[irow] * rad.Mm2()
                                + rad.Mm1_tr().dptr(idiag)[irow] * rad.S()
                                + rad.S().dptr(idiag)[irow] * rad.Mm1_tr();
                        
                        
                        // - update pivot
                        bcks[ill][irow+idiag] -= bck * (ibcks[irow] * bck);
                        
                    }
                }
                
                std::cout << "\t\t   [" << par.iproc() << "] ";
                std::cout << "factorization of the pivots\n";
                
                // convert and factorize the pivots
                size_t size = 0;
                for (int iblock = 0; iblock < Nspline; iblock++)
                {
                    
                    bcks_csr[ill][iblock] = bcks[ill][iblock].tocoo().tocsr();
                    bcks_lufts[ill][iblock] = bcks_csr[ill][iblock].factorize(droptol);
                    size += bcks_lufts[ill][iblock].size();
                }
                
                // log output
                std::chrono::duration<int> sec = std::chrono::duration_cast<std::chrono::duration<int>>(std::chrono::steady_clock::now()-start);
                std::cout << "\t\t   [" << par.iproc() << "] ";
                std::cout << "droptol " << droptol << ", ";
                std::cout << "time " << sec.count() / 60 << ":" << std::setfill('0') << std::setw(2) << sec.count() % 60 << ", ";
                std::cout << "average " << (sec.count() / Nspline) / 60 << ":" << std::setfill('0') << std::setw(2) << (sec.count() / Nspline) % 60 << ", ";
                std::cout << "mem " << size/1048576 << " MiB";
            }
            
            // drop-tolerance-incomplete single-electron LU factorization
            if (preconditioner == silu_prec)
            {
                // log output
                std::cout << "\n\t\t-> [" << par.iproc() << "] s-iLU factorization "
                          << ill << " of (" << l1 << "," << l2 << ") block started\n";
                
                // allocate space
                scsr_blocks[ill].resize(Nspline);
                siLU[ill].resize(Nspline);
                
                unsigned size = 0;
                
                // for all single-electron blocks
                for (int iblock = 0; iblock < Nspline; iblock++)
                {
                    // setup the single-electron block
                    scsr_blocks[ill][iblock] = (
                        E * rad.S().main_diagonal()[iblock] * rad.S()
                      - 0.5 * rad.D().main_diagonal()[iblock] * rad.S()
                      - 0.5 * rad.S().main_diagonal()[iblock] * rad.D()
                      - 0.5 * l1 * (l1 + 1.) * rad.Mm2().main_diagonal()[iblock] * rad.S()
                      - 0.5 * l2 * (l2 + 1.) * rad.S().main_diagonal()[iblock] * rad.Mm2()
                      + rad.Mm1_tr().main_diagonal()[iblock] * rad.S()
                      + rad.S().main_diagonal()[iblock] * rad.Mm1_tr()
                    ).tocoo().tocsr();
                    
//                     scsr_blocks[ill][iblock].plot(format("scsr-%d-%.03d.png", ill, iblock), 1.);
//                     scsr_blocks[ill][iblock].hdfsave(format("scsr-%d-%.03d.hdf", ill, iblock));
                    
                    // factorize the block
                    siLU[ill][iblock] = scsr_blocks[ill][iblock].factorize(droptol);
                    size += siLU[ill][iblock].size();
                }
                
                // log output
                std::chrono::duration<int> sec = std::chrono::duration_cast<std::chrono::duration<int>>(std::chrono::steady_clock::now()-start);
                std::cout << "\t\t   [" << par.iproc() << "] ";
                std::cout << "droptol " << droptol << ", ";
                std::cout << "time " << sec.count() / 60 << ":" << std::setfill('0') << std::setw(2) << sec.count() % 60 << ", ";
                std::cout << "average " << (sec.count() / Nspline) / 60 << ":" << std::setfill('0') << std::setw(2) << (sec.count() / Nspline) % 60 << ", ";
                std::cout << "mem " << size/1048576 << " MiB";
            }
            
            // sparse approximate inverse preconditioner
            if (preconditioner == spai_prec)
            {
                // use diagonal preconditioner
                spai[ill] = SymDiaMatrix (
                    Nspline * Nspline,
                    iArray(1),
                    cArray(Nspline * Nspline, 1.) / dia_blocks[ill].main_diagonal()
                );
            }
*/


/*
                // apply a block inversion preconditioner (parallel for ILU)
                # pragma omp parallel for if (preconditioner == ilu_prec)
                for (unsigned ill = 0; ill < l1_l2_.size(); ill++)
                {
                    // skip computation of unwanted blocks for this process
                    if (not par.isMyWork(ill))
                        continue;
                    
                    // create copy-to view of "z"
                    cArrayView zview(z, ill * Nspline * Nspline, Nspline * Nspline);
                    
                    // create copy-from view of "r"
                    cArrayView rview(r, ill * Nspline * Nspline, Nspline * Nspline);
                    
                    // preconditioner of the nested CG
                    auto apply_inner_preconditioner = [ & ](cArray const & r, cArray & z) -> void
                    {
                        // Incomplete LU factorization
                        if (preconditioner == ilu_prec)
                            z = iLU[ill].solve(r);
                        
                        // single-electron Incomplete LU factorization
                        // TODO Needs some modification to work!
                        if (preconditioner == silu_prec)
                        {
                            // for all single-electron blocks
                            # pragma omp parallel for
                            for (int iblock = 0; iblock < Nspline; iblock++)
                            {
                                // precondition by inverting a single diagonal block
                                cArrayView rview (r, iblock * Nspline, Nspline);
                                cArrayView zview (z, iblock * Nspline, Nspline);
                                zview = siLU[ill][iblock].solve(rview);
                            }
                        }
                        
                        // block D-ILU
                        if (preconditioner == bilu_prec)
                        {
                            // for all single-electron blocks
                            # pragma omp parallel for
                            for (int iblock = 0; iblock < Nspline; iblock++)
                            {
                                // precondition by inverting a single diagonal block
                                cArrayView rview (r, iblock * Nspline, Nspline);
                                cArrayView zview (z, iblock * Nspline, Nspline);
                                zview = bcks_lufts[ill][iblock].solve(rview);
                            }
                        }
                        
                        // Diagonal Incomplete Cholesky factorization
                        // TODO Needs to implement pivoting to work!
                        if (preconditioner == dic_prec)
                        {
                            z = DIC[ill].upperSolve( DIC[ill].dot( DIC[ill].lowerSolve(r), diagonal ) );
                            std::cout << z.norm()/r.norm() << "\n";
                            std::cout << cArray(DIC[ill].main_diagonal()).norm() << "\n";
                        }
                        
                        // Sparse Approximate Inverse preconditioner
                        if (preconditioner == spai_prec)
                            z = spai[ill].dot(r);
                        
                        // Symmetric Successive Over-Relaxation
                        // NOTE seems slower than Jacobi
                        if (preconditioner == ssor_prec)
                            z = SSOR[ill].upperSolve( SSOR[ill].dot( SSOR[ill].lowerSolve(r), diagonal ) );
                        
                        // Jacobi preconditioning
                        // NOTE seems the fastest
                        if (preconditioner == jacobi_prec)
                            z = SSOR[ill].dot(r,diagonal);
                        
                        // no preconditioning
                        // NOTE seems slower than Jacobi
                        if (preconditioner == no_prec)
                            z = r;
                    };
                    
                    // multiply by matrix block
                    auto inner_matrix_multiply = [ & ](cArray const & p, cArray & q) -> void
                    {
                        q = dia_blocks[ill].dot(p);
                    };
                    
                    // solve using the CG solver
                    cg_callbacks
                    (
                        rview,      // rhs
                        zview,      // solution
                        1e-11,      // tolerance
                        0,          // min. iterations
                        Nspline*Nspline,    // max. iteration
                        apply_inner_preconditioner,
                        inner_matrix_multiply,
                        true       // verbose
                    );
                }
                
*/
