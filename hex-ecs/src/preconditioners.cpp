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

#include <CL/cl.h>

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
    #define THRESHOLD 1e-5
    
    write_array(A.data(), "A.dat");
    
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
    
    // for all elements of the main diagonal
    for (register int irow = 0; irow < Nrows; irow++)
    {
        // pointer to the A's raw concatenated diagonal data (starting at first non-main diagonal)
        Complex const * restrict pA = A.data().data() + Nrows;
        
        // for all diagonals (except the main diagonal) constributing to DIC
        for (register int idiag = 1; idiag < Ndiag and pAd[idiag] <= irow; idiag++)
        {
            // update pivot
            pD[irow] -= pA[irow - pAd[idiag]] * pA[irow - pAd[idiag]] / pD[irow - pAd[idiag]];
            
            if (std::abs(pD[irow]) < THRESHOLD)
                throw exception ("[DIC] Pivot too small for irow = %d", irow);
            
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
    int const * const restrict pDICd = DIC.diag().data();
    
    // set the diagonal to the inverse of the DIC diagonal
    for (register int irow = 0; irow < Nrows; irow++)
        pDIC[irow] = 1. / pD[irow];
    
    // move on to the first non-main diagonal
    pDIC += Nrows;
    
    // for all non-main upper diagonals
    for (register int idiag = 1; idiag < Ndiag; idiag++)
    {
        // for all elements in the diagonal
        for (register int irow = 0; irow < Nrows - pDICd[idiag]; irow++)
        {
            // divide the off-diagonal element by the correct diagonal element
            pDIC[irow] *= pD[irow];
        }
        
        // move on to the next diagonal
        pDIC +=  Nrows - pDICd[idiag];
    }
    
    write_array(DIC.data(), "DIC.dat");
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
    # pragma omp parallel for
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
    cArray ji_expansion = s_rad_.S().tocoo().tocsr().solve(ji_overlaps, ji_overlaps.size() / Nspline);
    
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
    
    // multiply "p" by the matrix of the system
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
        cArrayView p_block (p, illp * Nspline * Nspline, Nspline * Nspline);
        
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
    par_.sync (q, Nspline * Nspline, l1_l2_.size());
}

void NoPreconditioner::precondition (const cArrayView r, cArrayView z) const
{
    z = r;
}

void CGPreconditioner::precondition (const cArrayView r, cArrayView z) const
{
    // shorthands
    int Nspline = s_rad_.bspline().Nspline();
    
    # pragma omp parallel for schedule (dynamic, 1)
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyWork(ill))
    {
        // create segment views
        cArrayView rview (r, ill * Nspline * Nspline, Nspline * Nspline);
        cArrayView zview (z, ill * Nspline * Nspline, Nspline * Nspline);
        
        // wrappers around the callbacks
        auto inner_mmul = [&](const cArrayView a, cArrayView b) { this->CG_mmul(ill, a, b); };
        auto inner_prec = [&](const cArrayView a, cArrayView b) { this->CG_prec(ill, a, b); };
        
        // solve using the CG solver
        cg_callbacks
        (
            rview,                  // rhs
            zview,                  // solution
            1e-11,                  // tolerance
            0,                      // min. iterations
            Nspline * Nspline,      // max. iteration
            inner_prec,             // preconditioner
            inner_mmul,             // matrix multiplication
            false
        );
    }
    
    // synchronize across processes
    par_.sync (z, Nspline * Nspline, l1_l2_.size());
}

void JacobiCGPreconditioner::setup ()
{
    NoPreconditioner::setup();
    
    // resize attributes
    invd_.resize(l1_l2_.size());
}

void JacobiCGPreconditioner::update (double E)
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

void SSORCGPreconditioner::setup ()
{
    // setup parent
    NoPreconditioner::setup();
    
    // resize array
    SSOR_.resize(l1_l2_.size());
}

void SSORCGPreconditioner::update (double E)
{
    // update parent
    NoPreconditioner::update(E);
    
    // compute preconditioner matrix
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyWork(ill))
        SSOR_[ill] = SSOR(dia_blocks_[ill]);
}

void ILUCGPreconditioner::setup ()
{
    NoPreconditioner::setup();
    
    // resize arrays
    csr_blocks_.resize(l1_l2_.size());
    lu_.resize(l1_l2_.size());
}

void ILUCGPreconditioner::update (double E)
{
    // update parent
    NoPreconditioner::update(E);
    
    std::cout << "\t[" << par_.iproc() << "] Update preconditioner\n";
    
    // for all diagonal blocks
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyWork(ill))
    {
        std::cout << "\t\t- block #" << ill << " (" << l1_l2_[ill].first << "," << l1_l2_[ill].second << ")..." << std::flush;
        
        // start timer
        Timer::timer().start();
        
        // create CSR block
        csr_blocks_[ill] = dia_blocks_[ill].tocoo().tocsr();
        
        // factorize the block
        lu_[ill] = csr_blocks_[ill].factorize(droptol_);
        
        // time usage
        int secs = Timer::timer().stop();
    
        // print info
        std::cout << "\b\b\b in " << secs / 60 << ":" << std::setw(2) << std::setfill('0') << secs % 60
                  << " (" << lu_[ill].size() / 1048576 << " MiB)\n";
    }
}

void DICCGPreconditioner::setup()
{
    CGPreconditioner::setup();
    DIC_.resize(l1_l2_.size());
}

void DICCGPreconditioner::update(double E)
{
    CGPreconditioner::update(E);
    
    // diagonal incomplete Cholesky factorization
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++)
    {
        DIC_[ill] = DIC(dia_blocks_[ill]);
        write_array(DIC_[ill].main_diagonal(), format("DIC-main_diagonal-%d.dat", ill));
        write_array(DIC_[ill].data(), format("DIC-al_data-%d.dat", ill));
    }
}


void SPAICGPreconditioner::setup()
{
    // setup parent
    NoPreconditioner::setup();
    
    // TODO
}

void SPAICGPreconditioner::update (double E)
{
    // update parent
    NoPreconditioner::update(E);
    
    // TODO
}


void TwoLevelPreconditioner::setup ()
{
    std::cout << "Setting up ML preconditioner.\n";
    std::cout << "\t- rknots = " << p_bspline_.rknots() << "\n";
    std::cout << "\t- cknots = " << p_bspline_.cknots() << "\n";
    std::cout << "\t- B-spline count = " << p_bspline_.Nspline() << "\n\n";
    
    // setup parent
    SSORCGPreconditioner::setup();
    
    // compute radial integrals
    p_rad_.setupOneElectronIntegrals();
    p_rad_.setupTwoElectronIntegrals(par_, s_rad_.maxlambda());
    
    // precompute kronecker products
    p_S_kron_S_ = p_rad_.S().kron(p_rad_.S());
    p_half_D_minus_Mm1_tr_kron_S_ = (0.5 * p_rad_.D() - p_rad_.Mm1_tr()).kron(p_rad_.S());
    p_S_kron_half_D_minus_Mm1_tr_ = p_rad_.S().kron(0.5 * p_rad_.D() - p_rad_.Mm1_tr());
    p_Mm2_kron_S_ = p_rad_.Mm2().kron(p_rad_.S());
    p_S_kron_Mm2_ = p_rad_.S().kron(p_rad_.Mm2());
    
    // compute the transition overlap matrix
    computeSigma_();
}

void TwoLevelPreconditioner::update (double E)
{
    // update parent
    SSORCGPreconditioner::update(E);
    
    // resize arrays
    p_csr_.resize(l1_l2_.size());
    p_lu_.resize(l1_l2_.size());
    
    // use total angular momentum from parent
    int L = NoPreconditioner::inp_.L;
    
    // construct hamiltonian blocks
    for (size_t ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyWork(ill))
    {
        int l1 = l1_l2_[ill].first;
        int l2 = l1_l2_[ill].second;
        
        std::cout << "\t[" << par_.iproc() << "] Update preconditioner for block ("
                  << l1 << "," << l2 << ")..." << std::flush;
        
        // start timer
        Timer::timer().start();
        
        // construct DIA block
        SymDiaMatrix p_dia = E * p_S_kron_S_
            - p_half_D_minus_Mm1_tr_kron_S_
            - p_S_kron_half_D_minus_Mm1_tr_
            - (0.5 * l1 * (l1 + 1.)) * p_Mm2_kron_S_
            - (0.5 * l2 * (l2 + 1.)) * p_S_kron_Mm2_;
        for (unsigned lambda = 0; lambda <= p_rad_.maxlambda(); lambda++)
        {
            Complex f = computef(lambda,l1,l2,l1,l2,L);
            if (f != 0.)
                p_dia -= f * p_rad_.R_tr_dia(lambda);
        }
        
        // convert DIA block to CSR and factorize
        p_csr_[ill] = p_dia.tocoo().tocsr();
        p_lu_[ill] = p_csr_[ill].factorize();
        
        // stop timer
        int secs = Timer::timer().stop();
        
        // print info
        std::cout << "\b\b\b in " << secs / 60 << ":" << std::setw(2) << std::setfill('0') << secs % 60
                  << " (" << p_lu_[ill].size() / 1048576 << " MiB)\n";
    }
}

Complex TwoLevelPreconditioner::computeSigma_iknot_ (int qord, int is, int iknots, int ip, int iknotp) const
{
    // get interval bounds
    Complex t1 = p_bspline_.t(iknotp);
    Complex t2 = p_bspline_.t(iknotp + 1);
    
    // get evaluation points and weights
    cArray x = g_.p_points(qord, t1, t2);
    cArray w = g_.p_weights(qord, t1, t2);
    
    // evaluate the B-splines
    cArray Bp(qord), Bs(qord);
    p_bspline_.B(ip, iknotp, qord, x.data(), Bp.data());
    s_bspline_.B(is, iknots, qord, x.data(), Bs.data());
    
    // integrate
    Complex sum = 0;
    for (int i = 0; i < qord; i++)
        sum += w[i] * Bp[i] * Bs[i];
    return sum;
}

void TwoLevelPreconditioner::computeSigma_()
{
    // sizes
    int Nss = s_bspline_.Nspline();
//     int Nks = s_bspline_.Nknot();
    int Nsp = p_bspline_.Nspline();
    int Nkp = p_bspline_.Nknot();
    int ords = s_bspline_.order();
    int ordp = p_bspline_.order();
    size_t N = Nss * Nsp;
    
    // quadrature order (use at least 2nd order rule)
    int qord = std::max (2, (ords + ordp + 1) / 2);
    
    // allocate memory
    spSigma = ColMatrix<Complex>(Nss, Nsp, cArray(N));
    SspSigma = ColMatrix<Complex>(Nss, Nsp, cArray(N));
    psSigma = ColMatrix<Complex>(Nsp, Nss, cArray(N));
    SpsSigma = ColMatrix<Complex>(Nsp, Nss, cArray(N));
    
    // for all B-splines
    for (int iss = 0; iss < Nss; iss++)
    for (int isp = 0; isp < Nsp; isp++)
    {
        // element of the Sigma matrix
        Complex elem = 0.;
        
        // for all preconditioner B-spline interknot intervals
        for (int iknotp = isp; iknotp <= isp + ordp and iknotp < Nkp - 1; iknotp++)
        {
            // find corresponding solution basis knot
            // WARNING : Assuming full multiplicity in zero and zero multiplicity elsewhere!
            int iknots = iknotp + ords - 1;
            
            // skip non-existent intervals
            if (iknots < 0)
                continue;
            
            // check that the knot indices correspond
            assert(p_bspline_.t(iknotp) == s_bspline_.t(iknots));
            
            // evaluate B-spline integral on this interval
            elem += computeSigma_iknot_(qord, iss, iknots, isp, iknotp);
        }
        
        // update the elements
        spSigma(iss, isp) = psSigma(isp, iss) = elem;
    }
    
    // compute inverse overlap matrices
    CsrMatrix csrSs = s_rad_.S().tocoo().tocsr();
    CsrMatrix csrSp = p_rad_.S().tocoo().tocsr();
    CsrMatrix::LUft luSs = csrSs.factorize();
    CsrMatrix::LUft luSp = csrSp.factorize();
    
    // multiply by S⁻¹(s)
    for (int icol = 0; icol < spSigma.cols(); icol++)
        SspSigma.col(icol) = luSs.solve(spSigma.col(icol));
    for (int icol = 0; icol < psSigma.cols(); icol++)
        SpsSigma.col(icol) = luSp.solve(psSigma.col(icol));
    
    // save the matrix
    spSigma_ = RowMatrix<Complex>(SspSigma);
    psSigma_ = RowMatrix<Complex>(SpsSigma);
}

void TwoLevelPreconditioner::rhs (const cArrayView chi, int ienergy, int instate) const
{
    NoPreconditioner::rhs(chi, ienergy, instate);
}

void TwoLevelPreconditioner::multiply (const cArrayView p, cArrayView q) const
{
    NoPreconditioner::multiply(p, q);
}

void TwoLevelPreconditioner::CG_prec (int iblock, const cArrayView rs, cArrayView zs) const
{
    // shorthands
    int Nss = s_bspline_.Nspline();
    int Nsp = p_bspline_.Nspline();
    size_t  N = Nss * Nss;
    
    // reinterpret rs as a column matrix
    ColMatrix<Complex> Rs (Nss, Nss, rs);
    
    // convert R^s to preconditioner basis
    ColMatrix<Complex> Rp = (psSigma_ * (psSigma_ * Rs).T()).T();
    
    // solve the low-order system of equations
    ColMatrix<Complex> Zp (Nsp, Nsp, p_lu_[iblock].solve(Rp.data()));
    
    // convert Z^p to solver basis
    ColMatrix<Complex> Zs = (spSigma_ * (spSigma_ * Zp).T()).T();
    
    // use the result
    zs = Zs.data();
}


MultiresPreconditioner::MultiresPreconditioner (
    Parallel const & par, InputFile const & inp, std::vector<std::pair<int,int>> const & ll, Bspline const & bspline
){
    // first order will be solved by ILU (-> ILUCGPreconditioner)
    p_.push_back
    (
        new ILUCGPreconditioner
        (
            par, inp, ll,
            Bspline
            (
                1,
                sorted_unique(bspline.rknots(), 1),
                bspline.ECStheta(),
                sorted_unique(bspline.cknots(), 1)
            )
        )
    );
    
    // all other levels need just common radial data (-> NoPreconditioner)
    for (int ord = 2; ord <= bspline.order(); ord++)
    {
        // construct preconditioner
        p_.push_back
        (
            new NoPreconditioner
            (
                par, inp, ll,
                Bspline
                (
                    ord,
                    sorted_unique(bspline.rknots(), ord),
                    bspline.ECStheta(),
                    sorted_unique(bspline.cknots(), ord)
                )
            )
        );
    }
}

void MultiresPreconditioner::precondition (const cArrayView r, cArrayView z) const
{
    // HOA WHOW... TODO
    
}
