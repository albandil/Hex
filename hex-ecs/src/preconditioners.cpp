/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2014                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>
#include <set>

#include "arrays.h"
#include "gauss.h"
#include "itersolve.h"
#include "misc.h"
#include "preconditioners.h"
#include "radial.h"

#ifndef NO_OPENCL
    #include "opencl.h"
#endif

cArray IC (cArrayView const & A, lArrayView const & I, lArrayView const & P)
{
    // this will be returned
    cArray LD(A.size());
    
    // check lengths
    assert(A.size() == (std::size_t)P.back());
    assert(I.size() == (std::size_t)P.back());
    
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

#ifndef NO_LAPACK
extern "C" void zgelsd_
(
    int * M,
    int * N,
    int * NRHS,
    Complex * A,
    int * LDA,
    Complex * B,
    int * LDB,
    double * S,
    double * RCOND,
    int * RANK,
    Complex * WORK,
    int * LWORK,
    double * RWORK,
    int * IWORK,
    int * INFO
);

void LapackLeastSquares (ColMatrix<Complex> const & A, cArrayView b)
{
    if ((std::size_t)A.rows() != b.size())
        throw exception ("[LapackLeastSquares] Matrix row count (%d) is different from RHS size (%d).", A.rows(), b.size());
    
    // set data to pass to Lapack
    int m = A.rows();
    int n = A.cols();
    int nrhs = 1;
    rArray svd (n);
    double rcond = -1;
    int rank = n;
    int info = 0;
    cArray work(1);
    rArray rwork(1);
    iArray iwork(1);
    
    // query for optimal dimensions of work arrays
    int lwork = -1;
    zgelsd_(&m, &n, &nrhs, const_cast<Complex*>(A.data().begin()), &m, b.data(), &m, svd.data(), &rcond, &rank, work.data(), &lwork, rwork.data(), iwork.data(), &info);
    
    if (info != 0)
        throw exception ("[SPAI] The Lapack routine ZGELSS returned INFO = %d.", info);
    
    // resize work arrays
    lwork = work.resize(work[0].real());
    rwork.resize(rwork[0]);
    iwork.resize(iwork[0]);
    
    // compute the least squares
    zgelsd_(&m, &n, &nrhs, const_cast<Complex*>(A.data().begin()), &m, b.data(), &m, svd.data(), &rcond, &rank, work.data(), &lwork, rwork.data(), iwork.data(), &info);
    
    if (info != 0)
        throw exception ("[SPAI] The Lapack routine ZGELSS returned INFO = %d.", info);
}

CooMatrix SPAI (SymDiaMatrix const & A, iArrayView diagonals)
{
    // dimensions
    int Asize = A.size();
    
    // preconditioner
    CooMatrix M (Asize, Asize);
    
    // for all columns of M
    for (int icol = 0; icol < Asize; icol++)
    {
        std::cout << icol << "\n";
        
        // assemble non-zero pattern of the column m[icol] (will be same as nonzero pattern of A[:,icol])
        //
        std::set<int> mnz;
        for (int id : diagonals)
        {
            if (id <= icol) // upper (or main) diagonal contributes to nz
                mnz.insert(icol - id);
            if (icol + id < Asize) // lower diagonal contributes to nz
                mnz.insert(icol + id);
        }
        
        // assemble non-zero pattern of the right hand side e_icol (defined by nonzero rows of A[:,{mnz}] submatrix)
        //
        std::set<int> enz;
        for (int acol : mnz)
        {
            // loop over all A'd diagonals
            for (int id : A.diag())
            {
                if (id <= acol) // upper (or main) diagonal contributes to nz
                    enz.insert(acol - id);
                if (acol + id < Asize) // lower diagonal contributes to nz
                    enz.insert(acol + id);
            }
        }
        
        // construct the dense sub-matrix of A
        //
        ColMatrix<Complex> submA (enz.size(), mnz.size()); // FORTRAN column-wise ordering necessary for Lapack.
        Complex * psubmA = submA.data().begin();
        for (int orgA_col : mnz) // original A's columns
        {
            for (int orgA_row : enz) // original A's rows
            {
                // which diagonal is this?
                int d = std::abs(orgA_col - orgA_row);
                
                // does it contain any data?
                iArray::const_iterator ptrd = std::find(A.diag().begin(),A.diag().end(),d);
                if (ptrd != A.diag().end())
                    *psubmA = A.dptr(ptrd-A.diag().begin())[std::min(orgA_row,orgA_col)];
                
                psubmA++;
            }
        }
        
        // construct the right-hand side
        //
        cArray e (enz.size());
        Complex * pe = e.data();
        for (int erow : enz)
        {
            // icol-th vector "e" has 'one' as the icol-th element (if not excluded by pattern)
            if (erow == icol)
                *pe = 1.;
            pe++;
        }
        
        // solve the overdetermined system Am = e using Lapack routine
        //
        LapackLeastSquares(submA, e); // overwrites 'e'
        
        // store computed data
        //
        Complex const * px = e.data();
        for (int mrow : mnz)
            M.add(mrow, icol, *(px++));
    }
    
    return M;
}
#endif

const std::string NoPreconditioner::name = "none";
const std::string NoPreconditioner::description = "Preconditioning by the identity matrix.";

void NoPreconditioner::setup ()
{
    // TODO : Determine which lambdas are needed by this process.
    // NOTE : At the moment each process holds in memory radial integrals
    //        for all lambdas, which needlessly raises memory requirements.
    Array<bool> lambdas (inp_.L + 2 * inp_.levels + 1, true);
    
    // compute large radial integrals
    s_rad_.setupOneElectronIntegrals();
    s_rad_.setupTwoElectronIntegrals(par_, cmd_, lambdas);
    
    std::cout << "Creating Kronecker products... ";
    
    // Kronecker products
    # pragma omp parallel sections if (!cmd_.outofcore)
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
    # pragma omp parallel sections if (!cmd_.outofcore)
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
    E_ = E;
    
    // setup diagonal blocks
    # pragma omp parallel for if (cmd_.parallel_block)
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyWork(ill))
    {
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
            Complex f = special::computef(lambda,l1,l2,l1,l2,inp_.L);
            
            // check that the "f" coefficient is valid (no factorial overflow etc.)
            if (not Complex_finite(f))
                throw exception ("Overflow in computation of f[%d](%d,%d,%d,%d).", inp_.L, l1, l2, l1, l2);
            
            // check that the "f" coefficient is nonzero
            if (f == 0.)
                continue;
            
            // load two-electron integrals if necessary
//             SymDiaMatrix & ref_R_tr_dia = const_cast<SymDiaMatrix&>(s_rad_.R_tr_dia(lambda));
//             if (cmd_.outofcore)
//             {
//                 # pragma omp critical
//                 ref_R_tr_dia.hdfload();
//             }
            
            // add two-electron contributions
//             Hdiag += f * ref_R_tr_dia;
            Hdiag += f * s_rad_.R_tr_dia(lambda);
            
            // unload the two-electron integrals if necessary
//             if (cmd_.outofcore)
//             {
//                 #pragma omp critical
//                 ref_R_tr_dia.drop();
//             }
        }
        
        // finalize the matrix
        dia_blocks_[ill] = E * S_kron_S_ - Hdiag;
        
        // if out-of-core is enabled, dump the matrix to disk
        if (cmd_.outofcore)
        {
            // link diagonal block to a disk file
            dia_blocks_[ill].link(format("dblk-%d.ooc", ill));
            
            // save diagonal block to disk
            # pragma omp critical   
            dia_blocks_[ill].hdfsave();
            
            // release memory
            dia_blocks_[ill].drop();
        }
    }
    
    par_.wait();
    std::cout << "ok\n";
}

void NoPreconditioner::rhs (cArrayView chi, int ie, int instate, int Spin) const
{
    // shorthands
    int li = std::get<1>(inp_.instates[instate]);
    int mi = std::get<2>(inp_.instates[instate]);
    
    // shorthands
    int Nspline = s_rad_.bspline().Nspline();
    
    // j-overlaps of shape [Nangmom × Nspline]
    cArray ji_overlaps = s_rad_.overlapj
    (
        inp_.maxell,
        inp_.ki.slice(ie, ie + 1), // use only one ki
        weightEdgeDamp(s_rad_.bspline())
    );
    ji_overlaps.hdfsave("ji_overlaps.hdf");
    if (not std::isfinite(ji_overlaps.norm()))
        throw exception ("Unable to compute Riccati-Bessel function B-spline overlaps!");
    
    // j-expansions
    cArray ji_expansion = s_rad_.S().tocoo().tocsr().solve(ji_overlaps, ji_overlaps.size() / Nspline);
    if (not std::isfinite(ji_expansion.norm()))
        throw exception ("Unable to expand Riccati-Bessel function in B-splines!");
    
    // compute P-overlaps and P-expansion
    cArray Pi_overlaps, Pi_expansion;
    Pi_overlaps = s_rad_.overlapP(inp_.ni, li, weightEndDamp(s_rad_.bspline()));
    Pi_expansion = s_rad_.S().tocoo().tocsr().solve(Pi_overlaps);
    if (not std::isfinite(Pi_expansion.norm()))
        throw exception ("Unable to expand hydrogen bound orbital in B-splines!");
    
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
        for (int l = std::abs(li - inp_.L); l <= li + inp_.L; l++)
        {
            // skip wrong parity
            if ((inp_.L + li + l) % 2 != inp_.Pi)
                continue;
            
            // (anti)symmetrization
            int Sign = ((Spin + inp_.Pi) % 2 == 0) ? 1 : -1;
            
            // compute energy- and angular momentum-dependent prefactor
            Complex prefactor = std::pow(Complex(0.,1.),l)
                              * std::sqrt(special::constant::two_pi * (2 * l + 1))
                              * special::ClebschGordan(li,mi, l,0, inp_.L,mi) / inp_.ki[ie];
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
                double f1 = special::computef(lambda, l1, l2, li, l, inp_.L);
                double f2 = special::computef(lambda, l1, l2, l, li, inp_.L);
                
                // abort if any of the coefficients is non-number (factorial overflow etc.)
                if (not std::isfinite(f1))
                    throw exception ("Invalid result of computef(%d,%d,%d,%d,%d,%d)\n", lambda,l1,l2,li,l,inp_.L);
                if (not std::isfinite(f2))
                    throw exception ("Invalid result of computef(%d,%d,%d,%d,%d,%d)\n", lambda,l1,l2,l,li,inp_.L);
                
                // skip contribution if both coefficients are zero
                if (f1 == 0. and f2 == 0.)
                    continue;
                
                // load two-electron integrals if necessary
//                 SymDiaMatrix & ref_R_tr_dia = const_cast<SymDiaMatrix&>(s_rad_.R_tr_dia(lambda));
//                 if (cmd_.outofcore)
//                 {
//                     # pragma omp critical
//                     ref_R_tr_dia.hdfload();
//                 }
                
                // add the contributions
                if (f1 != 0.)
                {
//                     chi_block += (prefactor * f1) * ref_R_tr_dia.dot(Pj1);
                    chi_block += (prefactor * f1) * s_rad_.R_tr_dia(lambda).dot(Pj1);
                }
                if (f2 != 0.)
                {
                    if (Sign > 0)
//                         chi_block += (prefactor * f2) * ref_R_tr_dia.dot(Pj2);
                        chi_block += (prefactor * f2) * s_rad_.R_tr_dia(lambda).dot(Pj2);
                    else
//                         chi_block -= (prefactor * f2) * ref_R_tr_dia.dot(Pj2);
                        chi_block -= (prefactor * f2) * s_rad_.R_tr_dia(lambda).dot(Pj2);
                }
                
                // unload two-electron integrals if necessary
//                 if (cmd_.outofcore)
//                 {
//                     ref_R_tr_dia.drop();
//                 }
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
    # pragma omp parallel for schedule (dynamic,1) collapse(2) if (cmd_.parallel_block)
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
            if (cmd_.outofcore)
            {
                // read diagonal block from a linked file
                # pragma omp critical
                dia_blocks_[ill].hdfload();
            }
            
            // use the diagonal block for multiplication
            q_contrib += dia_blocks_[ill].dot(p_block, both, cmd_.parallel_dot);
            
            if (cmd_.outofcore)
            {
                // release the memory
                dia_blocks_[ill].drop();
            }
        }
        else
        {
            // compute the offdiagonal block
            for (unsigned lambda = 0; lambda <= s_rad_.maxlambda(); lambda++)
            {
                double f = special::computef(lambda, l1, l2, l1p, l2p, inp_.L);
                
                // check finiteness
                if (not std::isfinite(f))
                    throw exception ("Invalid result of computef(%d,%d,%d,%d,%d,%d).", lambda, l1, l2, l1p, l2p, inp_.L);
                
                // check non-zero
                if (f == 0.)
                    continue;
                
                // load two-electron integrals if necessary
//                 SymDiaMatrix & ref_R_tr_dia = const_cast<SymDiaMatrix&>(s_rad_.R_tr_dia(lambda));
//                 if (cmd_.outofcore)
//                 {
//                     # pragma omp critical
//                     ref_R_tr_dia.hdfload();
//                 }
                
//                 q_contrib -= Complex(f) * ref_R_tr_dia.dot(p_block);
                q_contrib -= Complex(f) * s_rad_.R_tr_dia(lambda).dot(p_block, both, cmd_.parallel_dot);
                
                // unload two-electron integrals if necessary
//                 if (cmd_.outofcore)
//                 {
//                     ref_R_tr_dia.drop();
//                 }
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

const std::string CGPreconditioner::name = "cg";
const std::string CGPreconditioner::description = "Block inversion using plain conjugate gradients.";

void CGPreconditioner::precondition (const cArrayView r, cArrayView z) const
{
    // shorthands
    int Nspline = s_rad_.bspline().Nspline();
    
    // iterations
    iArray n (l1_l2_.size());
    
    # pragma omp parallel for schedule (dynamic, 1) if (cmd_.parallel_block)
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyWork(ill))
    {
        // create segment views
        cArrayView rview (r, ill * Nspline * Nspline, Nspline * Nspline);
        cArrayView zview (z, ill * Nspline * Nspline, Nspline * Nspline);
        
        // wrappers around the callbacks
        auto inner_mmul = [&](const cArrayView a, cArrayView b) { this->CG_mmul(ill, a, b); };
        auto inner_prec = [&](const cArrayView a, cArrayView b) { this->CG_prec(ill, a, b); };
        
        // solve using the CG solver
        n[ill] = cg_callbacks < cArray, cArrayView >
        (
            rview,                  // rhs
            zview,                  // solution
            cmd_.prec_itertol,      // preconditioner tolerance
            0,                      // min. iterations
            Nspline * Nspline,      // max. iteration
            inner_prec,             // preconditioner
            inner_mmul,             // matrix multiplication
            false                   // verbose output
        );
    }
    
    // inner preconditioner info (max and avg number of iterations)
    std::cout << " | ";
    std::cout << std::setw(4) << (*std::min_element(n.begin(), n.end()));
    std::cout << std::setw(4) << (*std::max_element(n.begin(), n.end()));
    std::cout << std::setw(4) << format("%g", std::accumulate(n.begin(), n.end(), 0) / float(n.size()));
    
    // synchronize across processes
    par_.sync (z, Nspline * Nspline, l1_l2_.size());
}

void CGPreconditioner::CG_mmul (int iblock, const cArrayView p, cArrayView q) const
{
    if (cmd_.outofcore)
    {
        // load matrix from disk file
        # pragma omp critical
        dia_blocks_[iblock].hdfload();
    }
    
    // multiply
    q = dia_blocks_[iblock].dot(p, both, cmd_.parallel_dot);
    
    if (cmd_.outofcore)
    {
        // release memory
        dia_blocks_[iblock].drop();
    }
}

void CGPreconditioner::CG_prec (int iblock, const cArrayView r, cArrayView z) const
{
    z = r;
}

#ifndef NO_OPENCL

const std::string GPUCGPreconditioner::name = "gpu";
const std::string GPUCGPreconditioner::description = "Block inversion using preconditioned conjugate gradients (GPU variant).";

// kernels' source as byte array, generated by "xxd" from the CL source
char kernels_cl [] = {
    #include "kernels_cl.c"
    , 0x00 // terminate the string by zero
};

// pointer to the source; to be used in setup
char * source = &kernels_cl[0];

void GPUCGPreconditioner::setup ()
{
    // initialize parent
    NoPreconditioner::setup();
    
    // reserve space for diagonal blocks
    block_.resize(l1_l2_.size());
    
    // compute memory requirements of the preconditioner
    int order = s_bspline_.order();
    int Nspline = s_bspline_.Nspline();
    int maxell = inp_.maxell;
    unsigned long req = ((maxell+1) * (Nspline + 2*Nspline*Nspline) + (2*order+1)*(2*order+1)*Nspline*Nspline + 3*Nspline*Nspline) * 16;
    
    std::cout << "Setting up OpenCL environment\n";
    char text [1000];
    
    // use platform 0
    clGetPlatformIDs (1, &platform_, nullptr);
    clGetPlatformInfo (platform_, CL_PLATFORM_NAME, sizeof(text), text, nullptr);
    std::cout << "\tplatform: " << text << " ";
    clGetPlatformInfo (platform_, CL_PLATFORM_VENDOR, sizeof(text), text, nullptr);
    std::cout << "(" << text << ")\n";
    clGetPlatformInfo (platform_, CL_PLATFORM_VERSION, sizeof(text), text, nullptr);
    std::cout << "\tavailable version: " << text << std::endl;
    
    // use device 0
    clGetDeviceIDs (platform_, CL_DEVICE_TYPE_GPU, 1, &device_, nullptr);
//     clGetDeviceIDs (platform_, CL_DEVICE_TYPE_CPU, 1, &device_, nullptr);
    clGetDeviceInfo (device_, CL_DEVICE_NAME, sizeof(text), text, nullptr);
    std::cout << "\tdevice: " << text << " ";
    clGetDeviceInfo (device_, CL_DEVICE_VENDOR, sizeof(text), text, nullptr);
    std::cout << "(" << text << ")\n";
    cl_ulong size;
    clGetDeviceInfo (device_, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &size, 0);
    std::cout << "\tlocal memory size: " << size/1024 << " kiB\n";
    clGetDeviceInfo (device_, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &size, 0);
    std::cout << "\tglobal memory size: " << std::setprecision(3) << size/pow(1024,3) << " GiB ";
    std::cout << "(appx. " << req * 100. / size << " % will be used)\n";
    clGetDeviceInfo (device_, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_ulong), &size, 0);
    std::cout << "\tmax compute units: " << size << "\n";
    clGetDeviceInfo (device_, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(cl_ulong), &size, 0);
    std::cout << "\tmax work group size: " << size << "\n\n";
    
    // choose (e.g.) the largest workgroup
    // NOTE : This may not be the most efficient choice.
    Nlocal_ = size;
    
    // create context and command queue
    context_ = clCreateContext (nullptr, 1, &device_, nullptr, nullptr, nullptr);
    queue_ = clCreateCommandQueue (context_, device_, 0, nullptr);
    
    // setup the structure of the matrices
    iArray diags;
    for (int i = -order; i <= order; i++)
        for (int j = -order; j <= order; j++)
            diags.push_back(i * Nspline + j);
    std::string diagonals = to_string(diags, ',');
    
    // setup compile flags
    std::ostringstream flags;
    flags << " -cl-fast-relaxed-math ";
    flags << " -D ORDER="     << order     << " ";
    flags << " -D NSPLINE="   << Nspline   << " ";
    flags << " -D DIAGONALS=" << diagonals << " ";
    flags << " -D NLOCAL="    << Nlocal_   << " ";
    
    // build program
    program_ = clCreateProgramWithSource (context_, 1, const_cast<const char**>(&source), nullptr, nullptr);
    clBuildProgram (program_, 1, &device_, flags.str().c_str(), nullptr, nullptr);
    
    cl_build_status status;
    clGetProgramBuildInfo(program_, device_, CL_PROGRAM_BUILD_STATUS, sizeof(status), &status, nullptr);
    if (status != CL_SUCCESS)
    {
        std::cout << "\nSource:\n" << source << "\n";
        std::cout << "\nCommand line:\n" << flags.str() << "\n\n";
        
        char log [100000];
        clGetProgramBuildInfo(program_, device_, CL_PROGRAM_BUILD_LOG, sizeof(log), log, nullptr);
        std::cout << "clGetProgramBuildInfo: log \n" << log << "\n";
        
        throw exception ("Failed to initialize OpenCL.");
    }
    
    // set program entry points
    mmul_ = clCreateKernel(program_, "DIA_dot_vec", nullptr);
    amul_ = clCreateKernel(program_, "vec_mul_vec", nullptr);
    axby_ = clCreateKernel(program_, "a_vec_b_vec", nullptr);
    vnrm_ = clCreateKernel(program_, "vec_norm",    nullptr);
    norm_ = clCreateKernel(program_, "norm",        nullptr);
    spro_ = clCreateKernel(program_, "scalar_product", nullptr);
    krd1_ = clCreateKernel(program_, "kron_dot1",   nullptr);
    krd2_ = clCreateKernel(program_, "kron_dot2",   nullptr);
    krdv_ = clCreateKernel(program_, "kron_div",    nullptr);
    
    // resize arrays
    invsqrtS_Cl_.resize(maxell + 1);
    invCl_invsqrtS_.resize(maxell + 1);
    Dl_.resize(maxell + 1);
    
    // diagonalize overlap matrix
    ColMatrix<Complex> S = s_rad_.S().torow().T();
    ColMatrix<Complex> CL, CR; cArray DS;
    std::tie(DS,CL,CR) = S.diagonalize();
    ColMatrix<Complex> invCR = CR.invert();
    
    // convert eigenvalues to diagonal matrix
    ColMatrix<Complex> DSmat(DS.size()), invDSmat(DS.size());
    ColMatrix<Complex> DSsqrtmat(DS.size()), invDSsqrtmat(DS.size());
    for (unsigned i = 0; i < DS.size(); i++)
    {
        DSmat(i,i) = DS[i];
        DSsqrtmat(i,i) = std::sqrt(DS[i]);
        invDSmat(i,i) = 1.0 / DS[i];
        invDSsqrtmat(i,i) = 1.0 / std::sqrt(DS[i]);
    }
    
    // NOTE : Now S = CR * DSmat * CR⁻¹
    std::cout << "Setting up GPU preconditioner" << std::endl;
    std::cout << "\tS factorization residual: " << cArray((RowMatrix<Complex>(S) - RowMatrix<Complex>(CR) * DSmat * invCR).data()).norm() << std::endl;
    
    // compute √S and √S⁻¹
    RowMatrix<Complex> sqrtS = RowMatrix<Complex>(CR) * DSsqrtmat * invCR;
    RowMatrix<Complex> invsqrtS = RowMatrix<Complex>(CR) * invDSsqrtmat * invCR;
    
    // diagonalize one-electron hamiltonians for all angular momenta
    for (int l = 0; l <= maxell; l++)
    {
        std::cout << "\tH(l=" << l << ") " << std::flush;
        
        // compose the one-electron hamiltonian
        ColMatrix<Complex> Hl ( (half_D_minus_Mm1_tr_ + (0.5*l*(l+1)) * rad().Mm2()).torow() );
        
        // symmetrically transform by inverse square root of the overlap matrix
        RowMatrix<Complex> tHl = invsqrtS * Hl * ColMatrix<Complex>(invsqrtS);
        
        // diagonalize the transformed matrix
        ColMatrix<Complex> ClL, ClR; cArray Dl;
        std::tie(Dl,ClL,ClR) = ColMatrix<Complex>(tHl).diagonalize();
        
        // store the data
        Dl_[l] = Dl;
        invsqrtS_Cl_[l] = (invsqrtS * ClR).data();
        invCl_invsqrtS_[l] = (RowMatrix<Complex>(ClR.invert()) * ColMatrix<Complex>(invsqrtS)).data();
        
        // transfer data to GPU
        Dl_[l].connect(context_, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR);
        invsqrtS_Cl_[l].connect(context_, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR);
        invCl_invsqrtS_[l].connect(context_, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR);
        
        // convert Dl to matrix form and print verification
        ColMatrix<Complex> Dlmat(Dl.size()), invDlmat(Dl.size());
        for (unsigned i = 0; i < Dl.size(); i++)
        {
            Dlmat(i,i) = Dl[i];
            invDlmat(i,i) = 1.0 / Dl[i];
        }
        std::cout << "factorization residual: " << cArray((tHl - RowMatrix<Complex>(ClR) * Dlmat * ClR.invert()).data()).norm() << std::endl;
    }
    
    std::cout << std::endl;
}

void GPUCGPreconditioner::update (double E)
{
    NoPreconditioner::update(E);
    E_ = E;
    
    std::cout << "\tUpdate preconditioner..." << std::flush;
    
    // zero-pad diagonals of the diagonal DIA blocks
    for (std::size_t ill = 0; ill < l1_l2_.size(); ill++)
    {
        if (cmd_.outofcore)
        {
            // load dia_blocks_[ill] from linked file
            # pragma omp critical
            dia_blocks_[ill].hdfload();
        }
        
        block_[ill] = dia_blocks_[ill].toPaddedCols();
        
        if (cmd_.outofcore)
        {
            // link this array to a HDF file
            block_[ill].link(format("pdblk-%d.ooc", ill));
            
            // save this array to a HDF file
            # pragma omp critical
            block_[ill].hdfsave();
            
            // free memory
            block_[ill].drop();
            dia_blocks_[ill].drop();
        }
    }
    
    std::cout << "ok\n";
}

void GPUCGPreconditioner::precondition (const cArrayView r, cArrayView z) const
{
    // shorthands
    std::size_t Nspline = s_rad_.bspline().Nspline();
    std::size_t Nsegsiz = Nspline * Nspline;
    
    // iterations
    iArray n (l1_l2_.size());
    
    // some OpenCL auxiliary storage arrays (used by kernels for temporary data)
    CLArray<Complex> tmp ((Nsegsiz + Nlocal_ - 1) / Nlocal_);  tmp.connect (context_, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR);
    CLArray<double>  nrm ((Nsegsiz + Nlocal_ - 1) / Nlocal_);  nrm.connect (context_, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR);
    CLArray<Complex> tmA (Nspline * Nspline);                  tmA.connect (context_, CL_MEM_READ_WRITE);
    
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyWork(ill))
    {
        int l1 = l1_l2_[ill].first;
        int l2 = l1_l2_[ill].second;
        
        if (cmd_.outofcore)
        {
            // load linked files from disk
            # pragma omp critical
            block_[ill].hdfload();
        }
        
        // create OpenCL representation of segment views + transfer data to GPU memory
        CLArrayView<Complex> rsegment (r, ill * Nsegsiz, Nsegsiz);  rsegment.connect (context_, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR);
        CLArrayView<Complex> zsegment (z, ill * Nsegsiz, Nsegsiz);  zsegment.connect (context_, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR);
        
        // create OpenCL representation of the matrix block + transfer data to GPU memory
        CLArrayView<Complex> A (block_[ill]); A.connect(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
        
        // allocation (and upload) of an OpenCL array
        auto new_opencl_array = [&](std::size_t n) -> CLArray<Complex>
        {
            // create array
            CLArray<Complex> a(n);
            
            // connect the array to GPU
            a.connect(context_, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR);
            
            // use this array
            return a;
        };
        
        // OpenCL implementation of basic operations
        auto axby_operation = [&](Complex a, CLArrayView<Complex> x, Complex b, const CLArrayView<Complex> y) -> void
        {
            // multiply
            //     a * x + b * y -> z
            clSetKernelArg (axby_, 0, sizeof(a),      &a);
            clSetKernelArg (axby_, 1, sizeof(cl_mem), &x.handle());
            clSetKernelArg (axby_, 2, sizeof(b),      &b);
            clSetKernelArg (axby_, 3, sizeof(cl_mem), &y.handle());
            clEnqueueNDRangeKernel (queue_, axby_, 1, nullptr, &Nsegsiz, nullptr, 0, nullptr, nullptr);
            clFinish (queue_);
        };
        auto scalar_product = [&](const CLArrayView<Complex> x, const CLArrayView<Complex> y) ->  Complex
        {
            // multiply
            //     x * y -> tmp
            clSetKernelArg (spro_, 0, sizeof(cl_mem), &x.handle());
            clSetKernelArg (spro_, 1, sizeof(cl_mem), &y.handle());
            clSetKernelArg (spro_, 2, sizeof(cl_mem), &tmp.handle());
            clEnqueueNDRangeKernel (queue_, spro_, 1, nullptr, &Nsegsiz, &Nlocal_, 0, nullptr, nullptr);
            tmp.EnqueueDownload(queue_);
            clFinish (queue_);
            
            // sum the product of the arrays
            return sum(tmp);
        };
        auto compute_norm = [&](const CLArrayView<Complex> x) -> double
        {
            // multiply
            //     |x|² -> nrm
            clSetKernelArg (norm_, 0, sizeof(cl_mem), &x.handle());
            clSetKernelArg (norm_, 1, sizeof(cl_mem), &nrm.handle());
            clEnqueueNDRangeKernel (queue_, norm_, 1, nullptr, &Nsegsiz, &Nlocal_, 0, nullptr, nullptr);
            nrm.EnqueueDownload(queue_);
            clFinish (queue_);
            
            // return square root of the sum
            return std::sqrt(sum(nrm));
        };
        
        // wrappers around the CG callbacks
        auto inner_mmul = [&](const CLArrayView<Complex> a, CLArrayView<Complex> b) -> void
        {
            // multiply
            //      b = A · a
            clSetKernelArg (mmul_, 0, sizeof(cl_mem), &A.handle());
            clSetKernelArg (mmul_, 1, sizeof(cl_mem), &a.handle());
            clSetKernelArg (mmul_, 2, sizeof(cl_mem), &b.handle());
            clEnqueueNDRangeKernel (queue_, mmul_, 1, nullptr, &Nsegsiz, &Nlocal_, 0, nullptr, nullptr);
            clFinish (queue_);
        };
        auto inner_prec = [&](const CLArrayView<Complex> x, CLArrayView<Complex> y) -> void
        {
            // multiply by approximate inverse block
            
            clSetKernelArg (krd1_, 0, sizeof(cl_mem), &invCl_invsqrtS_[l1].handle());
            clSetKernelArg (krd1_, 1, sizeof(cl_mem), &invCl_invsqrtS_[l2].handle());
            clSetKernelArg (krd1_, 2, sizeof(cl_mem), &x.handle());
            clSetKernelArg (krd1_, 3, sizeof(cl_mem), &tmA.handle());
            clEnqueueNDRangeKernel (queue_, krd1_, 1, nullptr, &Nspline, nullptr, 0, nullptr, nullptr);
            clFinish (queue_);
            
            clSetKernelArg (krd2_, 0, sizeof(cl_mem), &invCl_invsqrtS_[l1].handle());
            clSetKernelArg (krd2_, 1, sizeof(cl_mem), &invCl_invsqrtS_[l2].handle());
            clSetKernelArg (krd2_, 2, sizeof(cl_mem), &y.handle());
            clSetKernelArg (krd2_, 3, sizeof(cl_mem), &tmA.handle());
            clEnqueueNDRangeKernel (queue_, krd2_, 1, nullptr, &Nspline, nullptr, 0, nullptr, nullptr);
            clFinish (queue_);
            
            clSetKernelArg (krdv_, 0, sizeof(Complex), &E_);
            clSetKernelArg (krdv_, 1, sizeof(cl_mem),  &Dl_[l1].handle());
            clSetKernelArg (krdv_, 2, sizeof(cl_mem),  &Dl_[l2].handle());
            clSetKernelArg (krdv_, 3, sizeof(cl_mem),  &y.handle());
            clEnqueueNDRangeKernel (queue_, krdv_, 1, nullptr, &Nspline, nullptr, 0, nullptr, nullptr);
            clFinish (queue_);
            
            clSetKernelArg (krd1_, 0, sizeof(cl_mem), &invsqrtS_Cl_[l1].handle());
            clSetKernelArg (krd1_, 1, sizeof(cl_mem), &invsqrtS_Cl_[l2].handle());
            clSetKernelArg (krd1_, 2, sizeof(cl_mem), &y.handle());
            clSetKernelArg (krd1_, 3, sizeof(cl_mem), &tmA.handle());
            clEnqueueNDRangeKernel (queue_, krd1_, 1, nullptr, &Nspline, nullptr, 0, nullptr, nullptr);
            clFinish (queue_);
            
            clSetKernelArg (krd2_, 0, sizeof(cl_mem), &invsqrtS_Cl_[l1].handle());
            clSetKernelArg (krd2_, 1, sizeof(cl_mem), &invsqrtS_Cl_[l2].handle());
            clSetKernelArg (krd2_, 2, sizeof(cl_mem), &y.handle());
            clSetKernelArg (krd2_, 3, sizeof(cl_mem), &tmA.handle());
            clEnqueueNDRangeKernel (queue_, krd2_, 1, nullptr, &Nspline, nullptr, 0, nullptr, nullptr);
            clFinish (queue_);
        };
        
        // solve using the CG solver
        n[ill] = cg_callbacks < CLArray<Complex>, CLArrayView<Complex> >
        (
            rsegment,               // rhs
            zsegment,               // solution to be filled
            cmd_.prec_itertol,      // tolerance
            0,                      // min. iterations
            Nsegsiz,                // max. iteration
            inner_prec,             // preconditioner
            inner_mmul,             // matrix multiplication
            false,                  // verbose output
            new_opencl_array,       // return array that is initialized and connected to GPU
            axby_operation,         // a*x+b*y -> z operation
            scalar_product,         // scalar product of two CL arrays
            compute_norm            // evaluate norm of the vector
        );
        
        // download data arrays from the GPU
        zsegment.EnqueueDownload(queue_);
        clFinish(queue_);
        
        // free GPU memory
        rsegment.disconnect();
        zsegment.disconnect();
        A.disconnect();
        
        if (cmd_.outofcore)
        {
            // release memory
            block_[ill].drop();
        }
    }
    
    // free GPU memory
    tmp.disconnect();
    nrm.disconnect();
    tmA.disconnect();
    
    // inner preconditioner info (max and avg number of iterations)
    std::cout << " | ";
    std::cout << std::setw(4) << (*std::min_element(n.begin(), n.end()));
    std::cout << std::setw(4) << (*std::max_element(n.begin(), n.end()));
    std::cout << std::setw(4) << format("%g", std::accumulate(n.begin(), n.end(), 0) / float(n.size()));
    
    // synchronize across processes
    par_.sync (z, Nspline * Nspline, l1_l2_.size());
}

#endif

const std::string JacobiCGPreconditioner::name = "Jacobi";
const std::string JacobiCGPreconditioner::description = "Block inversion using Jacobi-preconditioned conjugate gradients.";

void JacobiCGPreconditioner::update (double E)
{
    // update parent
    NoPreconditioner::update(E);
    
    // compute inverse diagonals
    # pragma omp parallel for
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyWork(ill))
    {
        if (cmd_.outofcore)
        {
            // load DIA block from linked disk file
            # pragma omp critical
            dia_blocks_[ill].hdfload();
        }
        
        invd_[ill] = 1. / dia_blocks_[ill].main_diagonal();
        
        if (cmd_.outofcore)
        {
            // release memory
            dia_blocks_[ill].drop();
        }
    }
    
    par_.wait();
}

const std::string SSORCGPreconditioner::name = "SSOR";
const std::string SSORCGPreconditioner::description = "Block inversion using SSOR-preconditioned conjugate gradients.";

void SSORCGPreconditioner::update (double E)
{
    // update parent
    NoPreconditioner::update(E);
    
    // compute preconditioner matrix
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyWork(ill))
    {
        SSOR_[ill] = SSOR(dia_blocks_[ill]);
        
        if (cmd_.outofcore)
        {
            // link to a disk file
            SSOR_[ill].link(format("ssor-%d.ooc"));
            
            // save to the file
            # pragma omp critical
            SSOR_[ill].hdfsave();
            
            // release memory
            SSOR_[ill].drop();
        }
    }
}

void SSORCGPreconditioner::CG_prec (int iblock, const cArrayView r, cArrayView z) const
{
    if (cmd_.outofcore)
    {
        // load data from a linked disk file
        # pragma omp critical
        SSOR_[iblock].hdfload();
    }
    
    z = SSOR_[iblock].upperSolve
    (
            SSOR_[iblock].dot
            (
                SSOR_[iblock].lowerSolve(r),
                diagonal
            )
    );
    
    if (cmd_.outofcore)
    {
        // release memory
        SSOR_[iblock].drop();
    }
}

const std::string SepCGPreconditioner::name = "sep";
const std::string SepCGPreconditioner::description = "Block inversion using conjugate gradients preconditioned by separated electrons hamiltonian.";

void SepCGPreconditioner::setup ()
{
    NoPreconditioner::setup();
    
    // resize arrays
    invsqrtS_Cl_.resize(inp_.maxell + 1);
    invCl_invsqrtS_.resize(inp_.maxell + 1);
    Dl_.resize(inp_.maxell + 1);
    
    // diagonalize overlap matrix
    ColMatrix<Complex> S = s_rad_.S().torow().T();
    ColMatrix<Complex> CL, CR; cArray DS;
    std::tie(DS,CL,CR) = S.diagonalize();
    ColMatrix<Complex> invCR = CR.invert();
    
    // convert eigenvalues to diagonal matrix
    ColMatrix<Complex> DSmat(DS.size()), invDSmat(DS.size());
    ColMatrix<Complex> DSsqrtmat(DS.size()), invDSsqrtmat(DS.size());
    for (unsigned i = 0; i < DS.size(); i++)
    {
        DSmat(i,i) = DS[i];
        DSsqrtmat(i,i) = std::sqrt(DS[i]);
        invDSmat(i,i) = 1.0 / DS[i];
        invDSsqrtmat(i,i) = 1.0 / std::sqrt(DS[i]);
    }
    
    // NOTE : Now S = CR * DSmat * CR⁻¹
    std::cout << "Set up SEP preconditioner" << std::endl;
    std::cout << "\tS factorization residual: " << cArray((RowMatrix<Complex>(S) - RowMatrix<Complex>(CR) * DSmat * invCR).data()).norm() << std::endl;
    
    // compute √S and √S⁻¹
    RowMatrix<Complex> sqrtS = RowMatrix<Complex>(CR) * DSsqrtmat * invCR;
    RowMatrix<Complex> invsqrtS = RowMatrix<Complex>(CR) * invDSsqrtmat * invCR;
    
    // diagonalize one-electron hamiltonians for all angular momenta
    for (int l = 0; l <= inp_.maxell; l++)
    {
        std::cout << "\tH(l=" << l << ") " << std::flush;
        
        // compose the one-electron hamiltonian
        ColMatrix<Complex> Hl ( (half_D_minus_Mm1_tr_ + (0.5*l*(l+1)) * rad().Mm2()).torow() );
        
        // symmetrically transform by inverse square root of the overlap matrix
        RowMatrix<Complex> tHl = invsqrtS * Hl * ColMatrix<Complex>(invsqrtS);
        
        // diagonalize the transformed matrix
        ColMatrix<Complex> ClL, ClR;
        std::tie(Dl_[l],ClL,ClR) = ColMatrix<Complex>(tHl).diagonalize();
        
        // store the data
        invsqrtS_Cl_[l] = invsqrtS * ClR;
        invCl_invsqrtS_[l] = RowMatrix<Complex>(ClR.invert()) * ColMatrix<Complex>(invsqrtS);
        
        // covert Dl to matrix form and print verification
        ColMatrix<Complex> Dlmat(Dl_[l].size()), invDlmat(Dl_[l].size());
        for (unsigned i = 0; i < Dl_[l].size(); i++)
        {
            Dlmat(i,i) = Dl_[l][i];
            invDlmat(i,i) = 1.0 / Dl_[l][i];
        }
        std::cout << "factorization residual: " << cArray((tHl - RowMatrix<Complex>(ClR) * Dlmat * ClR.invert()).data()).norm() << std::endl;
    }
    
    std::cout << std::endl;
}

void SepCGPreconditioner::CG_prec (int iblock, const cArrayView r, cArrayView z) const
{
    // get angular momenta of this block
    int l1 = l1_l2_[iblock].first;
    int l2 = l1_l2_[iblock].second;
    
    // get the diagonalization matrices
    RowMatrix<Complex> const & Cl1S = invCl_invsqrtS_[l1];
    RowMatrix<Complex> const & Cl2S = invCl_invsqrtS_[l2];
    RowMatrix<Complex> const & SCl1 = invsqrtS_Cl_[l1];
    RowMatrix<Complex> const & SCl2 = invsqrtS_Cl_[l2];
    
    // get one-electron eigenvalues
    cArray const & Dl1 = Dl_[l1];
    cArray const & Dl2 = Dl_[l2];
    
    // construct common diagonal
    cArray diag (Dl1.size() * Dl2.size(), E_);
    diag -= outer_product(Dl1, cArray(Dl2.size(),1.));
    diag -= outer_product(cArray(Dl1.size(),1.), Dl2);
    
    // precondition
    z = kron_dot(SCl1, SCl2, kron_dot(Cl1S, Cl2S, r) / diag);
}

const std::string ILUCGPreconditioner::name = "ILU";
const std::string ILUCGPreconditioner::description = "Block inversion using ILU-preconditioned conjugate gradients. The drop tolerance can be given as the --droptol parameter.";

void ILUCGPreconditioner::update (double E)
{
    if (E != E_)
    {
        // release outdated LU factorizations
        for (auto & lu : lu_)
            lu.drop();
        
        // release outdated CSR diagonal blocks
        for (auto & csr : csr_blocks_)
            csr.drop();
    }
    
    // update parent
    CGPreconditioner::update(E);
}

void ILUCGPreconditioner::CG_prec (int iblock, const cArrayView r, cArrayView z) const
{
    // load data from linked disk files
    if (cmd_.outofcore)
    # pragma omp critical
    {
        csr_blocks_[iblock].hdfload();
        lu_[iblock].load();
    }
    
    // check that the factorization is loaded
    if (lu_[iblock].size() == 0)
    {
        // load diagonal block from disk (if OOC)
        if (cmd_.outofcore)
        # pragma omp critical
            dia_blocks_[iblock].hdfload();
        
        // create CSR block
        csr_blocks_[iblock] = dia_blocks_[iblock].tocoo().tocsr();
        
        // release memory
        if (cmd_.outofcore)
            dia_blocks_[iblock].drop();
        
        // start timer
        Timer timer;
        
        // factorize the block and store it
        lu_[iblock].transfer(csr_blocks_[iblock].factorize(droptol_));
        
        // print time and memory info for this block (one thread at a time)
        # pragma omp critical
        std::cout << std::endl << std::setw(37) << format
        (
            "\tLU #%d (%d,%d) in %d:%02d (%d MiB)",
            iblock, l1_l2_[iblock].first, l1_l2_[iblock].second,    // block identification (id, ℓ₁, ℓ₂)
            timer.seconds() / 60, timer.seconds() % 60,             // factorization time
            lu_[iblock].size() / 1048576                            // final memory size
        );
        
        // release memory
        if (cmd_.outofcore)
            csr_blocks_[iblock].drop();
    }
    
    // precondition by LU
    z = lu_[iblock].solve(r);
    
    // release memory
    if (cmd_.outofcore)
    {
        // link to a disk file
        lu_[iblock].link(format("lu-%d.ooc", iblock));
        # pragma omp critical
        lu_[iblock].save();
        lu_[iblock].drop();
    }
}
