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

#ifndef NO_OPENBLAS
    #include <omp.h>
    extern "C" void openblas_set_num_threads (int);
#endif

#include "arrays.h"
#include "gauss.h"
#include "input.h"
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
    if ((size_t)A.rows() != b.size())
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
        for (int l = std::abs(li - inp_.L); l <= li + inp_.L; l++)
        {
            // skip wrong parity
            if ((inp_.L + li + l) % 2 != inp_.Pi)
                continue;
            
            // (anti)symmetrization
            int Sign = ((Spin + inp_.Pi) % 2 == 0) ? 1. : -1.;
            
            // compute energy- and angular momentum-dependent prefactor
            Complex prefactor = std::pow(Complex(0.,1.),l) * std::sqrt(2*special::constant::pi*(2*l+1)) / Complex(inp_.ki[ie]); 
            prefactor *= special::ClebschGordan(li,mi,l,0,inp_.L,mi);
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
                    throw exception ("Invalid result of computef(%d,%d,%d,%d,%d,%d)\n", lambda, l1, l2, li, l, inp_.L);
                if (not std::isfinite(f2))
                    throw exception ("Invalid result of computef(%d,%d,%d,%d,%d,%d)\n", lambda, l1, l2, l, li, inp_.L);
                
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
        cg_callbacks < cArray, cArrayView >
        (
            rview,                  // rhs
            zview,                  // solution
            cmd_.itertol,        // tolerance
            0,                      // min. iterations
            Nspline * Nspline,      // max. iteration
            inner_prec,             // preconditioner
            inner_mmul,             // matrix multiplication
            false                   // verbose output
        );
    }
    
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

const std::string GPUCGPreconditioner::name = "gpuJacobi";
const std::string GPUCGPreconditioner::description = "Block inversion using Jacobi-preconditioned conjugate gradients (GPU variant).";

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
    csr_blocks_.resize(l1_l2_.size());
    block_.resize(l1_l2_.size());
    
    std::cout << "Setting up OpenCL environment\n";
    char text [1000];
    
    // use platform 0
    clGetPlatformIDs (1, &platform_, nullptr);
    clGetPlatformInfo (platform_, CL_PLATFORM_NAME, sizeof(text), text, nullptr);
    std::cout << "\tplatform: " << text << " ";
    clGetPlatformInfo (platform_, CL_PLATFORM_VENDOR, sizeof(text), text, nullptr);
    std::cout << "(" << text << ")\n";
    
    // use device 0
    clGetDeviceIDs (platform_, CL_DEVICE_TYPE_GPU, 1, &device_, nullptr);
    clGetDeviceInfo (device_, CL_DEVICE_NAME, sizeof(text), text, nullptr);
    std::cout << "\tdevice: " << text << " ";
    clGetDeviceInfo (device_, CL_DEVICE_VENDOR, sizeof(text), text, nullptr);
    std::cout << "(" << text << ")\n";
    cl_ulong size;
    clGetDeviceInfo (device_, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &size, 0);
    std::cout << "\tlocal memory size: " << size/1024 << " kiB\n";
    clGetDeviceInfo (device_, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &size, 0);
    std::cout << "\tglobal memory size: " << std::setprecision(3) << size/pow(1024,3) << " GiB\n";
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
    int order = s_bspline_.order();
    int Nspline = s_bspline_.Nspline();
    iArray diags;
    for (int i = -order; i <= order; i++)
        for (int j = -order; j <= order; j++)
            diags.push_back(i * Nspline + j);
    std::string diagonals = to_string(diags, ',');
    
    // setup compile flags
    std::ostringstream flags;
    flags << "-cl-strict-aliasing -cl-fast-relaxed-math ";
    flags << "-D ORDER="     << order     << " ";
    flags << "-D NSPLINE="   << Nspline   << " ";
    flags << "-D DIAGONALS=" << diagonals << " ";
    flags << "-D NLOCAL="    << Nlocal_   << " ";
    
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
}

void GPUCGPreconditioner::update (double E)
{
    NoPreconditioner::update(E);
    
    std::cout << "\tUpdate preconditioner..." << std::flush;
    
    // zero-pad diagonals of the diagonal DIA blocks
    for (size_t ill = 0; ill < l1_l2_.size(); ill++)
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
    size_t Nspline = s_rad_.bspline().Nspline();
    size_t Nsegsiz = Nspline * Nspline;
    
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyWork(ill))
    {
        std::cout << "\t\t\t -precondition block " << ill << "\n";
        
        if (cmd_.outofcore)
        {
            // load linked files from disk
            # pragma omp critical
            dia_blocks_[ill].hdfload();
            # pragma omp critical
            block_[ill].hdfload();
        }
        
        // create OpenCL representation of segment views + transfer data to GPU memory
        CLArrayView<Complex> rsegment (r, ill * Nsegsiz, Nsegsiz);  rsegment.connect (context_, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR);
        CLArrayView<Complex> zsegment (z, ill * Nsegsiz, Nsegsiz);  zsegment.connect (context_, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR);
        CLArray<Complex> tmp ((Nsegsiz + Nlocal_ - 1) / Nlocal_);  tmp.connect (context_, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR);
        CLArray<double>  nrm ((Nsegsiz + Nlocal_ - 1) / Nlocal_);  nrm.connect (context_, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR);
        
        // create OpenCL representation of the matrix block + transfer data to GPU memory
        CLArrayView<Complex> A (block_[ill]); A.connect(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
        
        // create OpenCL representation of the inverse diagonal + transfer data to GPU memory
        CLArray<Complex> invd (1. / dia_blocks_[ill].main_diagonal());  invd.connect(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
        
        // allocation (and upload) of an OpenCL array
        auto new_opencl_array = [&](size_t n) -> CLArray<Complex>
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
            clSetKernelArg (axby_, 0, sizeof(a),          &a);
            clSetKernelArg (axby_, 1, sizeof(x.handle()), &x.handle());
            clSetKernelArg (axby_, 2, sizeof(b),          &b);
            clSetKernelArg (axby_, 3, sizeof(y.handle()), &y.handle());
            clEnqueueNDRangeKernel (queue_, axby_, 1, nullptr, &Nsegsiz, nullptr, 0, nullptr, nullptr);
            clFinish (queue_);
        };
        auto scalar_product = [&](const CLArrayView<Complex> x, const CLArrayView<Complex> y) ->  Complex
        {
            // multiply
            //     x * y -> tmp
            clSetKernelArg (spro_, 0, sizeof(x.handle()), &x.handle());
            clSetKernelArg (spro_, 1, sizeof(y.handle()), &y.handle());
            clSetKernelArg (spro_, 2, sizeof(tmp.handle()), &tmp.handle());
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
            clSetKernelArg (norm_, 0, sizeof(x.handle()), &x.handle());
            clSetKernelArg (norm_, 1, sizeof(nrm.handle()), &nrm.handle());
            clEnqueueNDRangeKernel (queue_, norm_, 1, nullptr, &Nsegsiz, &Nlocal_, 0, nullptr, nullptr);
            nrm.EnqueueDownload(queue_);
            clFinish (queue_);
            
            // return square root of the sum
            return sqrt(sum(nrm));
        };
        
        // wrappers around the CG callbacks
        auto inner_mmul = [&](const CLArrayView<Complex> a, CLArrayView<Complex> b) -> void
        {
            // multiply
            //      b = A · a
            clSetKernelArg (mmul_, 0, sizeof(A.handle()),  &A.handle());
            clSetKernelArg (mmul_, 1, sizeof(a.handle()),  &a.handle());
            clSetKernelArg (mmul_, 2, sizeof(b.handle()),  &b.handle());
            clEnqueueNDRangeKernel (queue_, mmul_, 1, nullptr, &Nsegsiz, &Nlocal_, 0, nullptr, nullptr);
            clFinish (queue_);
        };
        auto inner_prec = [&](const CLArrayView<Complex> x, CLArrayView<Complex> y) -> void
        {
            // multiply by the inverse diagonal
            //     y = invd * x
            clSetKernelArg (amul_, 0, sizeof(invd.handle()), &invd.handle());
            clSetKernelArg (amul_, 1, sizeof(x.handle()),    &x.handle());
            clSetKernelArg (amul_, 2, sizeof(y.handle()),    &y.handle());
            clEnqueueNDRangeKernel (queue_, amul_, 1, nullptr, &Nsegsiz, nullptr, 0, nullptr, nullptr);
            clFinish (queue_);
        };
        
        // solve using the CG solver
        cg_callbacks < CLArray<Complex>, CLArrayView<Complex> >
        (
            rsegment,               // rhs
            zsegment,               // solution to be filled
            cmd_.itertol,           // tolerance
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
        tmp.disconnect();
        nrm.disconnect();
        A.disconnect();
        
        if (cmd_.outofcore)
        {
            // release memory
            block_[ill].drop();
            dia_blocks_[ill].drop();
        }
    }
    
    // synchronize across processes
    par_.sync (z, Nspline * Nspline, l1_l2_.size());
}

#endif

const std::string JacobiCGPreconditioner::name = "Jacobi";
const std::string JacobiCGPreconditioner::description = "Block inversion using Jacobi-preconditioned conjugate gradients.";

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

const std::string ILUCGPreconditioner::name = "ILU";
const std::string ILUCGPreconditioner::description = "Block inversion using ILU-preconditioned conjugate gradients. The drop tolerance can be given as the --droptol parameter.";

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
    
    // write info
    std::cout << "\t[" << par_.iproc() << "] Update preconditioner";
    if (cmd_.concurrent_factorizations > 1)
        std::cout << " using " << cmd_.concurrent_factorizations << " concurrent 1-thread parallelizations";
    else
        std::cout << " sequentially, using all available threads for each factorization";
    std::cout << std::endl;
    
    // for all diagonal blocks
    # pragma omp parallel for if (cmd_.concurrent_factorizations > 1) num_threads (cmd_.concurrent_factorizations)
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyWork(ill))
    {
        // delete old factorization first to maximally spare memory
        lu_[ill].drop();
        
        // load diagonal block from disk (if OOC)
        if (cmd_.outofcore)
        {
            // load DIA block from a linked disk file
            # pragma omp critical
            dia_blocks_[ill].hdfload();
        }
        
        // start timer
        Timer timer;
        
        // create CSR block
        csr_blocks_[ill] = dia_blocks_[ill].tocoo().tocsr();
        
        // factorize the block and store it in lu_[ill]
        lu_[ill].transfer
        (
            // transfer data (use &&-constructor)
            std::move
            (
                csr_blocks_[ill].factorize(droptol_)
            )
        );
        
        // print time and memory info for this block (one thread at a time)
        # pragma omp critical
        std::cout << format
        (
            "\t\t- block #%d (%d,%d) in %d:%02d (%d MiB)",
            ill, l1_l2_[ill].first, l1_l2_[ill].second, // block identification (id, ℓ₁, ℓ₂)
            timer.seconds() / 60,                       // factorization time: minutes
            timer.seconds() % 60,                       // factorization time: seconds
            lu_[ill].size() / 1048576                   // final memory size
        ) << std::endl;
        
        // save factorization to disk (if OOC)
        if (cmd_.outofcore)
        {
            // link CSR block to a disk file
            csr_blocks_[ill].link(format("csr-%d.ooc", ill));
            # pragma omp critical
            csr_blocks_[ill].hdfsave();
            
            // release memory
            csr_blocks_[ill].drop();
            dia_blocks_[ill].drop();
            
            // link to a disk file
            lu_[ill].link(format("lu-%d.ooc", ill));
            # pragma omp critical
            lu_[ill].save();
            lu_[ill].drop();
        }
    }
}

void ILUCGPreconditioner::CG_prec (int iblock, const cArrayView r, cArrayView z) const
{
    if (cmd_.outofcore)
    {
        // load data from linked disk files
        # pragma omp critical
        csr_blocks_[iblock].hdfload();
        # pragma omp critical
        lu_[iblock].load();
    }
    
    z = lu_[iblock].solve(r);
    
    if (cmd_.outofcore)
    {
        // release memory
        csr_blocks_[iblock].drop();
        lu_[iblock].drop();
    }
}

#if 0
const std::string DICCGPreconditioner::name = "DIC";
const std::string DICCGPreconditioner::description = "Block inversion using DIC-preconditioned conjugate gradients [not working].";

void DICCGPreconditioner::setup()
{
    CGPreconditioner::setup();
    DIC_.resize(l1_l2_.size());
}

void DICCGPreconditioner::update(double E)
{
    CGPreconditioner::update(E);
    
    // diagonal incomplete Cholesky factorization
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyWork(ill))
        DIC_[ill] = DIC(dia_blocks_[ill]);
}

#ifndef NO_LAPACK
const std::string SPAICGPreconditioner::name = "SPAI";
const std::string SPAICGPreconditioner::description = "Block inversion using SPAI-preconditioned conjugate gradients [not working].";

void SPAICGPreconditioner::setup()
{
    // setup parent
    NoPreconditioner::setup();
    
    // reserve space for preconditioner
    spai_.resize(l1_l2_.size());
}

void SPAICGPreconditioner::update (double E)
{
    // update parent
    NoPreconditioner::update(E);
    
    std::cout << "\tCompute SPAI preconditioner... " << std::flush;
    Timer::timer().start();
    
    // set diagonal pattern for the inverse
    iArray diagset = dia_blocks_[0].diag();
    
    // for all diagonal blocks sitting on this CPU compute SPAI
    for (unsigned ill = 0; ill < l1_l2_.size(); ill++) if (par_.isMyWork(ill))
        spai_[ill] = SPAI(dia_blocks_[ill], diagset).tocsr();
    
    int secs = Timer::timer().stop();
    std::cout << "done in " << (secs / 60) << std::setw(2) << std::setfill('0') << (secs % 60) << "\n";
}
#endif

const std::string TwoLevelPreconditioner::name = "two";
const std::string TwoLevelPreconditioner::description = "Block inversion using conjugate gradients preconditioned by solution of coarse system [not working].";

void TwoLevelPreconditioner::setup ()
{
    std::cout << "Setting up ML preconditioner.\n";
    std::cout << "\t- rknots = " << p_bspline_.rknots() << "\n";
    std::cout << "\t- cknots = " << p_bspline_.cknots() << "\n";
    std::cout << "\t- B-spline count = " << p_bspline_.Nspline() << "\n\n";
    
    // setup parent
    SSORCGPreconditioner::setup();
    
    // which lambdas to keep on this CPU TODO
    Array<bool> lambdas (s_rad_.maxlambda() + 1, true);
    
    // compute radial integrals
    p_rad_.setupOneElectronIntegrals();
    p_rad_.setupTwoElectronIntegrals(par_, lambdas);
    
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

const std::string MultiresPreconditioner::name = "res";
const std::string MultiresPreconditioner::description = "Multi-resolution preconditioner [implementation not finished].";

MultiresPreconditioner::MultiresPreconditioner (
    Parallel const & par, InputFile const & inp, std::vector<std::pair<int,int>> const & ll, Bspline const & bspline, CommandLine const & cmd
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
            ),
            cmd
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
                ),
                cmd
            )
        );
    }
}

void MultiresPreconditioner::precondition (const cArrayView r, cArrayView z) const
{
    // HOA WHOW... TODO
    
}

#endif
