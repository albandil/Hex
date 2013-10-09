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

#include "arrays.h"
#include "itersolve.h"
#include "matrix.h"

cArray IC_preconditioner(cArrayView const & A, lArrayView const & I, lArrayView const & P)
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

SymDiaMatrix DIC_preconditioner(SymDiaMatrix const & A)
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

SymDiaMatrix SSOR_preconditioner(SymDiaMatrix const & A)
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
