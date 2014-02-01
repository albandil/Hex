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

#ifndef HEX_CCC_LAPACK_SUBSET
#define HEX_CCC_LAPACK_SUBSET

//
// Lapack prototypes
//

/**
 * @brief Lapack DSYEV.
 * 
 * Compute eigenvalues and optionally also eigenvectors of a symmetric real
 * matrix. The eigenvectors will overwrite the contents of the matrix. It is
 * possible to call the functions with "lwork = -1" to query for the optimal
 * workspace size. That will be stored in the first element of "work".
 * 
 * @param jobz Compute eigenvalues ("N") or also eigenvectors ("V").
 * @param uplo Use upper triangle ("U") or lower triangle ("L").
 * @param n Order of the matrix.
 * @param a Matrix elements.
 * @param lda Leading dimension of the matrix "a".
 * @param w Eigenvalues.
 * @param work Workspace.
 * @param lwork Workspace size.
 * @param info "0" on success.
 */
extern "C" void dsyev_
(
    char* jobz,
    char* uplo,
    int* n,
    double* a,
    int* lda,
    double* w,
    double* work,
    int* lwork,
    int* info
);

#endif /* HEX_CCC_LAPACK_SUBSET */
