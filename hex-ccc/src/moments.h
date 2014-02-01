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

#ifndef HEX_CCC_MOMENTS
#define HEX_CCC_MOMENTS

#include "matrix.h"

namespace IntegralMoments
{
    /**
     * @brief Return @f$ I^{(-1)} @f$ integrals.
     * 
     * todo...
     */
    RowMatrix<double> Im1 (int size, int alpha);
    
    /**
     * @brief Return @f$ I^{(-2)} @f$ integrals.
     * 
     * todo...
     */
    RowMatrix<double> Im2 (int size, int alpha);
};

#endif /* HEX_CCC_MOMENTS */
