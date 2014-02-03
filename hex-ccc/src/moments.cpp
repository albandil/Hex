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

#include "matrix.h"
#include "moments.h"

RowMatrix<double> IntegralMoments::Im1 (int size, int alpha)
{
    // create a new empty matrix
    RowMatrix<double> I (size);
    
    // fill elements
    I.populate
    (
        [ = ] (int a, int b) -> double
        {
            int k = a + 1;
            int kp = b + 1;
            double sum = 0.;
            
            for (int i = 1; i <= std::min(k,kp); i++)
            {
                // compute fraction of factorials
                double prod = 1.;
                for (int j = i; j <= i + alpha - 2; j++)
                    prod *= j;
                
                // add the new term to sum
                sum += prod;
            }
            return sum;
        }
    );
    
    // return
    return I;
}

RowMatrix<double> IntegralMoments::Im2 (int size, int alpha)
{
    // create a new empty matrix
    RowMatrix<double> I (size);
    
    // fill elements
    I.populate
    (
        [ = ] (int a, int b) -> double
        {
            int k = a + 1;
            int kp = b + 1;
            double sum = 0.;
            
            for (int i = 1; i <= std::min(k,kp); i++)
            {
                // compute fraction of factorials
                double prod = 1.;
                for (int j = i; j <= i + alpha - 3; j++)
                    prod *= j;
                
                // add the new term to sum
                sum += prod * (k + 1 - i) * (kp + 1 - i);
            }
            
            return sum;
        }
    );
    
    // return
    return I;
}
