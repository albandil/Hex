//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2015, Jakub Benda, Charles University in Prague                    //
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
