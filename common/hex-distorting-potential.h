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

#ifndef HEX_DISTORTING_POTENTIAL
#define HEX_DISTORTING_POTENTIAL

#include "hex-special.h"

class DistortingPotentialBase : public special::RadialFunction<double>
{
    public:
        
        /// Evaluate the potential.
        virtual double operator() (double x) const = 0;
        
        /// Get far distance, where the participating hydrogen orbitals are small.
        virtual double getFarRadius () const
        {
            // if the rmax has been overriden, use the supplied value
            if (rmax_ > 0.)
                return rmax_;
            
            //
            // otherwise compute rmax using hunt & bisect run
            //
            
            // determine from parameters of the potential
            #define U_THRESHOLD    1e-100
            #define U_MAX_ITERS    1000
            
            // hunt for low value
            double far = 1., far_value;
            while ((far_value = fabs((*this)(far *= 2))) > U_THRESHOLD)
                /* continue hunting */;
            
            // bisect for exact value
            double near = 1.;
            for (int i = 0; i < U_MAX_ITERS and near != far; i++)
            {
                double middle = (near + far) * 0.5;
                double middle_val = fabs((*this)(middle));
                
                if (U_THRESHOLD > middle_val)
                    far = middle;
                    
                if (middle_val > U_THRESHOLD)
                    near = middle;
            }
            
            return (near + far) * 0.5;
        }
    
    protected:
        
        double rmax_;
};

#endif
