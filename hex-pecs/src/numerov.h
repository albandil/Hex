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

#ifndef HEX_PECS_NUMEROV_H
#define HEX_PECS_NUMEROV_H

#include "hex-arrays.h"

#include "ang.h"
#include "io.h"
#include "radial.h"

class Numerov2d
{
    public:
        
        Numerov2d (InputFile const & inp, AngularBasis const & ang, RadialBasis const & rad);
        
        void A (std::size_t i, Complex * M) const;
        void B (std::size_t i, Complex * M) const;
        void C (std::size_t i, Complex * M) const;
        
        void F (std::size_t i, Complex * v, unsigned istate) const;
    
    private:
        
        void calc_mat (std::size_t i, int term, Complex * M) const;
        
        template <class T> T coef_A (int i, int k, T h, T a, unsigned l) const
        {
//             std::cout << "coef_A " << i << " " << k << " " << h << " " << a << " " << l << std::endl;
            
            if (i > 1)
            {
                if (k  < i) return  12.0 * a;
                if (k == i) return -12.0 * (a + 1.0);
                if (k  > i) return  12.0;
            }
            if (l == 0)
            {
                T u = 2. * inp_.Z;
                
                if (k  < i) return 0;
                if (k == i) return (a + 1.) * (3.*u*u*a*a*h*h + 4.*u*u*h*h*a - 30.*u*a*h + u*u*h*h - 24.*u*h + 72.);
                if (k  > i) return  -u*u*h*h - 6.*u*a*h + 2.*u*u*h*h*a + 24.*u*h - 72.;
            }
            else
            {
                T u = 2. * inp_.Z;
                T v = l*(l+1.);
                
                if (k  < i) return 0.;
                if (k == i) return special::pow_int(1.+a,2+l) * (a*(2.+l)*(l*l+l+v+u*h) - 6. - 9.*l + v - 3.*l*l + u*h);
                if (k  > i) return a*a*u*h*(l+1.) + a*(l*l*l + 6.*l*l + (u*h + v + 11.)*l + v + 6.) + 9*l + 6. - u*h - v + 3.*l*l;
            }
        }
        
        template <class T> T coef_B (int i, int k, T h, T a, unsigned l) const
        {
//             std::cout << "coef_B " << i << " " << k << " " << h << " " << a << " " << l << std::endl;
            
            if (i > 1)
            {
                if (k  < i) return -a*a*a + a*a + a;
                if (k == i) return a*a*a + 4.*a*a + 4.*a + 1.;
                if (k  > i) return a*a + a - 1.;
            }
            if (l == 0)
            {
                T u = 2. * inp_.Z;
                
                if (k  < i) return 6.*a*(a*a - a - 1.);
                if (k == i) return (a + 1.)*(3.*u*a*a*h - 6.*a*a - 18.*a + 4.*u*a*h + u*h - 6.);
                if (k  > i) return 2.*u*a*a*h - 6.*a*a + u*a*h - u*h - 6.*a + 6.;
            }
            else
            {
                if (k  < i) return 0.;
                if (k == i) return special::pow_int(1.+a,2+l) * (l*a + 2.*a + 1.);
                if (k  > i) return special::pow_int(1.+a,2) * (l*a + a - 1.);
            }
        }
        
        InputFile const & inp_;
        AngularBasis const & ang_;
        RadialBasis const & rad_;
};

#endif
