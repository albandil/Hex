//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2014, Jakub Benda, Charles University in Prague                    //
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

#include <iostream>
#include <cstdlib>
#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#include "hex-special.h"
#include "hex-version.h"

int func (double t, const double y[12], double dydt[12], void * params)
{
    double p1x = y[0], p1y = y[1], p1z = y[2];
    double p2x = y[3], p2y = y[4], p2z = y[5];
    double r1x = y[6], r1y = y[7], r1z = y[8];
    double r2x = y[9], r2y = y[10], r2z = y[11];
    double r12x = r2x-r1x, r12y = r2y-r1y, r12z = r2z-r1z;
    double r21x = r1x-r2x, r21y = r1y-r2y, r21z = r1z-r2z;
    
    double r1 = std::sqrt(r1x*r1x + r1y*r1y + r1z*r1z);
    double r2 = std::sqrt(r2x*r2x + r2y*r2y + r2z*r2z);
    double r12 = std::sqrt(r12x*r12x + r12y*r12y + r12z*r12z);
    double r21 = r12;
    
    dydt[0] = -r1x/(r1*r1*r1) + r21x/(r21*r21*r21);
    dydt[1] = -r1y/(r1*r1*r1) + r21y/(r21*r21*r21);
    dydt[2] = -r1z/(r1*r1*r1) + r21z/(r21*r21*r21);
    dydt[3] = -r2x/(r2*r2*r2) + r12x/(r12*r12*r12);
    dydt[4] = -r2y/(r2*r2*r2) + r12y/(r12*r12*r12);
    dydt[5] = -r2z/(r2*r2*r2) + r12z/(r12*r12*r12);
    dydt[6] = p1x; dydt[7] = p1y; dydt[8] = p1z;
    dydt[9] = p2x; dydt[10] = p2y; dydt[11] = p2z;
    
    return GSL_SUCCESS;
}

int jac (double t, const double y[12], double * dfdy, double dfdt[12], void * params)
{
    // TODO
    
    return GSL_SUCCESS;
}

int main (int argc, char * argv[])
{
    // write the program header
    std::cout << logo(" ");
    
    // turn off GSL and HDF exceptions
    gsl_set_error_handler_off();
    
    // disable buffering of the standard output (-> immediate logging)
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    
    // Impact energy (Rydberg)
    double Ei = 4;
    
    // Projectile distance and impact parameter.
    double R = 100;
    
    // Initialize random number generator
    std::srand(std::time(nullptr));
    
    // Initialize atomic electron
    double r1[3], p1[3];
    double var, ppt, theta, phi;
    bool bound = false;
    
    do
    {
        //- Initialize position.
        do { var = std::rand() / (RAND_MAX + 1.); var /= 1. - var; ppt = std::rand() / (RAND_MAX + 1.); } while (var > R or 4*var*var*std::exp(-2*var) > ppt);
        theta = special::constant::pi * std::rand() / (RAND_MAX + 1.);
        phi = special::constant::two_pi * std::rand() / (RAND_MAX + 1.);
        r1[0] = var * std::sin(theta) * std::cos(phi);
        r1[1] = var * std::sin(theta) * std::sin(phi);
        r1[2] = var * std::cos(theta);
        
        //- Initialize momentum.
        do { var = std::rand() / (RAND_MAX + 1.); var /= 1. - var; ppt = std::rand() / (RAND_MAX + 1.); } while (var > R or 32*special::constant::inv_pi*var*var/gsl_sf_pow_int(1.+var*var,4) > ppt);
        theta = special::constant::pi * std::rand() / (RAND_MAX + 1.);
        phi = special::constant::two_pi * std::rand() / (RAND_MAX + 1.);
        p1[0] = var * std::sin(theta) * std::cos(phi);
        p1[1] = var * std::sin(theta) * std::sin(phi);
        p1[2] = var * std::cos(theta);
        
        //- Check if it is bound.
        bound = (0.5*(p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2])*std::sqrt(r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2]) < 1);
    
    }
    while (not bound);
    
    // Initialize projectile electron
    double b = R * std::rand() / (RAND_MAX + 1.);
    double r2[3] = { b, 0., -R }, p2[3] = { 0., 0., std::sqrt(Ei) };
    
    // Simulation time.
    double t = 0;
    
    gsl_odeiv2_system sys;
        sys.dimension = 12;
        sys.function = &func;
        sys.jacobian = nullptr;
        sys.params = nullptr;
    
    gsl_odeiv2_driver * driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, 1e-3, 1e-8, 1e-5);
    
    double y[12] = { p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], r1[0], r1[1], r1[2], r2[0], r2[1], r2[2] };
    std::cout << "0 ";
    for (int i = 0; i < 12; i++)
        std::cout << y[i] << " ";
    std::cout << std::endl;
    
    int err = 0;
    double lr1 = 0, lr2 = 0;
    
    while (err == 0 and lr1 < 2*R and lr2 < 2*R)
    {
        gsl_odeiv2_driver_apply(driver, &t, t + 1, y);
    
        std::cout << t << " ";
        for (int i = 0; i < 12; i++)
            std::cout << y[i] << " ";
        std::cout << std::endl;
        
        lr1 = std::sqrt(y[6]*y[6]+y[7]*y[7]+y[8]*y[8]);
        lr2 = std::sqrt(y[9]*y[9]+y[10]*y[10]+y[11]*y[11]);
    }
    
    gsl_odeiv2_driver_free(driver);
    
    return EXIT_SUCCESS;
}
