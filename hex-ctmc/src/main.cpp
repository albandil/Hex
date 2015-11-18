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
    // names for the twelve parameters
    double const *p1 = y, *p2 = y + 3, *r1 = y + 6, *r2 = y + 9;
    
    // relative particle position
    double r12[3] = { r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2] };
    double r21[3] = { r1[0] - r2[0], r1[1] - r2[1], r1[2] - r2[2] };
    
    // distance between a particle and the nucleus
    double d1 = std::sqrt(r1[0]*r1[0] + r1[1]*r1[1] + r1[2]*r1[2]);
    double d2 = std::sqrt(r2[0]*r2[0] + r2[1]*r2[1] + r2[2]*r2[2]);
    
    // distance between the particles
    double d12 = std::sqrt(r12[0]*r12[0] + r12[1]*r12[1] + r12[2]*r12[2]);
    double d21 = d12;
    
    // acceleration of the first particle
    dydt[ 0] = -r1[0]/(d1*d1*d1) + r21[0]/(d21*d21*d21);
    dydt[ 1] = -r1[1]/(d1*d1*d1) + r21[1]/(d21*d21*d21);
    dydt[ 2] = -r1[2]/(d1*d1*d1) + r21[2]/(d21*d21*d21);
    
    // acceleration of the second particle
    dydt[ 3] = -r2[0]/(d2*d2*d2) + r12[0]/(d12*d12*d12);
    dydt[ 4] = -r2[1]/(d2*d2*d2) + r12[1]/(d12*d12*d12);
    dydt[ 5] = -r2[2]/(d2*d2*d2) + r12[2]/(d12*d12*d12);
    
    // velocity of the first particle
    dydt[ 6] = p1[0];
    dydt[ 7] = p1[1];
    dydt[ 8] = p1[2];
    
    // velocity of the second particle
    dydt[ 9] = p2[0];
    dydt[10] = p2[1];
    dydt[11] = p2[2];
    
    return GSL_SUCCESS;
}

int jac (double t, const double y[12], double * dfdy, double dfdt[12], void * params)
{
    // no explicit time dependence
    std::memset(dfdt, 0, 12 * sizeof(double));
    
    // names for the twelve parameters
    double const *r1 = y + 6, *r2 = y + 9;
    
    // relative particle position
    double r12[3] = { r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2] };
    double r21[3] = { r1[0] - r2[0], r1[1] - r2[1], r1[2] - r2[2] };
    
    // distance between a particle and the nucleus
    double d1 = std::sqrt(r1[0]*r1[0] + r1[1]*r1[1] + r1[2]*r1[2]);
    double d2 = std::sqrt(r2[0]*r2[0] + r2[1]*r2[1] + r2[2]*r2[2]);
    
    // distance between the particles
    double d12 = std::sqrt(r12[0]*r12[0] + r12[1]*r12[1] + r12[2]*r12[2]);
    double d21 = d12;
    
    // --------------------------------------------------------------- //
    
    // da1[0]/dp1
    dfdy[ 0 * 12 +  0] = 0;
    dfdy[ 0 * 12 +  1] = 0;
    dfdy[ 0 * 12 +  2] = 0;
    
    // da1[0]/dp2
    dfdy[ 0 * 12 +  3] = 0;
    dfdy[ 0 * 12 +  4] = 0;
    dfdy[ 0 * 12 +  5] = 0;
    
    // da1[0]/dr1
    dfdy[ 0 * 12 +  6] = -(d1*d1 - 3*r1[0]*r1[0])/(d1*d1*d1*d1*d1) + (d21*d21 - 3*r21[0]*r21[0])/(d21*d21*d21*d21*d21);
    dfdy[ 0 * 12 +  7] = -(d1*d1 - 3*r1[0]*r1[1])/(d1*d1*d1*d1*d1) + (        - 3*r21[0]*r21[1])/(d21*d21*d21*d21*d21);
    dfdy[ 0 * 12 +  8] = -(d1*d1 - 3*r1[0]*r1[2])/(d1*d1*d1*d1*d1) + (        - 3*r21[0]*r21[2])/(d21*d21*d21*d21*d21);
    
    // da1[0]/dr2
    dfdy[ 0 * 12 +  9] = -(d21*d21 - 3*r21[0]*r21[0])/(d21*d21*d21*d21*d21);
    dfdy[ 0 * 12 + 10] = -(        - 3*r21[0]*r21[1])/(d21*d21*d21*d21*d21);
    dfdy[ 0 * 12 + 11] = -(        - 3*r21[0]*r21[2])/(d21*d21*d21*d21*d21);
    
    // --------------------------------------------------------------- //
    
    // da1[1]/dp1
    dfdy[ 1 * 12 +  0] = 0;
    dfdy[ 1 * 12 +  1] = 0;
    dfdy[ 1 * 12 +  2] = 0;
    
    // da1[1]/dp2
    dfdy[ 1 * 12 +  3] = 0;
    dfdy[ 1 * 12 +  4] = 0;
    dfdy[ 1 * 12 +  5] = 0;
    
    // da1[1]/dr1
    dfdy[ 1 * 12 +  6] = -(d1*d1 - 3*r1[1]*r1[0])/(d1*d1*d1*d1*d1) + (        - 3*r21[1]*r21[0])/(d21*d21*d21*d21*d21);
    dfdy[ 1 * 12 +  7] = -(d1*d1 - 3*r1[1]*r1[1])/(d1*d1*d1*d1*d1) + (d21*d21 - 3*r21[1]*r21[1])/(d21*d21*d21*d21*d21);
    dfdy[ 1 * 12 +  8] = -(d1*d1 - 3*r1[1]*r1[2])/(d1*d1*d1*d1*d1) + (        - 3*r21[1]*r21[2])/(d21*d21*d21*d21*d21);
    
    // da1[1]/dr2
    dfdy[ 1 * 12 +  9] = -(        - 3*r21[1]*r21[0])/(d21*d21*d21*d21*d21);
    dfdy[ 1 * 12 + 10] = -(d21*d21 - 3*r21[1]*r21[1])/(d21*d21*d21*d21*d21);
    dfdy[ 1 * 12 + 11] = -(        - 3*r21[1]*r21[2])/(d21*d21*d21*d21*d21);
    
    // --------------------------------------------------------------- //
    
    // da1[2]/dp1
    dfdy[ 2 * 12 +  0] = 0;
    dfdy[ 2 * 12 +  1] = 0;
    dfdy[ 2 * 12 +  2] = 0;
    
    // da1[2]/dp2
    dfdy[ 2 * 12 +  3] = 0;
    dfdy[ 2 * 12 +  4] = 0;
    dfdy[ 2 * 12 +  5] = 0;
    
    // da1[2]/dr1
    dfdy[ 2 * 12 +  6] = -(d1*d1 - 3*r1[2]*r1[0])/(d1*d1*d1*d1*d1) + (        - 3*r21[2]*r21[0])/(d21*d21*d21*d21*d21);
    dfdy[ 2 * 12 +  7] = -(d1*d1 - 3*r1[2]*r1[1])/(d1*d1*d1*d1*d1) + (        - 3*r21[2]*r21[1])/(d21*d21*d21*d21*d21);
    dfdy[ 2 * 12 +  8] = -(d1*d1 - 3*r1[2]*r1[2])/(d1*d1*d1*d1*d1) + (d21*d21 - 3*r21[2]*r21[2])/(d21*d21*d21*d21*d21);
    
    // da1[2]/dr2
    dfdy[ 2 * 12 +  9] = -(        - 3*r21[2]*r21[0])/(d21*d21*d21*d21*d21);
    dfdy[ 2 * 12 + 10] = -(        - 3*r21[2]*r21[1])/(d21*d21*d21*d21*d21);
    dfdy[ 2 * 12 + 11] = -(d21*d21 - 3*r21[2]*r21[2])/(d21*d21*d21*d21*d21);
    
    // --------------------------------------------------------------- //
    
    // da2[0]/dp1
    dfdy[ 3 * 12 +  0] = 0;
    dfdy[ 3 * 12 +  1] = 0;
    dfdy[ 3 * 12 +  2] = 0;
    
    // da2[0]/dp2
    dfdy[ 3 * 12 +  3] = 0;
    dfdy[ 3 * 12 +  4] = 0;
    dfdy[ 3 * 12 +  5] = 0;
    
    // da2[0]/dr1
    dfdy[ 3 * 12 +  6] = -(d12*d12 - 3*r12[0]*r12[0])/(d12*d12*d12*d12*d12);
    dfdy[ 3 * 12 +  7] = -(        - 3*r12[0]*r12[1])/(d12*d12*d12*d12*d12);
    dfdy[ 3 * 12 +  8] = -(        - 3*r12[0]*r12[2])/(d12*d12*d12*d12*d12);
    
    // da2[0]/dr2
    dfdy[ 3 * 12 +  9] = -(d2*d2 - 3*r2[0]*r2[0])/(d2*d2*d2*d2*d2) + (d12*d12 - 3*r12[0]*r12[0])/(d12*d12*d12*d12*d12);
    dfdy[ 3 * 12 + 10] = -(d2*d2 - 3*r2[0]*r2[1])/(d2*d2*d2*d2*d2) + (        - 3*r12[0]*r12[1])/(d12*d12*d12*d12*d12);
    dfdy[ 3 * 12 + 11] = -(d2*d2 - 3*r2[0]*r2[2])/(d2*d2*d2*d2*d2) + (        - 3*r12[0]*r12[2])/(d12*d12*d12*d12*d12);
    
    // --------------------------------------------------------------- //
    
    // da2[1]/dp1
    dfdy[ 4 * 12 +  0] = 0;
    dfdy[ 4 * 12 +  1] = 0;
    dfdy[ 4 * 12 +  2] = 0;
    
    // da2[1]/dp2
    dfdy[ 4 * 12 +  3] = 0;
    dfdy[ 4 * 12 +  4] = 0;
    dfdy[ 4 * 12 +  5] = 0;
    
    // da2[1]/dr1
    dfdy[ 4 * 12 +  6] = -(        - 3*r12[1]*r12[0])/(d12*d12*d12*d12*d12);
    dfdy[ 4 * 12 +  7] = -(d12*d12 - 3*r12[1]*r12[1])/(d12*d12*d12*d12*d12);
    dfdy[ 4 * 12 +  8] = -(        - 3*r12[1]*r12[2])/(d12*d12*d12*d12*d12);
    
    // da2[1]/dr2
    dfdy[ 4 * 12 +  9] = -(d2*d2 - 3*r2[1]*r2[0])/(d2*d2*d2*d2*d2) + (        - 3*r12[1]*r12[0])/(d12*d12*d12*d12*d12);
    dfdy[ 4 * 12 + 10] = -(d2*d2 - 3*r2[1]*r2[1])/(d2*d2*d2*d2*d2) + (d12*d12 - 3*r12[1]*r12[1])/(d12*d12*d12*d12*d12);
    dfdy[ 4 * 12 + 11] = -(d2*d2 - 3*r2[1]*r2[2])/(d2*d2*d2*d2*d2) + (        - 3*r12[1]*r12[2])/(d12*d12*d12*d12*d12);
    
    // --------------------------------------------------------------- //
    
    // da2[2]/dp1
    dfdy[ 5 * 12 +  0] = 0;
    dfdy[ 5 * 12 +  1] = 0;
    dfdy[ 5 * 12 +  2] = 0;
    
    // da2[2]/dp2
    dfdy[ 5 * 12 +  3] = 0;
    dfdy[ 5 * 12 +  4] = 0;
    dfdy[ 5 * 12 +  5] = 0;
    
    // da2[2]/dr1
    dfdy[ 5 * 12 +  6] = -(        - 3*r12[2]*r12[0])/(d12*d12*d12*d12*d12);
    dfdy[ 5 * 12 +  7] = -(        - 3*r12[2]*r12[1])/(d12*d12*d12*d12*d12);
    dfdy[ 5 * 12 +  8] = -(d12*d12 - 3*r12[2]*r12[2])/(d12*d12*d12*d12*d12);
    
    // da2[2]/dr2
    dfdy[ 5 * 12 +  9] = -(d2*d2 - 3*r2[2]*r2[0])/(d2*d2*d2*d2*d2) + (        - 3*r12[2]*r12[0])/(d12*d12*d12*d12*d12);
    dfdy[ 5 * 12 + 10] = -(d2*d2 - 3*r2[2]*r2[1])/(d2*d2*d2*d2*d2) + (        - 3*r12[2]*r12[1])/(d12*d12*d12*d12*d12);
    dfdy[ 5 * 12 + 11] = -(d2*d2 - 3*r2[2]*r2[2])/(d2*d2*d2*d2*d2) + (d12*d12 - 3*r12[2]*r12[2])/(d12*d12*d12*d12*d12);
    
    // --------------------------------------------------------------- //
    
    // dv1[0]/dp1
    dfdy[ 6 * 12 +  0] = 1;
    dfdy[ 6 * 12 +  1] = 0;
    dfdy[ 6 * 12 +  2] = 0;
    
    // dv1[0]/dp2
    dfdy[ 6 * 12 +  3] = 0;
    dfdy[ 6 * 12 +  4] = 0;
    dfdy[ 6 * 12 +  5] = 0;
    
    // dv1[0]/dr1
    dfdy[ 6 * 12 +  6] = 0;
    dfdy[ 6 * 12 +  7] = 0;
    dfdy[ 6 * 12 +  8] = 0;
    
    // dv1[0]/dr2
    dfdy[ 6 * 12 +  9] = 0;
    dfdy[ 6 * 12 + 10] = 0;
    dfdy[ 6 * 12 + 11] = 0;
    
    // --------------------------------------------------------------- //
    
    // dv1[1]/dp1
    dfdy[ 7 * 12 +  0] = 0;
    dfdy[ 7 * 12 +  1] = 1;
    dfdy[ 7 * 12 +  2] = 0;
    
    // dv1[1]/dp2
    dfdy[ 7 * 12 +  3] = 0;
    dfdy[ 7 * 12 +  4] = 0;
    dfdy[ 7 * 12 +  5] = 0;
    
    // dv1[1]/dr1
    dfdy[ 7 * 12 +  6] = 0;
    dfdy[ 7 * 12 +  7] = 0;
    dfdy[ 7 * 12 +  8] = 0;
    
    // dv1[1]/dr2
    dfdy[ 7 * 12 +  9] = 0;
    dfdy[ 7 * 12 + 10] = 0;
    dfdy[ 7 * 12 + 11] = 0;
    
    // --------------------------------------------------------------- //
    
    // dv1[2]/dp1
    dfdy[ 8 * 12 +  0] = 0;
    dfdy[ 8 * 12 +  1] = 0;
    dfdy[ 8 * 12 +  2] = 1;
    
    // dv1[2]/dp2
    dfdy[ 8 * 12 +  3] = 0;
    dfdy[ 8 * 12 +  4] = 0;
    dfdy[ 8 * 12 +  5] = 0;
    
    // dv1[2]/dr1
    dfdy[ 8 * 12 +  6] = 0;
    dfdy[ 8 * 12 +  7] = 0;
    dfdy[ 8 * 12 +  8] = 0;
    
    // dv1[2]/dr2
    dfdy[ 8 * 12 +  9] = 0;
    dfdy[ 8 * 12 + 10] = 0;
    dfdy[ 8 * 12 + 11] = 0;
    
    // --------------------------------------------------------------- //
    
    // dv2[0]/dp1
    dfdy[ 9 * 12 +  0] = 0;
    dfdy[ 9 * 12 +  1] = 0;
    dfdy[ 9 * 12 +  2] = 0;
    
    // dv2[0]/dp2
    dfdy[ 9 * 12 +  3] = 1;
    dfdy[ 9 * 12 +  4] = 0;
    dfdy[ 9 * 12 +  5] = 0;
    
    // dv2[0]/dr1
    dfdy[ 9 * 12 +  6] = 0;
    dfdy[ 9 * 12 +  7] = 0;
    dfdy[ 9 * 12 +  8] = 0;
    
    // dv2[0]/dr2
    dfdy[ 9 * 12 +  9] = 0;
    dfdy[ 9 * 12 + 10] = 0;
    dfdy[ 9 * 12 + 11] = 0;
    
    // --------------------------------------------------------------- //
    
    // dv2[1]/dp1
    dfdy[10 * 12 +  0] = 0;
    dfdy[10 * 12 +  1] = 0;
    dfdy[10 * 12 +  2] = 0;
    
    // dv2[1]/dp2
    dfdy[10 * 12 +  3] = 0;
    dfdy[10 * 12 +  4] = 1;
    dfdy[10 * 12 +  5] = 0;
    
    // dv2[1]/dr1
    dfdy[10 * 12 +  6] = 0;
    dfdy[10 * 12 +  7] = 0;
    dfdy[10 * 12 +  8] = 0;
    
    // dv2[1]/dr2
    dfdy[10 * 12 +  9] = 0;
    dfdy[10 * 12 + 10] = 0;
    dfdy[10 * 12 + 11] = 0;
    
    // --------------------------------------------------------------- //
    
    // dv2[2]/dp1
    dfdy[10 * 12 +  0] = 0;
    dfdy[10 * 12 +  1] = 0;
    dfdy[10 * 12 +  2] = 0;
    
    // dv2[2]/dp2
    dfdy[10 * 12 +  3] = 0;
    dfdy[10 * 12 +  4] = 0;
    dfdy[10 * 12 +  5] = 1;
    
    // dv2[2]/dr1
    dfdy[10 * 12 +  6] = 0;
    dfdy[10 * 12 +  7] = 0;
    dfdy[10 * 12 +  8] = 0;
    
    // dv2[2]/dr2
    dfdy[10 * 12 +  9] = 0;
    dfdy[10 * 12 + 10] = 0;
    dfdy[10 * 12 + 11] = 0;
    
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
        sys.jacobian = &jac;
        sys.params = nullptr;
    
    gsl_odeiv2_driver * driver = gsl_odeiv2_driver_alloc_y_new(&sys, /*gsl_odeiv2_step_rk4*/gsl_odeiv2_step_bsimp, 1e-3, 1e-8, 1e-5); // gsl_odeiv2_step_bsimp
    
    std::ofstream trajectory ("trajectory.log");
    
    double y[12] = { p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], r1[0], r1[1], r1[2], r2[0], r2[1], r2[2] };
    trajectory << "0 ";
    for (int i = 0; i < 12; i++)
        trajectory << y[i] << " ";
    trajectory << std::endl;
    
    int err = 0;
    double lr1 = 0, lr2 = 0;
    
    while (err == 0 and lr1 < 2*R and lr2 < 2*R)
    {
        gsl_odeiv2_driver_apply(driver, &t, t + 1, y);
    
        trajectory << t << " ";
        for (int i = 0; i < 12; i++)
            trajectory << y[i] << " ";
        trajectory << std::endl;
        
        lr1 = std::sqrt(y[6]*y[6]+y[7]*y[7]+y[8]*y[8]);
        lr2 = std::sqrt(y[9]*y[9]+y[10]*y[10]+y[11]*y[11]);
    }
    
    gsl_odeiv2_driver_free(driver);
    
    // pre-run summary
    {
        double L1[3] = { r1[1]*p1[2]-r1[2]*r1[1], r1[2]*p1[0]-r1[0]*p1[2], r1[0]*p1[1]-r1[1]*p1[2] };
        double L2[3] = { r2[1]*p2[2]-r2[2]*r2[1], r2[2]*p2[0]-r2[0]*p2[2], r2[0]*p2[1]-r2[1]*p2[2] };
        
        double d1 = std::sqrt(r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2]);
        double d2 = std::sqrt(r2[0]*r2[0]+r2[1]*r2[1]+r2[2]*r2[2]);
        
        double T1 = (p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2]), V1 = - 2./d1;
        double T2 = (p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2]), V2 = - 2./d2;
        
        std::cout << "Initial state" << std::endl;
        std::cout << "    atomic electron" << std::endl;
        std::cout << "        distance [a0]:    " << d1 << std::endl;;
        std::cout << "        kin. energy [Ry]: " << T1 << std::endl;
        std::cout << "        tot. energy [Ry]: " << T1 + V1 << std::endl;
        std::cout << "        ang. mom. [a.u.]: " << std::sqrt(L1[0]*L1[0]+L1[1]*L1[1]+L1[2]*L1[2]) << std::endl;
        std::cout << "    projectile electron" << std::endl;
        std::cout << "        distance [a0]:    " << d2<< std::endl;;
        std::cout << "        kin. energy [Ry]: " << T1 << std::endl;
        std::cout << "        tot. energy [Ry]: " << T1 + V1 << std::endl;
        std::cout << "        ang. mom. [a.u.]: " << std::sqrt(L2[0]*L2[0]+L2[1]*L2[1]+L2[2]*L2[2]) << std::endl;
        std::cout << "    total energy [Ry]:    " << T1 + V1 + T2 + V2 << std::endl;
        std::cout << std::endl;
    }
    
    // post-run summary
    {
        double *p1 = y, *p2 = y + 3, *r1 = y + 6, *r2 = y + 9;
        
        double L1[3] = { r1[1]*p1[2]-r1[2]*r1[1], r1[2]*p1[0]-r1[0]*p1[2], r1[0]*p1[1]-r1[1]*p1[2] };
        double L2[3] = { r2[1]*p2[2]-r2[2]*r2[1], r2[2]*p2[0]-r2[0]*p2[2], r2[0]*p2[1]-r2[1]*p2[2] };
        
        double d1 = std::sqrt(r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2]);
        double d2 = std::sqrt(r2[0]*r2[0]+r2[1]*r2[1]+r2[2]*r2[2]);
        
        double T1 = (p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2]), V1 = - 2./d1;
        double T2 = (p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2]), V2 = - 2./d2;
        
        std::cout << "Initial state" << std::endl;
        std::cout << "    atomic electron" << std::endl;
        std::cout << "        distance [a0]:    " << d1 << std::endl;;
        std::cout << "        kin. energy [Ry]: " << T1 << std::endl;
        std::cout << "        tot. energy [Ry]: " << T1 + V1 << std::endl;
        std::cout << "        ang. mom. [a.u.]: " << std::sqrt(L1[0]*L1[0]+L1[1]*L1[1]+L1[2]*L1[2]) << std::endl;
        std::cout << "    projectile electron" << std::endl;
        std::cout << "        distance [a0]:    " << d2<< std::endl;;
        std::cout << "        kin. energy [Ry]: " << T1 << std::endl;
        std::cout << "        tot. energy [Ry]: " << T1 + V1 << std::endl;
        std::cout << "        ang. mom. [a.u.]: " << std::sqrt(L2[0]*L2[0]+L2[1]*L2[1]+L2[2]*L2[2]) << std::endl;
        std::cout << "    total energy [Ry]:    " << T1 + V1 + T2 + V2 << std::endl;
        std::cout << std::endl;
    }
    
    return EXIT_SUCCESS;
}
