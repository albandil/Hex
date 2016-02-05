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

#include "hex-hydrogen.h"

#include "numerov.h"

Complex operator+ (Complex z, unsigned i) { return Complex(z.real() + i, z.imag()); }
Complex operator+ (unsigned i, Complex z) { return Complex(z.real() + i, z.imag()); }
Complex operator* (unsigned i, Complex z) { return Complex(z.real() * i, z.imag() * i); }
Complex operator* (Complex z, unsigned i) { return Complex(z.real() * i, z.imag() * i); }

Numerov2d::Numerov2d (InputFile const & inp, AngularBasis const & ang, RadialBasis const & rad)
    : inp_(inp), ang_(ang), rad_(rad) {}

void Numerov2d::A (std::size_t i, Complex * M) const { calc_mat(i, -1, M); }
void Numerov2d::B (std::size_t i, Complex * M) const { calc_mat(i,  0, M); }
void Numerov2d::C (std::size_t i, Complex * M) const { calc_mat(i, +1, M); }

void Numerov2d::F (std::size_t i, Complex * v, unsigned istate) const
{
    // 'i' is fixed, 'j' is the free index, 'k' and 'l' are summed
    
    // skip empty vectors
    if (i == 0)
        return;
    
    // for all angular momentum states
    for (std::size_t m = 0; m < ang_.size(); m++)
    {
        // get angular momenta
        unsigned l1 = ang_.state(m).first;
        unsigned l2 = ang_.state(m).second;
        
        // for all elements of the right-hand side
        for (std::size_t j = 1; j <= i; j++)
        {
            // clear the element
            v[m * i + (j - 1)] = 0.;
            
            // get mesh parameters
            double h = rad_.rgrid[i] - rad_.rgrid[i-1];
            double t = rad_.rgrid[j] - rad_.rgrid[j-1];
            double alpha = (rad_.rgrid[i+1] - rad_.rgrid[i]) / h;
            double beta  = (rad_.rgrid[j+1] - rad_.rgrid[j]) / t;
            
            // for all neighbours to sum
            for (std::size_t k = std::max<std::size_t>(1, i - 1); k <= i + 1 and k < rad_.rgrid.size(); k++)
            for (std::size_t l = std::max<std::size_t>(1, j - 1); l <= j + 1 and l < rad_.rgrid.size(); l++)
            {
                // get Numerov discretization coefficients
                double Bik = coef_B(i, k, h, alpha, l1);
                double Djl = coef_B(j, l, t, beta,  l2);
                
                // radial coordinates
                double r1 = rad_.rgrid[k];
                double r2 = rad_.rgrid[l];
                
                // smaller and larger radial coordinate (sorted by real part)
                double rmin = std::min(r1, r2, Complex_realpart_less);
                double rmax = std::max(r1, r2, Complex_realpart_less);
                
                // evaluate the right-hand side at (r1,r2)
                for (unsigned ell = 0; ell <= ang_.maxell(); ell++)
                {
                    // initial state quantum numbers
                    unsigned ni = inp_.istates[istate].ni;
                    unsigned li = inp_.istates[istate].li;
                    int      mi = inp_.istates[istate].mi;
                    
                    // impact momentum
                    double ki = std::sqrt(inp_.Etot + 1./(ni*ni));
                    
                    // total angular parameters
                    int S  = inp_.S;
                    int L  = inp_.L;
                    int Pi = inp_.Pi;
                    
                    // symmetry sign
                    double sign = ((S + Pi) % 2 == 0 ? +1. : -1.);
                    
                    // contribution to the right-hand side
                    Complex chikl = 0;
                    
                    // for all multipoles
                    for (unsigned lambda = 0; lambda <= ang_.maxlambda(); lambda++)
                    {
                        // monopole
                        if (lambda == 0)
                        {
                            // direct
                            if (ang_.f(lambda, l1, l2, li, ell) != 0. and r1 > r2)
                                chikl += ang_.f(lambda, l1, l2, li, ell) * (1./r1 - 1./r2) * Hydrogen::P(ni,li,r1) * special::ric_j(ell,ki*r2);
                            // exchange
                            if (ang_.f(lambda, l1, l2, ell, li) != 0. and r2 > r1)
                                chikl += ang_.f(lambda, l1, l2, ell, li) * (1./r2 - 1./r1) * special::ric_j(ell,ki*r1) * Hydrogen::P(ni,li,r2) * sign;
                        }
                        // higher multipoles
                        else
                        {
                            // direct
                            if (ang_.f(lambda, l1, l2, li, ell) != 0.)
                                chikl += ang_.f(lambda, l1, l2, li, ell) * special::pow_int(rmin/rmax, lambda)/rmax * Hydrogen::P(ni,li,r1) * special::ric_j(ell,ki*r2);
                            // exchange
                            if (ang_.f(lambda, l1, l2, ell, li) != 0.)
                                chikl += ang_.f(lambda, l1, l2, ell, li) * special::pow_int(rmin/rmax, lambda)/rmax * special::ric_j(ell,ki*r1) * Hydrogen::P(ni,li,r2) * sign;
                        }
                    }
                    
                    // add prefactor
                    chikl *= std::sqrt(special::constant::two_pi * (2*ell+1)) / ki * special::ClebschGordan(li,mi,ell,0,L,mi) * special::pow_int(1.0_i, ell);
                    
                    // update the vector
                    v[m * i + (j - 1)] += h*h*t*t*Bik*Djl*chikl;
                }
            }
        }
    }
}

void Numerov2d::calc_mat (std::size_t i, int term, Complex * M) const
{
    // skip empty matrices
    if (i == 0 and term < 0)
        return;
    
    // number of blocks
    std::size_t Nang = ang_.size();
    
    // number of rows and columns of a block in this matrix
    std::size_t Nrows = i;
    std::size_t Ncols = i + term;
    
    // expand the other index
    std::size_t k = i + term;
    
    // number of elements in every block
    std::size_t block_vol = Nrows * 3;
    
    // for all blocks
    for (std::size_t m = 0; m < Nang; m++)
    for (std::size_t n = 0; n < Nang; n++)
    {
        // get block pointer
        Complex * pM = M + (m * Nang + n) * block_vol;
        
        // get row angular momenta
        unsigned l1 = ang_.state(m).first;
        unsigned l2 = ang_.state(m).second;
        
        // for all elements in the block
        for (std::size_t j = 1; j <= Nrows and j < rad_.Npts; j++)
        {
        for (std::size_t l = std::max<std::size_t>(1, j - 1); l <= j + 1 and l <= Ncols and l < rad_.Npts; l++)
        {
            // erase any previous elements
            pM[(j - 1) * 3 + (l + 1 - j)] = 0;
            
            // skip the rest if out of grid
            if (k >= rad_.Npts)
                continue;
            
            // get mesh parameters
            Complex h = rad_.grid[i] - rad_.grid[i-1];
            Complex t = rad_.grid[j] - rad_.grid[j-1];
            Complex alpha = (rad_.grid[i+1] - rad_.grid[i]) / h;
            Complex beta  = (rad_.grid[j+1] - rad_.grid[j]) / t;
            
            // get Numerov discretization coefficients
            Complex Aik = coef_A(i, k, h, alpha, l1);
            Complex Bik = coef_B(i, k, h, alpha, l1);
            Complex Cjl = coef_A(j, l, t, beta,  l2);
            Complex Djl = coef_B(j, l, t, beta,  l2);
            
//             std::cout << "i = " << i << ", j = " << j << ", k = " << k << ", l = " << l << std::endl;
//             std::cout << "    r[" << i << "] = " << rad_.grid[i] << ", r[" << j << "] = " << rad_.grid[j] << ", r[" << k << "] = " << rad_.grid[k] << ", r[" << l << "] = " << rad_.grid[l] << std::endl;
//             std::cout << "    h = " << h << ", alpha = " << alpha << ", t = " << t << ", beta = " << beta << std::endl;
//             std::cout << "    Aik = " << Aik << ", Bik = " << Bik << ", Cjl = " << Cjl << ", Djl = " << Djl << std::endl;
            
            // calculate two-electron part
            for (unsigned lambda = 0; lambda <= ang_.maxlambda(); lambda++) if (ang_.f(lambda, m, n) != 0.)
            {
                double rmin = std::min(rad_.grid[k].real(), rad_.grid[l].real());
                double rmax = std::max(rad_.grid[k].real(), rad_.grid[l].real());
                pM[(j - 1) * 3 + (l + 1 - j)] -= h*h*t*t*Bik*Djl * ang_.f(lambda, m, n) * special::pow_int(rmin/rmax, lambda) / rmax;
            }
            
//             std::cout << "    E = " << inp_.Etot << std::endl;
//             std::cout << "    tD1: " << 0.5*h*h*Bik*Cjl << std::endl;
//             std::cout << "    tD2: " << 0.5*t*t*Aik*Djl << std::endl;
//             std::cout << "    tE : " << h*h*t*t*Bik*Djl*inp_.Etot << std::endl;
            
            // calculate one-electron (diagoanl) part, if any
            if (m == n)
            {
                // calculate potential
                Complex Vkl = 0.;
                if (k > 0) Vkl += 0.5 * l1 * (l1 + 1) / (rad_.grid[k].real() * rad_.grid[k].real()) - inp_.Z / rad_.grid[k].real();
                if (l > 0) Vkl += 0.5 * l2 * (l2 + 1) / (rad_.grid[l].real() * rad_.grid[l].real()) - inp_.Z / rad_.grid[l].real();
                
                // update element
                pM[(j - 1) * 3 + (l + 1 - j)] += 0.5*h*h*Bik*Cjl + 0.5*t*t*Aik*Djl + h*h*t*t*Bik*Djl*(inp_.Etot - Vkl);
            }
            
        }
//             std::cout << pM[0] << " " << pM[1] << " " << pM[2] << std::endl;
//             std::exit(0);
        }
    }
}
