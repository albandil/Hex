//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 20156 Jakub Benda, Charles University in Prague                    //
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
            Complex h = rad_.grid[i] - rad_.grid[i-1];
            Complex t = rad_.grid[j] - rad_.grid[j-1];
            Complex alpha = (rad_.grid[i+1] - rad_.grid[i]) / h;
            Complex beta  = (rad_.grid[j+1] - rad_.grid[j]) / t;
            
            // for all neighbours to sum
            for (std::size_t k = std::max<std::size_t>(1, i - 1); k <= i + 1 and k + 1 < rad_.rgrid.size(); k++)
            for (std::size_t l = std::max<std::size_t>(1, j - 1); l <= j + 1 and l + 1 < rad_.rgrid.size(); l++)
            {
                // radial coordinates
                double r1 = rad_.rgrid[k];
                double r2 = rad_.rgrid[l];
                
                // smaller and larger radial coordinate (sorted by real part)
                double rmin = std::min(r1, r2, Complex_realpart_less);
                double rmax = std::max(r1, r2, Complex_realpart_less);
                
                // get Numerov discretization coefficients
                Complex Bik = coef_B(i, k, h, alpha, l1);
                Complex Djl = coef_B(j, l, t, beta,  l2);
                
                // evaluate the right-hand side at (r1,r2)
                for (unsigned ell = 0; ell <= ang_.maxell(); ell++)
                {
                    // initial state quantum numbers
                    unsigned ni = inp_.istates[istate].n;
                    unsigned li = inp_.istates[istate].l;
                    int      mi = inp_.istates[istate].m;
                    
                    // impact momentum
                    double ki = std::sqrt(inp_.Etot + 1./(ni*ni));
                    
                    // total angular parameters
                    int L  = inp_.L;
                    int S  = inp_.S;
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
                    v[m * i + (j - 1)] += -h*h*t*t*Bik*Djl*chikl;
                }
            }
        }
    }
}

void Numerov2d::mask
(
    int term, int s,
    unsigned int l1, unsigned int l2, std::size_t i, std::size_t j,
    std::array<std::tuple<int,unsigned,unsigned,std::size_t,std::size_t>, 9> & mask
) const
{
    // erase the mask
    std::fill(mask.begin(), mask.end(), std::make_tuple(0,0,0,0,0));
    
    // are we far enough from the i = j diagonal?
    if (i > j + 1)
    {
        //   A  B  C
        //   A  B  C
        //   A  B  C
        
        if (term == -1) // A
        {
            mask[0] = std::make_tuple(1,l1,l2,i-1,j-1);
            mask[1] = std::make_tuple(1,l1,l2,i-1,j  );
            mask[2] = std::make_tuple(1,l1,l2,i-1,j+1);
        }
        else if (term == 0) // B
        {
            mask[3] = std::make_tuple(1,l1,l2,i  ,j-1);
            mask[4] = std::make_tuple(1,l1,l2,i  ,j  );
            mask[5] = std::make_tuple(1,l1,l2,i  ,j+1);
        }
        else if (term == +1) // C
        {
            mask[6] = std::make_tuple(1,l1,l2,i+1,j-1);
            mask[7] = std::make_tuple(1,l1,l2,i+1,j  );
            mask[8] = std::make_tuple(1,l1,l2,i+1,j+1);
        }
    }
    
    // are just below the i = j diagonal?
    else if (i == j + 1)
    {
        //   B* B  C
        //   A  B  C
        //   A  B  C
        
        if (term == -1) // A
        {
            mask[0] = std::make_tuple(1,l1,l2,i-1,j-1);
            mask[1] = std::make_tuple(1,l1,l2,i-1,j  );
        }
        else if (term == 0) // B
        {
            mask[2] = std::make_tuple(s,l2,l1,i  ,j  ); // symmetry
            
            mask[3] = std::make_tuple(1,l1,l2,i  ,j-1);
            mask[4] = std::make_tuple(1,l1,l2,i  ,j  );
            mask[5] = std::make_tuple(1,l1,l2,i  ,j+1);
        }
        else if (term == +1) // C
        {
            mask[6] = std::make_tuple(1,l1,l2,i+1,j-1);
            mask[7] = std::make_tuple(1,l1,l2,i+1,j  );
            mask[8] = std::make_tuple(1,l1,l2,i+1,j+1);
        }
    }
    
    // are we right at the i = j diagonal?
    else if (i == j)
    {
        //   C* C* C
        //   B* B  C
        //   A  B  C
        
        if (term == -1) // A
        {
            mask[0] = std::make_tuple(1,l1,l2,i-1,j-1);
        }
        else if (term == 0) // B
        {
            mask[1] = std::make_tuple(s,l2,l1,i  ,j-1); // symmetry
            
            mask[3] = std::make_tuple(1,l1,l2,i  ,j-1);
            mask[4] = std::make_tuple(1,l1,l2,i  ,j  );
        }
        else if (term == +1) // C
        {
            mask[2] = std::make_tuple(s,l2,l1,i+1,j-1); // symmetry
            
            mask[5] = std::make_tuple(s,l2,l1,i+1,j  ); // symmetry
            
            mask[6] = std::make_tuple(1,l1,l2,i+1,j-1);
            mask[7] = std::make_tuple(1,l1,l2,i+1,j  );
            mask[8] = std::make_tuple(1,l1,l2,i+1,j+1);
        }
    }
    
    // there is no other option...
    else
    {
        HexException("Runtime error.");
    }
}

void Numerov2d::calc_mat (std::size_t i, int term, Complex * M) const
{
    // skip empty matrices
    if (i == 0 and term < 0)
        return;
    
    // number of blocks and their volume
    std::size_t Nang = ang_.size();
    std::size_t block_vol = i * 3;
    
    // erase output array
    std::memset(M, 0, Nang * Nang * block_vol * sizeof(Complex));
    
    // discretization scheme mask
    std::array<std::tuple<int,unsigned,unsigned,std::size_t,std::size_t>, 9> msk;
    
    // symmetry factor
    int s = ((inp_.S + inp_.Pi) % 2 == 0 ? 1 : -1);
    
    // for all angular momentum blocks
    for (std::size_t m = 0; m < Nang; m++)
    for (std::size_t n = 0; n < Nang; n++)
    {
        // get angular momenta
        unsigned l1 = ang_.state(m).first;
        unsigned l2 = ang_.state(m).second;
        unsigned l1p= ang_.state(n).first;
        unsigned l2p= ang_.state(n).second;
        
        // for all grid points in the current (i-th) grid column, skip j = 0 (zero bc)
        for (std::size_t j = 1; j <= i; j++)
        {
            // get mesh parameters for this central point (i,j)
            Complex h = rad_.grid[i] - rad_.grid[i-1];
            Complex t = rad_.grid[j] - rad_.grid[j-1];
            Complex alpha = (rad_.grid[i+1] - rad_.grid[i]) / h;
            Complex beta  = (rad_.grid[j+1] - rad_.grid[j]) / t;
            
            // get discretization scheme mask for this central point (i,j)
            mask(term, s, l1p, l2p, i, j, msk);
            
            // for all nine points in neighbourhood
            for (int ip = -1; ip <= 1; ip++)
            for (int jp = -1; jp <= 1; jp++)
            {
                // get point position
                std::size_t k = i + ip;
                std::size_t l = j + jp;
                
                // get (k,l) neighbour symmetry info ('ks' is always just 'i + term')
                int signs; unsigned l1s, l2s; std::size_t ks, ls;
                std::tie(signs,l1s,l2s,ks,ls) = msk[(ip + 1) * 3 + (jp + 1)];
                unsigned ns = ang_.index(l1p,l2p);
                
                // skip this neighbour point (k,l) if it is not coupled to the current central element (i,j) through the matrix selected by 'term'
                if (signs == 0)
                    continue;
                
                // skip this neighbour point if it lies on the (zero) boundary
                if (k == 0 or l == 0 or /* k == rad_.Npts - 1 or */ l == rad_.Npts - 1)
                    continue;
                
                // get Numerov discretization coefficients for the position (i,j,k,l)
                Complex Aik = coef_A(i, k, h, alpha, l1);
                Complex Bik = coef_B(i, k, h, alpha, l1);
                Complex Cjl = coef_A(j, l, t, beta,  l2);
                Complex Djl = coef_B(j, l, t, beta,  l2);
                
//                 std::cout << "i = " << i << ", j = " << j << ", k = " << k << ", l = " << l << "; ks = " << ks << ", ls = " << ls << std::endl;
//                 std::cout << "\th = " << h << ", t = " << t << ", alpha = " << alpha << ", beta = " << beta << std::endl;
//                 std::cout << "\tAik = " << Aik << ", Bik = " << Bik << ", Cjl = " << Cjl << ", Djl = " << Djl << std::endl;
                
                // matrix element
                Complex el = 0;
                
                // calculate two-electron part
                for (unsigned lambda = 0; lambda <= ang_.maxlambda(); lambda++) if (ang_.f(lambda, m, n) != 0.)
                {
                    double rmin = std::min(rad_.grid[k].real(), rad_.grid[l].real());
                    double rmax = std::max(rad_.grid[k].real(), rad_.grid[l].real());
                    el -= h*h*t*t*Bik*Djl * ang_.f(lambda, m, n) * special::pow_int(rmin/rmax, lambda) / rmax;
                }
                
                // calculate one-electron (diagoanl) part, if any
                if (m == n)
                {
                    // calculate potential
                    Complex Vkl = 0.;
                    if (k > 0) Vkl += 0.5 * l1 * (l1 + 1) / (rad_.grid[k].real() * rad_.grid[k].real()) - inp_.Z / rad_.grid[k].real();
                    if (l > 0) Vkl += 0.5 * l2 * (l2 + 1) / (rad_.grid[l].real() * rad_.grid[l].real()) - inp_.Z / rad_.grid[l].real();
                    
                    // update element
                    el += 0.5*h*h*Bik*Cjl + 0.5*t*t*Aik*Djl + h*h*t*t*Bik*Djl*(0.5*inp_.Etot - Vkl);
                }
                
                // update matrix element
                M[(m * Nang + ns) * block_vol + (j - 1) * 3 + (ls + 1 - j)] += double(signs) * el;
            }
        }
    }
}
