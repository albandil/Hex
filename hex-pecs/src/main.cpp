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

#include <cstdlib>
#include <iostream>

#include "hex-arrays.h"
#include "hex-cmdline.h"
#include "hex-special.h"
#include "hex-version.h"

#include "ang.h"
#include "io.h"
#include "matops.h"
#include "numerov.h"
#include "os.h"
#include "radial.h"

namespace debug
{
    const bool enabled = true;
    
    void write_banded (std::ostream & os, std::size_t rows, std::size_t cols, std::size_t Ndiag, Complex const * M)
    {
        if (not enabled)
            return;
        
        for (unsigned i = 0; i < rows; i++)
        {
            os << "    [ ";
            for (unsigned j = 0; j < cols; j++)
            {
                if (std::max(i,j) - std::min(i,j) <= Ndiag)
                    os << M[i * (2*Ndiag + 1) + j + Ndiag - i].real() << " ";
                else
                    os << "0 ";
            }
            os << "]" << std::endl;
        }
    }
    
    void write_dense_rows (std::ostream & os, std::size_t rows, std::size_t cols, Complex const * M)
    {
        if (not enabled)
            return;
        
        for (unsigned i = 0; i < rows; i++)
        {
            os << "    [ ";
            for (unsigned j = 0; j < cols; j++)
                os << M[i * cols + j].real() << " ";
            os << "]" << std::endl;
        }
    }
    
    void write_dense_cols (std::ostream & os, std::size_t rows, std::size_t cols, Complex const * M)
    {
        if (not enabled)
            return;
        
        for (unsigned i = 0; i < rows; i++)
        {
            os << "    [ ";
            for (unsigned j = 0; j < cols; j++)
                os << M[i + j * rows].real() << " ";
            os << "]" << std::endl;
        }
    }
}

int main (int argc, char * argv[])
{
    //
    // Program initialization
    //
    
        // display logo
        std::cout << logo(" ") << std::endl;
        std::cout << "=== Propagating exterior complex scaling ===" << std::endl << std::endl;
        
        // echo command line
        std::cout << "Command line used" << std::endl;
        std::cout << "\t";
        for (int iarg = 0; iarg < argc; iarg++)
            std::cout << argv[iarg] << " ";
        std::cout << std::endl << std::endl;
        
        // turn off GSL exceptions
        gsl_set_error_handler_off();
        
        // disable buffering of the standard output (-> immediate logging)
        std::setvbuf(stdout, nullptr, _IONBF, 0);
    
    //
    // Preparations.
    //
    
        // analyze command line parameters
        CommandLine cmd (argc, argv);
        
        // read input file
        InputFile inp (cmd);
        
        // construct angular basis
        AngularBasis ang (inp);
        
        // construct radial basis
        RadialBasis rad (inp);
        
        // construct the discretization class
        Numerov2d num (inp, ang, rad);
        
        // shorthands
        std::size_t Nang = ang.size();
        std::size_t Npts = rad.Npts;
        
        // info
        std::cout << "Number of grid points: " << Npts << std::endl << std::endl;
        std::cout << "Matrix size: up to " << Nang * Npts - 2 << std::endl << std::endl;
        
        // allocate all banded matrices (row-major blocks, row-major bands)
        Complex * A = new Complex [Nang * Npts * 3 * Nang];
        Complex * B = new Complex [Nang * Npts * 3 * Nang];
        Complex * C = new Complex [Nang * Npts * 3 * Nang];
        
        // allocate all dense matrices (column-major)
        Complex * D   = new Complex [Nang * Npts * Nang * Npts];
        Complex * BLU = new Complex [Nang * Npts * Nang * Npts];
        Complex * psi = new Complex [Nang * Npts *        Npts];
        
        // allocate all vectors
        Complex * E = new Complex [Nang * Npts];
        Complex * F = new Complex [Nang * Npts];
        Complex * V = new Complex [Nang * Npts];
        
        // allocate pivots
        int * pivots = new int [Nang * Npts];
    
    //
    // Forward pass : construction of propagation matrices.
    //
    
    if (cmd.forward_grid)
    {
        for (std::size_t icol = 1; icol < Npts - 1; icol++)
        {
            std::cout << "Forward grid: column " << icol << " / " << Npts - 2 << std::endl;
            
            // create output directory
            os::mkdirp(format("%d", icol));
            
            // calculate the Numerov matrices
            std::cout << "  - calculate A-matrix" << std::endl;
            num.A(icol, A); // dim: Nang*i x Nang*(i-1)
//             debug::write_banded(std::cout, icol, icol - 1, 1, A);
            
            std::cout << "  - calculate B-matrix" << std::endl;
            num.B(icol, B); // dim: Nang*i x Nang*i
//             debug::write_banded(std::cout, icol, icol, 1, B);
            
            std::cout << "  - calculate C-matrix" << std::endl;
            num.C(icol, C); // dim: Nang*i x Nang*(i+1)
//             debug::write_banded(std::cout, icol, icol + 1, 1, C);
            
            // check that [iknot]/BLU.bin exists
            if (not matops::load(BLU, Nang * icol * Nang * icol, format("%d/BLU.bin", icol)) or
                not matops::load(pivots, Nang * icol, format("%d/pivots.bin", icol)))
            {
                // A * D -> M, dim: Nang*i x Nang*i
                std::cout << "  - multiply A * D -> M" << std::endl;
                matops::blockband_mul_dense(Nang, icol, icol - 1, icol, 1, A, D, BLU);
//                 debug::write_dense_rows(std::cout, icol, icol, BLU);
                
                // M + B -> M
                std::cout << "  - add (A * D) + B -> M" << std::endl;
                matops::dense_add_blockband(Nang, icol, 1, BLU, B);
//                 debug::write_dense_rows(std::cout, icol, icol, BLU);
                
                // LU-factorize M -> BLU
                std::cout << "  - LU-factorization of M" << std::endl;
                matops::dense_LU_factor(Nang * icol, BLU, pivots);
//                 debug::write_dense_rows(std::cout, icol, icol, BLU);
                
                // save factorization to disk
                std::cout << "  - save LU to disk" << std::endl;
                matops::save(BLU, Nang * icol * Nang * icol, format("%d/BLU.bin", icol));
                matops::save(pivots, Nang * icol, format("%d/pivots.bin", icol));
            }
            
            // update D, = -B^{-1} C
            std::cout << "  - solve B D = C for D" << std::endl;
            matops::dense_LU_solve_blockband(Nang, icol, icol + 1, 1, BLU, pivots, C, D, V);
//             debug::write_dense_cols(std::cout, icol, icol + 1, D);
            
            std::cout << "  - flip sign of D" << std::endl;
            matops::flip_sign(Nang * icol * Nang * (icol + 1), D);
//             debug::write_dense_cols(std::cout, icol, icol + 1, D);
            
            std::cout << std::endl;
            
//             if (icol == 10)
//                 std::exit(0);
        }
    }
    
    //
    // Forward pass : construction of constant vectors.
    //
    
    if (cmd.forward_states)
    {
        for (std::size_t icol = 1; icol < Npts - 1; icol++)
        {
            std::cout << "Forward states: column " << icol << " / " << Npts - 2 << std::endl;
            
            // calculate the Numerov matrices
            std::cout << "  - calculate A-matrix" << std::endl;
            num.A(icol, A);
//             debug::write_banded(std::cout, icol, icol - 1, 1, A);
            
            // load inverted B
            if (not matops::load(BLU, Nang * icol * Nang * icol, format("%d/BLU.bin", icol)) or
                not matops::load(pivots, Nang * icol, format("%d/pivots.bin", icol)))
                HexException("Missing precomputed LU-factorization for grid point %ld.", icol);
            
            // for all initial states
            std::cout << "  - prepare constant vectors" << std::endl;
            for (std::size_t istate = 0; istate < inp.istates.size(); istate++)
            {
                // calculate the constant vector
                std::cout << "  - calculate F-vector" << std::endl;
                num.F(icol, F, istate);
//                 debug::write_dense_rows(std::cout, 1, Nang * icol, F);
                matops::save(F, Nang * icol, format("%d/F-%d.bin", icol, istate));
                
                // A * E -> V
                std::cout << "  - multiply A * E -> V" << std::endl;
                matops::blockband_mul_vector(Nang, icol, icol - 1, 1, A, E, V);
//                 debug::write_dense_rows(std::cout, 1, Nang * icol, V);
                
                // F - (A * E)  ->  V
                std::cout << "  - subtract F - V -> V" << std::endl;
                matops::subtract(Nang * icol, F, V, V);
//                 debug::write_dense_rows(std::cout, 1, Nang * icol, V);
                
                // B^{-1} (F - (A * E))  ->  E
                std::cout << "  - solve B^{-1} V -> E" << std::endl;
                matops::dense_LU_solve_vector(Nang * icol, BLU, pivots, V, E);
//                 debug::write_dense_rows(std::cout, 1, Nang * icol, E);
                
                // save to disk
                matops::save(E, Nang * icol, format("%d/E-%d.bin", icol, istate));
            }
            
            std::cout << std::endl;
            
//             if (icol == 10)
//                 std::exit(0);
        }
    }
    
    //
    // Backward pass : solution of constant vectors.
    //
    
    if (cmd.backward)
    {
        for (std::size_t icol = Npts - 2; icol >= 1; icol--)
        {
            std::cout << "Backward states: column " << icol << " / " << Npts - 2 << std::endl;
            
            // calculate the Numerov matrices
            num.C(icol, C);
            
            // load inverted B
            if (not matops::load(BLU, Nang * icol * Nang * icol, format("%d/BLU.bin", icol)) or
                not matops::load(pivots, Nang * icol, format("%d/pivots.bin", icol)))
                HexException("Missing precomputed LU-factorization for grid point %ld.", icol);
            
            // update D, = -B^{-1} C
            std::cout << "  - solve B D = C for D" << std::endl;
            matops::dense_LU_solve_blockband(Nang, icol, icol + 1, 1, BLU, pivots, C, D, V);
//             debug::write_dense_cols(std::cout, Nang * icol, Nang * (icol + 1), D);
            
            std::cout << "  - flip sign of D" << std::endl;
            matops::flip_sign(Nang * icol * Nang * (icol + 1), D);
            
            // for all initial states
            std::cout << "  - propagate states" << std::endl;
            for (std::size_t istate = 0; istate < inp.istates.size(); istate++)
            {
                // load E
                matops::load(E, Nang * icol, format("%d/E-%d.bin", icol, istate));
                
                // load or initialize psi
                if (icol == rad.Npts - 1)
                    std::memset(psi, 0, Nang * (icol + 1) * sizeof(Complex));
                else
                    matops::load(psi, Nang * (icol + 1), format("%d/psi-%d.bin", icol + 1, istate));
                
                // D * psi  ->  V
                std::cout << "  - multiply D * psi -> V" << std::endl;
                matops::dense_mul_vector(Nang * icol, Nang * (icol + 1), D, psi, V);
//                 std::cout << "V = " << cArrayView(Nang * (icol + 1), V) << std::endl;
                
                // (D * psi) + E  ->  psi
                std::cout << "  - add V + E -> psi" << std::endl;
                matops::sum(Nang * icol, V, E, psi);
//                 std::cout << "psi = " << cArrayView(Nang * icol, psi) << std::endl;
                
                // save psi
                matops::save(psi, Nang * icol, format("%d/psi-%d.bin", icol, istate));
            }
            
            std::cout << std::endl;
        }
    }
    
    //
    // Extract the amplitudes.
    //
    
        // TODO
    
    
    return EXIT_SUCCESS;
}
