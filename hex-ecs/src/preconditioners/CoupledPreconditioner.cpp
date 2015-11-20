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

#ifdef WITH_MUMPS

#include "preconditioners.h"

#include <mumps/zmumps_c.h>

#define ICNTL(x) icntl[(x)-1]

const std::string CoupledPreconditioner::prec_name = "coupled";
const std::string CoupledPreconditioner::prec_description = "Coupled solver that uses MUMPS OOC.";

void CoupledPreconditioner::setup ()
{
    // update parent
    NoPreconditioner::setup();
    
    // some useful constants
    std::size_t order = bspline_atom_.order();
    std::size_t Nspline_atom = bspline_atom_.Nspline();
    std::size_t Nspline_proj = bspline_proj_.Nspline();
    std::size_t Nang = ang_.states().size(), Nchunk = Nspline_atom * Nspline_proj;
    
    // number of nonzero elements of atom/projectile basis overlap matrix
    std::size_t nz_atom = Nspline_atom * (1 + 2 * order) - order * (order + 1);
    std::size_t nz_proj = Nspline_proj * (1 + 2 * order) - order * (order + 1);
    
    // number of nonzero elements of single angular block
    std::size_t nz_block = nz_atom * nz_proj;
    
    // calculate participating blocks
    std::size_t nblocks = 0;
    for (unsigned ill  = 0; ill  < Nang; ill ++)
    for (unsigned illp = 0; illp < Nang; illp++)
    {
        for (int lambda = 0; lambda <= std::min(cmd_.coupling_limit, rad_.maxlambda()); lambda++)
        if (ang_.f(ill, illp, lambda) != 0)
        {
            nblocks++;
            break;
        }
    }
    
    // number of nonzero elements of the full matrix
    std::size_t nz_matrix = nblocks * nz_block;
    
    // number of elements on the main diagonal of the full matrix
    std::size_t nz_diagonal = Nang * Nspline_atom * Nspline_proj;
    
    // number of nonzero elements on and above the diagonal of the full matrix
    std::size_t nz = (nz_matrix - nz_diagonal) / 2 + nz_diagonal;
    
    // print some diagnostic information
    std::cout << "Setting up coupled preconditioner" << std::endl;
    std::cout << "\t- coupling: " << std::min(cmd_.coupling_limit,rad_.maxlambda()) << " of " << rad_.maxlambda() << " multipole moments" << std::endl;
    std::cout << "\t- blocks considered: " << nblocks << " of " << Nang * Nang << std::endl;
    std::cout << "\t- non-zero elements in upper triangle: " << nz << std::endl;
    std::cout << "\t- precomputing chosen subset of the hamiltonian matrix ... " << std::flush;
    
    // the block super-matrix in a single COO matrix
    I.resize(nz);
    J.resize(nz);
    A.resize(nz);
    std::size_t pos = 0;
    
    // fill the matrix
    for (unsigned ill  = 0; ill  < Nang; ill ++)
    for (unsigned illp = ill; illp < Nang; illp++)
    for (int i = 0; i < (int)Nspline_atom; i++)
    for (int j = 0; j < (int)Nspline_proj; j++)
    for (int k = std::max<int>(0, i - order); k <= std::min<int>(i + order, (int)Nspline_atom - 1); k++)
    for (int l = std::max<int>(0, j - order); l <= std::min<int>(j + order, (int)Nspline_proj - 1); l++)
    {
        // compose full-matrix row and column index
        std::size_t irow = (ill  * Nspline_atom + i) * Nspline_proj + j;
        std::size_t icol = (illp * Nspline_atom + k) * Nspline_proj + l;
        
        // skip lower triangle
        if (irow > icol)
            continue;
        
        Complex elem = 0;
        
        // add one-electron contribution
        if (ill == illp)
        {
            int l1 = ang_.states()[ill].first;
            int l2 = ang_.states()[ill].second;
            elem = E_ * rad_.S_atom()(i,k) * rad_.S_proj()(j,l)
                    - (0.5 * rad_.D_atom()(i,k) + 0.5 * l1 * (l1 + 1) * rad_.Mm2_atom()(i,k) - rad_.Mm1_tr_atom()(i,k)) * rad_.S_proj()(j,l)
                    - (0.5 * rad_.D_proj()(j,l) + 0.5 * l2 * (l2 + 1) * rad_.Mm2_proj()(j,l) - rad_.Mm1_tr_proj()(j,l)) * rad_.S_atom()(i,k);
        }
        
        // add two-electron contribution
        for (int lambda = 0; lambda <= std::min(cmd_.coupling_limit, rad_.maxlambda()); lambda++)
        {
            double f = ang_.f(ill, illp, lambda);
            if (f != 0)
                elem -= f * rad_.computeR(lambda,i,j,k,l);
        }
        
        // insert non-zero element
        if (elem != 0.)
        {
            I[pos] = irow + 1;
            J[pos] = icol + 1;
            A[pos] = elem;
            pos++;
        }
    }
    std::cout << "ok" << std::endl;
    
    // initialize
    settings.sym = 2;
    settings.par = 1;
    settings.job = -1;
    zmumps_c(&settings);
    
    // analyze
    settings.job = 1;
    settings.ICNTL(1) = 0; // errors to STDOUT (default: 6)
    settings.ICNTL(2) = 0; // diagnostics to /dev/null
    settings.ICNTL(3) = 0; // global info to STDOUT (default: 6)
    settings.ICNTL(4) = 0; // verbosity level (default: 2)
    settings.ICNTL(5) = 0; // COO format
    settings.ICNTL(22) = 1; // OOC factorization
    std::strcpy(settings.ooc_tmpdir, ".");
    std::strcpy(settings.ooc_prefix, "ooc-");
    settings.n = Nang * Nchunk;
    settings.nz = nz;
    settings.irn = I.data();
    settings.jcn = J.data();
    settings.a = reinterpret_cast<mumps_double_complex*>(A.data());
    std::cout << "\t- analyze the hamiltonian using MUMPS ... " << std::flush;
    zmumps_c(&settings);
    std::cout << "ok" << std::endl << std::endl;
    
    // resize also work vector
    X.resize(Nang * Nchunk);
    settings.nrhs = 1;
    settings.lrhs = Nang * Nchunk;
    settings.rhs = reinterpret_cast<mumps_double_complex*>(X.data());
}

void CoupledPreconditioner::update (double E)
{
    // store energy change
    double dE = E - E_;
    
    // update parent
    NoPreconditioner::update(E);
    
    // useful variables
    std::size_t Nspline_atom = bspline_atom_.Nspline();
    std::size_t Nspline_proj = bspline_proj_.Nspline();
    
    // update factorization
    if (settings.job < 2 or not cmd_.noluupdate)
    {
        // update matrix
        std::cout << "\tUpdate matrix (energy change from " << 2 * (E_ - dE) << " to " << 2 * E_ << " Ry) ... " << std::flush;
        for (std::size_t pos = 0; pos < A.size(); pos++)
        {
            // split row multi-index
            std::size_t ill = I[pos] - 1;
            std::size_t j = ill % Nspline_proj; ill /= Nspline_proj;
            std::size_t i = ill % Nspline_atom; ill /= Nspline_atom;
            
            // split column multi-index
            std::size_t illp = J[pos] - 1;
            std::size_t l = illp % Nspline_proj; illp /= Nspline_proj;
            std::size_t k = illp % Nspline_atom; illp /= Nspline_atom;
            
            // update diagonal blocks to use the new energy
            if (ill == illp)
                A[pos] += dE * rad_.S_atom()(i,k) * rad_.S_proj()(j,l);
        }
        std::cout << "ok" << std::endl;
        
        // factorize updated matrix
        std::cout << "\tFactorize the hamiltonian using MUMPS (out of core) ... " << std::flush;
        Timer t;
        settings.job = 2;
        zmumps_c(&settings);
        std::cout << "finished after " << t.nice_time() << std::endl;
    }
}

void CoupledPreconditioner::precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const
{
    // some useful constants
    std::size_t Nang = r.size(), Nchunk = r[0].size();
    
    // convert block array to monolithic array
    for (unsigned ill = 0; ill < Nang; ill++)
    for (unsigned i = 0; i < Nchunk; i++)
        X[ill * Nchunk + i] = r[ill][i];
    
    // solve
    settings.job = 3;
    zmumps_c(&settings);
    
    // copy solution to result
    for (unsigned ill = 0; ill < Nang; ill++)
    for (unsigned i = 0; i < Nchunk; i++)
        z[ill][i] = X[ill * Nchunk + i];
}

void CoupledPreconditioner::finish ()
{
    // destroy MUMPS instance
    settings.job = -2;
    zmumps_c(&settings);
    
    // release the arrays
    I.drop();
    J.drop();
    A.drop();
}

#endif // WITH_MUMPS
