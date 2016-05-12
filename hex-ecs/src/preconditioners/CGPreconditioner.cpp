//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2016, Jakub Benda, Charles University in Prague                    //
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
#include <cstdio>

#include "hex-arrays.h"
#include "hex-itersolve.h"
#include "hex-misc.h"

#include "preconditioners.h"

const std::string CGPreconditioner::prec_name = "CG";
const std::string CGPreconditioner::prec_description = 
    "Block inversion using plain conjugate gradients. "
    "Use --tolerance option to set the termination tolerance.";

int CGPreconditioner::solve_block (int ill, const cArrayView r, cArrayView z) const
{
    // shorthands
    int Nspline_atom = rad_.bspline_atom().Nspline();
    int Nspline_proj = rad_.bspline_proj().Nspline();
    
    // prepare the block-preconditioner for run
    this->CG_init(ill);
    
    // solve using the CG solver
    ConjugateGradients < Complex, cArray, cArrayView > CG;
    CG.reset();
    CG.verbose              = false;
    CG.apply_preconditioner = [&](const cArrayView a, cArrayView b) { this->CG_prec(ill, a, b); };
    CG.matrix_multiply      = [&](const cArrayView a, cArrayView b) { this->CG_mmul(ill, a, b); };
    CG.scalar_product       = [&](const cArrayView a, const cArrayView b) { return this->CG_scalar_product(a, b); };
    CG.compute_norm         = [&](const cArrayView a) { return this->CG_compute_norm(a); };
    CG.axby                 = [&](Complex a, cArrayView x, Complex b, const cArrayView y) { this->CG_axby_operation(a, x, b, y); };
    CG.new_array            = [&](std::size_t n, std::string name) { return cArray(n); };
    CG.constrain            = [&](cArrayView r) { this->CG_constrain(r); };
    int n = CG.solve(r, z, cmd_.prec_itertol, 0, Nspline_atom * Nspline_proj);
    
    // release block-preconditioner block-specific data
    this->CG_exit(ill);
    
    return n;
}

void CGPreconditioner::precondition (BlockArray<Complex> const & r, BlockArray<Complex> & z) const
{
    if (cmd_.ssor > 0)
    {
        // working arrays
        cArrays y (r.size()), x (r.size());
        
        // forward SOR
        for (int ill = 0; ill < (int)ang_.states().size(); ill++)
        {
            // start with right-hand side
            y[ill] = r[ill];
            
            // subtract lower block diagonals for ill-th block row
            for (int illp = 0; illp < ill; illp++)
            for (int lambda = 0; lambda <= rad_.maxlambda(); lambda++)
            if (ang_.f(ill, illp, lambda) != 0)
                rad_.R_tr_dia(lambda).dot(-ang_.f(ill, illp, lambda), y[illp], 1., y[ill], true);
            
            // use (preconditioned) conjugate gradients to invert a diagonal block
            x[ill].resize(y[ill].size());
            n_[ill] = solve_block(ill, cmd_.ssor * y[ill], x[ill]);
        }
        
        // normalize
        for (int ill = 0; ill < (int)ang_.states().size(); ill++)
            this->CG_mmul(ill, (2. - cmd_.ssor) / cmd_.ssor * x[ill], y[ill]);
        
        // backward SOR
        for (int ill = (int)ang_.states().size() - 1; ill >= 0; ill--)
        {
            // subtract upper block diagonals for ill-th block row
            for (int illp = ill + 1; illp < (int)ang_.states().size(); illp++)
            for (int lambda = 0; lambda <= rad_.maxlambda(); lambda++)
            if (ang_.f(ill, illp, lambda) != 0)
                rad_.R_tr_dia(lambda).dot(-ang_.f(ill, illp, lambda), y[illp], 1., y[ill], true);
            
            // use (preconditioned) conjugate gradients to invert a diagonal block
            n_[ill] += solve_block(ill, cmd_.ssor * y[ill], z[ill]);
        }
    }
    else
    {
#ifndef DISABLE_PARALLEL_PRECONDITION
        // NOTE : If the BLAS is multi-threaded, this will result in nested parallelism.
        //        Some combinations of system kernel and OpenMP implementation lead to system
        //        freeze. It can be avoided by using a serial BLAS.
        # pragma omp parallel for schedule (dynamic, 1) if (cmd_.parallel_precondition && cmd_.groupsize == 1)
#endif
        for (int ill = 0; ill < (int)ang_.states().size(); ill++) if (par_.isMyGroupWork(ill))
        {
            // load segment, if necessary
            if (cmd_.outofcore)
            {
                const_cast<BlockArray<Complex>&>(r).hdfload(ill);
                z.hdfload(ill);
            }

            // ivert diagonal block
            n_[ill] = solve_block(ill, r[ill], z[ill]);
            
            // unload segment
            if (cmd_.outofcore)
            {
                const_cast<BlockArray<Complex>&>(r)[ill].drop();
                
                z.hdfsave(ill);
                z[ill].drop();
            }
        }
    }
    
    // broadcast inner preconditioner iterations
    par_.sync_m(n_.data(), 1, ang_.states().size());
    par_.bcast_g(par_.igroup(), 0, n_.data(), ang_.states().size());
    
    // inner preconditioner info (max and avg number of iterations)
    std::cout << " | ";
    std::cout << std::setw(5) << (*std::min_element(n_.begin(), n_.end()));
    std::cout << std::setw(5) << (*std::max_element(n_.begin(), n_.end()));
    std::cout << std::setw(5) << format("%g", std::accumulate(n_.begin(), n_.end(), 0) / float(n_.size()));
}

void CGPreconditioner::CG_init (int iblock) const
{
    if (cmd_.outofcore and cmd_.wholematrix)
        dia_blocks_[iblock].hdfload();
}

void CGPreconditioner::CG_mmul (int iblock, const cArrayView p, cArrayView q) const
{
    if (cmd_.lightweight_full)
        HexException("Preconditioner %s is not compatible with the option --lightweight-full.", this->name().c_str());
        
    dia_blocks_[iblock].dot(1., p, 0., q, true);
}

void CGPreconditioner::CG_prec (int iblock, const cArrayView r, cArrayView z) const
{
    z = r;
}

void CGPreconditioner::CG_exit (int iblock) const
{
    if (cmd_.outofcore and cmd_.wholematrix)
        dia_blocks_[iblock].drop();
}

void CGPreconditioner::finish ()
{
    dia_blocks_.resize(0);
    n_.fill(-1);
    NoPreconditioner::finish();
}

double CGPreconditioner::CG_compute_norm (const cArrayView a) const
{
    // compute norm (part will be computed by every process in group, result will be synchronized)
        
    // calculate number of elements every group's process will sum
    std::size_t N = (a.size() + par_.groupsize() - 1) / par_.groupsize();
    
    // calculate part of the norm
    double norm2 = 0;
    for (std::size_t i = par_.igroupproc() * N; i < (par_.igroupproc() + 1) * N and i < a.size(); i++)
        norm2 += sqrabs(a[i]);
    
    // sum across group and return result
    par_.syncsum_g(&norm2, 1);
    return std::sqrt(norm2);
}
    
Complex CGPreconditioner::CG_scalar_product (const cArrayView a, const cArrayView b) const
{
    // compute scalar product (part will be computed by every process in group, result will be synchronized)
    
    assert(a.size() == b.size());
    
    // calculate number of elements every group's process will sum
    std::size_t N = (a.size() + par_.groupsize() - 1) / par_.groupsize();
    
    // calculate part of the norm
    Complex prod = 0;
    for (std::size_t i = par_.igroupproc() * N; i < (par_.igroupproc() + 1) * N and i < a.size(); i++)
        prod += a[i] * b[i];
    
    // sum across group and return result
    par_.syncsum_g(&prod, 1);
    return prod;
}

void CGPreconditioner::CG_axby_operation (Complex a, cArrayView x, Complex b, const cArrayView y) const
{
    // a*x + b*y operation
    
    assert(x.size() == y.size());
    
    // calculate number of elements each node should process
    std::size_t N = (x.size() + par_.groupsize() - 1) / par_.groupsize();
    
    // calculate linear combination
    for (std::size_t i = par_.igroupproc() * N; i < (par_.igroupproc() + 1) * N and i < x.size(); i++)
        x[i] = a * x[i] + b * y[i];
    
    // synchronize the segments
    for (int inode = 0; inode < par_.groupsize(); inode++)
    {
        // calculate 'begin' and 'end' of inode's data
        std::size_t begin = inode * N;
        std::size_t end = std::min(x.size(), (inode + 1) * N);
        
        // broadcast the data
        par_.bcast_g
        (
            par_.igroup(),          // broadcast within this node's group
            inode,                  // use inode's data
            &x[0] + begin,          // data pointer
            end - begin             // number of elements to synchronize
        );
    }
}

void CGPreconditioner::CG_constrain (cArrayView r) const
{
    // leave the resudual as it is
}
