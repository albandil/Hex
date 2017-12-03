#include <iostream>

#include "hex-blas.h"
#include "hex-densematrix.h"
#include "hex-itersolve.h"
#include "hex-vtkfile.h"

#include "bspline.h"
#include "ft.h"
#include "luft.h"
#include "parallel.h"
#include "radial.h"

// #define CG_JACOBI
// #define CG_SSOR
// #define CG_ILUT
// #define CG_DILU
// #define CG_ILUP
// #define CG_BLOCK_JACOBI
// #define CG_CIRCULANT
#define CG_KPA

int main (int argc, char* argv[])
{
    Complex Etot = -0.1;
    
    rArray rknots = concatenate(rArray{0,0,0}, linspace(0., 60., 61));
    rArray cknots = linspace(60., 100., 41);
    /*rArray rknots = concatenate(rArray{0,0,0}, linspace(0., 10., 21));
    rArray cknots = linspace(10., 20., 21);*/
    
    #include "createCmdParallel.h"
    #include "createBspline.h"
    #include "createRadialIntegrals.h"
    #include "createMainMatrix.h"
    #include "diagonalize1D.h"
    //#include "diagonalize2D.h"
    #include "createFields.h"
    #include "createLU.h"
    
    bool stationary = false;
    
    Real droptol = 1e-3;
    int level = 5;
    bool modified = false;
    
    #include "createILUT.h"
    #include "createDILU.h"
    #include "createILUP.h"
    #include "createBlockLU.h"
    
    #include "setupCG.h"
    {
        #if defined ( CG_JACOBI )
        
            // diagonal scaling
            z = r / A.diag();
        
        #elif defined ( CG_SSOR )
        
            // symmetric successive over-relaxation
            Real omega = 0.8;
            
            z = A.lowerSolve(r, omega);
            z = omega * (2 - omega) * z * A.diag();
            z = A.upperSolve(z, omega);
        
        #elif defined ( CG_ILUT )
        
            // incomplete LU decomposition (drop tolerance)
            ilut->solve(r, z, 1);
        
        #elif defined ( CG_DILU )
            
            // diagonal incomplete LU decomposition
            z = dilu.lowerSolve(r);
            z = z * A.diag();
            z = dilu.upperSolve(z);
        
        #elif defined ( CG_ILUP )
        
            // incomplete LU decomposition (pattern)
            z = ilupL.lowerSolve(r);
            z = z / ilupd;
            z = ilupU.upperSolve(z);
        
        #elif defined ( CG_BLOCK_JACOBI )
            
            // block-Jacobi preconditioner
            ColMatrixView<Complex> Z (N, N, z);
            ColMatrixView<Complex> R (N, N, r);
            
            for (int i = 0; i < N; i++)
                blocklu[i]->solve(R.col(i), Z.col(i), 1);
        
        #elif defined ( CG_CIRCULANT )
            
            // circulant matrix approximation
            z = r;
            
            cArray s (N), st (N), d (N), dt (N), w (N*N);
            
            for (int i = 0; i <= order; i++)  d[i] = d[(N - i) % N] = bD(20,20+i);
            for (int i = 0; i <= order; i++)  s[i] = s[(N - i) % N] = bS(20,20+i);
            
            DFT(+1, d, dt, 1, N);
            DFT(+1, s, st, 1, N);
            
            DFT(+1, z, w, N, N);  transpose(w, z, N, N);
            DFT(+1, z, w, N, N);  transpose(w, z, N, N);
            
            z /= (dt ^ st) + (st ^ dt);
            
            DFT(-1, z, w, N, N);  transpose(w, z, N, N);
            DFT(-1, z, w, N, N);  transpose(w, z, N, N);
        
        #elif defined ( CG_KPA )
        
            // Kronecker product approximation
            z = r;
            
            cArray w (N*N);
            
            RowMatrixView<Complex> Z (N, N, z);
            RowMatrixView<Complex> W (N, N, w);
            
            for (int turn = 1; turn <= 2; turn++)
            {
                blas::gemm(1., Z, V, 0., W);
                transpose(w, z, N, N);
            }
            
            for (int i = 0; i < N; i++) 
            for (int j = 0; j < N; j++)
                Z(i,j) /= Etot - e[i] - e[j];
            
            for (int turn = 1; turn <= 2; turn++)
            {
                blas::gemm(1., Z, Vt, 0., W);
                transpose(w, z, N, N);
            }
        
        #else
        
            // no preconditioner
            z = r;
        
        #endif
    };
    
    #include "runCG.h"
    #include "writeResults.h"
}
