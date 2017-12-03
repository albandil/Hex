#include <iostream>

#include "hex-blas.h"
#include "hex-densematrix.h"
#include "hex-itersolve.h"
#include "hex-vtkfile.h"

#include "bspline.h"
#include "luft.h"
#include "parallel.h"
#include "radial.h"

// #define CG_JACOBI
// #define CG_SSOR
// #define CG_ILUT
// #define CG_DILU
// #define CG_ILUP
#define CG_KPA

int main (int argc, char* argv[])
{
    Complex Etot = 1;
    
    rArray rknots = concatenate(rArray{0,0,0}, linspace(0., 20., 41));
    rArray cknots = linspace(20., 30., 21);
    
    #include "createCmdParallel.h"
    #include "createBspline.h"
    #include "createRadialIntegrals.h"
    #include "createMainMatrix.h"
    #include "diagonalize1D.h"
    //#include "diagonalize2D.h"
    #include "createFields.h"
    //#include "createLU.h"
    
    bool stationary = false;
    
    Real droptol = 1e-3;
    int level = 1;
    bool modified = false;
    
    //#include "createILUT.h"
    //#include "createDILU.h"
    #include "createILUP.h"
    
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
        
        #elif defined ( CG_KPA )
        
            // Kronecker product approximation
            z = r;
            
            cArray w (N*N*N);
            
            RowMatrixView<Complex> A (N, N, KPAA.data());
            RowMatrixView<Complex> B (N, N, KPAB.data());
            
            for (int turn = 1; turn <= 3; turn++)
            {
                for (int i = 0; i < N; i++)
                {
                    RowMatrixView<Complex> Z (N, N, z.slice(i*N*N, (i+1)*N*N));
                    RowMatrixView<Complex> W (N, N, w.slice(i*N*N, (i+1)*N*N));
                    
                    blas::gemm(1., Z, A, 0., W);
                }
                
                transpose(w, z, N, N, N);
            }
            
            for (int i = 0; i < N; i++) 
            for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                z[(i * N + j) * N + k] /= Etot - e[i] - e[j] - e[k];
            
            for (int turn = 1; turn <= 3; turn++)
            {
                for (int i = 0; i < N; i++)
                {
                    RowMatrixView<Complex> Z (N, N, z.slice(i*N*N, (i+1)*N*N));
                    RowMatrixView<Complex> W (N, N, w.slice(i*N*N, (i+1)*N*N));
                    
                    blas::gemm(1., Z, B, 0., W);
                }
                
                transpose(w, z, N, N, N);
            }
        
        #else
        
            // no preconditioner
            z = r;
        
        #endif
    };
    
    #include "runCG.h"
    #include "writeResults.h"
}
