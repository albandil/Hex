/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2014                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cstdlib>
#include <iostream>

#include "arrays.h"
#include "lapack_subset.h"
#include "matrix.h"
#include "moments.h"
#include "specf.h"
#include "version.h"

int main (int argc, char *argv[])
{
    // write the program header
    std::cout << logo_raw();
    
    //
    // initialization
    //
    
    // single-electron angular momentum limit
    int maxell = 2;
    
    // Laguerre basis sizes
    iArray Nl (maxell + 1);
    Nl.fill(50);
    
    // screening constants
    rArray lambda (maxell + 1);
    lambda.fill(2.);
    
    //
    // construct Hamiltonian matricex in the bases of Laguerre functions
    //
    
    // create matrix array
    Array<RowMatrix<double>> H (maxell + 1);
    
    // for all single-electron angular momenta
    for (int ell = 0; ell <= maxell; ell++)
    {
        // construct I-matrices
        RowMatrix<double> Im2 = IntegralMoments::Im2 (Nl[ell], 2*ell+2);
        RowMatrix<double> Im1 = IntegralMoments::Im1 (Nl[ell], 2*ell+2);
        RowMatrix<double> Im1_rowsums = Im1;
        RowMatrix<double> Im1_colsums = Im1;
        
        // replace elements by partial sums
        for (int irow = 0; irow < Nl[ell]; irow++)
        {
            double sum = 0, val = 0;
            for (int icol = 0; icol < Nl[ell]; icol++)
            {
                val = Im1_rowsums (irow,icol);
                Im1_rowsums (irow,icol) = sum;
                sum += val;
            }
        }
        for (int icol = 0; icol < Nl[ell]; icol++)
        {
            double sum = 0, val = 0;
            for (int irow = 0; irow < Nl[ell]; irow++)
            {
                val = Im1_colsums (irow,icol);
                Im1_colsums (irow,icol) = sum;
                sum += val;
            }
        }
        
        // compute normalization factors
        rArray Nkl(Nl[ell] + 1), Nkl_invsqrsum(Nl[ell] + 1);
        for (int k = 1; k <= Nl[ell]; k++)
        {
            // compute norm
            Nkl[k] = sqrt(1./pochhammer_up(k, 2*ell + 2));
            
            // partial sum of norm inverse squares
            Nkl_invsqrsum[k] = (k == 1) ? pow(Nkl[k],-2) : (Nkl_invsqrsum[k-1] + pow(Nkl[k],-2));
        }
        
        // calculate matrix elements of the second derivative
        RowMatrix<double> D (Nl[ell]);
        D = (ell+1.) * ((ell+1.)*Im2 - Im1 - Im1_rowsums - Im1_colsums);
        D.transform
        (
            [ & ] (int a, int b, double & x) -> void
            {
                int k = a + 1;
                int kp = b + 1;
                
                int kmin = std::min(k,kp);
                int kmax = std::max(k,kp);
                
                x += Nkl_invsqrsum[kmin-1];
                
                x *= Nkl[k] * Nkl[kp];
                
                if (k == kp)
                    x += -0.25;
                x += 0.5 * Nkl[kmax] / Nkl[kmin];
            }
        );
        D *= lambda[ell] * lambda[ell];
        
        // compose Hamiltonian
        H[ell] = lambda[ell] * (lambda[ell]*0.5*(ell*(ell+1.))*Im2 - Im1);
        
        // multiply elements by normalization factors
        H[ell].transform
        (
            [ & ](int a, int b, double & x) -> void
            {
                int k = a + 1;
                int kp = b + 1;
                
                x *= Nkl[k] * Nkl[kp];
            }
        );
        
        // add derivative matrix
        H[ell] += 0.5 * D;
    }
    
    //
    // diagonalize the Hamiltonian
    //
    
    // energies from diagonalizations
    rArrays El (maxell + 1);
    
    std::cout << "Diagonalization of Hamiltonian in Laguerre basis:\n\n";
    
    // for all single-electron angular momenta
    for (int ell = 0; ell <= maxell; ell++)
    {
        // allocate space for eigenenergies
        El[ell].resize(Nl[ell]);
        
        // parameters for diagonalization
        char jobz = 'V'; // also eigenvectors
        char uplo = 'U'; // use (e.g.) upper triangle
        
        // prepare workspace
        rArray workspace (1);
        int worksize = -1, info = 0;
        
        // get optimal workspace size
        dsyev_ (
            &jobz,
            &uplo,
            &Nl[ell],
            H[ell].data().data(),
            &Nl[ell],
            El[ell].data(),
            workspace.data(),
            &worksize,
            &info
        );
        
        // allocate workspace
        worksize = workspace[0];
        workspace.resize(worksize);
        
        // diagonalize
        dsyev_ (
            &jobz,
            &uplo,
            &Nl[ell],
            H[ell].data().data(),
            &Nl[ell],
            El[ell].data(),
            workspace.data(),
            &worksize,
            &info
        );
        
        // write energies
        std::cout << "\tOrbital momentum   ℓ = " << ell << "\n";
        std::cout << "\tScreening constant λ = " << lambda[ell] << "\n";
        std::cout << "\tBasis size         N = " << Nl[ell] << "\n";
        std::cout << "\tEigen-energies compared with hydrogen energies:\n";
        for (int n = 1; n < Nl[ell]; n++)
        {
            if (El[ell][n-1] >= 0)
            {
                std::cout << "\t\t ... and " << Nl[ell] - n + 1
                          << " higher states with ε > 0 (max: " << El[ell].back() << ")\n";
                break;
            }
            
            std::cout << "\t\t" << format("%11.7f", El[ell][n-1]) << "\t" << -0.5/((ell+n)*(ell+n)) << "\n";
        }
//         std::cout << "Eigenvectors:\n"; H[ell].write(std::cout, "\t", "");
        std::cout << "\n";
    }
    
    return EXIT_SUCCESS;
}
