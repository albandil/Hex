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

#include <iostream>

#include <gsl/gsl_sf.h>

#include "basis.h"
#include "matrix.h"
#include "misc.h"
#include "moments.h"
#include "specf.h"
#include "symbolic.h"

LaguerreBasis::LaguerreBasis (int maxell, const iArrayView Nl, const ArrayView<symbolic::rational> rat_lambda)
    : maxell_(maxell), rat_lambda_(rat_lambda)
{
    // sanity check
    assert (rat_lambda_.size() == (unsigned)(maxell + 1));
    
    // evaluate lambdas
    lambda_.resize(maxell + 1);
    for (int ell = 0; ell <= maxell; ell++)
        lambda_[ell] = symbolic::double_approx(rat_lambda_[ell]);
    
    //
    // construct Hamiltonian matricex in the bases of Laguerre functions
    //
    
    // resize array of matrices
    expansions_.resize(maxell + 1);
    
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
            Nkl[k] = sqrt(1./gsl_sf_poch(k, 2*ell + 2));
            
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
        D *= lambda_[ell] * lambda_[ell];
        
        // compose Hamiltonian
        expansions_[ell] = lambda_[ell] * (lambda_[ell]*0.5*(ell*(ell+1.))*Im2 - Im1);
        
        // multiply elements by normalization factors
        expansions_[ell].transform
        (
            [ & ](int a, int b, double & x) -> void
            {
                int k = a + 1;
                int kp = b + 1;
                
                x *= Nkl[k] * Nkl[kp];
            }
        );
        
        // add derivative matrix
        expansions_[ell] += 0.5 * D;
    }
    
    //
    // diagonalize the Hamiltonian
    //
    
    // resize energies from diagonalizations
    energies_.resize(maxell + 1);
    
    std::cout << "Diagonalization of Hamiltonian in Laguerre basis:\n\n";
    
    // for all single-electron angular momenta
    for (int ell = 0; ell <= maxell; ell++)
    {
        // allocate space for eigenenergies
        energies_[ell].resize(Nl[ell]);
        
        // parameters for diagonalization
        char jobz = 'V'; // also eigenvectors
        char uplo = 'U'; // use (e.g.) upper triangle
        
        // prepare workspace
        rArray workspace (1);
        int worksize = -1, info = 0, n = Nl[ell];
        
        // get optimal workspace size
        dsyev_ (
            &jobz,
            &uplo,
            &n,
            expansions_[ell].data().data(),
            &n,
            energies_[ell].data(),
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
            &n,
            expansions_[ell].data().data(),
            &n,
            energies_[ell].data(),
            workspace.data(),
            &worksize,
            &info
        );
        
        // write energies
        std::cout << "\tOrbital momentum   ℓ = " << ell << "\n";
        std::cout << "\tScreening constant λ = " << lambda_[ell] << "\n";
        std::cout << "\tBasis size         N = " << Nl[ell] << "\n";
        std::cout << "\tEigen-energies compared with hydrogen energies:\n";
        for (int n = 1; n < Nl[ell]; n++)
        {
            if (energies_[ell][n-1] >= 0)
            {
                std::cout << "\t\t ... and " << Nl[ell] - n + 1
                          << " higher states with ε > 0 (max: " << energies_[ell].back() << ")\n";
                break;
            }
            
            std::cout << "\t\t" << format("%11.7f", energies_[ell][n-1]) << "\t" << -0.5/((ell+n)*(ell+n)) << "\n";
        }
//         std::cout << "Eigenvectors:\n"; expansions_[ell].write(std::cout, "\t", "");
        std::cout << "\n";
    }
}

double LaguerreBasis::basestate (int ell, int k, double r) const
{
    if (ell < 0 or ell > maxell_)
        throw exception ("[LaguerreBasis::size] Can't access basis %d of %d.", ell, maxell_ + 1);
        
    if (k < 1 or k > expansions_[ell].cols())
        throw exception ("[LaguerreBasis::size] Can't access orbital %d of %d.", k, expansions_[ell].rows());
    
    return sqrt(lambda_[ell]/gsl_sf_poch(k, 2 * ell + 2)) * pow(lambda_[ell] * r, ell + 1)
           * exp(-lambda_[ell] * r / 2) * gsl_sf_laguerre_n(k - 1, 2 * ell + 1, lambda_[ell] * r);
}

int LaguerreBasis::size (int ell) const
{
    if (ell < 0)
        return maxell_ + 1;
    
    if (ell > maxell_)
        throw exception ("[LaguerreBasis::size] Can't access basis %d of %d.", ell, maxell_ + 1);
    
    return expansions_[ell].rows();
}

const RowMatrix< double >& LaguerreBasis::matrix (int ell) const
{
    if (ell < 0 or ell > maxell_)
        throw exception ("[LaguerreBasis::size] Can't access basis %d of %d.", ell, maxell_ + 1);
    
    return expansions_[ell];
}

const rArrayView LaguerreBasis::orbital (int ell, int n) const
{
    if (ell < 0 or ell > maxell_)
        throw exception ("[LaguerreBasis::size] Can't access basis %d of %d.", ell, maxell_ + 1);
    
    if (n < 1 or n > expansions_[ell].rows())
        throw exception ("[LaguerreBasis::size] Can't access orbital %d of %d.", n, expansions_[ell].rows());
    
    return expansions_[ell].row(n-1);
}

double LaguerreBasis::energy (int ell, int n) const
{
    if (ell < 0 or ell > maxell_)
        throw exception ("[LaguerreBasis::size] Can't access basis %d of %d.", ell, maxell_ + 1);
    
    if (n < 1 or n > expansions_[ell].cols())
        throw exception ("[LaguerreBasis::size] Can't access orbital %d of %d.", n, expansions_[ell].rows());
    
    return energies_[ell][n-1];
}
