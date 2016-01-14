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

#include "ProjPreconditioner.h"

const std::string ProjCGPreconditioner::prec_name = "proj";
const std::string ProjCGPreconditioner::prec_description = 
    "Block inversion using conjugate gradients preconditioned by channel projection.";

void ProjCGPreconditioner::setup ()
{
    CGPreconditioner::setup();
    
    std::cout << "Set up channel projection preconditioner" << std::endl;
    
    unsigned Nspline = bspline_atom_.Nspline();
    Timer timer;
    cArray D;
    ColMatrix<Complex> S = rad_.S_atom().torow().T(), CR, invCR, invsqrtS;
    
    std::cout << "\t\t- overlap matrix factorization" << std::endl;
    
    // diagonalize the overlap matrix
    S.diagonalize(D, nullptr, &CR);
    CR.invert(invCR);
    
    // Now S = CR * (D * CR⁻¹)
    std::cout << "\t\t\ttime: " << timer.nice_time() << std::endl;
    for (std::size_t i = 0; i < Nspline * Nspline; i++)
        invCR.data()[i] *= D[i % Nspline];
    std::cout << "\t\t\tresidual: " << (S - CR * invCR).data().norm() << std::endl;
    S = ColMatrix<Complex>();
    
    // compute √S⁻¹
    for (std::size_t i = 0; i < Nspline * Nspline; i++)
        invCR.data()[i] /= std::pow(D.data()[i % Nspline], 1.5);
    invsqrtS = std::move(CR * invCR);
    
    bound_states_.resize(inp_.maxell + 1);
    for (int l = 0; l <= inp_.maxell; l++)
    {
        // reset timer
        std::cout << "\t\t- one-electron Hamiltonian factorization (l = " << l << ")" << std::endl;
        timer.reset();
        
        // compose the symmetrical one-electron hamiltonian
        ColMatrix<Complex> tHl = (Complex(0.5) * rad_.D_atom() - rad_.Mm1_tr_atom() + Complex(0.5*l*(l+1)) * rad_.Mm2_atom()).torow().T();
        
        // symmetrically transform by inverse square root of the overlap matrix
        tHl = std::move(invsqrtS * tHl * invsqrtS);
        
        // diagonalize the transformed matrix
        tHl.diagonalize(D, nullptr, &CR);
        CR.invert(invCR);
        
        // analyze bound states for this angular momentum
        int maxn = 0;
        for (unsigned i = 0; i < Nspline; i++)
        {
            Complex E = D[i];
            if (E.real() < 0)
            {
                // get closest principal quantum number
                int n = std::floor(0.5 + 1.0 / std::sqrt(-2.0 * E.real()));
                
                if (n > 0)
                {
                    // get theoretical energy
                    double E0 = -0.5 / (n * n);
                    
                    // add this bound state to the data
                    if ((int)bound_states_[l].size() < n - l)
                        bound_states_[l].resize(n - l);
                    bound_states_[l][n - l - 1] = invsqrtS * CR.col(i);
                    
                    // if the energy is close enough to the theoretical value, raise number of well-represented bound states
                    if (std::abs(E0 - E) < 1e-3 * std::abs(E0))
                        maxn = std::max(maxn, n);
                }
            }
        }
        
        // Now Hl = ClR * D * ClR⁻¹
        std::cout << "\t\t\t- time: " << timer.nice_time() << std::endl;
        for (std::size_t i = 0; i < Nspline * Nspline; i++)
            invCR.data()[i] *= D[i % Nspline];
        std::cout << "\t\t\t- residual: " << (tHl - CR * invCR).data().norm() << std::endl;
        std::cout << "\t\t\t- bound states with energy within 0.1 % from exact value: " << l + 1 << " <= n <= " << maxn << std::endl;
    }
}

void ProjCGPreconditioner::CG_init (int iblock) const
{
    CGPreconditioner::CG_init(iblock);
}

void ProjCGPreconditioner::CG_mmul (int iblock, const cArrayView r, cArrayView z) const
{
    CGPreconditioner::CG_mmul(iblock, r, z);
}

void ProjCGPreconditioner::CG_prec (int iblock, const cArrayView r, cArrayView z) const
{
    int Nspline_atom = rad_.bspline_atom().Nspline();
    int Nspline_proj = rad_.bspline_proj().Nspline();
    int order = rad_.bspline_atom().order();
    
    // number of opened bound state channels
    int m = 0;
    for (unsigned n = 1; n <= bound_states_[0].size(); n++)
    {
        if (-0.5 / (n*n) < E_)
            m++;
    }
    
    /// DEBUG : write bound states
    /*for (int i = 0; i < m; i++)
    {
        write_array(bound_states_[0][i], format("bound-%ds", i + 1).c_str());
        rArray grid = linspace(0., bspline_atom_.Rmax(), 500);
        write_array(grid, bspline_atom_.zip(bound_states_[0][i], grid), format("bound-%ds-zip", i + 1).c_str());
    }*/
    /// ----
    
    // 1. Project the residual 'r' onto the bound state channels,
    //        r = sum_(l,i) a_(l,i) phi_(l,i) + sum_(l',j) phi_(l',j) b_(l',j) 
    
    RowMatrixView<Complex> R (Nspline_atom, Nspline_proj, r);
    cArray ab (2 * m * Nspline_atom);
    for (int i = 0; i < m; i++)
    {
        // a = r S psi
        cArrayView ai (ab, i * Nspline_atom, Nspline_atom);
        ai = R * rad_.S_atom().dot(bound_states_[0][i]);
        
        // b = psi^T S r
        cArrayView bi (ab, (m + i) * Nspline_atom, Nspline_atom);
        bi = rad_.S_proj().dot(bound_states_[0][i]) * R;
        
        // orthogonalize projections to the bound orbitals
        for (int j = 0; j < m; j++)
        {
            // a = a - (a S psi) psi
            ai -= (ai | rad_.S_atom() | bound_states_[0][i]) * bound_states_[0][i];
            
            // b = b - ( b S psi ) psi
            bi -= (bi | rad_.S_proj() | bound_states_[0][i]) * bound_states_[0][i];
        }
    }
    
    /// DEBUG : verify
    cArray v (Nspline_atom * Nspline_proj);
    static int run = 0; run++;
    for (int i = 0; i < m; i++)
    {
        v += outer_product(ab.slice(i * Nspline_atom, (i + 1) * Nspline_atom), bound_states_[0][i]);
        v += outer_product(bound_states_[0][i], ab.slice((i + m) * Nspline_atom, (i + m + 1) * Nspline_atom));
    }
    write_2D_data(Nspline_atom, Nspline_proj, format("r-%d-orig", run).c_str(),  [&](int i, int j) { return r[i * Nspline_proj + j].real(); });
    write_2D_data(Nspline_atom, Nspline_proj, format("r-%d-verif", run).c_str(), [&](int i, int j) { return v[i * Nspline_proj + j].real(); });
    /// ----
    
    // 2. Construct the bound state channel matrix 'A'.
    SymBandMatrix<Complex> U (m * Nspline_atom, order + 1);
    U.populate
    (
        [&](int irow, int icol) -> Complex
        {
            unsigned i = irow / Nspline_atom, a = irow % Nspline_atom;
            unsigned j = icol / Nspline_atom, b = icol % Nspline_atom;
            
            Complex elem = 0;
            
            if (i == j)
            {
                unsigned n2i = i + 1;
                double E2i = -0.5 / (n2i * n2i);
                elem = (E_ - E2i) * rad_.S_atom()(a,b) - 0.5 * rad_.D_atom()(a,b) /*- 0.5 * l1 * (l1 + 1) * rad_.Mm2_atom()(a,b)*/ + rad_.Mm1_atom()(a,b);
            }
            
            for (int lambda = 0; lambda <= rad_.maxlambda(); lambda++) if (ang_.f(iblock,iblock,lambda) != 0)
                elem -= ang_.f(iblock,iblock,lambda) * (bound_states_[0][i] | rad_.R_tr_dia(lambda)(a,b) | bound_states_[0][j]);
            
            return elem;
        }
    );
    CooMatrix<LU_int_t,Complex> U_coo = U.tocoo<LU_int_t>();
    SymBandMatrix<Complex> V (m * Nspline_atom, order + 1);
    V.populate
    (
        [&](int irow, int icol) -> Complex
        {
            unsigned i = irow / Nspline_atom, a = irow % Nspline_atom;
            unsigned j = icol / Nspline_atom, b = icol % Nspline_atom;
            
            Complex elem = 0;
            
            if (i == j)
            {
                unsigned n1i = i + 1;
                double E1i = -0.5 / (n1i * n1i);
                elem = (E_ - E1i) * rad_.S_atom()(a,b) - 0.5 * rad_.D_atom()(a,b) /*- 0.5 * l2 * (l2 + 1) * rad_.Mm2_atom()(a,b)*/ + rad_.Mm1_atom()(a,b);
            }
            
            for (int lambda = 0; lambda <= rad_.maxlambda(); lambda++) if (ang_.f(iblock,iblock,lambda) != 0)
                elem -= ang_.f(iblock,iblock,lambda) * (bound_states_[0][i] | rad_.R_tr_dia(lambda)(a,b) | bound_states_[0][j]);
            
            return elem;
        }
    );
    CooMatrix<LU_int_t,Complex> V_coo = V.tocoo<LU_int_t>();
    SymBandMatrix<Complex> W (m * Nspline_atom, order + 1);
    W.populate
    (
        [&](int irow, int icol) -> Complex
        {
            unsigned i = irow / Nspline_atom, a = irow % Nspline_atom;
            unsigned j = icol / Nspline_atom, b = icol % Nspline_atom;
            
            Complex elem = 0;
            
            for (int lambda = 0; lambda <= rad_.maxlambda(); lambda++) if (ang_.f(iblock,iblock,lambda) != 0)
            {
                Complex el = 0;
                
                for (int u = 0; u < Nspline_atom; u++)
                for (int v = 0; v < Nspline_proj; v++)
                    el += bound_states_[0][i][u] * rad_.R_tr_dia(lambda)(a,b,u,v) * bound_states_[0][j][v];
                
                elem -= ang_.f(iblock,iblock,lambda) * el;
            }
            
            return elem;
        }
    );
    CooMatrix<LU_int_t,Complex> W_coo = W.tocoo<LU_int_t>();
    NumberArray<LU_int_t> I_arr = concatenate(U_coo.i(), W_coo.i(), W_coo.i() + m * Nspline_atom, V_coo.i() + m * Nspline_atom);
    NumberArray<LU_int_t> J_arr = concatenate(U_coo.j(), W_coo.j() + m * Nspline_atom, W_coo.j(), V_coo.j() + m * Nspline_atom);
    NumberArray<Complex>  V_arr = concatenate(U_coo.v(), W_coo.v(), W_coo.v(), V_coo.v());
    CooMatrix<LU_int_t,Complex> A_coo
    (
        2 * m * Nspline_atom,
        2 * m * Nspline_proj,
        I_arr, J_arr, V_arr
    );
    CsrMatrix<LU_int_t,Complex> A_csr = A_coo.tocsr();
    
    // 3. Solve the matrix.
    std::shared_ptr<LUft<LU_int_t,Complex>> lu = A_csr.factorize();
    cArray uv = lu->solve(ab);
    
    /// DEBUG : check solution
    
    /// ----
    
    // 4. Reconstruct the error vector 'z'.
    z.fill(0.);
    for (int i = 0; i < m; i++)
    {
        cArrayView ui (uv, i * Nspline_atom, Nspline_proj);
        cArrayView vi (uv, (m + i) * Nspline_atom, Nspline_proj);
        z += outer_product(ui, bound_states_[0][i]);
        z += outer_product(bound_states_[0][i], vi);
    }
    
    std::exit(0);
}

void ProjCGPreconditioner::CG_exit (int iblock) const
{
    CGPreconditioner::CG_exit(iblock);
}

void ProjCGPreconditioner::finish ()
{
    CGPreconditioner::finish();
}
