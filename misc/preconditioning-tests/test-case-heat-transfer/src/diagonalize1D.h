
    cArray e (N);
    ColMatrix<Complex> dD = bD.torow().T();
    ColMatrix<Complex> dS = bS.torow().T();
    ColMatrix<Complex> C (N,N), invC (N,N);
    
    // diagonalize the overlap matrix
    dS.diagonalize(e, nullptr, &C);
    C.invert(invC);
    
    // calculate inverse square root of S
    ColMatrix<Complex> invsqrtS (N,N);
    for (unsigned i = 0; i < N; i++)
        C.col(i) *= 1.0 / std::sqrt(e[i]);
    blas::gemm(1.0, C, invC, 0.0, invsqrtS);
   
    // transform hamiltonian symmetrically by √S⁻¹
    blas::gemm(1., invsqrtS, dD, 0., dS);
    blas::gemm(1., dS, invsqrtS, 0., dD);
    
    // diagonalize the one-electron hamiltonian
    dD.diagonalize(e, nullptr, &C);
    C.invert(invC);
    
    // write eigenvalues
    rArray eigs = realpart(e);
    std::sort(eigs.begin(), eigs.end());
    write_array(eigs, "1D-eigenvalues.txt");
    
    // set up KPA preconditioner matrices KPAA = C⁻¹√S⁻¹ and KPAB = √S⁻¹C
    ColMatrix<Complex> KPAA (N,N), KPAB (N,N);
    blas::gemm(1.0, invC,     invsqrtS, 0.0, KPAA);
    blas::gemm(1.0, invsqrtS, C,        0.0, KPAB);
    
    // orthonormalize C's columns
    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
            C.col(j) -= (C.col(i) | bS | C.col(j)) / (C.col(i) | bS | C.col(i)) * C.col(i);
        
        C.col(i) /= C.col(i).norm();
        C.col(i) = bS.dot(C.col(i));
    }
    
    
