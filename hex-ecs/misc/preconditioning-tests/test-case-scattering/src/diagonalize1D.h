
    cArray e (N);
    ColMatrix<Complex> dH = (0.5_z * bD - bM).torow().T();
    ColMatrix<Complex> dS = bS.torow().T();
    ColMatrix<Complex> V (N,N);
    
    // diagonalize the one-electron hamiltonian
    dH.diagonalize_g(dS, e, nullptr, &V);
    
    // write eigenvalues
    rArray eigs = realpart(e);
    std::sort(eigs.begin(), eigs.end());
    write_array(eigs, "1D-eigenvalues.txt");
    
    // normalize eigenvectors
    for (int i = 0; i < N; i++)
        V.col(i) /= std::sqrt(V.col(i) | bS | V.col(i));
    
    // also setup transpose
    RowMatrixView<Complex> Vt (N, N, V.data());
